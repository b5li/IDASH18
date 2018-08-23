/*
 * @file       CipherPvalues.cpp, cpp file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#ifdef USE_NTL
#include <NTL/BasicThreadPool.h>
#endif

#include "threadpool.h"


#include "../src/Ciphertext.h"
#include "../src/Context.h"
#include "../src/Scheme.h"
#include "../src/SecretKey.h"
#include "../src/TimeUtils.h"
#include "../src/EvaluatorUtils.h"
#include "../src/SchemeAlgo.h"

#include "Database.h"
#include "BasicTest.h"
#include "TestPvals.h"
#include "CipherPvals.h"




//! return the covariance matrix of size d * d



//!@ encrypt an input value and generate a fully-packed ciphertext
//!@ Output: E(data, ..., data)
void CipherPvals::encValue(Ciphertext& encData, double data,  long nslots, long L){
    complex<double>* cmsg = new complex<double>[nslots];
    TP_EXEC_RANGE(nslots, first, last);
    for(int l = first; l < last; ++l){
        cmsg[l].real(data);
    }
    TP_EXEC_RANGE_END;
    encData = scheme.encrypt(cmsg, nslots, L);
    delete[] cmsg;
}

void CipherPvals::encFullyPackedVec(Ciphertext& encData, double* data, long nslots, long L){
    complex<double>* cmsg = new complex<double>[nslots];
    TP_EXEC_RANGE(nslots, first, last);
    for(int l = first; l < last; ++l){
        cmsg[l].real(data[l]);
    }
    TP_EXEC_RANGE_END;
    encData = scheme.encrypt(cmsg, nslots, L);
    delete[] cmsg;
}

//!@ Input: data with a length "len"
//! encrypt an input vector

void CipherPvals::encSparselyPackedVec(Ciphertext& encData, double* data, long len, long nslots, long L){
    complex<double>* cmsg = new complex<double>[nslots];
    TP_EXEC_RANGE(len, first, last);
    for(int l = first; l < last; ++l){
        cmsg[l].real(data[l]);
    }
    TP_EXEC_RANGE_END;
    encData = scheme.encrypt(cmsg, nslots, L);
    delete[] cmsg;
}



void CipherPvals::decSingleData(double& Data, Ciphertext encData){
    complex<double>* dcw = scheme.decrypt(secretKey, encData);
    Data = dcw[0].real();
}

void CipherPvals::decVector(double*& Data, Ciphertext encData, long len){
    complex<double>* dcw = scheme.decrypt(secretKey, encData);
    Data = new double[len];
    for (long j = 0; j < len; ++j) {
        Data[j] = dcw[j].real();
    }
    delete[] dcw;
}

//! nterms = 10
void CipherPvals::computeCov(double*& covariance, double* data, long dim, long nterms, long scalefactor){
    long k = 0;
    for(long i = 0; i < dim; ++i){
        for(long j = i; j < dim; ++j){
            covariance[k] = (data[i] * data[j])/scalefactor;
            k++;
        }
    }
}

 

//! encYXData[i][k]: = Enc(y[i] * xDta[i][k]), 0 <= i <n
//! encCovData[i][k] = Enc(cov_i[k]) (10s)

void CipherPvals::encryptXData(Ciphertext**& encYXData, Ciphertext**& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long dim, long nterms, long scalefactor, long nslots, long L) {

    double** cov = new double*[sampleDim];
    double** xiData = new double*[sampleDim];
    TP_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        if(yData[i] == 1){
            for(long k = 0; k < factorDim; ++k){
                encValue(encYXData[i][k], xData[i][k], nslots, L);
            }
        }
        else{
            for(long k = 0; k < factorDim; ++k){
                encValue(encYXData[i][k], - xData[i][k], nslots, L);
            }
        }
        
        xiData[i] = new double[dim];
        for(long k = 0; k < factorDim; ++k){
            xiData[i][k] = xData[i][k];
        }
        for(long k = factorDim; k < dim; ++k){
            xiData[i][k] = 0.0;
        }
        cov[i] = new double[nterms];
        computeCov(cov[i], xiData[i], dim, nterms, scalefactor);  //! dim -> dim * (dim-1)/2
   
        for(long k = 0; k < nterms; ++k){
            encValue(enccovData[i][k], cov[i][k], nslots, L);
        }
    }
    TP_EXEC_RANGE_END;

    for(long i = 0; i < sampleDim; i++) {
       delete[] xiData[i];
       delete[] cov[i];
    }
    
    delete[] cov;
    delete[] xiData;

}


//!@ EncSData[i] = Enc(s[i][0], ... , s[i][p-1], - 0 -)
//!@ EncYSData[i] = Enc(y[i] * s[i][0], ... , y[i] * s[i][p-1], - 0 -)
//!@ EncSXData[i][k] = Enc(s[i][0] * X[i][k], ..., s[i][p-1] * X[i][k], - 0 -)

void CipherPvals::encryptSData(Ciphertext**& encSData, Ciphertext**& encYSData, Ciphertext***& encSXData,  double* yData, double** xData, double** sData,  long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots, long L) {
    long nslots1 = nsnp - (nencsnp-1) * nslots;   // number of slots in the final ctxts
    
    double** scaled_sData = new double*[sampleDim];
    double*** sxData = new double**[sampleDim];
    double** fullvec = new double*[sampleDim];
    double** sparsevec = new double*[sampleDim];
    
    TP_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        fullvec[i] = new double[nslots];
        sparsevec[i] = new double[nslots1];
        
        //! 1. encryption of sData (because s = 0/1)
        long j1 = 0;
        //! full slot
        for(long j = 0; j < nencsnp - 1; ++j){
            for(long l = 0; l < nslots; ++l){
                fullvec[i][l] = sData[i][j1] ;
                j1++;
            }
            encFullyPackedVec(encSData[i][j], fullvec[i], nslots, L);
        }
        //! not full slot
        for(long l = 0; l < nslots1; ++l){
            sparsevec[i][l] = sData[i][j1] ;
            j1++;
        }
        encSparselyPackedVec(encSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, L);
        
        
        //! 1. encryption of YSData = Y * S
        j1 = 0;
        if(yData[i] == 1){
            //! full slot
            for(long j = 0; j < nencsnp - 1; ++j){
                for(long l = 0; l < nslots; ++l){
                    fullvec[i][l] = sData[i][j1] ;
                    j1++;
                }
                encFullyPackedVec(encYSData[i][j], fullvec[i], nslots, L);
            }
            //! not full slot
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = sData[i][j1] ;
                j1++;
            }
            encSparselyPackedVec(encYSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, L);
        }
        else{
            for(long j = 0; j < nencsnp - 1; ++j){
                for(long l = 0; l < nslots; ++l){
                    fullvec[i][l] = - sData[i][j1];
                    j1++;
                }
                encFullyPackedVec(encYSData[i][j], fullvec[i], nslots, L);
            }
            //! not full slot
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = - sData[i][j1];
                j1++;
            }
            encSparselyPackedVec(encYSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, L);
        }
        
        //! 2. encryption of sxData
        sxData[i] = new double*[factorDim];
        for(long k = 0; k < factorDim; ++k){
            sxData[i][k] = new double[nsnp];
            for(long j = 0; j < nsnp; ++j){
                sxData[i][k][j] = sData[i][j] * xData[i][k];
            }
            j1 = 0;
            //! full slot
            for(long j = 0; j < nencsnp - 1; ++j){
                for(long l = 0; l < nslots; ++l){
                    fullvec[i][l] = sxData[i][k][j1];
                    j1++;
                }
                encFullyPackedVec(encSXData[i][k][j], fullvec[i], nslots, L);
            }
            //! not full slot
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = sxData[i][k][j1];
                j1++;
            }
            encSparselyPackedVec(encSXData[i][k][nencsnp-1], sparsevec[i], nslots1, nslots, L);
        }
    }
    TP_EXEC_RANGE_END;
    
    assert(!IDASH::isThreadPoolActive()&&"has active thread!");

    for(long i = 0; i < sampleDim; i++) {
       for(long k = 0; k < factorDim; k++) {
          delete[] sxData[i][k];
       }
       delete[] sxData[i];
       delete[] fullvec[i];
       delete[] sparsevec[i];
    }

    delete[] scaled_sData;
    delete[] sxData;
    // delete[] fullvec;
    // delete[] sparsevec;
}


//! sum_{i,j} encMatrix[i][j] * (encData1[i] * encData2[j]) : 2*k^2 HM (2 levels)
void CipherPvals::QuadForm(Ciphertext& res, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim){
    Ciphertext** tensoring = new Ciphertext*[factorDim];
    
    // 0: 0 1 2 3
    // 1: 1 4 5 6
    // 2: 2 5 7 8
    // 3: 3 6 8 9
    
    TP_EXEC_RANGE(factorDim, first, last);
    for(long i = first; i< last; ++i){
        tensoring[i] = new Ciphertext[factorDim];
        
        //! (encData1[i] * encData2[j])
        for(long j = 0; j < factorDim; ++j){
            if(encData2[j].l < encData1[i].l){
                tensoring[i][j] = scheme.modDownTo(encData1[i], encData2[j].l);
                scheme.multAndEqual(tensoring[i][j], encData2[j]);
            } else{
                tensoring[i][j] = scheme.modDownTo(encData2[j], encData1[i].l);
                scheme.multAndEqual(tensoring[i][j], encData1[i]);
            }
            scheme.reScaleByAndEqual(tensoring[i][j], 1);
        }
        
        //! (encData1[i] * encData2[j]) * encMatrix[i][j]
        //! lower
        for(long j = 0; j < i; ++j){
            long l1 = j* factorDim - (j * (j-1))/2 + (i - j);   // l = j* factorDim - j*(j-1)/2 + i - j;
            Ciphertext ctemp = encMatrix[l1];
            if(ctemp.l < tensoring[i][j].l){
                scheme.modDownToAndEqual(tensoring[i][j], ctemp.l);
            } else{
                scheme.modDownToAndEqual(ctemp, tensoring[i][j].l);
            }
            scheme.multAndEqual(tensoring[i][j], ctemp);
            //scheme.reScaleByAndEqual(tensoring[i][j], 1);
        }
        
        //! upper
        long diag_index = i* factorDim - (i * (i-1))/2;
        for(long j = i; j < factorDim; ++j){
            long l1 = diag_index + (j - i);   // l = i* factorDim - i*(i-1)/2 + j - i ;
            Ciphertext ctemp = encMatrix[l1];
            if(ctemp.l < tensoring[i][j].l){
                scheme.modDownToAndEqual(tensoring[i][j], ctemp.l);
            } else{
                scheme.modDownToAndEqual(ctemp, tensoring[i][j].l);
            }
            scheme.multAndEqual(tensoring[i][j], ctemp);
            //scheme.reScaleByAndEqual(tensoring[i][j], 1);
        }
        for(long j = 1; j <factorDim; ++j){
            scheme.addAndEqual(tensoring[i][0], tensoring[i][j]);
        }
    }
    TP_EXEC_RANGE_END;
    
    res = tensoring[0][0];
    for(long i = 1; i < factorDim; ++i){
        scheme.addAndEqual(res, tensoring[i][0]);
    }
    scheme.reScaleByAndEqual(res, 1);
}


//! sum_{i,j} encMatrix[i][j] * (encData1[i] * encData2[j]) : 2* (k + (k2 - k)/2) HM
//! lvl: pBits (1 level)
void CipherPvals::SqrQuadForm(Ciphertext& res, Ciphertext* encData, Ciphertext* encMatrix, long factorDim){
    Ciphertext** tensoring = new Ciphertext*[factorDim];
    
    TP_EXEC_RANGE(factorDim, first, last);
    for(int i = first; i< last; ++i){
        tensoring[i] = new Ciphertext[factorDim];
        
        //! Firstly, compute the diagonal (using squaring)
        tensoring[i][i] = scheme.square(encData[i]);
        scheme.reScaleByAndEqual(tensoring[i][i], 1);
    
        long diag_index = i* factorDim - (i * (i-1))/2;
        Ciphertext ctemp = encMatrix[diag_index];
        
        if(ctemp.l < tensoring[i][i].l){
            scheme.modDownToAndEqual(tensoring[i][i], ctemp.l);
        } else{
            scheme.modDownToAndEqual(ctemp, tensoring[i][i].l);
        }
    
        scheme.multAndEqual(tensoring[i][i], ctemp);
        //scheme.reScaleByAndEqual(tensoring[i][i], 1); //

        for(long j = i + 1; j < factorDim; ++j){
            long l1 = diag_index + (j - i);
            tensoring[i][j] = encData[i];
            scheme.multAndEqual(tensoring[i][j], encData[j]);
            scheme.reScaleByAndEqual(tensoring[i][j], 1);
            
            Ciphertext ctemp = encMatrix[l1];
            if(ctemp.l < tensoring[i][j].l){
                scheme.modDownToAndEqual(tensoring[i][j], ctemp.l);
            } else{
                scheme.modDownToAndEqual(ctemp, tensoring[i][j].l);
            }
            scheme.multAndEqual(tensoring[i][j], ctemp);
            //scheme.reScaleByAndEqual(tensoring[i][j], 1);
        }
        
#if defined(__DEBUG_)
        cout << i << ": " ;
        for(long j = i; j < factorDim; ++j){
            double res1;
            decSingleData(res1, tensoring[i][j]);
            cout <<  res1 << ", " ;
        }
        cout << endl;
#endif
        
        for(long j = i + 2; j <factorDim; ++j){
            scheme.addAndEqual(tensoring[i][i+1], tensoring[i][j]);
        }
    }
    TP_EXEC_RANGE_END;
    
    TP_EXEC_RANGE(factorDim - 1, first, last);
    for(int i = first; i< last; ++i){
        scheme.addAndEqual(tensoring[i][i], tensoring[i][i+1]);
        scheme.addAndEqual(tensoring[i][i], tensoring[i][i+1]);
    }
    TP_EXEC_RANGE_END;
    
    res = tensoring[0][0];
    for(long i = 1; i < factorDim; ++i){
        scheme.addAndEqual(res, tensoring[i][i]);
    }
    
    scheme.reScaleByAndEqual(res, 1);
    
}




void CipherPvals::HesInverse(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext** enccovData, long dim, long nslots, long L){
    double adjoint[4][4]= {
        {1741.836911,    -1305.78448,    -779.6744328,    -1808.451025},
        {-1305.755823,    3375.054075,    -8.241993338,    106.974673},
        {-779.6744328,    -8.241993338,    5899.023691,    -3124.853411},
        {-1808.451025,    106.974673,    -3124.853411,    6171.985412},
    };
    double det = 19152.84924;
    double scalefactor = 32768;  // power(32,3)
    
    encValue(encDet, det/(scalefactor*32), nslots, L - 2);
    
    long k = 0;
    for(long i = 0; i < dim; ++i){
        for(long j = i; j < dim; ++j){
            encValue(encAdj[k], adjoint[i][j]/scalefactor, nslots, L - 2);
            k++;
        }
    }
}

//!@ Intput encData[nterms], encryption of covariant matrix
//!@ Output: endAdj[nterms], encrpytion of adjoint matrix
//!@         encDet, encryption of determinant

void CipherPvals::encAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* encData){
    double table[20][2] = {
        {1, 5},{1, 7},{1, 8},{1, 9},
        {2, 4},{2, 5},{2, 6},{2, 9},
        {3, 5},{3, 6},{3, 8},
        {4, 7},{4, 8},{4, 9},
        {5, 6},{5, 8},{5, 9},
        {6, 7},{6, 8},{7, 9},
    };
    
    //! 1. squaring and multiplication
    Ciphertext* sqrtemp = new Ciphertext[4];
    Ciphertext* temp = new Ciphertext[20];
    
    sqrtemp[0] = scheme.square(encData[3]); //! {{0,3} -> (x[0][3])^2
    sqrtemp[1] = scheme.square(encData[5]); //! {1,2},
    sqrtemp[2] = scheme.square(encData[6]); //! {1,3}
    sqrtemp[3] = scheme.square(encData[8]); //! {2,3}
    
    TP_EXEC_RANGE(20, first, last);
    for(long i = first; i < last; ++i){
        long j0= table[i][0];
        long j1= table[i][1];    
        temp[i] = scheme.mult(encData[j0], encData[j1]);
    }
    TP_EXEC_RANGE_END;
    
    //! 2. two polynomials mult
    Ciphertext* adj = new Ciphertext[30];
    
    adj[0] = scheme.sub(temp[19], sqrtemp[3]);
    adj[1] = scheme.sub(temp[18], temp[16]);
    adj[2] = sqrtemp[2];
    
    adj[3] = scheme.sub(sqrtemp[3], temp[19]);
    adj[4] = scheme.negate(adj[1]);
    adj[5] = scheme.sub(temp[17], temp[15]);
    
    adj[6] = adj[4];
    adj[7] = scheme.sub(sqrtemp[2], temp[13]);
    adj[8] = scheme.sub(temp[12], temp[14]);
    
    adj[9]  = adj[5];
    adj[10] = adj[8];
    adj[11] = scheme.sub(sqrtemp[1], temp[11]);
    
    adj[12] = scheme.negate(adj[3]);
    adj[13] = scheme.sub(temp[10], temp[7]);
    adj[14] = sqrtemp[0];
    
    adj[15] = adj[1];
    adj[16] = scheme.sub(temp[3], temp[9]);
    adj[17] = scheme.sub(temp[8], temp[2]);
    
    adj[18] = scheme.negate(adj[5]);
    adj[19] = scheme.sub(temp[6], temp[2]);
    adj[20] = scheme.sub(temp[1], temp[5]);
    
    adj[21] = scheme.negate(adj[7]);
    adj[22] = scheme.negate(adj[16]);
    adj[23] = sqrtemp[0];
    
    adj[24] = scheme.negate(adj[8]);
    adj[25] = scheme.negate(adj[19]);
    adj[26] = scheme.sub(temp[4], temp[0]);
    
    adj[27] = scheme.sub(temp[11], sqrtemp[1]);
    adj[28] = scheme.negate(adj[20]);
    adj[29] = temp[4];
    
    Ciphertext* xtemp = new Ciphertext[8];
    for(long i = 0; i < 8; ++i){
        xtemp[i] = encData[i];
        //xtemp[i] = scheme.modDownTo(encData[i], adj[0].l);
    }
    
    
    double diag[4][5] = {
        {0, 4, 5, 18, 7}, // 0th
        {4, 0, 2, 10, 7}, // 4th
        {7, 0, 1, 9, 4},  // 7th
        {9, 0, 1, 5, 2},  // 9th
    };
    
    double nondiag[6][3] = {
        {1, 1, 2},
        {2, 1, 2},
        {3, 1, 2},
        {5, 0, 2},
        {6, 0, 2},
        {8, 0, 1},
    };
    
    //! diagonal computation
    TP_EXEC_RANGE(4, first, last);
    for(long i = first; i < last; ++i){
        long k = 3 * diag[i][0];
        long k0 = diag[i][0];
        long k1 = diag[i][1];
        long k2 = diag[i][2];
        long k3 = diag[i][3];
        long k4 = diag[i][4];
        
        encAdj[k0] = scheme.mult(xtemp[k1], adj[k]);
        Ciphertext ctemp = scheme.add(adj[k + 1], temp[k3]);
        scheme.multAndEqual(ctemp, xtemp[k2]);
        scheme.addAndEqual(encAdj[k0], ctemp);
        ctemp = scheme.mult(xtemp[k4], adj[k + 2]);
        scheme.subAndEqual(encAdj[k0], ctemp);
    }
    TP_EXEC_RANGE_END;
    
    //! non-diagonal computation: x[k1] * adj[k] + x[k2] * adj[k + 1] + x[3] * adj[k + 2]
    TP_EXEC_RANGE(6, first, last);
    for(long i = first; i < last; ++i){
        long k = 3 * nondiag[i][0];
        long k0 = nondiag[i][0];
        long k1 = nondiag[i][1];
        long k2 = nondiag[i][2];
        
        encAdj[k0] = scheme.mult(xtemp[k1], adj[k]);
        Ciphertext ctemp = scheme.mult(xtemp[k2], adj[k + 1]);
        scheme.addAndEqual(encAdj[k0], ctemp);
        ctemp = scheme.mult(xtemp[3], adj[k + 2]);
        scheme.addAndEqual(encAdj[k0], ctemp);
    }
    TP_EXEC_RANGE_END;
    
    // "+------------------------------------+"
    //               Deternimant
    // "+------------------------------------+"
    
    
    Ciphertext* encDetemp = new Ciphertext[4];
    TP_EXEC_RANGE(4, first, last);
    for(long i = first; i < last; ++i){
        encDetemp[i] = scheme.mult(xtemp[i], encAdj[i]);
    }
    TP_EXEC_RANGE_END;
    encDet = encDetemp[0];
    for(long i = 1; i < 4; ++i){
        scheme.addAndEqual(encDet, encDetemp[i]);
    }
    
    
    //! scale by 2 * logp
    TP_EXEC_RANGE(10, first, last);
    for(long i = first; i < last; ++i){
        scheme.reScaleByAndEqual(encAdj[i], 2);
    }
    TP_EXEC_RANGE_END;
        
    scheme.reScaleByAndEqual(encDet, 3);
    
#if  1
    //defined(__DEBUG_)
    for(long i = 0; i < 10; ++i){
        double data;
        decSingleData(data, encAdj[i]);
        cout << i<< ":" << data << endl;
    }
    
    double data;
    decSingleData(data, encDet);
    cout << "0.0182656 ?= " << data << endl;
#endif
    
    
}

