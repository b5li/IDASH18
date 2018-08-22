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

void CipherPvals::computeCov(double**& covariance, double* data, long dim){
    covariance = new double*[dim];
    for(long i = 0; i < dim; ++i){
        covariance[i] = new double[dim];
        for(long j = 0; j < dim; ++j){
            covariance[i][j] = data[i] * data[j];
        }
    }
}

 

//! encYXData[i][k]: = Enc(y[i] * xDta[i][k]), 0 <= i <n

void CipherPvals::encryptXData(Ciphertext**& encYXData, Ciphertext**& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long dim, long nslots, long L) {

    double*** temp = new double**[sampleDim];
    double** xiData = new double*[sampleDim];

    TP_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        if(yData[i] == 1){
            for(long k = 0; k < factorDim; ++k){
                encValue(encYXData[i][k], xData[i][k], nslots, L);
                //scaled_xData1[i][k] = xData[i][k];
            }
        }
        else{
            for(long k = 0; k < factorDim; ++k){
                encValue(encYXData[i][k], - xData[i][k], nslots, L);
                //scaled_xData1[i][k] = - xData[i][k];
            }
        }
        
        xiData[i] = new double[dim];
        for(long k = 0; k < factorDim; ++k){
            xiData[i][k] = xData[i][k];
        }
        for(long k = factorDim; k < dim; ++k){
            xiData[i][k] = 0.0;
        }
        computeCov(temp[i], xiData[i], dim);  //! dim -> dim * dim
   
        
    }
    TP_EXEC_RANGE_END;
    
    delete[] temp;
    delete[] xiData;
    //delete[] scaled_xData1;
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
                fullvec[i][l] = sData[i][j1];
                j1++;
            }
            encFullyPackedVec(encSData[i][j], fullvec[i], nslots, L);
        }
        //! not full slot
        for(long l = 0; l < nslots1; ++l){
            sparsevec[i][l] = sData[i][j1];
            j1++;
        }
        encSparselyPackedVec(encSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, L);
        
        
        //! 1. encryption of YSData = Y * S
        j1 = 0;
        if(yData[i] == 1){
            //! full slot
            for(long j = 0; j < nencsnp - 1; ++j){
                for(long l = 0; l < nslots; ++l){
                    fullvec[i][l] = sData[i][j1];
                    j1++;
                }
                encFullyPackedVec(encYSData[i][j], fullvec[i], nslots, L);
            }
            //! not full slot
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = sData[i][j1];
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
    
    delete[] scaled_sData;
    delete[] sxData;
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
    double scalefactor = 16384;
    
    encValue(encDet, det/scalefactor, nslots, L - 2);
    
    long k = 0;
    for(long i = 0; i < dim; ++i){
        for(long j = i; j < dim; ++j){
            encValue(encAdj[k], adjoint[i][j]/scalefactor, nslots, L - 2);
            k++;
        }
    }
}


/*
//! encXData[i] = Enc(x[i][0], x[i][1], ... , x[i][k-1])
void CipherPvals::encryptXData(Ciphertext*& encXData, Ciphertext*& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long dim, long nbatching){
    double*** temp = new double**[sampleDim];
    double** scaled_xData = new double*[sampleDim];
    
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        scaled_xData[i] = new double[dim];
        for(long j = 0; j < factorDim; ++j){
            scaled_xData[i][j] = xData[i][j] * (yData[i]);
        }
        for(long j = factorDim; j < dim; ++j){
            scaled_xData[i][j] = 0.0;
        }
        
        matrixform(temp[i], scaled_xData[i], dim);  //! dim -> dim * dim
        encSingleData(encXData[i], temp[i], dim, nbatching);
        
        computeCov(temp[i], scaled_xData[i], dim);  //! dim -> dim * dim
        encSingleData(enccovData[i], temp[i], dim, HEmatpar.nbatching);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] temp;
    delete[] scaled_xData;
}

//! encXData[i] = Enc(x[i][0], x[i][1], ... , x[i][k-1])
void CipherPvals::encryptSData(Ciphertext**& encSData, Ciphertext**& encSXData,  double* yData, double** xData, double** sData,  long factorDim, long sampleDim, long dim, long nsnp, long nencsnp, long nbatching) {
    double**** temp1 = new double***[sampleDim];
    double**** temp2 = new double***[sampleDim];
    
    double*** sxData = new double**[sampleDim];
    double*** scaled_sData = new double**[sampleDim];
    
    long j1 = (nencsnp - 1) * nbatching;
    long nbatching1 = nsnp - j1; // number of slots in the final ctxt
    
    cout << j1 << endl;
    cout << nbatching1 << endl;
    
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        //! encSData
        scaled_sData[i] = new double*[nbatching];
        temp1[i] = new double**[nbatching];
        
        for(long j = 0; j < nencsnp - 1; ++j){
            long j1 = j * nbatching;
            for(long l = 0; l < nbatching; ++l){
                scaled_sData[i][l] = new double[dim];
                double syData = sData[i][l + j1] * (yData[i]);  // scaled by "y"
                for(long k = 0; k < factorDim; ++k){
                    scaled_sData[i][l][k] = syData;
                }
                for(long k = factorDim; k < dim; ++k){
                    scaled_sData[i][l][k] = 0.0;
                }
                matrixform(temp1[i][l], scaled_sData[i][l], dim); //! dim -> dim * dim
            }
            encMultipleData(encSData[i][j], temp1[i], dim, nbatching);
        }

        for(long l = 0; l < nbatching1; ++l){
            scaled_sData[i][l] = new double[dim];
            double syData = sData[i][l + j1] * (yData[i]);
            for(long k = 0; k < factorDim; ++k){
                scaled_sData[i][l][k] = syData;
            }
            for(long k = factorDim; k < dim; ++k){
                scaled_sData[i][l][k] = 0.0;
            }
            matrixform(temp1[i][l], scaled_sData[i][l], dim); //! dim -> dim * dim
        }
        encMultipleData(encSData[i][nencsnp - 1], temp1[i], dim, nbatching1);

        //! encSX
        sxData[i] = new double*[nbatching];
        temp2[i] = new double**[nbatching];
        
        for(long j = 0; j < nencsnp - 1; ++j){
            long j1 = j * nbatching;
            for(long l = 0; l < nbatching; ++l){
                sxData[i][l] = new double[dim];
                for(long k = 0; k < factorDim; ++k){
                    sxData[i][l][k] = sData[i][l + j1] * xData[i][k];
                }
                for(long k = factorDim; k < dim; ++k){
                    sxData[i][l][k] = 0.0;
                }
                matrixform(temp2[i][l], sxData[i][l], dim);
            }
            encMultipleData(encSXData[i][j], temp2[i], dim, nbatching);
        }

        for(long l = 0; l < nbatching1; ++l){
            sxData[i][l] = new double[dim];
            for(long k = 0; k < factorDim; ++k){
                sxData[i][l][k] = sData[i][l + j1] * xData[i][k];
            }
            for(long k = factorDim; k < dim; ++k){
                sxData[i][l][k] = 0.0;
            }
            matrixform(temp2[i][l], sxData[i][l], dim);
        }
        encMultipleData(encSXData[i][nencsnp - 1], temp2[i], dim, nbatching1);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] temp1;
    delete[] temp2;
    delete[] sxData;
    delete[] scaled_sData;
}




//! Input: encData1, encData2, encMatrix, logdim2 = log2(dim^2)
//! Complexity: 2 HM + logdim2 Rot
void CipherPvals::QuadForm(Ciphertext& res, Ciphertext encData1, Ciphertext encMatrix, Ciphertext encData2, long logdim2){
    
    //! res = encData1 * encData2 (tensoring for the corresponding plaintext vectors)
    if(encData2.logq < encData1.logq){
        res = encData1;
        scheme.modDownToAndEqual(res, encData2.logq);
        scheme.multAndEqual(res, encData2);
    } else{
        res = encData2;
        scheme.modDownToAndEqual(res, encData1.logq);
        scheme.multAndEqual(res, encData1);
    }
    scheme.reScaleByAndEqual(res, res.logp);
    
    //! res = encData1 * encMatrix * encData2
    Ciphertext ctemp = encMatrix;
    if(res.logq > ctemp.logq){
        scheme.modDownToAndEqual(res, ctemp.logq);
    } else{
        scheme.modDownToAndEqual(ctemp, res.logq);
    }
    scheme.multAndEqual(res, ctemp);
  
    scheme.reScaleByAndEqual(res, res.logp);
    
    //! Allsum
    for(long i = 0; i < logdim2; ++i){
        long l = (1<<i);
        Ciphertext ctemp = scheme.leftRotate(res, l*HEmatpar.nbatching);
        scheme.addAndEqual(res, ctemp);
    }
}



//! Input: encMatrix (symmetric matrix)
//! complexity: (k^2+K)/2 * (SM + 2logk Rot), lvl: cBits
void CipherPvals::replicate(Ciphertext**& res, Ciphertext encMatrix, long dim, long logdim2, long nbatching){
    res = new Ciphertext*[dim];
    ZZX** poly = new ZZX*[dim];
    
    NTL_EXEC_RANGE(dim, first, last);
    for(int i = first; i< last; ++i){
        res[i] = new Ciphertext[dim];
        poly[i] = new ZZX[dim];
        for(long j = i; j < dim; ++j){
            long start = (dim * i + j) * nbatching;
            long end = start + nbatching;
            complex<double>* pvals = new complex<double>[HEmatpar.nslots];
            for(long l = start; l < end; ++l){
                pvals[l].real(1.0);
            }
            poly[i][j] = scheme.context.encode(pvals, HEmatpar.nslots, HEmatpar.cBits);
            delete[] pvals;
            res[i][j] = scheme.multByPoly(encMatrix, poly[i][j], HEmatpar.cBits);
            scheme.reScaleByAndEqual(res[i][j], HEmatpar.cBits);
            
            for(long k = 0; k < logdim2; ++k){
                long l = (1 << k);
                Ciphertext ctemp = scheme.rightRotate(res[i][j], l * nbatching);
                scheme.addAndEqual(res[i][j], ctemp);
            }
        }
    }
    NTL_EXEC_RANGE_END;
}
 
 //! Input: data of size (dim * dim)
 //! Return Enc(x) = (- x[0][0] -|  - x[0][1] - | ... | - x[d-1][d-1] -)
 //!  each of values is replicated with nbatching times
 
 void CipherPvals::encSingleData(Ciphertext& encData, double** data, long dim, long PoTdim, long nbatching, long nslots, long L){
 complex<double>* cmsg = new complex<double>[nslots];
 NTL_EXEC_RANGE(nbatching, first, last);
 for(int l = first; l < last; ++l){
 for(long i = 0; i < dim; ++i){
 for(long j = 0; j < dim; ++j){
 cmsg[(i*PoTdim + j)* nbatching + l].real(data[i][j]);
 }
 }
 }
 NTL_EXEC_RANGE_END;
 encData = scheme.encrypt(cmsg, nslots, L);
 delete[] cmsg;
 }
*/

