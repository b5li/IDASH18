/*
 * @file       CipherLinRegPvalues.cpp, cpp file
 * @brief      Homomorphic Evaluation of Linear Regression
 *
 * @author     Miran Kim
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#define USE_NTL 1

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

#include "../src/ExtCiphertext.h"
#include "../src/ExtScheme.h"

#include "Database.h"
#include "BasicTest.h"
#include "TestLinRegPvals.h"
#include "CipherLinRegPvals.h"


//!@ encrypt an input value and generate a fully-packed ciphertext
//!@ Output: E(data, ..., data)
void CipherPvals::encValue(Ciphertext& encData, double data,  long nslots, long L){
    complex<double>* cmsg = new complex<double>[nslots];
    NTL_EXEC_RANGE(nslots, first, last);
    for(int l = first; l < last; ++l){
        cmsg[l].real(data);
    }
    NTL_EXEC_RANGE_END;
    encData = scheme.encrypt(cmsg, nslots, L);
    delete[] cmsg;
}

void CipherPvals::encFullyPackedVec(Ciphertext& encData, double* data, long nslots, long L){
    complex<double>* cmsg = new complex<double>[nslots];
    NTL_EXEC_RANGE(nslots, first, last);
    for(int l = first; l < last; ++l){
        cmsg[l].real(data[l]);
    }
    NTL_EXEC_RANGE_END;
    encData = scheme.encrypt(cmsg, nslots, L);
    delete[] cmsg;
}

//!@ Input: data with a length "len"
//! encrypt an input vector

void CipherPvals::encSparselyPackedVec(Ciphertext& encData, double* data, long len, long nslots, long L){
    complex<double>* cmsg = new complex<double>[nslots];
    NTL_EXEC_RANGE(len, first, last);
    for(int l = first; l < last; ++l){
        cmsg[l].real(data[l]);
    }
    NTL_EXEC_RANGE_END;
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
void CipherPvals::computeCov(double*& covariance, double* data, long dim, long scaleBits){
    long k = 0;
    long scalefactor = (1 << scaleBits);
    for(long i = 0; i < dim; ++i){
        for(long j = i; j < dim; ++j){
            covariance[k] = (data[i] * data[j]) / scalefactor;
            k++;
        }
    }
}

/********************************************************************/

//!@ EncSData[i] = Enc(s[i][0], ... , s[i][p-1], - 0 -)
//!@ EncYSData[i] = Enc(y[i] * s[i][0], ... , y[i] * s[i][p-1], - 0 -)
//!@ EncSXData[i][k] = Enc(s[i][0] * X[i][k], ..., s[i][p-1] * X[i][k], - 0 -)

void CipherPvals::encryptSData(Ciphertext**& encSData, Ciphertext**& encYSData, Ciphertext***& encSXData,  double* yData, double** xData, double** sData,  long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots, long L) {
    long nslots1 = nsnp - (nencsnp-1) * nslots;   // number of slots in the final ctxts
    
    double** scaled_sData = new double*[sampleDim];
    double*** sxData = new double**[sampleDim];
    double** fullvec = new double*[sampleDim];
    double** sparsevec = new double*[sampleDim];
    
    NTL_EXEC_RANGE(sampleDim, first, last);
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
            encFullyPackedVec(encSData[i][j], fullvec[i], nslots, L - 2);
        }
        //! not full slot
        for(long l = 0; l < nslots1; ++l){
            sparsevec[i][l] = sData[i][j1] ;
            j1++;
        }
        encSparselyPackedVec(encSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, L - 2);
        
        //! 1. encryption of YSData = Y * S
        j1 = 0;
        if(yData[i] == 1){
            //! full slot
            for(long j = 0; j < nencsnp - 1; ++j){
                for(long l = 0; l < nslots; ++l){
                    fullvec[i][l] = sData[i][j1] ;
                    j1++;
                }
                encFullyPackedVec(encYSData[i][j], fullvec[i], nslots, L - 2);
            }
            //! not full slot
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = sData[i][j1] ;
                j1++;
            }
            encSparselyPackedVec(encYSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, L - 2);
        }
        else{
            for(long j = 0; j < nencsnp - 1; ++j){
                for(long l = 0; l < nslots; ++l){
                    fullvec[i][l] = - sData[i][j1];
                    j1++;
                }
                encFullyPackedVec(encYSData[i][j], fullvec[i], nslots, L - 2);
            }
            //! not full slot
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = - sData[i][j1];
                j1++;
            }
            encSparselyPackedVec(encYSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, L - 2);
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

            delete [] sxData[i][k];
        }
        delete [] sxData[i];
        delete [] fullvec[i];
        delete [] sparsevec[i];
    }
    NTL_EXEC_RANGE_END;
    
    delete[] scaled_sData;
    delete[] sxData;
    delete[] fullvec;
    delete[] sparsevec;
}

/********************************************************************/

void CipherPvals::encryptSIMDXData(Ciphertext& encYXData, Ciphertext*& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long sampleDim2, long nXbatching, long nCovbatching, long nterms, long scaleBits, long nslots, long L) {
    
    long nslots1 = factorDim * nXbatching;  // total number of slots of a xData[i] = 16
    double* temp = new double[nslots];
    
    // "+------------------------------------+"
    //!  encryption of YXData
    // "+------------------------------------+"
    
    //! encoding of YXData
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        long start = nslots1 * i;
        if(yData[i] == 1){
            for(long j = 0; j < nXbatching; ++j){
                long start1 = start + j * factorDim;
                for(long k = 0; k < factorDim; ++k){
                    temp[start1 + k] = xData[i][k];
                }
            }
        }
        else{
            for(long j = 0; j < nXbatching; ++j){
                long start1 = start + j * factorDim;
                for(long k = 0; k < factorDim; ++k){
                    temp[start1 + k] = - xData[i][k];
                }
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    long nslots2 = nslots1 * sampleDim;  // number of used for encoding data
    for(long i = nslots2; i < nslots; ++i){
        temp[i] = 0.0;
    }
    
    //! encryption
    encFullyPackedVec(encYXData, temp, nslots, L);
    
    // "+------------------------------------+"
    //!  encryption of covariance
    // "+------------------------------------+"
   
    double** cov = new double*[sampleDim];
    double** xiData = new double*[sampleDim];
    
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        cov[i] = new double[nterms];
        computeCov(cov[i], xData[i], factorDim, scaleBits);  //! dim -> dim * (dim-1)/2
    }
    NTL_EXEC_RANGE_END;
    
    //! encoding rule for computing adjoint
    double Index[10][3][6]= {
        {{4,6,5,5,6,4},{7,7,8,5,5,8},{9,6,6,9,8,8}}, //00
        {{3,2,2,3,1,1},{7,8,5,5,8,7},{6,6,9,8,8,9}}, //01
        {{1,3,2,2,3,1},{5,5,6,4,4,6},{9,6,6,9,8,8}}, //02
        {{3,2,2,3,1,1},{5,6,4,4,6,5},{5,5,8,7,7,8}}, //03
        {{0,3,2,2,3,0},{7,7,8,2,2,8},{9,3,3,9,8,8}}, //11
        {{3,2,2,3,0,0},{5,6,1,1,6,5},{3,3,9,8,8,9}}, //12
        {{0,3,2,2,3,0},{5,5,6,1,1,6},{8,2,2,8,7,7}}, //13
        {{0,3,1,1,3,0},{4,4,6,1,1,6},{9,3,3,9,6,6}}, //22
        {{3,1,1,3,0,0},{4,6,1,1,6,4},{2,2,8,5,5,8}}, //23
        {{0,2,1,1,2,0},{4,4,5,1,1,5},{7,2,2,7,5,5}}, //33
    };
    
    long nslots3 = 8 * nCovbatching;        //! total number of slots of a xData[i] = 16
    long nslots4 = nslots3 * sampleDim;     //! number of used for encoding data
    double** xData2 = new double*[30];
    
    //! encoding of adjoint
    NTL_EXEC_RANGE(30, first, last);
    for (long l = first; l < last; ++l){
        xData2[l] = new double[nslots];
        long lrow = (long) l/3;
        long lcol = (l % 3);
        
        for (long i = 0; i < sampleDim; ++i) {
            long start = nslots3 * i;
            for(long j = 0; j < nCovbatching; ++j){
                long start1 = start + j * 8;
                for(long k = 0; k < 6; k += 2){
                    long ind = Index[lrow][lcol][k];
                    xData2[l][start1 + k] = cov[i][ind];
                    ind = Index[lrow][lcol][k + 1];
                    xData2[l][start1 + k + 1] = - cov[i][ind];
                }
                xData2[l][start1 + 6] = 0.0;
                xData2[l][start1 + 7] = 0.0;
            }
        }
        for(long i = nslots4; i < nslots; ++i){
            xData2[l][i] = 0.0;
        }
        encFullyPackedVec(enccovData[l], xData2[l], nslots, L);
    }
    NTL_EXEC_RANGE_END;
    
    //! encoding of determinant
    long scalefactor = (1 << scaleBits);
    NTL_EXEC_RANGE(4, first, last);
    for (long l = first; l < last; ++l){
        for (long i = 0; i <sampleDim; ++i) {
            double cov1 = (cov[i][l] * scalefactor); //! original covariant component
            long start = nslots3 * i;
            for(long j = 0; j < nCovbatching; ++j){
                long start1 = start + j * 8;
                for(long k = 0; k < 8; ++k){
                    xData2[l][start1 + k] = cov1;
                }
            }
        }
        for(long i = nslots4; i < nslots; ++i){
            xData2[l][i] = 0.0;
        }
        encFullyPackedVec(enccovData[30 + l], xData2[l], nslots, L);
    }
    NTL_EXEC_RANGE_END;
}

void CipherPvals::decryptResult(double& Ynorm, double**& YSnorm, double**& Snorm, Ciphertext encYnorm, Ciphertext* encYSnorm, Ciphertext* encSnorm, long nencsnp, long nslots){
    
    decSingleData(Ynorm, encYnorm);
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        YSnorm[j] = new double[nslots];
        Snorm[j] = new double[nslots];
        
        decVector(YSnorm[j], encYSnorm[j], nslots);
        decVector(Snorm[j], encSnorm[j], nslots);
    }
    NTL_EXEC_RANGE_END;
}
/********************************************************************/
//! encYX[k] = encryption of YXvec = y^T * X[k] = [-29, -6.506140997, -6.368385672, -12.88024343] (0 <= k < 4)
//! encYX[k] = E((y^T * X)[k])
//! 8 Rot + 16 Rot
void CipherPvals::aggYXData(Ciphertext*& encYX, Ciphertext encYXData, long sdimBits, long nbatching, long factorDim, long nslots){
    long nslots1 = nbatching * factorDim;  // nbatching = replicated number of a user's data in a single ciphertext
    Ciphertext res = encYXData;
    
    //! Allsum over "sampleDim"
    for(long i = 0; i < sdimBits; ++i){
        Ciphertext ctemp = scheme.leftRotateFast(res, (1 << i) * nslots1); // by 2^i * 16
        scheme.addAndEqual(res, ctemp);
    }
    fullReplicate4(encYX, res, nslots);
}

//! (34 * 8) Rot 
void CipherPvals::aggCovData(Ciphertext*& encCov, Ciphertext* enccovData,  long sdimBits, long nbatching){
    long nslots1 = nbatching * 8;
    NTL_EXEC_RANGE(34, first, last);
    for (long l = first; l < last; ++l){
        encCov[l] = enccovData[l];
        //! Allsum
        for(long i = 0; i < sdimBits; ++i){
            Ciphertext ctemp = scheme.leftRotateFast(encCov[l], (1 << i) * nslots1); // by 2^i * 16
            scheme.addAndEqual(encCov[l], ctemp);
        }
    }
    NTL_EXEC_RANGE_END;
}

//!@ Full replication of vector of size 4  (v0,v1,v2,v3) -> (v0), (v1), (v2), (v3)
//!@ Return ciphertexts of L-1
//!@ (2*8) Rot

void CipherPvals::fullReplicate4(Ciphertext*& res, Ciphertext Data, long nslots){
    //! full replication of size "4"
    complex<double>** pvals = new complex<double>*[4];
    uint64_t** poly = new uint64_t*[4];
    
    NTL_EXEC_RANGE(4, first, last);
    for (long i = first; i < last; ++i) {
        pvals[i] = new complex<double>[nslots];
        for(long j = 0; j < nslots ; j += 4){
            pvals[i][j + i].real(1.0);
        }
        poly[i] = new uint64_t[scheme.context.L << scheme.context.logN];
        scheme.context.encode(poly[i], pvals[i], nslots, scheme.context.L);
        
        res[i] = Data;
        scheme.multByPolyAndEqual(res[i], poly[i]); // (1,0,0,0,|1,0,0,0|....)
        
        for(long j = 0; j < 2; ++j){
            Ciphertext ctemp = scheme.leftRotateFast(res[i], (1 << j)); // by 2^i * 16
            scheme.addAndEqual(res[i], ctemp);
        }
        scheme.reScaleByAndEqual(res[i], 1);
    }
    NTL_EXEC_RANGE_END;
}

//! Output:
//! 1) encAdj[10], lvl = L - 2 / comp: 10 * (2raw mult + 1 KS +  2Rot )
//! 2) encDet = sum Data[i] * encAdj[i] , 0 <= i < 3, lvl = L - 3/ comp: 4 raw mul + 1 KS

void CipherPvals::encSIMDAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* encData){
    NTL_EXEC_RANGE(10, first, last);
    for (long l = first; l < last; ++l){
        ExtCiphertext temp = extscheme.rawmult(encData[3 * l], encData[3 * l + 1]);
        temp = extscheme.rawmult(temp, encData[3 * l + 2]);
        encAdj[l] = extscheme.ModRaiseKeySwitch(temp);
        //encAdj[l] = scheme.mult(encData[3 * l], encData[3 * l + 1]);
        //scheme.multAndEqual(encAdj[l], encData[3 * l + 2]);
        
        for(long j = 0; j < 3; ++j){
            Ciphertext ctemp = scheme.leftRotateFast(encAdj[l], (1 << j));
            scheme.addAndEqual(encAdj[l], ctemp);
        }
        scheme.reScaleByAndEqual(encAdj[l], 2);
    }
    NTL_EXEC_RANGE_END;
    
    // "+------------------------------------+"
    //!  determinant
    // "+------------------------------------+"
  
    ExtCiphertext* temp = new ExtCiphertext[4];
    
    NTL_EXEC_RANGE(4, first, last);
    for (long l = first; l < last; ++l){
        Ciphertext ctemp = scheme.modDownTo(encData[30 + l], encAdj[l].l);
        temp[l] = extscheme.rawmult(ctemp, encAdj[l]);
    }
    NTL_EXEC_RANGE_END;
    
    for (long l = 1; l < 4; ++l){
        extscheme.addAndEqual(temp[0], temp[l]);
    }
    encDet = extscheme.ModRaiseKeySwitch(temp[0]);
    scheme.reScaleByAndEqual(encDet, 1);
    
    cout << "encDet.l : " << encDet.l << endl;
    
    delete[] temp;
}


//! Output: encAdj[10] = E(adj(X^T * X)[0,0]), lvl = L - 2
//!         encDet = sum Data[i] * encAdj[i] , 0 <= i < 3, L - 2
// 28 raw mult + 22 KS + 28 Rot
void CipherPvals::encSIMDAdjoint_smalldep(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* encData){
    Ciphertext* encDet1 = new Ciphertext[4];
    
    //!  4* (4 HM + 4 Rot) = 4 * (4 raw mult + 4 KS + 4 Rot)
    NTL_EXEC_RANGE(4, first, last);
    for (long l = first; l < last; ++l){
        encAdj[l] = scheme.mult(encData[3 * l], encData[3 * l + 1]);
        scheme.reScaleByAndEqual(encAdj[l], 1);
        
        encDet1[l] = scheme.mult(encData[3 * l + 2], encData[30 + l]);  // required for det
        scheme.reScaleByAndEqual(encDet1[l], 1);
        
        scheme.multAndEqual(encDet1[l], encAdj[l]);
        
        Ciphertext tmp = scheme.modDownTo(encData[3 * l + 2], encAdj[l].l);
        scheme.multAndEqual(encAdj[l], tmp);
        
        for(long j = 0; j < 3; ++j){
            Ciphertext ctemp = scheme.leftRotateFast(encAdj[l], (1 << j));
            scheme.addAndEqual(encAdj[l], ctemp);
            
            Ciphertext ctemp1 = scheme.leftRotateFast(encDet1[l], (1 << j));
            scheme.addAndEqual(encDet1[l], ctemp1);
        }
        scheme.reScaleByAndEqual(encAdj[l], 1);
    }
    NTL_EXEC_RANGE_END;
    
    //! encAdj[4], ..., encAdj[9]
    //! 6 * (2raw mult + 1 KS +  2Rot )
    NTL_EXEC_RANGE(6, first, last);
    for (long l = first; l < last; ++l){
        ExtCiphertext temp = extscheme.rawmult(encData[3 * (l + 4)], encData[3 * (l + 4) + 1]);
        temp = extscheme.rawmult(temp, encData[3 * (l + 4) + 2]);
        encAdj[l + 4] = extscheme.ModRaiseKeySwitch(temp);
        
        for(long j = 0; j < 3; ++j){
            Ciphertext ctemp = scheme.leftRotateFast(encAdj[l + 4], (1 << j));
            scheme.addAndEqual(encAdj[l + 4], ctemp);
        }
        scheme.reScaleByAndEqual(encAdj[l + 4], 2);
    }
    NTL_EXEC_RANGE_END;
    
    
    encDet = encDet1[0];
    for (long l = 1; l < 4; ++l){
        scheme.addAndEqual(encDet, encDet1[l]);
    }
    scheme.reScaleByAndEqual(encDet, 1);
    
    //cout << "encAdj.l: " << encAdj[0].l << "," << encAdj[5].l << endl;
    //cout << "encDet.l : " << encDet.l << endl;
    
    delete[] encDet1;
}

//! encYX[k] = E((y^T * X)[k])
void CipherPvals::aggYXData_DecompKS(Ciphertext*& encYX, Ciphertext encYXData, long sdimBits, long nbatching, long factorDim, long nslots){
    long nslots1 = nbatching * factorDim;  // nbatching = replicated number of a user's data in a single ciphertext
    Ciphertext res = encYXData;
    
    //! Allsum over "sampleDim"
    for(long i = 0; i < sdimBits; ++i){
        Ciphertext ctemp = extscheme.leftRotateFast(res, (1 << i) * nslots1);
        scheme.addAndEqual(res, ctemp);
    }
    fullReplicate4_DecompKS(encYX, res, nslots);
}

//! (34 * 8) Rot
void CipherPvals::aggCovData_DecompKS(Ciphertext*& encCov, Ciphertext* enccovData,  long sdimBits, long nbatching){
    long nslots1 = nbatching * 8;
    NTL_EXEC_RANGE(34, first, last);
    for (long l = first; l < last; ++l){
        encCov[l] = enccovData[l];
        //! Allsum
        for(long i = 0; i < sdimBits; ++i){
            Ciphertext ctemp = extscheme.leftRotateFast(encCov[l], (1 << i) * nslots1); // by 2^i * 16
            scheme.addAndEqual(encCov[l], ctemp);
        }
    }
    NTL_EXEC_RANGE_END;
}

//!@ Full replication of vector of size 4  (v0,v1,v2,v3) -> (v0), (v1), (v2), (v3)
//!@ Return ciphertexts of L-1
//!@ (2*8) Rot

void CipherPvals::fullReplicate4_DecompKS(Ciphertext*& res, Ciphertext Data, long nslots){
    //! full replication of size "4"
    complex<double>** pvals = new complex<double>*[4];
    uint64_t** poly = new uint64_t*[4];
    
    NTL_EXEC_RANGE(4, first, last);
    for (long i = first; i < last; ++i) {
        pvals[i] = new complex<double>[nslots];
        for(long j = 0; j < nslots ; j += 4){
            pvals[i][j + i].real(1.0);
        }
        poly[i] = new uint64_t[scheme.context.L << scheme.context.logN];
        scheme.context.encode(poly[i], pvals[i], nslots, scheme.context.L);
        
        res[i] = Data;
        scheme.multByPolyAndEqual(res[i], poly[i]); // (1,0,0,0,|1,0,0,0|....)
        
        for(long j = 0; j < 2; ++j){
            Ciphertext ctemp = extscheme.leftRotateFast(res[i], (1 << j));
            scheme.addAndEqual(res[i], ctemp);
        }
        
        scheme.reScaleByAndEqual(res[i], 1);
    }
    NTL_EXEC_RANGE_END;
}


//! Output: encAdj[10], lvl = L - 2
//!         encDet = sum Data[i] * encAdj[i] , 0 <= i < 3, L - 2
// 28 raw mult + 22 KS + 28 Rot
void CipherPvals::encSIMDAdjoint_DecompKS(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* encData){
    Ciphertext* encDet1 = new Ciphertext[4];

    //!  4* (4 HM + 4 Rot) = 4 * (4 raw mult + 4 KS + 4 Rot)
    //!  encAdj[l] = \sum (encData[3 * l] * encData[3 * l + 1]) * encData[3 * l + 2]
    //!  encDet[l] = \sum (encData[3 * l] * encData[3 * l + 1]) * (encData[3 * l + 2] * encData[30 + l])
    
    NTL_EXEC_RANGE(4, first, last);
    for (long l = first; l < last; ++l){
        encAdj[l] = extscheme.mult(encData[3 * l], encData[3 * l + 1]);
        scheme.reScaleByAndEqual(encAdj[l], 1);
       
        encDet1[l] = extscheme.mult(encData[3 * l + 2], encData[30 + l]);  // required for det
        scheme.reScaleByAndEqual(encDet1[l], 1);
        
        encDet1[l] = extscheme.mult(encDet1[l], encAdj[l]);
        
        Ciphertext tmp = scheme.modDownTo(encData[3 * l + 2], encAdj[l].l);
        encAdj[l] = extscheme.mult(encAdj[l], tmp);
      
        for(long j = 0; j < 3; ++j){
            Ciphertext ctemp = extscheme.leftRotateFast(encAdj[l], (1 << j));
            scheme.addAndEqual(encAdj[l], ctemp);
            
            Ciphertext ctemp1 = extscheme.leftRotateFast(encDet1[l], (1 << j));
            scheme.addAndEqual(encDet1[l], ctemp1);
        }
        scheme.reScaleByAndEqual(encAdj[l], 1);
    }
    NTL_EXEC_RANGE_END;
    
    
    //! encAdj[4], ..., encAdj[9]
    //! 6 * (2raw mult + 1 KS +  2Rot )
    NTL_EXEC_RANGE(6, first, last);
    for (long l = first; l < last; ++l){
        ExtCiphertext tmp = extscheme.rawmult3(encData[3 * (l + 4)], encData[3 * (l + 4) + 1], encData[3 * (l + 4) + 2]);
        encAdj[l + 4] = extscheme.DecompKeySwitch(tmp);
        
        for(long j = 0; j < 3; ++j){
            Ciphertext tmp1 = extscheme.leftRotateFast(encAdj[l + 4], (1 << j));
            scheme.addAndEqual(encAdj[l + 4], tmp1);
        }
        scheme.reScaleByAndEqual(encAdj[l + 4], 2);
    }
    NTL_EXEC_RANGE_END;
    

    encDet = encDet1[0];
    for (long l = 1; l < 4; ++l){
        scheme.addAndEqual(encDet, encDet1[l]);
    }
    scheme.reScaleByAndEqual(encDet, 1);

    //cout << "encAdj.l: " << encAdj[0].l << "," << encAdj[5].l << endl;
    //cout << "encDet.l : " << encDet.l << endl;
    
    delete[] encDet1;
}

/********************************************************************/


//! sum_{i,j} encMatrix[i][j] * (encData1[i] * encData2[j]) :
//! 16 * (2 raw mult) + 1 KS
void CipherPvals::extQuadForm(Ciphertext& res, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim){
    ExtCiphertext** tensoring = new ExtCiphertext*[factorDim];

    //! (encData1[i] * encData2[j]) * encMatrix[i][j]: 1 lvl
    //! two raw mult (lower levels) + one mult
    NTL_EXEC_RANGE(factorDim, first, last);
    for(long i = first; i< last; ++i){
        tensoring[i] = new ExtCiphertext[factorDim];

        //! lower
        for(long j = 0; j < i; ++j){
            long l1 = j* factorDim - (j * (j-1))/2 + (i - j);   // l = j* factorDim - j*(j-1)/2 + i - j;
            tensoring[i][j] = extscheme.rawmult(encData1[i], encData2[j]);
            extscheme.reScaleByAndEqual(tensoring[i][j], 1);
            extscheme.modDownToAndEqual(tensoring[i][j], encMatrix[l1].l);
            tensoring[i][j] = extscheme.rawmult(tensoring[i][j], encMatrix[l1]);
        }
        
        //! upper
        long diag_index = i* factorDim - (i * (i-1))/2;
        for(long j = i; j < factorDim; ++j){
            long l1 = diag_index + (j - i);   // l = i* factorDim - i*(i-1)/2 + j - i ;
            tensoring[i][j] = extscheme.rawmult(encData1[i], encData2[j]);
            extscheme.reScaleByAndEqual(tensoring[i][j], 1);
            extscheme.modDownToAndEqual(tensoring[i][j], encMatrix[l1].l);
            tensoring[i][j] = extscheme.rawmult(tensoring[i][j], encMatrix[l1]);
        }
        
        for(long j = 1; j <factorDim; ++j){
            extscheme.addAndEqual(tensoring[i][0], tensoring[i][j]);
        }
        
    }
    NTL_EXEC_RANGE_END;
    
    ExtCiphertext extres = tensoring[0][0];
    for(long i = 1; i < factorDim; ++i){
        extscheme.addAndEqual(extres, tensoring[i][0]);
    }
    
    //! Final KS rescaling
    res = extscheme.ModRaiseKeySwitch(extres);
    scheme.reScaleByAndEqual(res, 1);   // rawmult -> rescaleBy "2"
    
     delete[] tensoring;
}

//! res = \sum Matrix[i][j] * Data[i] * Data[j], 0 <= i,j < 4
//! 10 * (2 raw mult) + 1 KS
void CipherPvals::extSqrQuadForm(Ciphertext& res, Ciphertext* encData, Ciphertext* encMatrix, long factorDim){
    ExtCiphertext** tensoring = new ExtCiphertext*[factorDim];
    
    NTL_EXEC_RANGE(factorDim, first, last);
    for(int i = first; i< last; ++i){
        tensoring[i] = new ExtCiphertext[factorDim];
        
        //! Firstly, compute the diagonal (using squaring)
        tensoring[i][i] = extscheme.rawsquare(encData[i]);
        extscheme.reScaleByAndEqual(tensoring[i][i], 1);
        
        long diag_index = i* factorDim - (i * (i-1))/2;
        extscheme.modDownToAndEqual(tensoring[i][i], encMatrix[diag_index].l);
        tensoring[i][i] = extscheme.rawmult(tensoring[i][i], encMatrix[diag_index]);  //!if we multiply together, then it requires 2 lvls
        
        for(long j = i + 1; j < factorDim; ++j){
            long l1 = diag_index + (j - i);
            tensoring[i][j] = extscheme.rawmult(encData[i], encData[j]);
            extscheme.reScaleByAndEqual(tensoring[i][j], 1);
            extscheme.modDownToAndEqual(tensoring[i][j], encMatrix[l1].l);
            tensoring[i][j] = extscheme.rawmult(tensoring[i][j], encMatrix[l1]);
        }
    
        for(long j = i + 2; j <factorDim; ++j){
            extscheme.addAndEqual(tensoring[i][i+1], tensoring[i][j]);
        }
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(factorDim - 1, first, last);
    for(int i = first; i< last; ++i){
        extscheme.addAndEqual(tensoring[i][i], tensoring[i][i+1]);
        extscheme.addAndEqual(tensoring[i][i], tensoring[i][i+1]);
    }
    NTL_EXEC_RANGE_END;
    
    ExtCiphertext extres = tensoring[0][0];
    for(long i = 1; i < factorDim; ++i){
        extscheme.addAndEqual(extres, tensoring[i][i]);
    }

    res = extscheme.ModRaiseKeySwitch(extres);
    scheme.reScaleByAndEqual(res, 1);
    
    delete[] tensoring;
}


//! sum_{i,j} encMatrix[i][j] * (encData1[i] * encData2[j]) :
//! 16 * (2 raw mult) + 1 KS
void CipherPvals::extQuadForm_DecompKS(Ciphertext& res, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim){
    ExtCiphertext** tensoring = new ExtCiphertext*[factorDim];
    
    //! (encData1[i] * encData2[j]) * encMatrix[i][j]: 1 lvl
    //! two raw mult (lower levels) + one mult
    NTL_EXEC_RANGE(factorDim, first, last);
    for(long i = first; i< last; ++i){
        tensoring[i] = new ExtCiphertext[factorDim];
        
        //! lower
        for(long j = 0; j < i; ++j){
            long l1 = j* factorDim - (j * (j-1))/2 + (i - j);   // l = j* factorDim - j*(j-1)/2 + i - j;
            tensoring[i][j] = extscheme.rawmult(encData1[i], encData2[j]);
            extscheme.reScaleByAndEqual(tensoring[i][j], 1);
            extscheme.modDownToAndEqual(tensoring[i][j], encMatrix[l1].l);
            tensoring[i][j] = extscheme.rawmult(tensoring[i][j], encMatrix[l1]);
        }
        
        //! upper
        long diag_index = i* factorDim - (i * (i-1))/2;
        for(long j = i; j < factorDim; ++j){
            long l1 = diag_index + (j - i);   // l = i* factorDim - i*(i-1)/2 + j - i ;
            tensoring[i][j] = extscheme.rawmult(encData1[i], encData2[j]);
            extscheme.reScaleByAndEqual(tensoring[i][j], 1);
            extscheme.modDownToAndEqual(tensoring[i][j], encMatrix[l1].l);
            tensoring[i][j] = extscheme.rawmult(tensoring[i][j], encMatrix[l1]);
        }
        
        for(long j = 1; j <factorDim; ++j){
            extscheme.addAndEqual(tensoring[i][0], tensoring[i][j]);
        }
        
    }
    NTL_EXEC_RANGE_END;
    
    ExtCiphertext extres = tensoring[0][0];
    for(long i = 1; i < factorDim; ++i){
        extscheme.addAndEqual(extres, tensoring[i][0]);
    }
    
    //! Final KS rescaling
    res = extscheme.DecompKeySwitch(extres);
    scheme.reScaleByAndEqual(res, 1);   // rawmult -> rescaleBy "2"
    
    delete[] tensoring;
}

//! res = \sum Matrix[i][j] * Data[i] * Data[j], 0 <= i,j < 4
//! 10 * (2 raw mult) + 1 KS
void CipherPvals::extSqrQuadForm_DecompKS(Ciphertext& res, Ciphertext* encData, Ciphertext* encMatrix, long factorDim){
    ExtCiphertext** tensoring = new ExtCiphertext*[factorDim];
    
    NTL_EXEC_RANGE(factorDim, first, last);
    for(int i = first; i< last; ++i){
        tensoring[i] = new ExtCiphertext[factorDim];
        
        //! Firstly, compute the diagonal (using squaring)
        tensoring[i][i] = extscheme.rawsquare(encData[i]);
        extscheme.reScaleByAndEqual(tensoring[i][i], 1);
        
        long diag_index = i* factorDim - (i * (i-1))/2;
        extscheme.modDownToAndEqual(tensoring[i][i], encMatrix[diag_index].l);
        tensoring[i][i] = extscheme.rawmult(tensoring[i][i], encMatrix[diag_index]);  //!if we multiply together, then it requires 2 lvls
        
        for(long j = i + 1; j < factorDim; ++j){
            long l1 = diag_index + (j - i);
            tensoring[i][j] = extscheme.rawmult(encData[i], encData[j]);
            extscheme.reScaleByAndEqual(tensoring[i][j], 1);
            extscheme.modDownToAndEqual(tensoring[i][j], encMatrix[l1].l);
            tensoring[i][j] = extscheme.rawmult(tensoring[i][j], encMatrix[l1]);
        }
        
        for(long j = i + 2; j <factorDim; ++j){
            extscheme.addAndEqual(tensoring[i][i+1], tensoring[i][j]);
        }
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(factorDim - 1, first, last);
    for(int i = first; i< last; ++i){
        extscheme.addAndEqual(tensoring[i][i], tensoring[i][i+1]);
        extscheme.addAndEqual(tensoring[i][i], tensoring[i][i+1]);
    }
    NTL_EXEC_RANGE_END;
    
    ExtCiphertext extres = tensoring[0][0];
    for(long i = 1; i < factorDim; ++i){
        extscheme.addAndEqual(extres, tensoring[i][i]);
    }
    
    res = extscheme.DecompKeySwitch(extres);
    scheme.reScaleByAndEqual(res, 1);
    
    delete[] tensoring;
}


/********************************************************************/

//! sum_{i,j} encMatrix[i][j] * (encData1[i] * encData2[j]) : 2*k^2 HM (2 levels)
void CipherPvals::QuadForm(Ciphertext& res, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim){
    Ciphertext** tensoring = new Ciphertext*[factorDim];
    
    // 0: 0 1 2 3
    // 1: 1 4 5 6
    // 2: 2 5 7 8
    // 3: 3 6 8 9
    
    NTL_EXEC_RANGE(factorDim, first, last);
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
    NTL_EXEC_RANGE_END;
    
    res = tensoring[0][0];
    for(long i = 1; i < factorDim; ++i){
        scheme.addAndEqual(res, tensoring[i][0]);
    }
    scheme.reScaleByAndEqual(res, 1);
    
    delete[] tensoring;
}


//! sum_{i,j} encMatrix[i][j] * (encData1[i] * encData2[j]) : 2* (k + (k2 - k)/2) HM
//! lvl: pBits (1 level)
void CipherPvals::SqrQuadForm(Ciphertext& res, Ciphertext* encData, Ciphertext* encMatrix, long factorDim){
    Ciphertext** tensoring = new Ciphertext*[factorDim];
    
    NTL_EXEC_RANGE(factorDim, first, last);
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
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(factorDim - 1, first, last);
    for(int i = first; i< last; ++i){
        scheme.addAndEqual(tensoring[i][i], tensoring[i][i+1]);
        scheme.addAndEqual(tensoring[i][i], tensoring[i][i+1]);
    }
    NTL_EXEC_RANGE_END;
    
    res = tensoring[0][0];
    for(long i = 1; i < factorDim; ++i){
        scheme.addAndEqual(res, tensoring[i][i]);
    }
    
    scheme.reScaleByAndEqual(res, 1);
    
     delete[] tensoring;
}



//! encYXData[i][k]: = Enc(y[i] * xDta[i][k]), 0 <= i <n
//! encCovData[i][k] = Enc(cov_i[k]) (10s)

void CipherPvals::encryptXData(Ciphertext**& encYXData, Ciphertext**& enccovData, double* yData, double** xData, long factorDim, long sampleDim,  long nterms, long scalefactor, long nslots, long L) {
    
    double** cov = new double*[sampleDim];
    double** xiData = new double*[sampleDim];
    
    NTL_EXEC_RANGE(sampleDim, first, last);
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
        
        xiData[i] = new double[factorDim];
        for(long k = 0; k < factorDim; ++k){
            xiData[i][k] = xData[i][k];
        }
    
        cov[i] = new double[nterms];
        computeCov(cov[i], xiData[i], factorDim, scalefactor);  //! dim -> dim * (dim-1)/2
        
        for(long k = 0; k < nterms; ++k){
            encValue(enccovData[i][k], cov[i][k], nslots, L);
        }
    }
    NTL_EXEC_RANGE_END;
    
    delete[] cov;
    delete[] xiData;
    
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
    
    NTL_EXEC_RANGE(20, first, last);
    for(long i = first; i < last; ++i){
        long j0= table[i][0];
        long j1= table[i][1];    
        temp[i] = scheme.mult(encData[j0], encData[j1]);
    }
    NTL_EXEC_RANGE_END;
    
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
        {1, 1, 2},{2, 1, 2}, {3, 1, 2},{5, 0, 2},{6, 0, 2},{8, 0, 1},
    };
    
    
    //! diagonal computation
    NTL_EXEC_RANGE(4, first, last);
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
    NTL_EXEC_RANGE_END;
    
    //! non-diagonal computation: x[k1] * adj[k] + x[k2] * adj[k + 1] + x[3] * adj[k + 2]
    NTL_EXEC_RANGE(6, first, last);
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
    NTL_EXEC_RANGE_END;
    
    // "+------------------------------------+"
    //               Deternimant
    // "+------------------------------------+"
    
    
    Ciphertext* encDetemp = new Ciphertext[4];
    NTL_EXEC_RANGE(4, first, last);
    for(long i = first; i < last; ++i){
        encDetemp[i] = scheme.mult(xtemp[i], encAdj[i]);
    }
    NTL_EXEC_RANGE_END;
    encDet = encDetemp[0];
    for(long i = 1; i < 4; ++i){
        scheme.addAndEqual(encDet, encDetemp[i]);
    }
    
    
    //! scale by 2 * logp
    NTL_EXEC_RANGE(10, first, last);
    for(long i = first; i < last; ++i){
        scheme.reScaleByAndEqual(encAdj[i], 2);
    }
    NTL_EXEC_RANGE_END;
        
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
    
    delete[] sqrtemp;
    delete[] temp;
    delete[] adj;
    delete[] xtemp;
    delete[] encDetemp;
}

//! Output: encAdj[10], 54 raw mul + 10 KS
//!         encDet, 4 raw mul + 1 KS

void CipherPvals::extencAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* encData){
    double table[20][2] = {
        {1, 5},{1, 7},{1, 8},{1, 9},
        {2, 4},{2, 5},{2, 6},{2, 9},
        {3, 5},{3, 6},{3, 8},
        {4, 7},{4, 8},{4, 9},
        {5, 6},{5, 8},{5, 9},
        {6, 7},{6, 8},{7, 9},
    };
    
    //! 1. squaring and multiplication
    ExtCiphertext* sqrtemp = new ExtCiphertext[4];
    ExtCiphertext* temp = new ExtCiphertext[20];
    
    NTL_EXEC_RANGE(4, first, last);
    for(long i = first; i < last; ++i){
        if(i == 0){
            sqrtemp[i] = extscheme.rawsquare(encData[3]); //! {{0,3} -> (x[0][3])^2
        }
        else if(i == 1){
             sqrtemp[i] = extscheme.rawsquare(encData[5]); //! {1,2},
        }
        else if(i == 2){
            sqrtemp[i] = extscheme.rawsquare(encData[6]);  //! {1,3}
        }
        else{
            sqrtemp[i] = extscheme.rawsquare(encData[8]); //! {2,3}
        }
    }
    NTL_EXEC_RANGE_END;
    
    
    NTL_EXEC_RANGE(20, first, last);
    for(long i = first; i < last; ++i){
        long j0= table[i][0];
        long j1= table[i][1];
        temp[i] = extscheme.rawmult(encData[j0], encData[j1]);
    }
    NTL_EXEC_RANGE_END;
   
    //! 2. two polynomials mult
    ExtCiphertext* adj = new ExtCiphertext[30];
    
    adj[0] = extscheme.sub(temp[19], sqrtemp[3]);
    adj[1] = extscheme.sub(temp[18], temp[16]);
    adj[2] = sqrtemp[2];
    
    
    adj[3] = extscheme.sub(sqrtemp[3], temp[19]);
    adj[4] = extscheme.negate(adj[1]);
    adj[5] = extscheme.sub(temp[17], temp[15]);
    
    adj[6] = adj[4];
    adj[7] = extscheme.sub(sqrtemp[2], temp[13]);
    adj[8] = extscheme.sub(temp[12], temp[14]);
    
    adj[9]  = adj[5];
    adj[10] = adj[8];
    adj[11] = extscheme.sub(sqrtemp[1], temp[11]);
    
    adj[12] = extscheme.negate(adj[3]);
    adj[13] = extscheme.sub(temp[10], temp[7]);
    adj[14] = sqrtemp[0];
    
    adj[15] = adj[1];
    adj[16] = extscheme.sub(temp[3], temp[9]);
    adj[17] = extscheme.sub(temp[8], temp[2]);
    
    adj[18] = extscheme.negate(adj[5]);
    adj[19] = extscheme.sub(temp[6], temp[2]);
    adj[20] = extscheme.sub(temp[1], temp[5]);
    
    adj[21] = extscheme.negate(adj[7]);
    adj[22] = extscheme.negate(adj[16]);
    adj[23] = sqrtemp[0];
    
    adj[24] = extscheme.negate(adj[8]);
    adj[25] = extscheme.negate(adj[19]);
    adj[26] = extscheme.sub(temp[4], temp[0]);
    
    adj[27] = extscheme.sub(temp[11], sqrtemp[1]);
    adj[28] = extscheme.negate(adj[20]);
    adj[29] = temp[4];
    

    Ciphertext* xtemp = new Ciphertext[8];
    for(long i = 0; i < 8; ++i){
        xtemp[i] = encData[i];
    }
    
    /********************************************************/
    double diag[4][5] = {{0, 4, 5, 18, 7},{4, 0, 2, 10, 7},{7, 0, 1, 9, 4},{9, 0, 1, 5, 2},};
    double nondiag[6][3] = {{1, 1, 2},{2, 1, 2}, {3, 1, 2},{5, 0, 2},{6, 0, 2},{8, 0, 1},};
    
    //! diagonal computation
    ExtCiphertext* extAdj = new ExtCiphertext[10];
    
    NTL_EXEC_RANGE(4, first, last);
    for(long i = first; i < last; ++i){
        long k = 3 * diag[i][0];
        long k0 = diag[i][0];
        long k1 = diag[i][1];
        long k2 = diag[i][2];
        long k3 = diag[i][3];
        long k4 = diag[i][4];
        
        extAdj[k0] = extscheme.rawmult(adj[k], xtemp[k1]);  // deg = 3
        ExtCiphertext ctemp = extscheme.add(adj[k + 1], temp[k3]);
        ctemp = extscheme.rawmult(ctemp, xtemp[k2]);
        extscheme.addAndEqual(extAdj[k0], ctemp);
        ctemp = extscheme.rawmult(adj[k + 2], xtemp[k4]);
        extscheme.subAndEqual(extAdj[k0], ctemp);
    }
    NTL_EXEC_RANGE_END;

    //! non-diagonal computation: x[k1] * adj[k] + x[k2] * adj[k + 1] + x[3] * adj[k + 2]
    NTL_EXEC_RANGE(6, first, last);
    for(long i = first; i < last; ++i){
        long k = 3 * nondiag[i][0];
        long k0 = nondiag[i][0];
        long k1 = nondiag[i][1];
        long k2 = nondiag[i][2];
        
        extAdj[k0] = extscheme.rawmult(adj[k], xtemp[k1]);
        ExtCiphertext ctemp = extscheme.rawmult(adj[k + 1], xtemp[k2]);
        extscheme.addAndEqual(extAdj[k0], ctemp);
        ctemp = extscheme.rawmult(adj[k + 2], xtemp[3]);
        extscheme.addAndEqual(extAdj[k0], ctemp);
    }
    NTL_EXEC_RANGE_END;
    
    //! KS + scale by 2 * logp
    NTL_EXEC_RANGE(10, first, last);
    for(long i = first; i < last; ++i){
        encAdj[i] = extscheme.ModRaiseKeySwitch(extAdj[i]);
        scheme.reScaleByAndEqual(encAdj[i], 2);
    }
    NTL_EXEC_RANGE_END;
   
    // "+------------------------------------+"
    //               Deternimant
    // "+------------------------------------+"
    
    ExtCiphertext* encDetemp = new ExtCiphertext[4];
    NTL_EXEC_RANGE(4, first, last);
    for(long i = first; i < last; ++i){
        scheme.modDownToAndEqual(xtemp[i], encAdj[i].l);
        encDetemp[i] = extscheme.rawmult(xtemp[i], encAdj[i]);
    }
    NTL_EXEC_RANGE_END;
   
    for(long i = 1; i < 4; ++i){
        extscheme.addAndEqual(encDetemp[0], encDetemp[i]);
    }
    
    encDet = extscheme.ModRaiseKeySwitch(encDetemp[0]);
    scheme.reScaleByAndEqual(encDet, 1);

    
    delete[] sqrtemp;
    delete[] temp;
    delete[] adj;
    delete[] xtemp;
    delete[] encDetemp;
}

