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
#include "Matrix.h"

//!@ encrypt an input value and generate a fully-packed ciphertext at level L
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

void CipherPvals::encryptSData(Ciphertext**& encSData, Ciphertext**& encYSData, Ciphertext***& encSXData,  double* yData, double** xData, double** sData,  long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots, long Slvl, long YSlvl, long SXlvl) {
    long nslots1 = nsnp - (nencsnp-1) * nslots;   // number of slots in the final ctxts
    
    double** scaled_sData = new double*[sampleDim];
    double*** sxData = new double**[sampleDim];
    double** fullvec = new double*[sampleDim];
    double** sparsevec = new double*[sampleDim];
    
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        fullvec[i] = new double[nslots];
        sparsevec[i] = new double[nslots1];
        /***************************************/
        //! 1. encryption of sData (because s = 0/1)
        long j1 = 0;
        //! full slot
        for(long j = 0; j < nencsnp - 1; ++j){
            for(long l = 0; l < nslots; ++l){
                fullvec[i][l] = sData[i][j1] ;
                j1++;
            }
            encFullyPackedVec(encSData[i][j], fullvec[i], nslots, Slvl);
        }
        //! not full slot
        for(long l = 0; l < nslots1; ++l){
            sparsevec[i][l] = sData[i][j1] ;
            j1++;
        }
        encSparselyPackedVec(encSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, Slvl);
        /***************************************/
        //! 1. encryption of YSData = Y * S
        j1 = 0;
        if(yData[i] == 1){
            //! full slot
            for(long j = 0; j < nencsnp - 1; ++j){
                for(long l = 0; l < nslots; ++l){
                    fullvec[i][l] = sData[i][j1] ;
                    j1++;
                }
                encFullyPackedVec(encYSData[i][j], fullvec[i], nslots, YSlvl);
            }
            //! not full slot
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = sData[i][j1] ;
                j1++;
            }
            encSparselyPackedVec(encYSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, YSlvl);
        }
        else{
            for(long j = 0; j < nencsnp - 1; ++j){
                for(long l = 0; l < nslots; ++l){
                    fullvec[i][l] = - sData[i][j1];
                    j1++;
                }
                encFullyPackedVec(encYSData[i][j], fullvec[i], nslots, YSlvl);
            }
            //! not full slot
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = - sData[i][j1];
                j1++;
            }
            encSparselyPackedVec(encYSData[i][nencsnp-1], sparsevec[i], nslots1, nslots, YSlvl);
        }
        /***************************************/
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
                encFullyPackedVec(encSXData[i][k][j], fullvec[i], nslots, SXlvl);
            }
            //! not full slot
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = sxData[i][k][j1];
                j1++;
            }
            encSparselyPackedVec(encSXData[i][k][nencsnp-1], sparsevec[i], nslots1, nslots, SXlvl);
            delete [] sxData[i][k];
        }
        delete[] sxData[i];
    }
    NTL_EXEC_RANGE_END;
    
    delete[] scaled_sData;
    delete[] sxData;
    delete[] fullvec;
    delete[] sparsevec;
}

/********************************************************************/

void CipherPvals::encryptSIMDXData(Ciphertext& encYXData, Ciphertext*& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long sampleDim2, long nXbatching, long nCovbatching, long nterms, long scaleBits, long nslots, long YXlvl, long Covlvl) {
    
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
    encFullyPackedVec(encYXData, temp, nslots, YXlvl);
    
    // "+------------------------------------+"
    //!  encryption of covariance
    // "+------------------------------------+"
   
    double** cov = new double*[sampleDim];
    
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
        encFullyPackedVec(enccovData[l], xData2[l], nslots, Covlvl);
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
        encFullyPackedVec(enccovData[30 + l], xData2[l], nslots, Covlvl);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] temp;
    delete[] cov;
    delete[] xData2;
}

// ! encSData: ciphertext for individual entries in S, packed
//   encSData[i][j] = E{S[i][j']} ... E{S[i][j'+nencsnp-1]} packed in a ciphertext
// ! encXSData: ciphertext for decomposed entries in X^T * S
//   Notice that (X^T * S)_i,j = \sum_{l=1}^n X_l,i S_l,i
//   So X^T * S = \sum_{l=1}^n XS_l, where (XS_l)_i,j = X_l,i S_l,j
//   encXSData[l][i][j] = E{ X_l,i S_l,j }, l in [n], i in [k], j in [nencsnp] = [m/slots]
//   where n = sampleDim
//         k = # feature = factorDim
//         m = # SNP
void CipherPvals::encryptTrivialSData(Ciphertext**& encSData, Ciphertext***& encXSData, Ciphertext**& encYSData, Matrix const& matX, Matrix const& matS, Matrix const& matY, long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots) {
    long nslots1 = nsnp - (nencsnp-1) * nslots;   // number of slots in the final ctxts
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long l = first; l < last; ++l) {
       /***************************************/
       //! 1. encryption of sData (because s = 0/1)
       long j = 0;
       //! full slot
       double * svec = matS.data[l];
       for(long j2 = 0; j2 < nencsnp - 1; ++j2){ // index over packed j's
          encFullyPackedVec(encSData[l][j2], svec, nslots, 1);
          svec += nslots;
       }
       //! not full slot
       encSparselyPackedVec(encSData[l][nencsnp-1], svec, nslots1, nslots, 1);
       /***************************************/
       //! 2. encryption of X^T * S
       double * xsData = new double[nsnp]; // i-th row of X_l^T \otimes S_l
       double * xsvec = 0;
       for(long i = 0; i < factorDim; ++i){
          for(long j0 = 0; j0 < nsnp; ++j0){
             xsData[j0] = matX.at(l,i) * matS.at(l,j0);
          }
          j = 0;             // index over xsData
          //! full slot
          xsvec = xsData;
          for(long j2 = 0; j2 < nencsnp - 1; ++j2){ // index over packed j's
             encFullyPackedVec(encXSData[l][i][j2], xsvec, nslots, 1);
             xsvec += nslots;
          }
          //! not full slot
          encSparselyPackedVec(encXSData[l][i][nencsnp-1], xsvec, nslots1, nslots, 1);
       }
       delete [] xsData;

       /***************************************/
       //! 3. encryption of Y[l] * S[l][j]
       double * fullvec = new double[nslots];
       double * sparsevec = new double[nslots1];
       j = 0;
       for(long j2 = 0; j2 < nencsnp - 1; ++j2){ // index over packed j's
          for(long j1 = 0; j1 < nslots; ++j1){  // index over slots
             if(matY.at(l,0) == 1) {
                fullvec[j1] = matS.at(l,j);      // y[l] * S[l][j]
             } else {
                fullvec[j1] = 0;                 // y[l] * S[l][j]
             }
             j++;
          }
          encFullyPackedVec(encYSData[l][j2], fullvec, nslots, 1);
       }
       //! not full slot
       for(long j1 = 0; j1 < nslots1; ++j1){
          if(matY.at(l,0) == 1) {
             sparsevec[j1] = matS.at(l,j);      // y[l] * S[l][j]
          } else {
             sparsevec[j1] = 0;                 // y[l] * S[l][j]
          }
          j++;
       }
       encSparselyPackedVec(encYSData[l][nencsnp-1], sparsevec, nslots1, nslots, 1);
       delete[] fullvec;
       delete[] sparsevec;
    }
    NTL_EXEC_RANGE_END;
    
}

void CipherPvals::encryptTrivialYData(Ciphertext& encYData, Ciphertext*& encXYData, Ciphertext*& enccovData, Matrix const& matX, Matrix const& matY,  long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots) {
   assert(matY.n <= nslots);

   double* tempY = new double[nslots];
   long const nslots1 = nsnp - (nencsnp-1) * nslots;   // number of slots in the final ctxts
   long const covSizeEachLevel = factorDim*factorDim;

   NTL_EXEC_RANGE(sampleDim, first, last);
   for(long l = first; l < last; ++l) {
      /**************************************************/
      //! encryption of Y and (X^T y)
      double * tempXY = new double[nslots];
      if(matY.at(l,0) == 1) {
         tempY[l] = 1;                   // column vector
         for(long i = 0; i < factorDim; i++) {
            tempXY[i] = matX.at(l,i);
         }
      } else {
         tempY[l] = 0;                   // column vector
         for(long i = 0; i < factorDim; i++) {
            tempXY[i] = 0;
         }
      }
      encSparselyPackedVec(encXYData[l], tempXY, factorDim, nslots, 1);
      delete [] tempXY;

      /**************************************************/
      //!  encryption of covariance
      double * tempCov = new double[covSizeEachLevel];
      for(size_t i = 0; i < factorDim; i++) {
         for(size_t h = 0; h < factorDim; h++) {
            tempCov[i*factorDim+h] = matX.at(l,i) * matX.at(l,h);
         }
      }
      encFullyPackedVec(enccovData[l], tempCov, covSizeEachLevel, 1);
      delete [] tempCov;
   }
   NTL_EXEC_RANGE_END;

   encSparselyPackedVec(encYData, tempY, sampleDim, nslots, 1);
   delete [] tempY;
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
//! encYX[k] = E((y^T * X)[k])
void CipherPvals::aggYXData(Ciphertext*& encYX, Ciphertext encYXData, long sdimBits, long nbatching, long factorDim, long nslots){
    long nslots1 = nbatching * factorDim;  // nbatching = replicated number of a user's data in a single ciphertext
    Ciphertext res = encYXData;
    
    //! Allsum over "sampleDim"
    for(long i = 0; i < sdimBits; ++i){
        Ciphertext ctemp = extscheme.leftRotateFastMT(res, (1 << i) * nslots1);
        scheme.addAndEqual(res, ctemp);
    }
    //fullReplicate4(encYX, res, nslots);
    
    //!@ Full replication of vector of size 4  (v0,v1,v2,v3) -> (v0), (v1), (v2), (v3), 1 level
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
        
        encYX[i] = res;
        scheme.multByPolyAndEqual(encYX[i], poly[i]); // (1,0,0,0,|1,0,0,0|....)
        scheme.reScaleByAndEqual(encYX[i], 1);
        
        for(long j = 0; j < 2; ++j){
            Ciphertext ctemp = extscheme.leftRotateFast(encYX[i], (1 << j));
            scheme.addAndEqual(encYX[i], ctemp);
        }
    }
    NTL_EXEC_RANGE_END;
}


//! Output: encAdj[10], lvl = L - 2
//!         encDet = sum Data[i] * encAdj[i] , 0 <= i < 3, L - 2

void CipherPvals::encAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* enccovData, long sdimBits, long nbatching){
    //! 1. aggregate
    long nslots1 = nbatching * 8;
    Ciphertext* encCov = new Ciphertext[34];
    NTL_EXEC_RANGE(34, first, last);
    for (long l = first; l < last; ++l){
        encCov[l] = enccovData[l];
        for(long i = 0; i < sdimBits; ++i){
            Ciphertext ctemp = extscheme.leftRotateFast(encCov[l], (1 << i) * nslots1); // by 2^i * 16
            scheme.addAndEqual(encCov[l], ctemp);
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! 2. Compute the adj and det
    //!  4* (4 HM + 4 Rot) = 4 * (4 raw mult + 4 KS + 4 Rot)
    //!  encAdj[l] = \sum (encData[3 * l] * encData[3 * l + 1]) * encData[3 * l + 2]
    //!  encDet[l] = \sum (encData[3 * l] * encData[3 * l + 1]) * (encData[3 * l + 2] * encData[30 + l])
    Ciphertext* encDet1 = new Ciphertext[4];
    
    NTL_EXEC_RANGE(4, first, last);
    for (long l = first; l < last; ++l){
        encAdj[l] = extscheme.mult(encCov[3 * l], encCov[3 * l + 1]);
        scheme.reScaleByAndEqual(encAdj[l], 1);
        
        encDet1[l] = extscheme.mult(encCov[3 * l + 2], encCov[30 + l]);  // required for det
        scheme.reScaleByAndEqual(encDet1[l], 1);
        
        encDet1[l] = extscheme.mult(encDet1[l], encAdj[l]);
        
        Ciphertext tmp = scheme.modDownTo(encCov[3 * l + 2], encAdj[l].l);
        encAdj[l] = extscheme.mult(encAdj[l], tmp);
        scheme.reScaleByAndEqual(encAdj[l], 1);
        
        for(long j = 0; j < 3; ++j){
            Ciphertext ctemp = extscheme.leftRotateFast(encAdj[l], (1 << j));
            scheme.addAndEqual(encAdj[l], ctemp);
            
            Ciphertext ctemp1 = extscheme.leftRotateFast(encDet1[l], (1 << j));
            scheme.addAndEqual(encDet1[l], ctemp1);
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! encAdj[4], ..., encAdj[9]
    //! 6 * (2raw mult + 1 KS +  2Rot )
    NTL_EXEC_RANGE(6, first, last);
    for (long l = first; l < last; ++l){
        ExtCiphertext tmp = extscheme.rawmult3(encCov[3 * (l + 4)], encCov[3 * (l + 4) + 1], encCov[3 * (l + 4) + 2]);
        extscheme.reScaleByAndEqual(tmp, 2);
        encAdj[l + 4] = extscheme.DecompKeySwitch(tmp);
        
        for(long j = 0; j < 3; ++j){
            Ciphertext tmp1 = extscheme.leftRotateFast(encAdj[l + 4], (1 << j));
            scheme.addAndEqual(encAdj[l + 4], tmp1);
        }
        //scheme.reScaleByAndEqual(encAdj[l + 4], 2);
    }
    NTL_EXEC_RANGE_END;
    
    encDet = encDet1[0];
    for (long l = 1; l < 4; ++l){
        scheme.addAndEqual(encDet, encDet1[l]);
    }
    scheme.reScaleByAndEqual(encDet, 1);
    
    delete[] encCov;
    delete[] encDet1;
}


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
    extscheme.reScaleByAndEqualMT(extres, 1);   // rawmult -> rescaleBy "2"
    res = extscheme.DecompKeySwitchMT(extres);
    //scheme.reScaleByAndEqual(res, 1);   // rawmult -> rescaleBy "2"
    
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
    
    extscheme.reScaleByAndEqualMT(extres, 1);
    res = extscheme.DecompKeySwitchMT(extres);
    //scheme.reScaleByAndEqual(res, 1);
    
    delete[] tensoring;
}

