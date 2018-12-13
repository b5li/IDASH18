#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/resource.h>   //! check the memory usage

#define USE_NTL 1

#ifdef USE_NTL
#include <NTL/BasicThreadPool.h>
#endif

 
#include "../src/Ciphertext.h"
#include "../src/EvaluatorUtils.h"

#include "Database.h"
#include "BasicTest.h"

 
#include "TestLogRegPvals.h"
  
#include "CipherLogRegPvals.h"


//! encYXData.l = L (fresh)
//! encXData.l = L - 9 (= beta.lvl)
//! encYData.l = L - 13 (= encPr.lvl)
//! enccovData.l = L - 14 = w.lvl

void CipherLRPvals::encXData(Ciphertext& encYXData, Ciphertext& encXData, Ciphertext& encYData, Ciphertext*& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long nXbatching, long nCovbatching, long nterms, long YXscaleBits, long covscaleBits, long nslots, long YXlvl, long Xlvl, long Ylvl, long Covlvl, long Covlvl1) {
    
    double* YXtmp = new double[nslots];
    double* Xtmp = new double[nslots];
    double* Ytmp = new double[nslots];
    long YXscalefactor = (1 << YXscaleBits);
    
    //! encoding for YXData, XData, YData
    long nslots1 = factorDim * nXbatching;  // total number of slots for xData[i], (1 << xBatchingBits)
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        long start = nslots1 * i;
        if(yData[i] == 1){
            for(long j = 0; j < nXbatching; ++j){
                long start1 = start + j * factorDim;
                for(long k = 0; k < factorDim; ++k){
                    YXtmp[start1 + k] = (xData[i][k]/YXscalefactor);
                    Xtmp[start1 + k] = xData[i][k];
                    Ytmp[start1 + k] = 1.0;
                }
            }
        }
        else{
            for(long j = 0; j < nXbatching; ++j){
                long start1 = start + j * factorDim;
                for(long k = 0; k < factorDim; ++k){
                    YXtmp[start1 + k] = - (xData[i][k]/YXscalefactor);
                    Xtmp[start1 + k] = xData[i][k];
                    Ytmp[start1 + k] = 0.0;
                }
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    long nslots2 = nslots1 * sampleDim;     //! number of packed msgs in the slots for encoding data
    for(long i = nslots2; i < nslots; ++i){
        YXtmp[i] = 0.0;
        Xtmp[i] = 0.0;
        Ytmp[i] = 0.0;
    }
    //! encryption
    cipherPvals.encFullyPackedVec(encYXData, YXtmp, nslots, YXlvl);
    cipherPvals.encFullyPackedVec(encXData, Xtmp, nslots, Xlvl);
    cipherPvals.encFullyPackedVec(encYData, Ytmp, nslots, Ylvl);
    
    // "+------------------------------------+"
    //!  encryption of covariance
    // "+------------------------------------+"
    double** cov = new double*[sampleDim];
    
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        cov[i] = new double[nterms];
        cipherPvals.computeCov(cov[i], xData[i], factorDim, covscaleBits);  //! dim -> dim * (dim-1)/2
    }
    NTL_EXEC_RANGE_END;
    
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
        cipherPvals.encFullyPackedVec(enccovData[l], xData2[l], nslots, Covlvl);
    }
    NTL_EXEC_RANGE_END;
    
    //! encoding of determinant
    long scalefactor = (1 << covscaleBits);
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
        cipherPvals.encFullyPackedVec(enccovData[30 + l], xData2[l], nslots, Covlvl1);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] YXtmp;
    delete[] Xtmp;
    delete[] Ytmp;
    delete[] cov;
    delete[] xData2;
}


//! encYXData.l = L (fresh)
//! encXData.l = L - 9 (= beta.lvl)
//! encYData.l = L - 13 (= encPr.lvl)
//! enccovData.l = L - 14 = w.lvl

void CipherLRPvals::new_encXData(Ciphertext& encYXData, Ciphertext& encXData, Ciphertext& encYData, Ciphertext*& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long nXbatching, long nCovbatching, long nterms, long YXscaleBits, long covscaleBits, long nslots, long YXlvl, long Xlvl, long Ylvl, long Covlvl, long Covlvl1) {
    
    double* YXtmp = new double[nslots];
    double* Xtmp = new double[nslots];
    double* Ytmp = new double[nslots];
    long YXscalefactor = (1 << YXscaleBits);
    
    //! encoding for YXData, XData, YData
    long nslots1 = factorDim * nXbatching;  // total number of slots for xData[i], (1 << xBatchingBits)
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        long start = nslots1 * i;
        if(yData[i] == 1){
            for(long j = 0; j < nXbatching; ++j){
                long start1 = start + j * factorDim;
                for(long k = 0; k < factorDim; ++k){
                    YXtmp[start1 + k] = (xData[i][k]/YXscalefactor);
                    Xtmp[start1 + k] = xData[i][k];
                    Ytmp[start1 + k] = 1.0;
                }
            }
        }
        else{
            for(long j = 0; j < nXbatching; ++j){
                long start1 = start + j * factorDim;
                for(long k = 0; k < factorDim; ++k){
                    YXtmp[start1 + k] = - (xData[i][k]/YXscalefactor);
                    Xtmp[start1 + k] = xData[i][k];
                    Ytmp[start1 + k] = 0.0;
                }
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    long nslots2 = nslots1 * sampleDim;     //! number of packed msgs in the slots for encoding data
    for(long i = nslots2; i < nslots; ++i){
        YXtmp[i] = 0.0;
        Xtmp[i] = 0.0;
        Ytmp[i] = 0.0;
    }
    //! encryption
    cipherPvals.encFullyPackedVec(encYXData, YXtmp, nslots, YXlvl);
    cipherPvals.encFullyPackedVec(encXData, Xtmp, nslots, Xlvl);
    cipherPvals.encFullyPackedVec(encYData, Ytmp, nslots, Ylvl);
    
    // "+------------------------------------+"
    //!  encryption of covariance
    // "+------------------------------------+"
    double** cov = new double*[sampleDim];
    
    NTL_EXEC_RANGE(sampleDim, first, last);
    for (long i = first; i < last; ++i) {
        cov[i] = new double[nterms];
        cipherPvals.computeCov(cov[i], xData[i], factorDim, covscaleBits);  //! dim -> dim * (dim-1)/2
    }
    NTL_EXEC_RANGE_END;
    
    
    long DiagIndex[4][3][6] = {
        {{4,6,5,5,6,4},{7,7,8,5,5,8},{9,6,6,9,8,8}}, //00: + - + - + -, gen: 0, 1, 2
        {{0,3,2,2,3,0},{7,7,8,2,2,8},{9,3,3,9,8,8}}, //11: + - + - + -, gen: 12, 13, 14
        {{0,3,1,1,3,0},{4,4,6,1,1,6},{9,3,3,9,6,6}}, //22: + - + - + -, gen: 21, 22, 23
        {{0,2,1,1,2,0},{4,4,5,1,1,5},{7,2,2,7,5,5}}, //33: + - + - + -, gen: 27, 28, 29
    };
    
    long offDiagIndex[6][6] = {
        {1,3,2,2,3,1}, {5,5,6,4,4,6}, {8,5,5,8,7,7},
        {5,5,6,1,1,6}, {8,2,2,8,7,7}, {8,2,2,8,5,5},
    };
    
    long nslots3 = 8 * nCovbatching;        //! total number of slots of a xData[i] = 16
    long nslots4 = nslots3 * sampleDim;     //! number of used for encoding data
    double** xData2 = new double*[18];
    
    //! 1. Diagonal components
    NTL_EXEC_RANGE(12, first, last);
    for (long l = first; l < last; ++l){
        xData2[l] = new double[nslots];
        long lrow = (long) l/3;
        long lcol = (l % 3);
        
        for (long i = 0; i < sampleDim; ++i) {
            long start = nslots3 * i;
            for(long j = 0; j < nCovbatching; ++j){
                long start1 = start + j * 8;
                for(long k = 0; k < 6; k += 2){
                    long ind = DiagIndex[lrow][lcol][k];
                    xData2[l][start1 + k] = cov[i][ind];
                    ind = DiagIndex[lrow][lcol][k + 1];
                    xData2[l][start1 + k + 1] = - cov[i][ind];
                }
                xData2[l][start1 + 6] = 0.0;
                xData2[l][start1 + 7] = 0.0;
            }
        }
        for(long i = nslots4; i < nslots; ++i){
            xData2[l][i] = 0.0;
        }
        cipherPvals.encFullyPackedVec(enccovData[l], xData2[l], nslots, Covlvl);
    }
    NTL_EXEC_RANGE_END;
    

    //! 2. Generation: 13 ~ 17
    NTL_EXEC_RANGE(6, first, last);
    for (long l = first; l < last; ++l){
        xData2[l + 12] = new double[nslots];
        for (long i = 0; i < sampleDim; ++i) {
            long start = nslots3 * i;
            for(long j = 0; j < nCovbatching; ++j){
                long start1 = start + j * 8;
                for(long k = 0; k < 6; k += 2){
                    long ind = offDiagIndex[l][k];
                    xData2[l + 12][start1 + k] = cov[i][ind];
                    ind = offDiagIndex[l][k + 1];
                    xData2[l + 12][start1 + k + 1] = - cov[i][ind];
                }
                xData2[l + 12][start1 + 6] = 0.0;
                xData2[l + 12][start1 + 7] = 0.0;
            }
        }
        for(long i = nslots4; i < nslots; ++i){
            xData2[l + 12][i] = 0.0;
        }
        cipherPvals.encFullyPackedVec(enccovData[l + 12], xData2[l + 12], nslots, Covlvl);
    }
    NTL_EXEC_RANGE_END;
    
    
//    NTL_EXEC_RANGE(30, first, last);
//    for (long l = first; l < last; ++l){
//        xData2[l] = new double[nslots];
//        long lrow = (long) l/3;
//        long lcol = (l % 3);
//
//        for (long i = 0; i < sampleDim; ++i) {
//            long start = nslots3 * i;
//            for(long j = 0; j < nCovbatching; ++j){
//                long start1 = start + j * 8;
//                for(long k = 0; k < 6; k += 2){
//                    long ind = Index[lrow][lcol][k];
//                    if((l == 0)||(l == 2)||(l == 4)||(l == 6)||(l == 7)||(l == 9)){
//                        xData2[l][start1 + k] = cov[i][ind];
//                        ind = Index[lrow][lcol][k + 1];
//                        xData2[l][start1 + k + 1] = - cov[i][ind];
//                    }
//                    else{
//                        xData2[l][start1 + k] = - cov[i][ind];
//                        ind = Index[lrow][lcol][k + 1];
//                        xData2[l][start1 + k + 1] = cov[i][ind];
//                    }
//                }
//                xData2[l][start1 + 6] = 0.0;
//                xData2[l][start1 + 7] = 0.0;
//            }
//        }
//        for(long i = nslots4; i < nslots; ++i){
//            xData2[l][i] = 0.0;
//        }
//        cipherPvals.encFullyPackedVec(enccovData[l], xData2[l], nslots, Covlvl);
//    }
//    NTL_EXEC_RANGE_END;
    

    //! encoding of determinant
    long scalefactor = (1 << covscaleBits);
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
        cipherPvals.encFullyPackedVec(enccovData[18 + l], xData2[l], nslots, Covlvl1);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] YXtmp;
    delete[] Xtmp;
    delete[] Ytmp;
    delete[] cov;
    delete[] xData2;
}


//!@ EncSData[i] = Enc(s[i][0], ... , s[i][p-1], - 0 -)
//!@ EncSXData[i][k] = Enc(s[i][0] * X[i][k], ..., s[i][p-1] * X[i][k], - 0 -)

void CipherLRPvals::encryptSData(Ciphertext**& encSData, Ciphertext***& encSXData, double** xData, double** sData,  long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots, long Slvl, long SXlvl) {
    long nslots1 = nsnp - (nencsnp-1) * nslots;   // number of slots in the final ctxts
    
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
        for(long j = 0; j < nencsnp - 1; ++j){
            for(long l = 0; l < nslots; ++l){
                fullvec[i][l] = sData[i][j1] ;
                j1++;
            }
            cipherPvals.encFullyPackedVec(encSData[i][j], fullvec[i], nslots, Slvl);
        }
        for(long l = 0; l < nslots1; ++l){
            sparsevec[i][l] = sData[i][j1];
            j1++;
        }
        cipherPvals.encSparselyPackedVec(encSData[i][nencsnp - 1], sparsevec[i], nslots1, nslots, Slvl);
        
        /***************************************/
        //! 2. encryption of sxData
        sxData[i] = new double*[factorDim];
        for(long k = 0; k < factorDim; ++k){
            sxData[i][k] = new double[nsnp];
            for(long j = 0; j < nsnp; ++j){
                sxData[i][k][j] = sData[i][j] * xData[i][k];
            }
            j1 = 0;
            for(long j = 0; j < nencsnp - 1; ++j){
                for(long l = 0; l < nslots; ++l){
                    fullvec[i][l] = sxData[i][k][j1];
                    j1++;
                }
                cipherPvals.encFullyPackedVec(encSXData[i][k][j], fullvec[i], nslots, SXlvl);
            }
            for(long l = 0; l < nslots1; ++l){
                sparsevec[i][l] = sxData[i][k][j1];
                j1++;
            }
            cipherPvals.encSparselyPackedVec(encSXData[i][k][nencsnp - 1], sparsevec[i], nslots1, nslots, SXlvl);
        }
    }
    NTL_EXEC_RANGE_END;
    
    delete[] sxData;
    delete[] fullvec;
    delete[] sparsevec;
}

void CipherLRPvals::decryptResult(double**& ZSnorm, double**& Snorm, double& det, Ciphertext* encZSnorm, Ciphertext* encSnorm, Ciphertext encDet, long nencsnp, long nslots){
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        ZSnorm[j] = new double[nslots];
        Snorm[j] = new double[nslots];
        
        cipherPvals.decVector(ZSnorm[j], encZSnorm[j], nslots);
        cipherPvals.decVector(Snorm[j], encSnorm[j], nslots);
    }
    NTL_EXEC_RANGE_END;
    cipherPvals.decSingleData(det, encDet);
}

/****************************************************************************************************************/
//! Output: encWData

void CipherLRPvals::encNLGDiteration(Ciphertext& encWData, Ciphertext encYXData, long factorDim,  long sampleDim, long fdimBits, long sdimBits, long xBatchingBits, uint64_t* poly, long numIter){
    
    //double alpha0[3] = {0.01, 1.0001, 1.61812};
    //double alpha1[3] = {1.0001, 1.61812, 2.19361};
    double eta[5] = {0.989901, 0, -0.281783, -0.434061, -0.531076};
    double gamma = 1.0 / sampleDim;
    
    //----------------------------------------
    //! 0th iteration (1 lvl)
    Ciphertext* encW0Data = new Ciphertext[2];   //! sum gamma * 0.5 * XData
    Ciphertext encGrad0;
    encNLGDiteration0(encW0Data, encGrad0, encYXData, sdimBits, xBatchingBits, gamma, eta[0]);
    
    //----------------------------------------
    //! 1th iteration (2 lvl)
    encWData = encW0Data[0];
    Ciphertext encVData = encW0Data[1];
    Ciphertext encGrad2 = extscheme.multByConstMT(encYXData, gamma * scaledsigmoid3[2]);       //! a[2] * x
    scheme.reScaleByAndEqual(encGrad2, 1);
    
//    encNLGDiteration1(encWData, encVData, encYXData, encW0Data[0], encGrad2, poly, fdimBits, sdimBits, xBatchingBits, gamma);

//    for (long iter = 2; iter < numIter - 1 ; ++iter) {
//        encNLGDiteration2(encWData, encVData, encYXData, encW0Data[0], encGrad0, encGrad2, poly, fdimBits, sdimBits, xBatchingBits, gamma, eta[iter]);
//    }
 
    //----------------------------------------
    //! final iteration
    encNLGDiteration_final(encWData, encVData, encYXData, encW0Data[0], encGrad2, poly, fdimBits, sdimBits, xBatchingBits, gamma);

    delete[] encW0Data;
}


//!@ Input: nslots1 (real number of used slots)
//!@        nslots 
uint64_t* CipherLRPvals::generateNLGDAuxPoly(long nslots, long nslots1, long batch) {
	complex<double>* pvals = new complex<double>[nslots];
	for (long j = 0; j < nslots1; j += batch) {
		pvals[j].real(1.0);
	}
	uint64_t* poly = new uint64_t[scheme.context.L << scheme.context.logN];
	scheme.context.encode(poly, pvals, nslots, scheme.context.L);
	delete[] pvals;
	return poly;
}

//!@ Function: 1st iteratin of NLGD (starting point: w = v = 0)
//!@ Output: encWData (L - 1),  encVData (L - 1)
//!@         encGrad = sum (YXData[i]/4)
//!@         encWData = gamma * (a[0] * 4) * sum (YXData[i]/4)
//!@         encVData = (1 - eta) * gamma * (a[0] * 4) * sum (YXData[i]/4)

void CipherLRPvals::encNLGDiteration0(Ciphertext*& encData, Ciphertext& encGrad, Ciphertext encYXData, long sdimBits, long xBatchingBits, double gamma, double eta){
    
    //! encGrad = sum_i YX[i]
    encGrad = encYXData;
    for (long l = xBatchingBits; l < xBatchingBits + sdimBits; ++l) {
        Ciphertext tmp = extscheme.leftRotateFastMT(encGrad, (1 << l));
        scheme.addAndEqual(encGrad, tmp);
    }
    
    NTL_EXEC_RANGE(2, first, last);
    for (long i = first; i < last; ++i){
        if(i == 0){
            encData[0] = extscheme.multByConstMT(encGrad, -(gamma * scaledsigmoid3[0])); //! encW = - gamma * a[0] * (sum_i YX[i])
        }
        else{
            encData[1] = extscheme.multByConstMT(encGrad, -((1.0 - eta) * gamma * scaledsigmoid3[0])); //! v = eta * 0 + (1-eta) * (- g) = - (1 - eta) * gamma * a[0] * (sum_i YX[i])
        }
        scheme.reScaleByAndEqual(encData[i], 1);
    }
    NTL_EXEC_RANGE_END
}

//!@ Function: NLGD
//!@ Output: encWData (w.lvl - 4), encVData (w.lvl - )
//! suppose that eta ~ 0
void CipherLRPvals::encNLGDiteration1(Ciphertext& encWData, Ciphertext& encVData, Ciphertext encYXData, Ciphertext encW0Data, Ciphertext encGrad2, uint64_t* poly, long fdimBits, long sdimBits, long xBatchingBits, double gamma){

    /*******************************************************/
    //! 1. encIP = 1/ 4 * (y^T * X) * W, Allsum over "factorDim"  (v.lvl - 1)
    //!  and  Replicate the inner product to the other slots (v.lvl - 2)
    
    Ciphertext encIP = scheme.modDownTo(encYXData, encVData.l);
    extscheme.multAndEqualMT(encIP, encVData);
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext rot = extscheme.leftRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, rot);
    }
    
    extscheme.multByPolyAndEqualMT(encIP, poly);
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext rot = extscheme.rightRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, rot);
    }
    
    /*******************************************************/
    //! 2. evaluate  "encGrad = sum sig(encIP) * Z[i]" for 0 <= i < n, (w.lvl - 4)
    //!                       = gamma * x * (a[0] + a[1] * IP + a[2] * IP^3)
    Ciphertext encIP2 = extscheme.squareMT(encIP);
    scheme.reScaleByAndEqual(encIP2, 1);
    scheme.addConstAndEqual(encIP2, scaledsigmoid3[1] / scaledsigmoid3[2]); //! IP^2 + a[1]/a[2]
    
    Ciphertext encGrad = scheme.modDownTo(encGrad2, encIP.l);
    extscheme.multAndEqualMT(encGrad, encIP);
    scheme.reScaleByAndEqual(encGrad, 1);           //! (gamma * a[2] * x) * IP
    
    extscheme.multAndEqualMT(encGrad, encIP2);      //! (gamma * a[2] * x * IP) * (IP^2 + a[1]/a[2])
    scheme.reScaleByAndEqual(encGrad, 1);
    
    //! 3. Aggregate over "sampleDim"
    for (long l = xBatchingBits; l < xBatchingBits + sdimBits; ++l) {
        Ciphertext tmp = extscheme.leftRotateFast(encGrad, (1 << l));
        scheme.addAndEqual(encGrad, tmp);
    }
    
    /*******************************************************/
    //! 4. Update
    // w = v - gamma * g, (-> w.lvl - 4)
    // v = (1 - eta) * (v - gamma * grad) +  (eta * w)
    
    Ciphertext encGrad0 = encW0Data; //! constant term, gamma * x * a[0]
    
    scheme.modDownToAndEqual(encVData, encGrad.l);
    Ciphertext Wtmp = scheme.sub(encVData, encGrad);   //! v - (gamma) * g
    
    scheme.modDownToAndEqual(encGrad0, Wtmp.l);
    scheme.addAndEqual(Wtmp, encGrad0);
    
    encVData = Wtmp;
    encWData = Wtmp;
}


//!@ Function: NLGD
//!@ Output: encWData (w.lvl - 4), encVData (w.lvl - 5)
//! suppose that eta ~ 0
void CipherLRPvals::encNLGDiteration2(Ciphertext& encWData, Ciphertext& encVData, Ciphertext encYXData, Ciphertext encW0Data, Ciphertext encGrad0, Ciphertext encGrad2, uint64_t* poly, long fdimBits, long sdimBits, long xBatchingBits, double gamma, double eta){

    /*******************************************************/
    //! 1. encIP = 1/ 4 * (y^T * X) * W, Allsum over "factorDim"  (v.lvl - 1)
    //! and Replicate the inner product to the other slots (v.lvl - 2)
    
    Ciphertext encIP = scheme.modDownTo(encYXData, encVData.l);
    extscheme.multAndEqualMT(encIP, encVData);
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext rot = extscheme.leftRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, rot);
    }
    
    extscheme.multByPolyAndEqualMT(encIP, poly);
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext rot = extscheme.rightRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, rot);
    }
    
    /*******************************************************/
    //! 3. evaluate  "encGrad = sum sig(encIP) * Z[i]" for 0 <= i < n, (w.lvl - 4)
    //!                       = gamma * x * (a[0] + a[1] * IP + a[2] * IP^3)
    
    Ciphertext encIP2 = extscheme.squareMT(encIP);
    scheme.reScaleByAndEqual(encIP2, 1);
    scheme.addConstAndEqual(encIP2, scaledsigmoid3[1] / scaledsigmoid3[2]); //! IP^2 + a[1]/a[2]
    
    Ciphertext* encGrad = new Ciphertext[2];
    encGrad[0] = encGrad2;  //! gamma * a[2] * YX
    encGrad[1] = extscheme.multByConstMT(encYXData, (1.0 - eta) * gamma * scaledsigmoid3[2]);
    scheme.reScaleByAndEqual(encGrad[1], 1);
    
    NTL_EXEC_RANGE(2, first, last);
    for(long i = first; i < last; ++i){
        scheme.modDownToAndEqual(encGrad[i], encIP.l);  //! (gamma * a[2] * x) * IP
        extscheme.multAndEqualMT(encGrad[i], encIP);
        scheme.reScaleByAndEqual(encGrad[i], 1);
        
        extscheme.multAndEqualMT(encGrad[i], encIP2);     //! (gamma * a[2] * x * IP) * (IP^2 + a[1]/a[2])
        scheme.reScaleByAndEqual(encGrad[i], 1);
    
        for (long l = xBatchingBits; l < xBatchingBits + sdimBits; ++l) {
            Ciphertext tmp = extscheme.leftRotateFastMT(encGrad[i], (1 << l));
            scheme.addAndEqual(encGrad[i], tmp);
        }
    }
    NTL_EXEC_RANGE_END

#if defined(__DEBUG_)
    double* res = new double[4];
    cipherPvals.decVector(res, encGrad[0], 4);
    for(long l = 0; l < 4; ++l) cout << res[l] << "," ;
    cout << endl;
    
    cipherPvals.decVector(res, encGrad[1], 4);
    for(long l = 0; l < 4; ++l) cout << res[l] << "," ;
    cout << endl;
#endif
    
    /*******************************************************/
    //! 4. Update
    //! w = v - gamma * g, (-> w.lvl - 4)
    Ciphertext Wtmp = scheme.modDownTo(encVData, encGrad[0].l);
    scheme.subAndEqual(Wtmp, encGrad[0]);
    
    Ciphertext tmp0 = encW0Data;
    scheme.modDownToAndEqual(tmp0, Wtmp.l);
    scheme.addAndEqual(Wtmp, tmp0);     //! v - (gamma) * g + (- gamma * g0)
    
    //  (1-eta) * v - (1- eta) * (gamma) * g + (eta * w)
    extscheme.multByConstAndEqualMT(encVData, (1.0 - eta));
    scheme.reScaleByAndEqual(encVData, 1);
    
    scheme.modDownToAndEqual(encVData, encGrad[1].l);
    scheme.subAndEqual(encVData, encGrad[1]);
    
    tmp0 = extscheme.multByConstMT(encGrad0,  -((1.0 - eta) * gamma * scaledsigmoid3[0]));
    scheme.reScaleByAndEqual(tmp0, 1); //! sum (1-eta) * gamma * a[0] * YX
    scheme.modDownToAndEqual(tmp0, encVData.l);
    scheme.addAndEqual(encVData, tmp0);
    
    tmp0 = extscheme.multByConstMT(encWData, eta);
    scheme.reScaleByAndEqual(tmp0, 1);
    scheme.modDownToAndEqual(tmp0, encVData.l);
    scheme.addAndEqual(encVData, tmp0);
    
    encWData = Wtmp;
    
    delete[] encGrad;
}


//!@ Function: NLGD
//!@ Output: encWData (w.lvl - 4), encVData (w.lvl - 5)
void CipherLRPvals::encNLGDiteration_final(Ciphertext& encWData, Ciphertext& encVData, Ciphertext encYXData, Ciphertext encW0Data, Ciphertext encGrad2, uint64_t* poly, long fdimBits, long sdimBits, long xBatchingBits, double gamma){
    
    /*******************************************************/
    //! 1. encIP = 1/ 4 * (y^T * X) * W, Allsum over "factorDim"  (v.lvl - 1)
    //! and Replicate the inner product to the other slots (v.lvl - 2)
    
    Ciphertext encIP = scheme.modDownTo(encYXData, encVData.l);
    extscheme.multAndEqualMT(encIP, encVData);
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext rot = extscheme.leftRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, rot);
    }
    
    extscheme.multByPolyAndEqualMT(encIP, poly);
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext rot = extscheme.rightRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, rot);
    }
    
    /*******************************************************/
    //! 2. evaluate  "encGrad = sum sig(encIP) * Z[i]" for 0 <= i < n, (w.lvl - 4)
    //!                       = gamma * x * (a[0] + a[1] * IP + a[2] * IP^3)
    Ciphertext encIP2 = extscheme.squareMT(encIP);
    scheme.reScaleByAndEqual(encIP2, 1);
    scheme.addConstAndEqual(encIP2, scaledsigmoid3[1] / scaledsigmoid3[2]); //! IP^2 + a[1]/a[2]
    
    Ciphertext encGrad = scheme.modDownTo(encGrad2, encIP.l);              //! (gamma * a[2] * x) * IP
    extscheme.multAndEqualMT(encGrad, encIP);
    scheme.reScaleByAndEqual(encGrad, 1);
    
    extscheme.multAndEqualMT(encGrad, encIP2);                    //! (gamma * a[2] * x * IP) * (IP^2 + a[1]/a[2])
    scheme.reScaleByAndEqual(encGrad, 1);
    
    //! 3. Aggregate over "sampleDim"
    for (long l = xBatchingBits; l < xBatchingBits + sdimBits; ++l) {
        Ciphertext tmp = extscheme.leftRotateFastMT(encGrad, (1 << l));
        scheme.addAndEqual(encGrad, tmp);
    }

    //! 4. Update:  w = v - gamma * g, (-> w.lvl - 4)
    Ciphertext encGrad0 = encW0Data; //! constant term, gamma * x * a[0]
    
    scheme.modDownToAndEqual(encVData, encGrad.l);
    encWData = scheme.sub(encVData, encGrad);
    
    scheme.modDownToAndEqual(encGrad0, encWData.l);
    scheme.addAndEqual(encWData, encGrad0);
}

void CipherLRPvals::encZWData(Ciphertext& encWData, Ciphertext& encZWData, Ciphertext encBeta, Ciphertext encXData, Ciphertext encYData, uint64_t* poly, long fdimBits, long sdeg, long scale){
    
    //! 1. encIP = XData * beta (b.lvl - 2 = 14 - 2 = 12)
    //! If sdeg = 5 or 7, replace XData into (XData / scale)
    //! Then the input of sigmoid function is (IP / scale)
    
    Ciphertext encIP;
    if(sdeg == 3){
        encIP = extscheme.multMT(encXData, encBeta);
    }
    else{
        encIP = extscheme.multByConstMT(encXData, (double) 1/scale);
        scheme.reScaleByAndEqual(encIP, 1);
        extscheme.multAndEqualMT(encIP, encBeta);
    }
    
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext rot = extscheme.leftRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, rot);
    }
    
    extscheme.multByPolyAndEqualMT(encIP, poly);
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext tmp = extscheme.rightRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, tmp);
    }

    /**************************************************/
    //! 2. evaluate  "encPr = sig(XData * beta)" for 0 <= i < n, (b.lvl - 4),
    Ciphertext encIP2 = extscheme.squareMT(encIP);
    scheme.reScaleByAndEqual(encIP2, 1);
    
    Ciphertext encPr;
    Ciphertext tmp;
    
    switch(sdeg){
        case 3: //! a0 + (a2 * IP) * (a1/a2 + IP^2)
            scheme.addConstAndEqual(encIP2, sigmoid3[1] / sigmoid3[2]);
            encPr = extscheme.multByConstMT(encIP, sigmoid3[2]);
            scheme.reScaleByAndEqual(encPr, 1);
            extscheme.multAndEqualMT(encPr, encIP2);
            scheme.reScaleByAndEqual(encPr, 1);
            scheme.addConstAndEqual(encPr, sigmoid3[0]);
            break;
            
        case 5: //! (a0 + a1 * IP) + (a3 * IP^2) * (IP^2 + a2/a3)
            encPr = scheme.addConst(encIP2, scaledsigmoid5[2] / scaledsigmoid5[3]);
            extscheme.multByConstAndEqualMT(encIP2, scaledsigmoid5[3]); //! 0.04
            scheme.reScaleByAndEqual(encIP2, 1);
            scheme.modDownToAndEqual(encPr, encIP2.l);
            extscheme.multAndEqualMT(encPr, encIP2);
            scheme.reScaleByAndEqual(encPr, 1);
            
            tmp = extscheme.multByConstMT(encIP, scaledsigmoid5[1]);    //! 0.765
            scheme.reScaleByAndEqual(tmp, 1);
            scheme.addConstAndEqual(tmp, scaledsigmoid5[0]);
            
            scheme.modDownToAndEqual(tmp, encPr.l);
            scheme.addAndEqual(encPr, tmp);
            break;
            
        case 7: //! (a0 + a1) * IP + (a4 * IP3) * (IP4 + a3/a4 * IP2 + a2/a4)
            
            encPr = extscheme.squareMT(encIP2);
            scheme.reScaleByAndEqual(encPr, 1);
            
            tmp = extscheme.multByConstMT(encIP2, (double) scaledsigmoid7[3] / scaledsigmoid7[4]);   //! -0.542
            scheme.reScaleByAndEqual(tmp, 1);
            
            scheme.addAndEqual(encPr, tmp);
            scheme.addConstAndEqual(encPr, scaledsigmoid7[2] / scaledsigmoid7[4]);  //! 0.10454
            
            tmp = extscheme.multByConstMT(encIP, (long) scaledsigmoid7[4] - 1);         //! -320.998 ~ -321
            scheme.modDownToAndEqual(tmp, encIP2.l);
            //scheme.reScaleByAndEqual(tmp, 1);
            extscheme.multAndEqualMT(tmp, encIP2);
            scheme.reScaleByAndEqual(tmp, 1);
            extscheme.multAndEqualMT(encPr, tmp);
            scheme.reScaleByAndEqual(encPr, 1);
            
            tmp = extscheme.multByConstMT(encIP, scaledsigmoid7[1]);                //! (p * 3.47)
            scheme.reScaleByAndEqual(tmp, 1);
            scheme.addConstAndEqual(tmp, scaledsigmoid7[0]);
            scheme.modDownToAndEqual(tmp, encPr.l);
            scheme.addAndEqual(encPr, tmp);
            break;
    }
    
    /**************************************************/
    //! 3. encW = p * (1 - p), (b.lvl - 5)
    //!    encW2 = enc(W^2),   (b.lvl - 6)
    
    tmp = scheme.subConst(encPr, 1.0);
    scheme.negateAndEqual(tmp);
    encWData = extscheme.multMT(encPr, tmp);
    scheme.reScaleByAndEqual(encWData, 1);
    
    /**************************************************/
    //! 4. encZW = z * p * (1 - p) = IP * W + (y - Pr),
    
    if(sdeg != 3){
        extscheme.multByConstAndEqualMT(encIP, (long) scale);   //! real IP
        //scheme.reScaleByAndEqual(encIP, 1);
    }
    
    encZWData = scheme.modDownTo(encIP, encWData.l);
    extscheme.multAndEqualMT(encZWData, encWData);
    scheme.reScaleByAndEqual(encZWData, 1);
    
    tmp = scheme.sub(encYData, encPr); //! encryption level
    scheme.modDownToAndEqual(tmp, encZWData.l);
    scheme.addAndEqual(encZWData, tmp);
    
#if  defined(__DEBUG_)
    double* res = new double[16 * 5];
    cipherPvals.decVector(res, encPr, 16 * 5);
    cout << "Pr : ";
    for(long l = 0; l < 16 * 5; l += 16) cout << res[l] << "," ;
    cout << endl;
    
    cipherPvals.decVector(res, encWData, 16 * 5);
    cout << "W : ";
    for(long l = 0; l < 16 * 5; l += 16) cout << res[l] << "," ;
    cout << endl;
    
    cipherPvals.decVector(res, encZWData, 16 * 5);
    cout << "ZW : ";
    for(long l = 0; l < 16 * 5; l += 16) cout << res[l] << "," ;
    cout << endl;
#endif
}



//**************************************************************************/
//!@ encCov[34] = sum_{0<= i < n} W[i] * cov[i]
//!@ and compute the adjoint/det
void CipherLRPvals::encAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext encWData, Ciphertext* enccovData, long sdimBits, long nCovbatching){
    
    long nslots1 = nCovbatching * 8;
    Ciphertext* encCov = new Ciphertext[30];
    
    NTL_EXEC_RANGE(30, first, last);
    for(long l = first; l < last; ++l){
        encCov[l] = scheme.modDownTo(enccovData[l], encWData.l);
        extscheme.multAndEqual(encCov[l], encWData);
        scheme.reScaleByAndEqual(encCov[l], 1);
        
        //! Aggregate over n ( = cipherPvals.aggCovData_DecompKS)
        for(long i = 0; i < sdimBits; ++i){
            Ciphertext ctemp = extscheme.leftRotateFast(encCov[l], (1 << i) * nslots1);     //! AllSum
            scheme.addAndEqual(encCov[l], ctemp);
        }
    }
    NTL_EXEC_RANGE_END;
   
    //---------------------------------------
    //! encAdj[0], ..., encAdj[9]: 5 -> 3
    
    NTL_EXEC_RANGE(10, first, last);
    for (long l = first; l < last; ++l){
        ExtCiphertext tmp = extscheme.rawmult3(encCov[3 * (l)], encCov[3 * (l) + 1], encCov[3 * (l) + 2]);
        extscheme.reScaleByAndEqual(tmp, 2);
        encAdj[l] = extscheme.DecompKeySwitch(tmp);
        
        for(long j = 0; j < 3; ++j){
            Ciphertext tmp1 = extscheme.leftRotateFast(encAdj[l], (1 << j));
            scheme.addAndEqual(encAdj[l], tmp1);
        }
    }
    NTL_EXEC_RANGE_END;
    
    //---------------------------------------
    //! determinant
    Ciphertext tmp = scheme.modDownTo(encWData, enccovData[30].l);  //! lvl = 4
    
    NTL_EXEC_RANGE(4, first, last);
    for(long l = first; l < last; ++l){
        encCov[l] = enccovData[30 + l];
        extscheme.multAndEqual(encCov[l], tmp);
        scheme.reScaleByAndEqual(encCov[l], 1);
        
        //! Aggregate over n ( = cipherPvals.aggCovData_DecompKS)
        for(long i = 0; i < sdimBits; ++i){
            Ciphertext ctemp = extscheme.leftRotateFast(encCov[l], (1 << i) * nslots1);     //! AllSum
            scheme.addAndEqual(encCov[l], ctemp);
        }
    }
    NTL_EXEC_RANGE_END;
    
    ExtCiphertext* encDet1 = new ExtCiphertext[4];
    NTL_EXEC_RANGE(4, first, last);
    for (long i = first; i < last; ++i){
        Ciphertext tmp = scheme.modDownTo(encCov[i], encAdj[i].l); //! lvl = 3
        encDet1[i] = extscheme.rawmult(tmp, encAdj[i]);
        extscheme.reScaleByAndEqual(encDet1[i], 1);
    }
    NTL_EXEC_RANGE_END;
    
    for (long i = 1; i < 4; ++i){
        extscheme.addAndEqual(encDet1[0], encDet1[i]);
    }
    encDet = extscheme.DecompKeySwitchMT(encDet1[0]);
    
    delete[] encCov;
    delete[] encDet1;
}


//!@ encCov[34] = sum_{0<= i < n} W[i] * cov[i],
//!@ and compute the adjoint/det
void CipherLRPvals::new_encAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext encWData, Ciphertext* enccovData, long sdimBits, long nCovbatching){
    
    long nslots1 = nCovbatching * 8;
    Ciphertext* encCov = new Ciphertext[18];
    
    NTL_EXEC_RANGE(18, first, last);
    for(long l = first; l < last; ++l){
        encCov[l] = scheme.modDownTo(enccovData[l], encWData.l);
        extscheme.multAndEqual(encCov[l], encWData);
        scheme.reScaleByAndEqual(encCov[l], 1);
        
        //! Aggregate over n ( = cipherPvals.aggCovData_DecompKS)
        for(long i = 0; i < sdimBits; ++i){
            Ciphertext ctemp = extscheme.leftRotateFast(encCov[l], (1 << i) * nslots1);     //! AllSum
            scheme.addAndEqual(encCov[l], ctemp);
        }
    }
    NTL_EXEC_RANGE_END;
    
    
    //---------------------------------------
    //! diagonal: encAdj[0],[4],[7],[9]: 5 -> 3
    
    ExtCiphertext* extmp = new ExtCiphertext[3];
    
    long index[3][3] = {{1,2,0}, {3,5,4}, {6,7,8}};
    
    NTL_EXEC_RANGE(4, first, last);
    for (long l = first; l < last; ++l){
        ExtCiphertext tmp;
        if(l == 3){
            tmp = extscheme.rawmult3(encCov[9], encCov[10], encCov[11]);
        }
        else{
            long index0 = index[l][0];
            long index1 = index[l][1];
            long index2 = index[l][2];
            extmp[l] = extscheme.rawmult(encCov[index0], encCov[index1]);
            tmp = extscheme.rawmult(extmp[l], encCov[index2]);
        }
        extscheme.reScaleByAndEqual(tmp, 2);
        
        long l1 = (l * (9 - l))/2;
        encAdj[l1] = extscheme.DecompKeySwitch(tmp);
        
        for(long j = 0; j < 3; ++j){
            Ciphertext tmp1 = extscheme.leftRotateFast(encAdj[l1], (1 << j));
            scheme.addAndEqual(encAdj[l1], tmp1);
        }
    }
    NTL_EXEC_RANGE_END;
   
    //! encAdj[l1], l1 = 1, 5, 8
    //! l1 = 1: 12 * (1 * 2)
    //! l2 = 5 : 3 * 15 * 5
    //! l3 = 8: (6 * 7) * 17
    NTL_EXEC_RANGE(3, first, last);
    for (long l = first; l < last; ++l){
        long ind = 12 + (l * (7 - l))/2;
        ExtCiphertext tmp = extscheme.rawmult(extmp[l], encCov[ind]);
        extscheme.reScaleByAndEqual(tmp, 2);
        
        long l1 = 1 + (l * (9 - l))/2;
        encAdj[l1] = extscheme.DecompKeySwitch(tmp);
        for(long j = 0; j < 3; ++j){
            Ciphertext tmp1 = extscheme.leftRotateFast(encAdj[l1], (1 << j));
            scheme.addAndEqual(encAdj[l1], tmp1);
        }
        scheme.negateAndEqual(encAdj[l1]);
    }
    NTL_EXEC_RANGE_END;
    
    //! encAdj[l1], l1 = 2, 3, 6
    extmp[0] = extscheme.rawmult(encCov[12], encCov[13]);
    
    NTL_EXEC_RANGE(3, first, last);
    for (long l = first; l < last; ++l){
        ExtCiphertext tmp;
        if(l == 0) {
            tmp = extscheme.rawmult(extmp[0], encCov[2]);
        }
        else if(l == 1) {
            tmp = extscheme.rawmult(extmp[0], encCov[14]);
        }
        else{
            tmp = extscheme.rawmult3(encCov[3], encCov[15], encCov[16]);
        }
        extscheme.reScaleByAndEqual(tmp, 2);
        
        long l1 = 2 + l * l;
        encAdj[l1] = extscheme.DecompKeySwitch(tmp);
        for(long j = 0; j < 3; ++j){
            Ciphertext tmp1 = extscheme.leftRotateFast(encAdj[l1], (1 << j));
            scheme.addAndEqual(encAdj[l1], tmp1);
        }
    }
    NTL_EXEC_RANGE_END;
    scheme.negateAndEqual(encAdj[3]);
    
    //---------------------------------------
    //! determinant
    Ciphertext tmp = scheme.modDownTo(encWData, enccovData[0].l);  //! lvl = 4
    
    NTL_EXEC_RANGE(4, first, last);
    for(long l = first; l < last; ++l){
        encCov[l] = enccovData[18 + l];
        extscheme.multAndEqual(encCov[l], tmp);
        scheme.reScaleByAndEqual(encCov[l], 1);
        
        //! Aggregate over n ( = cipherPvals.aggCovData_DecompKS)
        for(long i = 0; i < sdimBits; ++i){
            Ciphertext ctemp = extscheme.leftRotateFast(encCov[l], (1 << i) * nslots1);     //! AllSum
            scheme.addAndEqual(encCov[l], ctemp);
        }
    }
    NTL_EXEC_RANGE_END;
    
    ExtCiphertext* encDet1 = new ExtCiphertext[4];
    NTL_EXEC_RANGE(4, first, last);
    for (long i = first; i < last; ++i){
        Ciphertext tmp = scheme.modDownTo(encCov[i], encAdj[i].l); //! lvl = 3
        encDet1[i] = extscheme.rawmult(tmp, encAdj[i]);
        extscheme.reScaleByAndEqual(encDet1[i], 1);
    }
    NTL_EXEC_RANGE_END;
    
    for (long i = 1; i < 4; ++i){
        extscheme.addAndEqual(encDet1[0], encDet1[i]);
    }
    encDet = extscheme.DecompKeySwitchMT(encDet1[0]);
    
    delete[] encCov;
    delete[] encDet1;
}



//! encYX[k] = E((y^T * X)[k])
void CipherLRPvals::encZXData(Ciphertext*& encZX, Ciphertext encXData, Ciphertext encZData, long sdimBits, long nbatching, long factorDim, long nslots){
    long nslots1 = nbatching * factorDim;  // nbatching = replicated number of a user's data in a single ciphertext
    
    Ciphertext encZXData = scheme.modDownTo(encXData, encZData.l);
    extscheme.multAndEqualMT(encZXData, encZData);
    scheme.reScaleByAndEqual(encZXData, 1);
    
    //! Allsum over "sampleDim"
    for(long i = 0; i < sdimBits; ++i){
        Ciphertext ctemp = extscheme.leftRotateFastMT(encZXData, (1 << i) * nslots1);
        scheme.addAndEqual(encZXData, ctemp);
    }
    
    //! full replication of size "4" (fullReplicate4(encZX, encZXData, nslots);)
    complex<double>** pvals = new complex<double>*[4];
    long lvl = encZXData.l;
    //uint64_t** poly = new uint64_t*[4];
    
    NTL_EXEC_RANGE(4, first, last);
    for (long i = first; i < last; ++i) {
        pvals[i] = new complex<double>[nslots];
        for(long j = 0; j < nslots ; j += 4){
            pvals[i][j + i].real(1.0);
        }
        uint64_t* poly = new uint64_t[lvl << scheme.context.logN];
        scheme.context.encode(poly, pvals[i], nslots, lvl);
        
        encZX[i] = encZXData;
        scheme.multByPolyAndEqual(encZX[i], poly); // (1,0,0,0,|1,0,0,0|....)
        scheme.reScaleByAndEqual(encZX[i], 1);
        
        for(long j = 0; j < 2; ++j){
            Ciphertext ctemp = extscheme.leftRotateFast(encZX[i], (1 << j));
            scheme.addAndEqual(encZX[i], ctemp);
        }
    }
    NTL_EXEC_RANGE_END;
    
    delete[] pvals;
    //delete[] poly;
}

//!@ Firstly replicate, raw-multiply, RS and KS
//!@ suppose that encVec.l = encSData.l + 5
//!@ replicate(encVec.l) = encSData.l
//!@ Output: encZS[nencsnp]

void CipherLRPvals::encVecSData(Ciphertext*& encZS, Ciphertext encZData, Ciphertext** encSData, uint64_t** poly0, uint64_t** poly, long sampleDim, long nencsnp,  long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot){
    
    //cout << encSData[0][0].l << "," << (encVec.l - 1) << endl;
    // "+------------------------------------+"
    //! 1. first step, 1 lvl, (sampleDim/subblocksize) * logrot Rot
    //!    second step: Enc(z0,..., z15) -> Enc(z0),...,Enc(z15)
    // "+------------------------------------+"
    
    ExtCiphertext** res = new ExtCiphertext*[niter];
    
    long logrot = (long) ceil(log2(nslots/nblock));
    long niter1 = niter - 1;
    
    NTL_EXEC_RANGE(niter, first, last);
    for (long i = first; i < last; ++i) {
        //! first layer
        Ciphertext* tmp = new Ciphertext[subblocksize];
        tmp[0] = encZData;
        scheme.multByPolyAndEqual(tmp[0], poly0[i]);
        scheme.reScaleByAndEqual(tmp[0], 1);
        
        for(long j = 0; j < logrot; ++j){
            Ciphertext ctemp = extscheme.leftRotateFast(tmp[0], (1 << j) * nblock);
            scheme.addAndEqual(tmp[0], ctemp);
        }
        
        //! second layer: res[niter][nencsnp]
        res[i] = new ExtCiphertext[nencsnp];
        if(i != niter1){
            fullReplicate16(tmp, tmp[0], poly, rot);      //! subReplicate(tmp[i], tmp[i][0], poly, subblocksize, nstep, rot);
            long i1 = (i * subblocksize);
            for(long j = 0; j < nencsnp; ++j){
                res[i][j] = extscheme.rawmult(encSData[i1][j], tmp[0]);
                for(long l = 1; l < subblocksize; ++l){
                    Ciphertext tmp1 = scheme.modDownTo(encSData[i1 + l][j], tmp[l].l);
                    ExtCiphertext extmp = extscheme.rawmult(tmp1, tmp[l]);
                    extscheme.addAndEqual(res[i][j], extmp);
                }
            }
        }
        
        //! final sub-block
        else{
            long nvals = sampleDim - niter1 * subblocksize;     //! number of final slots
            sparseReplicate16(tmp, tmp[0], poly, nslots, nvals, rot);
            long i1 = (i * subblocksize);
            for(long j = 0; j < nencsnp; ++j){
                res[i][j] = extscheme.rawmult(encSData[i1][j], tmp[0]);
                for(long l = 1; l < nvals; ++l){
                    Ciphertext tmp1 = scheme.modDownTo(encSData[i1 + l][j], tmp[l].l);
                    ExtCiphertext extmp = extscheme.rawmult(tmp1, tmp[l]);
                    extscheme.addAndEqual(res[i][j], extmp);
                }
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! Aggregate sum res[i] * S[i]
    for(long j = 0; j < nencsnp; ++j){
        for(long i = 1; i < niter; ++i){
            extscheme.addAndEqual(res[0][j], res[i][j]);
        }
        extscheme.reScaleByAndEqual(res[0][j], 1);
        encZS[j] = extscheme.DecompKeySwitchMT(res[0][j]);
    }
    //delete[] tmp;
    delete[] res;
}



//! SWS[k] = encSWX[k]
void CipherLRPvals::encVecSData(Ciphertext** encSWX, Ciphertext encWData, Ciphertext*** encSXData, uint64_t** poly0, uint64_t** poly, long factorDim, long sampleDim, long nencsnp, long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot){
    
    //ExtCiphertext** res = new ExtCiphertext*[niter];
    ExtCiphertext*** res1 = new ExtCiphertext**[niter];
    
    long logrot = (long) ceil(log2(nslots/nblock));
    long niter1 = niter - 1;
    
    NTL_EXEC_RANGE(niter, first, last);
    for (long i = first; i < last; ++i) {
        //! first layer:
        Ciphertext* tmp = new Ciphertext[subblocksize];
        tmp[0] = encWData;
        scheme.multByPolyAndEqual(tmp[0], poly0[i]);
        scheme.reScaleByAndEqual(tmp[0], 1);
        
        for(long j = 0; j < logrot; ++j){
            Ciphertext ctemp = extscheme.leftRotateFast(tmp[0], (1 << j) * nblock);
            scheme.addAndEqual(tmp[0], ctemp);
        }
        
        //! second layer: res[niter][nencsnp]
        //! second layer: res[niter][nencsnp][factorDim]
        //res[i] = new ExtCiphertext[nencsnp];
        res1[i] = new ExtCiphertext*[nencsnp];
        if(i != niter1){
            fullReplicate16(tmp, tmp[0], poly, rot);      //! subReplicate(tmp[i], tmp[i][0], poly, subblocksize, nstep, rot); lvl(tmp) = 4
            
            long i1 = (i * subblocksize);
            for(long j = 0; j < nencsnp; ++j){
                res1[i][j] = new ExtCiphertext[factorDim];
                for(long k = 0; k < factorDim; ++k){
                    res1[i][j][k] = extscheme.rawmult(encSXData[i1][k][j], tmp[0]);
                    for(long l = 1; l < subblocksize; ++l){
                        Ciphertext tmp1 = scheme.modDownTo(encSXData[i1 + l][k][j], tmp[l].l);
                        ExtCiphertext extmp = extscheme.rawmult(tmp1, tmp[l]);
                        extscheme.addAndEqual(res1[i][j][k], extmp);
                    }
                }
            }
        }
        
        //! final sub-block
        else{
            long nvals = sampleDim - niter1 * subblocksize;                          //! number of final slots
            sparseReplicate16(tmp, tmp[0], poly, nslots, nvals, rot);
            
            long i1 = (i * subblocksize);
            for(long j = 0; j < nencsnp; ++j){
                res1[i][j] = new ExtCiphertext[factorDim];
                for(long k = 0; k < factorDim; ++k){
                    res1[i][j][k] = extscheme.rawmult(encSXData[i1][k][j], tmp[0]);
                    
                    for(long l = 1; l < nvals; ++l){
                        Ciphertext tmp1 = scheme.modDownTo(encSXData[i1 + l][k][j], tmp[l].l);
                        ExtCiphertext extmp = extscheme.rawmult(tmp1, tmp[l]);
                        extscheme.addAndEqual(res1[i][j][k], extmp);
                    }
                }
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! Aggregate sum res[i] * S[i]
    for(long j = 0; j < nencsnp; ++j){
        encSWX[j] = new Ciphertext[factorDim];
        NTL_EXEC_RANGE(factorDim, first, last);
        for (long k = first; k < last; ++k) {
            for(long i = 1; i < niter; ++i){
                extscheme.addAndEqual(res1[0][j][k], res1[i][j][k]);
            }
            extscheme.reScaleByAndEqual(res1[0][j][k], 1);
            encSWX[j][k] = extscheme.DecompKeySwitchMT(res1[0][j][k]);
        }
        NTL_EXEC_RANGE_END
    }
    //delete[] tmp;
    delete[] res1;
}



/****************************************************************************************************************/
//!@ Generate auxiliary polynomials for replicate
//! poly0: 1st sub-block
//! poly1: 2nd sub-blocks

void CipherLRPvals::generateRepAuxPoly(uint64_t**& poly0, uint64_t**& poly, long nslots, long niter, long nstep, long* nblock, long* rot, long sdeg){
    
    // poly0lvl = encWData.l = 9
    long nslots_block = nblock[0];
    long polylvl = 9;
    
    if(sdeg != 3){
        polylvl = 7;
    }
    complex<double>** pvals0 = new complex<double>*[niter];
    NTL_EXEC_RANGE(niter, first, last);
    for (long i = first; i < last; ++i) {
        pvals0[i] = new complex<double>[nslots];
        for(long j = 0; j < nslots_block; ++j){
            pvals0[i][i * nslots_block + j].real(1.0);
        }
        poly0[i] = new uint64_t[polylvl << scheme.context.logN];
        scheme.context.encode(poly0[i], pvals0[i], nslots, polylvl);
    }
    NTL_EXEC_RANGE_END;
    delete[] pvals0;
    
    //-----------------------------------------------
    polylvl -= 1;
    complex<double>** pvals = new complex<double>*[nstep];
    NTL_EXEC_RANGE(nstep, first, last);
    for (long l = first; l < last; ++l) {
        long niter1 = niter * (1 << l);
        pvals[l] = new complex<double>[nslots];
        for(long i = 0; i < niter1; ++i){
            long i1 = i * nblock[l];
            for(long j = 0; j < rot[l]; ++j){
                pvals[l][i1 + j].real(1.0);
            }
        }
        poly[l] = new uint64_t[polylvl << scheme.context.logN];
        scheme.context.encode(poly[l], pvals[l], nslots, polylvl);
    }
    NTL_EXEC_RANGE_END;
    delete[] pvals;
}

//!@ Function: Enc(z0,z1,...,z15) -> Enc(z0),Enc(z1),...,Enc(z15)
//!@ niter: iteration number of data in a ciphertext
//! 4 lvls
void CipherLRPvals::fullReplicate16(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long* rot){
    
    Ciphertext* tmp = new Ciphertext[8];
    
    // "+--------------------------------------------------------+"
    //!  0st layer : E(0,1,2,3,4,5,6,7) / E(8,..., 15)
    // "+--------------------------------------------------------+"
    tmp[0] = encData;
    scheme.multByPolyAndEqual(tmp[0], poly[0]);
    scheme.reScaleByAndEqual(tmp[0], 1);
    
    tmp[1] = scheme.modDownTo(encData, tmp[0].l);
    scheme.subAndEqual(tmp[1], tmp[0]);
    
    for(long i = 0; i < 2; ++i){    //! = 8 / 4
        Ciphertext ctemp = extscheme.leftRotateFast(tmp[i], rot[0]);
        scheme.addAndEqual(tmp[i], ctemp);
    }
    
    // "+--------------------------------------------------------+"
    //!  1st layer : E(0,1,2,3), E(4,5,6,7), E(8, -), E(12, -)
    // "+--------------------------------------------------------+"
    
    for(long i = 0; i < 2; ++i){
        long j = (i << 1);
        res[j] = tmp[i];
        scheme.multByPolyAndEqual(res[j], poly[1]);
        scheme.reScaleByAndEqual(res[j], 1);
        
        res[j + 1] = scheme.modDownTo(tmp[i], res[j].l); //! E(2,3)
        scheme.subAndEqual(res[j + 1], res[j]);
    }
    
    for(long i = 0; i < 4; ++i){    // 8/2
        Ciphertext ctemp = extscheme.leftRotateFast(res[i], rot[1]);
        scheme.addAndEqual(res[i], ctemp);
    }
    
    // "+--------------------------------------------------------+"
    //!  2nd layer : E(0,1) E(2,3), ..., E(12,13), E(14,15)
    // "+--------------------------------------------------------+"
    
    for(long i = 0; i < 4; ++i){
        long j = (i << 1);
        tmp[j] = res[i];
        scheme.multByPolyAndEqual(tmp[j], poly[2]);
        scheme.reScaleByAndEqual(tmp[j], 1);
        
        tmp[j + 1] = scheme.modDownTo(res[i], tmp[j].l);
        scheme.subAndEqual(tmp[j + 1], tmp[j]);
    }
    
    for(long i = 0; i < 8; ++i){
        Ciphertext ctemp = extscheme.leftRotateFast(tmp[i], rot[2]);
        scheme.addAndEqual(tmp[i], ctemp);
    }
    
    // "+--------------------------------------------------------+"
    //!  3rd layer : E(0),E(1) ...  E(6) E(7)
    // "+--------------------------------------------------------+"
    
    for(long i = 0; i < 8; ++i){
        long j = (i << 1);
        res[j] = tmp[i];
        scheme.multByPolyAndEqual(res[j], poly[3]);
        scheme.reScaleByAndEqual(res[j], 1);
        
        res[j + 1] = scheme.modDownTo(tmp[i], res[j].l);
        scheme.subAndEqual(res[j + 1], res[j]);
    }
    
    for(long i = 0; i < 16; ++i){
        Ciphertext ctemp = extscheme.leftRotateFast(res[i], rot[3]);
        scheme.addAndEqual(res[i], ctemp);
    }
    
    delete[] tmp;
}



void CipherLRPvals::fullReplicate4(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long* rot){
    
    Ciphertext* tmp = new Ciphertext[2];
    
    // "+--------------------------------------------------------+"
    //!  0st layer : E(0,1) / E(2,3)
    // "+--------------------------------------------------------+"
    tmp[0] = encData;
    scheme.multByPolyAndEqual(tmp[0], poly[0]);
    scheme.reScaleByAndEqual(tmp[0], 1);
    
    tmp[1] = scheme.modDownTo(encData, tmp[0].l);
    scheme.subAndEqual(tmp[1], tmp[0]);
    
    for(long i = 0; i < 2; ++i){    //! = 8 / 4
        Ciphertext ctemp = extscheme.leftRotateFast(tmp[i], rot[0]);
        scheme.addAndEqual(tmp[i], ctemp);
    }
    
    // "+--------------------------------------------------------+"
    //!  1st layer : E(0,1,2,3), E(4,5,6,7), E(8, -), E(12, -)
    // "+--------------------------------------------------------+"
    
    for(long i = 0; i < 2; ++i){
        long j = (i << 1);
        res[j] = tmp[i];
        scheme.multByPolyAndEqual(res[j], poly[1]);
        scheme.reScaleByAndEqual(res[j], 1);
        
        res[j + 1] = scheme.modDownTo(tmp[i], res[j].l); //! E(2,3)
        scheme.subAndEqual(res[j + 1], res[j]);
    }
    
    for(long i = 0; i < 4; ++i){    // 8/2
        Ciphertext ctemp = extscheme.leftRotateFast(res[i], rot[1]);
        scheme.addAndEqual(res[i], ctemp);
    }
    
    delete[] tmp;
}


void CipherLRPvals::sparseReplicate16(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots,  long nvals, long* rot){
    
    Ciphertext* tmp = new Ciphertext[8];
    
    // "+------------------------------------+"
    //!  0st layer : E(0,...,7), E(8,...,15)
    // "+------------------------------------+"
    
    tmp[0] = encData;
    scheme.multByPolyAndEqual(tmp[0], poly[0]);
    scheme.reScaleByAndEqual(tmp[0], 1);
    
    if(nvals > 8){
        tmp[1] = scheme.modDownTo(encData, tmp[0].l);
        scheme.subAndEqual(tmp[1], tmp[0]);
    }
    
    long nctxt = (long) ceil((double)nvals/8.0);
    for(long i = 0; i < nctxt; ++i){
        Ciphertext ctemp = extscheme.leftRotateFast(tmp[i], rot[0]);
        scheme.addAndEqual(tmp[i], ctemp);
    }
    
    // "+------------------------------------+"
    //!  1st layer : E(0,...,3), ,,,, E(12,...,15)
    // "+------------------------------------+"
    nctxt = (long) ceil((double)nvals/4.0);
    long nvals1 = (long) ceil((double)(nctxt - 1)/2.0);
    
    for(long i = 0; i < nvals1; ++i){
        long i1 = (i << 1);
        res[i1] = tmp[i];
        scheme.multByPolyAndEqual(res[i1], poly[1]);
        scheme.reScaleByAndEqual(res[i1], 1);
        
        res[i1 + 1] = scheme.modDownTo(tmp[i], res[i1].l);
        scheme.subAndEqual(res[i1 + 1], res[i1]);
    }
    
    if(nctxt % 2){
        res[nctxt - 1] = tmp[nvals1];
        scheme.multByPolyAndEqual(res[nctxt - 1], poly[1]);
        scheme.reScaleByAndEqual(res[nctxt - 1], 1);
    }
    
    for(long i = 0; i < nctxt; ++i){
        Ciphertext ctemp = extscheme.leftRotateFast(res[i], rot[1]);
        scheme.addAndEqual(res[i], ctemp);
    }
    
    // "+------------------------------------+"
    //!  2nd layer : E(0,1), E(2,3),...
    // "+------------------------------------+"
    nctxt = (long) ceil((double)nvals/2.0);
    nvals1 = (long) ceil((double)(nctxt - 1)/2.0);
    
    for(long i = 0; i < nvals1; ++i){
        long i1 = (i << 1);
        tmp[i1] = res[i];
        scheme.multByPolyAndEqual(tmp[i1], poly[2]);
        scheme.reScaleByAndEqual(tmp[i1], 1);
        
        tmp[i1 + 1] = scheme.modDownTo(res[i], tmp[i1].l);
        scheme.subAndEqual(tmp[i1 + 1], tmp[i1]);
    }
    
    //! odd
    if(nctxt % 2){
        tmp[nctxt - 1] = res[nvals1];
        scheme.multByPolyAndEqual(tmp[nctxt - 1], poly[2]);
        scheme.reScaleByAndEqual(tmp[nctxt - 1], 1);
    }
    
    for(long i = 0; i < nctxt; ++i){
        Ciphertext ctemp = extscheme.leftRotateFast(tmp[i], rot[2]);
        scheme.addAndEqual(tmp[i], ctemp);
    }
    
    // "+------------------------------------+"
    //!  3rd layer : E(0),E(1) ...  E(6) E(7)
    // "+------------------------------------+"
    nvals1 = (long) ceil((double)(nvals - 1)/2.0);
    
    for(long i = 0; i < nvals1; ++i){
        long i1 = (i << 1);
        res[i1] = tmp[i];
        scheme.multByPolyAndEqual(res[i1], poly[1]);
        scheme.reScaleByAndEqual(res[i1], 1);
        
        res[i1 + 1] = scheme.modDownTo(tmp[i], res[i1].l);
        scheme.subAndEqual(res[i1 + 1], res[i1]);
    }
    
    //! odd
    if(nctxt % 2){
        res[nvals - 1] = tmp[nvals1];
        scheme.multByPolyAndEqual(res[nvals - 1], poly[1]);
        scheme.reScaleByAndEqual(res[nvals - 1], 1);
    }
    
    for(long i = 0; i < nvals; ++i){
        Ciphertext ctemp = extscheme.leftRotateFast(res[i], rot[1]);
        scheme.addAndEqual(res[i], ctemp);
    }
    delete[] tmp;
}


void CipherLRPvals::sparseReplicate4(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots,  long nvals, long* rot){
    
    Ciphertext* tmp = new Ciphertext[2];
    
    // "+------------------------------------+"
    //!  0st layer : E(0,...,7), E(8,...,15)
    // "+------------------------------------+"
    
    tmp[0] = encData;
    scheme.multByPolyAndEqual(tmp[0], poly[0]);
    scheme.reScaleByAndEqual(tmp[0], 1);
    
    if(nvals > 2){
        tmp[1] = scheme.modDownTo(encData, tmp[0].l);
        scheme.subAndEqual(tmp[1], tmp[0]);
    }
    
    long nctxt = (long) ceil((double)nvals/2.0);
    for(long i = 0; i < nctxt; ++i){
        Ciphertext ctemp = extscheme.leftRotateFast(tmp[i], rot[0]);
        scheme.addAndEqual(tmp[i], ctemp);
    }
    
    // "+------------------------------------+"
    //!  3rd layer : E(0),E(1) ...  E(6) E(7)
    // "+------------------------------------+"
    long nvals1 = (long) ceil((double)(nvals - 1)/2.0);
    
    for(long i = 0; i < nvals1; ++i){
        long i1 = (i << 1);
        res[i1] = tmp[i];
        scheme.multByPolyAndEqual(res[i1], poly[1]);
        scheme.reScaleByAndEqual(res[i1], 1);
        
        res[i1 + 1] = scheme.modDownTo(tmp[i], res[i1].l);
        scheme.subAndEqual(res[i1 + 1], res[i1]);
    }
    
    if(nctxt % 2){
        res[nvals - 1] = tmp[nvals1];
        scheme.multByPolyAndEqual(res[nvals - 1], poly[1]);
        scheme.reScaleByAndEqual(res[nvals - 1], 1);
    }
    
    for(long i = 0; i < nvals; ++i){
        Ciphertext ctemp = extscheme.leftRotateFast(res[i], rot[1]);
        scheme.addAndEqual(res[i], ctemp);
    }
    delete[] tmp;
}


//! Data0 * mat * Data2
//! Data2 * mat * Data2
void CipherLRPvals::extQuadForm(Ciphertext& res0, Ciphertext& res1, Ciphertext* encData0, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim){
    ExtCiphertext** tensoring0 = new ExtCiphertext*[factorDim];
    ExtCiphertext** tensoring1 = new ExtCiphertext*[factorDim];
    //cout << encData0[0].l << "," << encMatrix[0].l << "," << encData2[0].l << endl;
    
    NTL_EXEC_RANGE(factorDim, first, last);
    for(long i = first; i< last; ++i){
        tensoring0[i] = new ExtCiphertext[factorDim];
        tensoring1[i] = new ExtCiphertext[factorDim];
        
        //! lower
        for(long j = 0; j < i; ++j){
            long l1 = j* factorDim - (j * (j-1))/2 + (i - j);
            
            ExtCiphertext tmp0 = extscheme.rawmult(encData2[j], encMatrix[l1]);
            extscheme.reScaleByAndEqual(tmp0, 1);
            
            Ciphertext tmp = scheme.modDownTo(encData0[i], tmp0.l);
            tensoring0[i][j] = extscheme.rawmult(tmp0, tmp);
            
            tmp = scheme.modDownTo(encData2[i], tmp0.l);
            tensoring1[i][j] = extscheme.rawmult(tmp0, tmp);
        }
        
        //! upper
        long diag_index = i* factorDim - (i * (i-1))/2;
        for(long j = i; j < factorDim; ++j){
            long l1 = diag_index + (j - i);
            
            ExtCiphertext tmp0 = extscheme.rawmult(encData2[j], encMatrix[l1]);
            extscheme.reScaleByAndEqual(tmp0, 1);
            
            Ciphertext tmp = scheme.modDownTo(encData0[i], tmp0.l);
            tensoring0[i][j] = extscheme.rawmult(tmp0, tmp);
            
            tmp = scheme.modDownTo(encData2[i], tmp0.l);
            tensoring1[i][j] = extscheme.rawmult(tmp0, tmp);
        }
        
        for(long j = 1; j <factorDim; ++j){
            extscheme.addAndEqual(tensoring0[i][0], tensoring0[i][j]);
            extscheme.addAndEqual(tensoring1[i][0], tensoring1[i][j]);
        }
    }
    NTL_EXEC_RANGE_END;
    
    ExtCiphertext extres0 = tensoring0[0][0];
    ExtCiphertext extres1 = tensoring1[0][0];
    for(long i = 1; i < factorDim; ++i){
        extscheme.addAndEqual(extres0, tensoring0[i][0]);
        extscheme.addAndEqual(extres1, tensoring1[i][0]);
    }
    
    //! Final KS rescaling
    res0 = extscheme.DecompKeySwitchMT(extres0);
    res1 = extscheme.DecompKeySwitchMT(extres1);
    
    scheme.reScaleByAndEqual(res0, 1);
    scheme.reScaleByAndEqual(res1, 1);
    
    delete[] tensoring0;
    delete[] tensoring1;
}







//-------------------------------------------------------------------

//! Output: encZData, encWData, encW2Data
//! Z[i] = X[i] * beta - (y[i] - Pr[i])/(Pr[i] * (1- Pr[i])),
//! W[i] = Pr[i] * (1- Pr[i]),
//! W2[i] = (W[i])^2,
//! ZW[i] = z[i] * w[i]

void CipherLRPvals::encZWData(Ciphertext& encZData, Ciphertext& encWData, Ciphertext& encW2Data, Ciphertext& encZWData, Ciphertext encBeta, Ciphertext encXData, Ciphertext encYData, uint64_t* poly, long fdimBits, long sdeg, long scale){
    
    //! 1. encIP = XData * beta (b.lvl - 2)
    //! If sdeg = 5 or 7, replace XData into (XData / scale)
    //! Then the input of sigmoid function is (IP / scale)
    
    Ciphertext encIP;
    if(sdeg == 3){
        encIP = extscheme.multMT(encXData, encBeta);
    }
    else{
        encIP = extscheme.multByConstMT(encXData, (double) 1/scale);
        scheme.reScaleByAndEqual(encIP, 1);
        extscheme.multAndEqualMT(encIP, encBeta);
    }
    
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext rot = extscheme.leftRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, rot);
    }
    
    extscheme.multByPolyAndEqualMT(encIP, poly);
    scheme.reScaleByAndEqual(encIP, 1);
    for (long l = 0; l < fdimBits; ++l) {
        Ciphertext tmp = extscheme.rightRotateFastMT(encIP, (1 << l));
        scheme.addAndEqual(encIP, tmp);
    }

    /**************************************************/
    //! 2. evaluate  "encPr = sig(XData * beta)" for 0 <= i < n, (b.lvl - 4),
    Ciphertext encIP2 = extscheme.squareMT(encIP);
    scheme.reScaleByAndEqual(encIP2, 1);
    
    Ciphertext encPr;
    Ciphertext tmp;
    
    switch(sdeg){
        case 3: //! a0 + (a2 * IP) * (a1/a2 + IP^2)
            scheme.addConstAndEqual(encIP2, sigmoid3[1] / sigmoid3[2]);
            encPr = extscheme.multByConstMT(encIP, sigmoid3[2]);
            scheme.reScaleByAndEqual(encPr, 1);
            extscheme.multAndEqualMT(encPr, encIP2);
            scheme.reScaleByAndEqual(encPr, 1);
            scheme.addConstAndEqual(encPr, sigmoid3[0]);
            break;
            
        case 5: //! (a0 + a1 * IP) + (a3 * IP^2) * (IP^2 + a2/a3)
            encPr = scheme.addConst(encIP2, scaledsigmoid5[2] / scaledsigmoid5[3]);
            extscheme.multByConstAndEqualMT(encIP2, scaledsigmoid5[3]); //! 0.04
            scheme.reScaleByAndEqual(encIP2, 1);
            scheme.modDownToAndEqual(encPr, encIP2.l);
            extscheme.multAndEqualMT(encPr, encIP2);
            scheme.reScaleByAndEqual(encPr, 1);
            
            tmp = extscheme.multByConstMT(encIP, scaledsigmoid5[1]);    //! 0.765
            scheme.reScaleByAndEqual(tmp, 1);
            scheme.addConstAndEqual(tmp, scaledsigmoid5[0]);
            
            scheme.modDownToAndEqual(tmp, encPr.l);
            scheme.addAndEqual(encPr, tmp);
            break;
            
        case 7: //! (a0 + a1) * IP + (a4 * IP3) * (IP4 + a3/a4 * IP2 + a2/a4)
             
            encPr = extscheme.squareMT(encIP2);
            scheme.reScaleByAndEqual(encPr, 1);
            
            tmp = extscheme.multByConstMT(encIP2, (double) scaledsigmoid7[3] / scaledsigmoid7[4]);   //! -0.542
            scheme.reScaleByAndEqual(tmp, 1);
            
            scheme.addAndEqual(encPr, tmp);
            scheme.addConstAndEqual(encPr, scaledsigmoid7[2] / scaledsigmoid7[4]);  //! 0.10454
            
            tmp = extscheme.multByConstMT(encIP, (long) scaledsigmoid7[4] - 1);         //! -320.998 ~ -321
            scheme.modDownToAndEqual(tmp, encIP2.l);
            //scheme.reScaleByAndEqual(tmp, 1);
            extscheme.multAndEqualMT(tmp, encIP2);
            scheme.reScaleByAndEqual(tmp, 1);
            extscheme.multAndEqualMT(encPr, tmp);
            scheme.reScaleByAndEqual(encPr, 1);
            
            tmp = extscheme.multByConstMT(encIP, scaledsigmoid7[1]);                //! (p * 3.47)
            scheme.reScaleByAndEqual(tmp, 1);
            scheme.addConstAndEqual(tmp, scaledsigmoid7[0]);
            scheme.modDownToAndEqual(tmp, encPr.l);
            scheme.addAndEqual(encPr, tmp);
            break;
    }
    
    /**************************************************/
    //! 3. encW = p * (1 - p), (b.lvl - 5)
    //!    encW2 = enc(W^2),   (b.lvl - 6)
    
    Ciphertext encPr1 = scheme.subConst(encPr, 1.0);
    scheme.negateAndEqual(encPr1);
    encWData = extscheme.multMT(encPr, encPr1);
    scheme.reScaleByAndEqual(encWData, 1);
    
    encW2Data = extscheme.squareMT(encWData);
    scheme.reScaleByAndEqual(encW2Data, 1);
   
    /**************************************************/
    //! 4. encZW = z * p * (1 - p) = IP * W + (y - Pr),
    
    if(sdeg != 3){
        extscheme.multByConstAndEqualMT(encIP, (long) scale);   //! real IP
        //scheme.reScaleByAndEqual(encIP, 1);
    }

    encZWData = scheme.modDownTo(encIP, encWData.l);
    extscheme.multAndEqualMT(encZWData, encWData);
    scheme.reScaleByAndEqual(encZWData, 1);
    
    encZData = scheme.sub(encYData, encPr); //! encryption level
    scheme.modDownToAndEqual(encZData, encZWData.l);
    scheme.addAndEqual(encZWData, encZData);
    
    /**************************************************/
    //! 5. encW^-1 = 1 / p * (1 - p),
    //! 6. encZ = (y - pr) * Winv + IP,
    //encTwoInverse(encWinv, encPr, encPr1, steps);
    
    Ciphertext encWinv;
    encWinverse(encWinv, encPr);

    long invscale = 16;
    extscheme.multByConstAndEqualMT(encZData, invscale);
    
    scheme.modDownToAndEqual(encZData, encWinv.l);
    extscheme.multAndEqualMT(encZData, encWinv);
    scheme.reScaleByAndEqual(encZData, 1);
    
    tmp = scheme.modDownTo(encIP, encZData.l);
    scheme.addAndEqual(encZData, tmp);
    
#if defined(__DEBUG_)
    double* res = new double[16*5];
    cipherPvals.decVector(res, encPr, 16 * 5);
    cout << "Pr : ";
    for(long l = 0; l < 16 * 5; l += 16) cout << res[l] << "," ;
    cout << endl;
    
    cipherPvals.decVector(res, encWinv, 16 * 5);
    cout << "Winv : ";
    for(long l = 0; l < 16 * 5; l += 16) cout << res[l] << "," ;
    cout << endl;
    
    cipherPvals.decVector(res, encZData, 16 * 5);
    cout << "Z : ";
    for(long l = 0; l < 16 * 5; l += 16) cout << res[l] << "," ;
    cout << endl;

    cipherPvals.decVector(res, encZWData, 16 * 5);
    cout << "ZW : ";
    for(long l = 0; l < 16 * 5; l += 16) cout << res[l] << "," ;
    cout << endl;
#endif
    
}


//! Output: Winv = 1/(pr)*(1-pr)
//! a[0] = 4, a[2] = 0.227599, a[4] = 254.447, a[6] = 1024
//! 1/(1/4 - y^2) = a[0] + a[2] * y2 + a[4] * y4 + a[6] * y6
//!               = a[0] + (a[2]/2^2) * (2*y)^2 + (a[4]/2^4) * (2*y)^4 + (a[6]/2^6) * (2*y)^6, a[6]/2^6 ~ 16
//!               = 16 * (a[0]/16 + (a[2]/2^2)/16 * (2x-1)^2 + (a[4]/2^4)/16 * (2x - 1)^4 + (2x - 1)^6)
//! computation: 3 lvls + 1 lvl (Mod-embeding)

void CipherLRPvals::encWinverse(Ciphertext& encWinv, Ciphertext encPr){
    double a0 = 0.24989289;     //! a[0]/16
    double a2 = 0.266938701;    //! (a[2]/2^2)/16
    double a4 = - 0.01339687;   //! (a[4]/2^4)/16
   
    //! y = (x - 0.5), 2y = 2x -  1
    Ciphertext encX = scheme.subConst(encPr, 0.5);
    extscheme.multByConstAndEqualMT(encX, (long) 2);
    scheme.modDownByAndEqual(encX, 1);
    
    Ciphertext encX2 = extscheme.squareMT(encX);
    scheme.reScaleByAndEqual(encX2, 1);
    
    Ciphertext encX4 = extscheme.squareMT(encX2);
    scheme.reScaleByAndEqual(encX4, 1);
    scheme.modDownToAndEqual(encX2, encX4.l);
    
    encWinv = extscheme.multMT(encX4, encX2);
    scheme.reScaleByAndEqual(encWinv, 1);
    
    extscheme.multByConstAndEqualMT(encX4, a4);
    extscheme.multByConstAndEqualMT(encX2, a2);
    scheme.addAndEqual(encX4, encX2);
    scheme.reScaleByAndEqual(encX4, 1);
    
    scheme.modDownToAndEqual(encX4, encWinv.l);
    scheme.addAndEqual(encWinv, encX4);
    scheme.addConstAndEqual(encWinv, a0);

}

//! encPr1 = p, encPr2 = (1 - p)
//!  Winv = (1 + p) * (1 + p^2) * (1 + p^4) * (1 + p^8) * (1 + p1) * (1 + p1^2) * (1 + p1^4) * (1 + p1^8)   where p1 = 1 - p
//!  consumed lvl: steps

void CipherLRPvals::encTwoInverse(Ciphertext& encWinv, Ciphertext encPr1, Ciphertext encPr2, long steps){
    Ciphertext* res = new Ciphertext[2];
    Ciphertext* cpow = new Ciphertext[2];
    Ciphertext* tmp = new Ciphertext[2];
    
    NTL_EXEC_RANGE(2, first, last);
    for (long i = first; i < last; ++i) {
        if(i == 0){
            cpow[i] = encPr1;
        }
        else{
            cpow[i] = encPr2;
        }
        
        tmp[i] = scheme.addConst(cpow[i], 1.0);
        scheme.modDownByAndEqual(tmp[i], 1);
        res[i] = tmp[i];
        
        for (long j = 1; j < steps; ++j) {
            extscheme.squareAndEqual(cpow[i]);
            scheme.reScaleByAndEqual(cpow[i], 1);
            tmp[i] = scheme.addConst(cpow[i], 1.0);
            extscheme.multAndEqual(tmp[i], res[i]);
            scheme.reScaleByAndEqual(tmp[i], 1);
            res[i] = tmp[i];
        }
    }
    NTL_EXEC_RANGE_END;
    
    encWinv = extscheme.multMT(res[0], res[1]);
    scheme.reScaleByAndEqual(encWinv, 1);
    
    delete[] res;
    delete[] cpow;
    delete[] tmp;
}




//! encSX[k][nencsnp] = sum_{0 <= i < n} SXData[i][k][nencsnp]
void CipherLRPvals::encSXData(Ciphertext**& encSX, Ciphertext*** encSXData, long factorDim, long sampleDim, long nencsnp){
    NTL_EXEC_RANGE(factorDim, first, last);
    for(long k = first; k < last; ++k){
        encSX[k] = new Ciphertext[nencsnp];
        for(long j = 0; j < nencsnp; ++j){
            encSX[k][j] = encSXData[0][k][j];
            for(long i = 1; i < sampleDim; ++i){
                scheme.addAndEqual(encSX[k][j], encSXData[i][k][j]);
            }
        }
    }
    NTL_EXEC_RANGE_END;
}


//! encData0: L - 21, encData1 = L - 21
//! encMatrix : L - 21, encData2 = L - 22

void CipherLRPvals::extQuadForm(Ciphertext& res0, Ciphertext& res1, Ciphertext* encData0, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim){
    ExtCiphertext** tensoring0 = new ExtCiphertext*[factorDim];
    ExtCiphertext** tensoring1 = new ExtCiphertext*[factorDim];
    
    NTL_EXEC_RANGE(factorDim, first, last);
    for(long i = first; i< last; ++i){
        tensoring0[i] = new ExtCiphertext[factorDim];
        tensoring1[i] = new ExtCiphertext[factorDim];
        
        //! lower
        for(long j = 0; j < i; ++j){
            long l1 = j* factorDim - (j * (j-1))/2 + (i - j);
            
            ExtCiphertext tmp0 = extscheme.rawmult(encData0[i], encMatrix[l1]);
            extscheme.reScaleByAndEqual(tmp0, 1);
            tensoring0[i][j] = extscheme.rawmult(tmp0, encData2[j]);
            
            ExtCiphertext tmp1 = extscheme.rawmult(encData1[i], encMatrix[l1]);
            extscheme.reScaleByAndEqual(tmp1, 1);
            tensoring1[i][j] = extscheme.rawmult(tmp1, encData2[j]);
        }
        
        //! upper
        long diag_index = i* factorDim - (i * (i-1))/2;
        for(long j = i; j < factorDim; ++j){
            long l1 = diag_index + (j - i);
            
            ExtCiphertext tmp0 = extscheme.rawmult(encData0[i], encMatrix[l1]);
            extscheme.reScaleByAndEqual(tmp0, 1);
            tensoring0[i][j] = extscheme.rawmult(tmp0, encData2[j]);
            
            ExtCiphertext tmp1 = extscheme.rawmult(encData1[i], encMatrix[l1]);
            extscheme.reScaleByAndEqual(tmp1, 1);
            tensoring1[i][j] = extscheme.rawmult(tmp1, encData2[j]);
        }
        
        for(long j = 1; j <factorDim; ++j){
            extscheme.addAndEqual(tensoring0[i][0], tensoring0[i][j]);
            extscheme.addAndEqual(tensoring1[i][0], tensoring1[i][j]);
        }
    }
    NTL_EXEC_RANGE_END;
    
    ExtCiphertext extres0 = tensoring0[0][0];
    ExtCiphertext extres1 = tensoring1[0][0];
    for(long i = 1; i < factorDim; ++i){
        extscheme.addAndEqual(extres0, tensoring0[i][0]);
        extscheme.addAndEqual(extres1, tensoring1[i][0]);
    }
    
    //! Final KS rescaling
    res0 = extscheme.DecompKeySwitchMT(extres0);
    res1 = extscheme.DecompKeySwitchMT(extres1);
    
    scheme.reScaleByAndEqual(res0, 1);
    scheme.reScaleByAndEqual(res1, 1);
    
    delete[] tensoring0;
    delete[] tensoring1;
}




//!@ Input: encVec, encSXData[n][k][j]
//!@ first replicate: encVec -> encVec[n]
//!@ Output: encZX[nencsnp][k]
void CipherLRPvals::encVecMultipleSData(Ciphertext**& encZS, Ciphertext encVec, Ciphertext*** encSData, uint64_t** poly0, uint64_t** poly, long factorDim, long sampleDim, long nencsnp,  long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot){

    //Ciphertext** tmp = new Ciphertext*[niter];
    ExtCiphertext*** res = new ExtCiphertext**[niter];
    
    long logrot = (long) ceil(log2(nslots/nblock));
    long niter1 = niter - 1;
    
    NTL_EXEC_RANGE(niter, first, last);
    for (long i = first; i < last; ++i) {
        //! first layer
        Ciphertext* tmp = new Ciphertext[subblocksize];
        tmp[0] = encVec;
        scheme.multByPolyAndEqual(tmp[0], poly0[i]);
        scheme.reScaleByAndEqual(tmp[0], 1);
        
        for(long j = 0; j < logrot; ++j){
            Ciphertext ctemp = extscheme.leftRotateFast(tmp[0], (1 << j) * nblock);
            scheme.addAndEqual(tmp[0], ctemp);
        }
        
        //! second layer: res[niter][nencsnp][factorDim]
        res[i] = new ExtCiphertext*[nencsnp];
        if(i != niter1){
            fullReplicate16(tmp, tmp[0], poly, rot);      //! subReplicate(tmp[i], tmp[i][0], poly, subblocksize, nstep, rot);
            
            long i1 = (i * subblocksize);
            for(long j = 0; j < nencsnp; ++j){
                res[i][j] = new ExtCiphertext[factorDim];
                for(long k = 0; k < factorDim; ++k){
                    //Ciphertext tmp1 = scheme.modDownTo(encSData[i1][k][j], tmp[i][0].l);         //! encryption lvl
                    res[i][j][k] = extscheme.rawmult(encSData[i1][k][j], tmp[0]);
                    
                    for(long l = 1; l < subblocksize; ++l){
                        Ciphertext tmp1 = scheme.modDownTo(encSData[i1 + l][k][j], tmp[l].l);
                        ExtCiphertext extmp = extscheme.rawmult(tmp1, tmp[l]);
                        extscheme.addAndEqual(res[i][j][k], extmp);
                    }
                }
            }
        }
        
        //! final sub-block
        else{
            long nvals = sampleDim - niter1 * subblocksize;                          //! number of final slots
            sparseReplicate16(tmp, tmp[0], poly, nslots, nvals, rot);
            
            long i1 = (i * subblocksize);
            for(long j = 0; j < nencsnp; ++j){
                res[i][j] = new ExtCiphertext[factorDim];
                for(long k = 0; k < factorDim; ++k){
                    //Ciphertext tmp1 = scheme.modDownTo(encSData[i1][k][j], tmp[i][0].l);         //! encryption lvl
                    res[i][j][k] = extscheme.rawmult(encSData[i1][k][j], tmp[0]);
                    
                    for(long l = 1; l < nvals; ++l){
                        Ciphertext tmp1 = scheme.modDownTo(encSData[i1 + l][k][j], tmp[l].l);
                        ExtCiphertext extmp = extscheme.rawmult(tmp1, tmp[l]);
                        extscheme.addAndEqual(res[i][j][k], extmp);
                    }
                }
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 0; j < nencsnp; ++j){
        encZS[j] = new Ciphertext[factorDim];
        NTL_EXEC_RANGE(factorDim, first, last);
        for (long k = first; k < last; ++k) {
            for(long i = 1; i < niter; ++i){
                extscheme.addAndEqual(res[0][j][k], res[i][j][k]);
            }
            extscheme.reScaleByAndEqual(res[0][j][k], 1);
            encZS[j][k] = extscheme.DecompKeySwitchMT(res[0][j][k]);
        }
        NTL_EXEC_RANGE_END
    }
    
    delete[] res;
}




