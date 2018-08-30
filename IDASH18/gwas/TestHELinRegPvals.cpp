/*
 * @file       TestHEPvalues.cpp, cpp file
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
#include "TestHELinRegPvals.h"

 
//!@ XData: fully packed in a single ciphertext
//!@ Function: using mod-raising KS -> logN = 14

void TestHEPvals::testHESIMDLinReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp,  string filename){
    
    //! Parameters for approx-HE
    long logp = 45;
    long logN = 14;
    long nslots = (1 << (logN-1));      //! total number of plaintext slots
    long L = 4;                         //! 2 (det,adj) + 1 (multiplied by det,adj) + 1
    long K = L + 1;                     //! need for "correct" rotation
    long h = 90;

    
    //! Parameters for GWAS
    long scaleBits = 4;            //! scale factor for covariance
    long sampleDim2 = (1 << (long)ceil(log2(sampleDim)));   //! closet PoT
    long sdimBits = (long)ceil(log2(sampleDim));            //! log2(sampleDim)
    long factorDim2 = factorDim * factorDim;
    
    long nXbatching = nslots/(sampleDim2 * factorDim);      //! replicated number of a user's Xdata in a single ciphertext
    long nCovbatching = nslots/(sampleDim2 * 8);            //! replicated number of a user's (X^T * X) in a single ciphertext
    long nencsnp = (long)ceil((double)nsnp/nslots);          //! number of ciphertexts for snp encryptions
    long nterms = factorDim * (factorDim + 1)/2;
    
    long logQ = ((Q0_BIT_SIZE) + logp * (L - 1)) + logp * K;
    cout << "(logN,logp,L,K,logQ,h) = ("  << logN << "," << logp << "," << L << "," << K << "," << logQ << ","  << h << ")" << endl;
    cout << "(dim,nslots,nencsnp) = ("  << factorDim << "," << nslots  << "," << nencsnp << ")" << endl;
 
    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
 
    auto start = chrono::steady_clock::now();
    
    Context context(logN, logp, L, K, h);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    ExtScheme extscheme(secretKey, context, scheme);
    
    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);
    
    extscheme.addThreeMultKey(secretKey);
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time= " << timeElapsed << " s" << endl;
    
    CipherPvals cipherPvals(scheme, secretKey, extscheme);
 
    cout << "+------------------------------------+" << endl;
    cout << "|             Encryption             |" << endl;
    cout << "+------------------------------------+" << endl;

    start= chrono::steady_clock::now();
    
    Ciphertext encYXData;   //! E( - y[1] * X[1] - , ..., - y[n] * X[n] -)
    Ciphertext* enccovData = new Ciphertext[34];
    
    Ciphertext** encSData = new Ciphertext*[sampleDim];     //! encSData[i] = E( - S[i] - ), 1 <= i < nsnp
    Ciphertext** encYSData = new Ciphertext*[sampleDim];    //! encYSData[i] = E( - y[i] * S[i] - ), 1 <= i < nsnp
    Ciphertext*** encSXData = new Ciphertext**[sampleDim];  //! encSXData[i][k] = E( - X[i][k] * S[i] - ), 1 <= i < nsnp, 1 <= k < 4
    
    for(long i = 0; i < sampleDim; ++i){
        encYSData[i] = new Ciphertext[nencsnp];
        encSData[i] = new Ciphertext[nencsnp];
        encSXData[i] = new Ciphertext*[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSXData[i][k] = new Ciphertext[nencsnp];
        }
    }
    //! encryption of XData
    cipherPvals.encryptSIMDXData(encYXData, enccovData, yData, xData, factorDim, sampleDim, sampleDim2, nXbatching, nCovbatching, nterms, scaleBits, nslots, L) ;
    
    //! encryption of snp
    cipherPvals.encryptSData(encSData, encYSData, encSXData, yData, xData, sData, factorDim, sampleDim, nsnp, nencsnp, nslots, L);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;

    double totalEvaltime;
    
    // "+------------------------------------+"
    //! 1. Aggregation : (18n ADD)
    // "+------------------------------------+"
   
    start= chrono::steady_clock::now();
    
    //! encSX[factorDim][nencsnp]
    Ciphertext* encYS = new Ciphertext[nencsnp];        //! encYS[j] = Y^T * S = sum_i encYSData[i][j]
    Ciphertext* encS = new Ciphertext[nencsnp];
    Ciphertext** encSX = new Ciphertext*[factorDim];
    
    //! (3n ADD) * 2
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        encYS[j] = encYSData[0][j];
        encS[j] = encSData[0][j];
        for(long i = 1; i < sampleDim; ++i){
            scheme.addAndEqual(encYS[j], encYSData[i][j]);
            scheme.addAndEqual(encS[j], encSData[i][j]);
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! 4*3 ADD
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
 
    //! encYX[k] = YXvec = y^T * X[k] = [-29, -6.506140997, -6.368385672, -12.88024343]
    //! 24 Rot
    Ciphertext* encYX = new Ciphertext[factorDim];
    cipherPvals.aggYXData(encYX, encYXData, sdimBits, nXbatching, factorDim, nslots);
    
    //! (34 * 8) Rot
    Ciphertext* encCov = new Ciphertext[34];
    cipherPvals.aggCovData(encCov, enccovData, sdimBits, nCovbatching);
    
#if  defined(__DEBUG_)
    for(long k = 0; k < nterms; ++k){
        double res;
        cipherPvals.decSingleData(res, encCov[k]);
        cout << res << endl;
    }
    for(long k = 0; k < 4; ++k){
        double res;
        cipherPvals.decSingleData(res, encYX[k]);
        cout << res << endl;
    }
#endif

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "1. Aggregation = " << timeElapsed << " s" << endl;
    totalEvaltime = timeElapsed;
    
    // "+------------------------------------+"
    //! 2. Hessian
    // "+------------------------------------+"
    
    start = chrono::steady_clock::now();
    
    Ciphertext* encAdj = new Ciphertext[nterms];
    Ciphertext encDet;
    
    //cipherPvals.encSIMDAdjoint(encDet, encAdj, encCov);  //! L = 5
    cipherPvals.encSIMDAdjoint_smalldep(encDet, encAdj, encCov);     //! L = 4

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "2. (X^T * X)^-1 = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
#if 1
    double dtemp;
    cout << "1/(s^3) * Adj: [" ;
    for(long l = 0; l < nterms; ++l){
        cipherPvals.decSingleData(dtemp, encAdj[l]);
        cout << dtemp << ", " ;
    }
    cout << "] " << endl ;
    
    cipherPvals.decSingleData(dtemp, encDet);
    cout << "1/(s^3) * det: 4.6759 ?= " << dtemp << endl;
#endif
 
    // "+------------------------------------+"
    //! 3. Norm
    // "+------------------------------------+"
    //! encwtYX = (y^T * X) * adj(X^T * X) * (X^T * y)
    //! Ynorm = |A| <Y,Y> - wtYX = |A| * sampleDim - wtYX

    start= chrono::steady_clock::now();
    
    Ciphertext encWtYX;
    cipherPvals.extSqrQuadForm(encWtYX, encYX, encAdj, factorDim);          //! L - 3, 12.7
    
    Ciphertext encYnorm = scheme.multByConst(encDet, (double) sampleDim);   //! L - 1, (0.58 * n * 32)
    scheme.reScaleByAndEqual(encYnorm, 1);
    
    scheme.modDownToAndEqual(encWtYX, encYnorm.l);
    scheme.subAndEqual(encYnorm, encWtYX);

#if defined(__DEBUG_)
    double dtemp
    cipherPvals.decSingleData(dtemp, encYnorm);
    cout << dtemp << endl;
#endif
    
    Ciphertext* encWtYS = new Ciphertext[nencsnp];
    Ciphertext* encYSnorm = new Ciphertext[nencsnp];
    Ciphertext* encWtSS  = new Ciphertext[nencsnp];
    Ciphertext* encSnorm = new Ciphertext[nencsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        Ciphertext* encSX1 = new Ciphertext[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSX1[k] = encSX[k][j];
            scheme.modDownToAndEqual(encSX1[k], encYX[0].l);    //! L - 1
        }
        //! YSnorm = <ystar, sstar> = encYS - encWtYS
        encYSnorm[j] = scheme.modDownTo(encYS[j], encDet.l);    //! L - 2
        scheme.multAndEqual(encYSnorm[j], encDet);
        scheme.reScaleByAndEqual(encYSnorm[j], 1);              //! L - 3
       
        cipherPvals.extQuadForm(encWtYS[j], encYX, encAdj, encSX1, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtYS[j], encYSnorm[j].l);   //! L - 3
        scheme.subAndEqual(encYSnorm[j], encWtYS[j]);
        
        //! Snorm = <sstar, sstar> = encS - encWtS
        encSnorm[j] = scheme.modDownTo(encS[j], encDet.l);      //! L - 2
        scheme.multAndEqual(encSnorm[j], encDet);
        scheme.reScaleByAndEqual(encSnorm[j], 1);               //! L - 3
        
        cipherPvals.extSqrQuadForm(encWtSS[j], encSX1, encAdj, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtSS[j], encSnorm[j].l);    //! L - 3
        scheme.subAndEqual(encSnorm[j], encWtSS[j]);
    }
    NTL_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3. Norm = " << timeElapsed << " s" << endl;
    cout << "logq: " << encYSnorm[0].l  << "," << encSnorm[0].l << endl;   // L - 3
    totalEvaltime += timeElapsed;
    cout << "Total Evaluation Timing = " << totalEvaltime << " s" << endl;
 
    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    ofstream outf(filename);
    outf.close();
    
    start = chrono::steady_clock::now();

    double Ynorm;
    double** YSnorm = new double*[nencsnp];
    double** Snorm = new double*[nencsnp];
    
    cipherPvals.decryptResult(Ynorm, YSnorm, Snorm, encYnorm, encYSnorm, encSnorm, nencsnp, nslots);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Decryption time = " << timeElapsed << " s" << endl;
    
    //! zscore^2 = (n-k-2) * YSnorm^2 / (Snorm * (Ynorm - YSnorm))
    double deg = (double)(sampleDim - factorDim - 1);
    long nslots1 = nsnp - (nencsnp-1) * nslots;
    
    zScore = new double[nsnp];   //! absolute value of beta/sigma
    pVals = new double[nsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        if(j != nencsnp - 1){
            long starting = j * nslots;
            for(long l = 0; l < nslots; ++l){
                zScore[starting + l] =  sqrt((pow(YSnorm[j][l],2) * deg)/ (Snorm[j][l] * (Ynorm - YSnorm[j][l])));
                pVals[starting + l] =  pnorm(zScore[starting + l]);
            }
        }
        else{
            long starting = (nencsnp - 1) * nslots;
            for(long l = 0; l < nslots1; ++l){
                zScore[starting + l] = sqrt((pow(YSnorm[j][l],2) * deg)/ (Snorm[j][l] * (Ynorm - YSnorm[j][l])));
                pVals[starting + l] =  pnorm(zScore[starting + l]);
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //printvectorToFile(zScore, filename, nsnp);
    printvectorToFile(pVals, filename, nsnp);
}



void TestHEPvals::testHESIMDLinReg_DecompKS(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp,  string filename){
    
    //! Parameters for approx-HE
    long logp = 45;
    long logN = 13;
    long nslots = (1<< (logN-1));      //! total number of plaintext slots
    long L = 4;                        //! 2 (det,adj) + 1 (multiplied by det,adj) + 1
    long K = 1;
    long h = 131;
    
    //! Parameters for GWAS
    long scaleBits = 4;            //! scale factor for covariance
    long sampleDim2 = (1 << (long)ceil(log2(sampleDim)));   //! closet PoT
    long sdimBits = (long)ceil(log2(sampleDim));            //! log2(sampleDim)
    long factorDim2 = factorDim * factorDim;
    
    long nXbatching = nslots/(sampleDim2 * factorDim);     //! replicated number of a user's Xdata in a single ciphertext
    long nCovbatching = nslots/(sampleDim2 * 8);           //! replicated number of a user's (X^T * X) in a single ciphertext
    long nencsnp = (long)ceil((double)nsnp/nslots);          //! number of ciphertexts for snp encryptions
    long nterms = factorDim * (factorDim + 1)/2;
    
    long logQ = (Q0_BIT_SIZE) + logp * (L - 1) + (P0_BIT_SIZE);
    cout << "(logN,logp,L,K,logQ) = ("  << logN << "," << logp << "," << L << "," << K << "," << logQ << ")" << endl;
    cout << "(dim,nslots,nencsnp) = ("  << factorDim << "," << nslots  << "," << nencsnp  << ")" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    auto start = chrono::steady_clock::now();
    
    Context context(logN, logp, L, K, h);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    ExtScheme extscheme(secretKey, context, scheme);
    
    //! evk for multiplication
    extscheme.addDecompTwoKeys(secretKey);
    extscheme.addDecompThreeKeys(secretKey);
    
    //! evk for rotation
    extscheme.addDecompLeftRotKeys(secretKey);
    extscheme.addDecompRightRotKeys(secretKey);
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time= " << timeElapsed << " s" << endl;
    
    CipherPvals cipherPvals(scheme, secretKey, extscheme);
    
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Encryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    start= chrono::steady_clock::now();
    
    Ciphertext encYXData;   //! E( - y[1] * X[1] - , ..., - y[n] * X[n] -)
    Ciphertext* enccovData = new Ciphertext[34];
    
    Ciphertext** encSData = new Ciphertext*[sampleDim];     //! encSData[i] = E( - S[i] - ), 1 <= i < nsnp
    Ciphertext** encYSData = new Ciphertext*[sampleDim];    //! encYSData[i] = E( - y[i] * S[i] - ), 1 <= i < nsnp
    Ciphertext*** encSXData = new Ciphertext**[sampleDim];  //! encSXData[i][k] = E( - X[i][k] * S[i] - ), 1 <= i < nsnp, 1 <= k < 4
    
    for(long i = 0; i < sampleDim; ++i){
        encYSData[i] = new Ciphertext[nencsnp];
        encSData[i] = new Ciphertext[nencsnp];
        encSXData[i] = new Ciphertext*[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSXData[i][k] = new Ciphertext[nencsnp];
        }
    }
    
    cipherPvals.encryptSIMDXData(encYXData, enccovData, yData, xData, factorDim, sampleDim, sampleDim2, nXbatching, nCovbatching, nterms, scaleBits, nslots, L) ;
    
    cipherPvals.encryptSData(encSData, encYSData, encSXData, yData, xData, sData, factorDim, sampleDim, nsnp, nencsnp,  nslots, L);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    double totalEvaltime;
    
    // "+------------------------------------+"
    //! 1. Aggregation : (18n ADD)
    // "+------------------------------------+"
    
    start= chrono::steady_clock::now();
    
    //! encSX[factorDim][nencsnp]
    Ciphertext* encYS = new Ciphertext[nencsnp];  //! encYS[j] = Y^T * S = sum_i encYSData[i][j]
    Ciphertext* encS = new Ciphertext[nencsnp];
    Ciphertext** encSX = new Ciphertext*[factorDim];
    
    //! (3n ADD) * 2
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        encYS[j] = encYSData[0][j];
        encS[j] = encSData[0][j];
        for(long i = 1; i < sampleDim; ++i){
            scheme.addAndEqual(encYS[j], encYSData[i][j]);
            scheme.addAndEqual(encS[j], encSData[i][j]);
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! 4*3 ADD
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
    
    //! encYX[k] = YXvec = y^T * X[k] = [-29, -6.506140997, -6.368385672, -12.88024343]
    //! 24 Rot
    Ciphertext* encYX = new Ciphertext[factorDim];
    cipherPvals.aggYXData_DecompKS(encYX, encYXData, sdimBits, nXbatching, factorDim, nslots);
    
    //! (34 * 8) Rot
    Ciphertext* encCov = new Ciphertext[34];
    cipherPvals.aggCovData_DecompKS(encCov, enccovData, sdimBits, nCovbatching);
    
#if defined(__DEBUG_)
    for(long k = 0; k < nterms; ++k){
        double res;
        cipherPvals.decSingleData(res, encCov[k]);
        cout << res << endl;
    }
    for(long k = 0; k < 4; ++k){
        double res;
        cipherPvals.decSingleData(res, encYX[k]);
        cout << res << endl;
    }
#endif
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "1. Aggregation = " << timeElapsed << " s" << endl;
    totalEvaltime = timeElapsed;
 
    // "+------------------------------------+"
    //! 2. Hessian
    // "+------------------------------------+"
   
    start = chrono::steady_clock::now();
    
    Ciphertext* encAdj = new Ciphertext[10];
    Ciphertext encDet;
    
    cipherPvals.encSIMDAdjoint_DecompKS(encDet, encAdj, encCov);     //! L = 4
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "2. (X^T * X)^-1 = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
#if 1
    double dtemp;
    cout << "1/(s^3) * Adj: [" ;
    for(long l = 0; l < nterms; ++l){
        cipherPvals.decSingleData(dtemp, encAdj[l]);
        cout << dtemp << ", " ;
    }
    cout << "] " << endl ;
    
    cipherPvals.decSingleData(dtemp, encDet);
    cout << "1/(s^3) * det: 4.6759 ?= " << dtemp << endl;
#endif
    
    // "+------------------------------------+"
    //! 3. Norm
    // "+------------------------------------+"
    //! encwtYX = (y^T * X) * adj(X^T * X) * (X^T * y)
    //! Ynorm = |A| <Y,Y> - wtYX = |A| * sampleDim - wtYX
    
    start= chrono::steady_clock::now();
    
    Ciphertext encWtYX;
    cipherPvals.extSqrQuadForm_DecompKS(encWtYX, encYX, encAdj, factorDim);     //! L - 3
    
    Ciphertext encYnorm = scheme.multByConst(encDet, (double) sampleDim);       //! L - 1
    scheme.reScaleByAndEqual(encYnorm, 1);
    
    scheme.modDownToAndEqual(encWtYX, encYnorm.l);
    scheme.subAndEqual(encYnorm, encWtYX);
    
#if defined(__DEBUG_)
    double dtemp
    cipherPvals.decSingleData(dtemp, encYnorm);
    cout << dtemp << endl;
#endif
    
    //! YSnorm = <ystar, sstar> = encYS - encWtYS
    Ciphertext* encWtYS = new Ciphertext[nencsnp];
    Ciphertext* encYSnorm = new Ciphertext[nencsnp];
    Ciphertext* encWtSS  = new Ciphertext[nencsnp];
    Ciphertext* encSnorm = new Ciphertext[nencsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        Ciphertext* encSX1 = new Ciphertext[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSX1[k] = encSX[k][j];
            scheme.modDownToAndEqual(encSX1[k], encYX[0].l);      // L - 1
        }
        
        //! YSnorm = <ystar, sstar> = encYS - encWtYS
        encYSnorm[j] = scheme.modDownTo(encYS[j], encDet.l);    //! L - 2
        encYSnorm[j] = extscheme.mult(encYSnorm[j], encDet);
        scheme.reScaleByAndEqual(encYSnorm[j], 1);              //! L - 3
        
        cipherPvals.extQuadForm_DecompKS(encWtYS[j], encYX, encAdj, encSX1, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtYS[j], encYSnorm[j].l);
        scheme.subAndEqual(encYSnorm[j], encWtYS[j]);
        
        //! Snorm = <sstar, sstar> = encS - encWtS
        encSnorm[j] = scheme.modDownTo(encS[j], encDet.l);      //! L - 2
        encSnorm[j] = extscheme.mult(encSnorm[j], encDet);
        scheme.reScaleByAndEqual(encSnorm[j], 1);               //! L - 3
        
        cipherPvals.extSqrQuadForm_DecompKS(encWtSS[j], encSX1, encAdj, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtSS[j], encSnorm[j].l);    //! L - 3
        scheme.subAndEqual(encSnorm[j], encWtSS[j]);
    }
    NTL_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3. Norm = " << timeElapsed << " s" << endl;
    cout << "logq: " << encYSnorm[0].l  << "," << encSnorm[0].l << endl;   // L - 3
    totalEvaltime += timeElapsed;
    cout << "Total Evaluation Timing = " << totalEvaltime << " s" << endl;
    
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    ofstream outf(filename);
    outf.close();
    
    start = chrono::steady_clock::now();
    
    double Ynorm;
    double** YSnorm = new double*[nencsnp];
    double** Snorm = new double*[nencsnp];
    
    cipherPvals.decryptResult(Ynorm, YSnorm, Snorm, encYnorm, encYSnorm, encSnorm, nencsnp, nslots);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Decryption time = " << timeElapsed << " s" << endl;
    
    //! zscore^2 = (n-k-2) * YSnorm^2 / (Snorm * (Ynorm - YSnorm))
    double deg = (double)(sampleDim - factorDim - 1);
    long nslots1 = nsnp - (nencsnp-1) * nslots;
    
    zScore = new double[nsnp];   //! absolute value of beta/sigma
    pVals = new double[nsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        if(j != nencsnp - 1){
            long starting = j * nslots;
            for(long l = 0; l < nslots; ++l){
                zScore[starting + l] =  sqrt((pow(YSnorm[j][l],2) * deg)/ (Snorm[j][l] * (Ynorm - YSnorm[j][l])));
                pVals[starting + l] =  pnorm(zScore[starting + l]);
            }
        }
        else{
            long starting = (nencsnp - 1) * nslots;
            for(long l = 0; l < nslots1; ++l){
                zScore[starting + l] = sqrt((pow(YSnorm[j][l],2) * deg)/ (Snorm[j][l] * (Ynorm - YSnorm[j][l])));
                pVals[starting + l] =  pnorm(zScore[starting + l]);
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //printvectorToFile(zScore, filename, nsnp);
    printvectorToFile(pVals, filename, nsnp);
 
}


void TestHEPvals::testHELinReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp,  string filename){
    
    //! Parameters for approx-HE
    long logp = 40;
    long logN = 13;
    long nslots = (1<< (logN-1)); //! total number of plaintext slots
    long L = 5;
    long K = L;
    
    //! Parameters for GWAS
    long factorDim2 = factorDim * factorDim;
    long scalefactor = 32.0;  //! scale factor for covariance
    long nencsnp = (long)ceil((double)nsnp/nslots);  // number of ciphertext for snp encryptions
    long nterms = factorDim * (factorDim + 1)/2;
    
    long logq = (Q0_BIT_SIZE) + logp * (L - 1);
    cout << "(logN,logp,L,K, logq) = ("  << logN << "," << logp << "," << L << "," << K << "," << logq << ")" << endl;
    cout << "(dim,nslots,nencsnp, nterms) = ("  << factorDim << "," << nslots  << "," << nencsnp<< "," << nterms  << ")" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    auto start = chrono::steady_clock::now();
    
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);   // binary
    Scheme scheme(secretKey, context);
    ExtScheme extscheme(secretKey, context, scheme);
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time= " << timeElapsed << " s" << endl;
    
    CipherPvals cipherPvals(scheme, secretKey, extscheme);
    
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Encryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    //! encXData[i][k]: encryption of y[i] * xDta[i][k] (scaled) = Enc(xDta[i][k].....xDta[i][k])
    //! encYSData[i]: encryption of y[i] * sDta[i]   = Enc(y[i] * sDta[i][0] , ... , y[i]* sDta[i][p-1])
    //! encSXData[i][k] = Enc(x[i][k] * sDta[i][0] , ... , (x[i][k] * sDta[i][p-1])
    //! encS2Data[i] = encSDatat[i] = Enc(sDta[i][0] , ... , sDta[i][p-1])
    
    start= chrono::steady_clock::now();
    
    Ciphertext** enccovData = new Ciphertext*[sampleDim];
    Ciphertext** encYXData = new Ciphertext*[sampleDim];
    
    Ciphertext** encSData = new Ciphertext*[sampleDim];
    Ciphertext** encYSData = new Ciphertext*[sampleDim];
    Ciphertext*** encSXData = new Ciphertext**[sampleDim];
    
    for(long i = 0; i < sampleDim; ++i){
        enccovData[i] = new Ciphertext[nterms];
        encYXData[i]  = new Ciphertext[factorDim];
        encYSData[i] = new Ciphertext[nencsnp];
        encSData[i] = new Ciphertext[nencsnp];
        encSXData[i] = new Ciphertext*[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSXData[i][k] = new Ciphertext[nencsnp];
        }
    }
    
    cipherPvals.encryptXData(encYXData, enccovData, yData, xData, factorDim, sampleDim, nterms, scalefactor, nslots, L);
    cipherPvals.encryptSData(encSData, encYSData, encSXData, yData, xData, sData, factorDim, sampleDim, nsnp, nencsnp,  nslots, L);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    double totalEvaltime;
    
    // "+------------------------------------+"
    //! 1. Aggregation
    // "+------------------------------------+"
    //! encYS[j] = Y^T * S = sum_i encYSData[i][j]
    //! encS[j]  = sum_i SData[i][j]^2
    //! (10 * n + 3 * 2 * n + 4 * 2 * n ) ADD
    
    start= chrono::steady_clock::now();
    
    Ciphertext* encCov = new Ciphertext[nterms];
    Ciphertext* encYS = new Ciphertext[nencsnp];
    Ciphertext* encS = new Ciphertext[nencsnp];
    
    //! 10n ADD
    NTL_EXEC_RANGE(nterms, first, last);
    for(long k = first; k < last; ++k){
        encCov[k] = enccovData[0][k];
        for(long i = 1; i < sampleDim; ++i){
            scheme.addAndEqual(encCov[k], enccovData[i][k]);
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! (3n ADD) * 2
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        encYS[j] = encYSData[0][j];
        encS[j] = encSData[0][j];
        for(long i = 1; i < sampleDim; ++i){
            scheme.addAndEqual(encYS[j], encYSData[i][j]);
            scheme.addAndEqual(encS[j], encSData[i][j]);
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! encYX[k] = YXvec = y^T * X[k] = [-29, -6.506140997, -6.368385672, -12.88024343]
    //! encSX[factorDim][nencsnp]
    
    Ciphertext* encYX = new Ciphertext[factorDim];
    Ciphertext** encSX = new Ciphertext*[factorDim];
    
    //! (4n + 4*3*n) ADD
    NTL_EXEC_RANGE(factorDim, first, last);
    for(long k = first; k < last; ++k){
        encYX[k] = encYXData[0][k];
        for(long i = 1; i < sampleDim; ++i){
            scheme.addAndEqual(encYX[k], encYXData[i][k]);
        }
        encSX[k] = new Ciphertext[nencsnp];
        for(long j = 0; j < nencsnp; ++j){
            encSX[k][j] = encSXData[0][k][j];
            for(long i = 1; i < sampleDim; ++i){
                scheme.addAndEqual(encSX[k][j], encSXData[i][k][j]);
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
#if  defined(__DEBUG_)
    for(long k = 0; k < nterms; ++k){
        double res;
        cipherPvals.decSingleData(res, encCov[k]);
        cout << res << endl;
    }
    
    for(long k = 0; k < 4; ++k){
        double res;
        cipherPvals.decSingleData(res, encYX[k]);
        cout << res << endl;
    }
#endif
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "1. Aggregation = " << timeElapsed << " s" << endl;
    totalEvaltime = timeElapsed;
    
    // "+------------------------------------+"
    //! 2. Hessian
    // "+------------------------------------+"
    
    start = chrono::steady_clock::now();
    Ciphertext encDet;
    Ciphertext* encAdj = new Ciphertext[nterms];
    //cipherPvals.HesInverse(encDet, encAdj, enccovData, dim, nslots, L);   //! trivial encryption
    //cipherPvals.encAdjoint(encDet, encAdj, encCov);  // usual mult
    cipherPvals.extencAdjoint(encDet, encAdj, encCov); // extended
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "2. (X^T * X)^-1 = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    double dtemp;
    cout << "Adj: [" ;
    for(long l = 0; l < 10; ++l){
        cipherPvals.decSingleData(dtemp, encAdj[l]);
        cout << dtemp << ", " ;
    }
    cout << "] " << endl ;
    
    cipherPvals.decSingleData(dtemp, encDet);
    cout << "det: 0.0182656 ?= " << dtemp << endl;
    
    
    // "+------------------------------------+"
    //! 3. Norm
    // "+------------------------------------+"
    //! encwtYX = (y^T * X) * adj(X^T * X) * (X^T * y)
    //! Ynorm = |A| <Y,Y> - wtYX = |A| * sampleDim - wtYX
    
    
    start= chrono::steady_clock::now();
    
    Ciphertext encWtYX;
    cipherPvals.extSqrQuadForm(encWtYX, encYX, encAdj, factorDim);  // L - 3, 12.7
    
    Ciphertext encYnorm = scheme.multByConst(encDet, (double) sampleDim * scalefactor);  // L - 4, (0.58 * n * 32)
    scheme.reScaleByAndEqual(encYnorm, 1);
    
    scheme.modDownToAndEqual(encWtYX, encYnorm.l);
    scheme.subAndEqual(encYnorm, encWtYX); // L - 4
    
#if defined(__DEBUG_)
    double res;
    cipherPvals.decSingleData(res, encDet);
    cout << res << endl;
    
    cipherPvals.decSingleData(res, encYnorm);
    cout << res << endl;
#endif
    
    //! YSnorm = <ystar, sstar> = encYS - encWtYS
    Ciphertext* encWtYS = new Ciphertext[nencsnp];
    Ciphertext* encYSnorm = new Ciphertext[nencsnp];
    Ciphertext* encWtSS  = new Ciphertext[nencsnp];
    Ciphertext* encSnorm = new Ciphertext[nencsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        Ciphertext* encSX1 = new Ciphertext[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSX1[k] = encSX[k][j];
        }
        
        //! encYSnorm = encDet * (encYS[j] * scalefactor) - encWtYS
        encYSnorm[j] = scheme.multByConst(encYS[j], (double) scalefactor);
        scheme.reScaleByAndEqual(encYSnorm[j], 1);              //! L - 1
        
        scheme.modDownToAndEqual(encYSnorm[j], encDet.l);       //! L - 3
        scheme.multAndEqual(encYSnorm[j], encDet);
        scheme.reScaleByAndEqual(encYSnorm[j], 1);              //! L - 4
        
        cipherPvals.extQuadForm(encWtYS[j], encYX, encAdj, encSX1, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtYS[j], encYSnorm[j].l);   //! L - 4
        scheme.subAndEqual(encYSnorm[j], encWtYS[j]);
        
        //! encSnorm
        encSnorm[j] = scheme.multByConst(encS[j], (double) scalefactor);
        scheme.reScaleByAndEqual(encSnorm[j], 1);               //! L - 1
        
        scheme.modDownToAndEqual(encSnorm[j], encDet.l);        //! L - 2
        scheme.multAndEqual(encSnorm[j], encDet);
        scheme.reScaleByAndEqual(encSnorm[j], 1);               //! L - 3
        
        cipherPvals.extSqrQuadForm(encWtSS[j], encSX1, encAdj, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtSS[j], encSnorm[j].l);    //! L - 4
        scheme.subAndEqual(encSnorm[j], encWtSS[j]);
    }
    NTL_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3. Norm = " << timeElapsed << " s" << endl;
    cout << "logq: " << encYSnorm[0].l  << "," << encSnorm[0].l << endl;   // L - 3
    totalEvaltime += timeElapsed;
    cout << "Total Evaluation Timing = " << totalEvaltime << " s" << endl;
    
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    ofstream outf(filename);
    outf.close();
    
    start = chrono::steady_clock::now();
    
    double Ynorm;
    cipherPvals.decSingleData(Ynorm, encYnorm);
    
    double deg = (double)(sampleDim - factorDim - 1);
    long nslots1 = nsnp - (nencsnp-1) * nslots;
    
    // zscore^2 = (n-k-2) * YSnorm^2 / (Snorm * (Ynorm - YSnorm))
    zScore = new double[nsnp];   //! absolute value of beta/sigma
    pVals = new double[nsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        double* YSnorm;
        double* Snorm;
        
        cipherPvals.decVector(YSnorm, encYSnorm[j], nslots);
        cipherPvals.decVector(Snorm, encSnorm[j], nslots);
        
        if(j != nencsnp - 1){
            long starting = j * nslots;
            for(long l = 0; l < nslots; ++l){
                zScore[starting + l] =  sqrt((pow(YSnorm[l],2) * deg)/ (Snorm[l] * (Ynorm - YSnorm[l])));
                pVals[starting + l] =  pnorm(zScore[starting + l]);
            }
        }
        else{
            long starting = (nencsnp - 1) * nslots;
            for(long l = 0; l < nslots1; ++l){
                zScore[starting + l] = sqrt((pow(YSnorm[l],2) * deg)/ (Snorm[l] * (Ynorm - YSnorm[l])));
                pVals[starting + l] =  pnorm(zScore[starting + l]);
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Decryption time = " << timeElapsed << " s" << endl;
    
    //printvectorToFile(zScore, filename, nsnp);
    printvectorToFile(pVals, filename, nsnp);
    
}
