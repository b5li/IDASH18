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
#include <sys/resource.h>   //! check the memory usage

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
//!@ Function: using decomposition KS -> logN = 13

void TestHEPvals::testFastHELinReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, vector<string> snptag, string filename){
    
    struct rusage usage;
    long memoryscale = (1 << 20);
    
    //! Parameters for approx-HE
    long logp = 43;                    //! all the msg are scaled by "p", logq0 - logp = (final bits of precision)
    long logq0 = 55;
    long logp0 = 61;
    
    long logN = 13;
    long nslots = (1<< (logN-1));      //! total number of plaintext slots
    long L = 4;                        //! 2 (det,adj) + 1 (multiplied by det,adj) + 1 (decoding)
    long K = 1;
    long h = 130;                      //! Hamming weight of sk
    
    //! encryption level
    long Slvl = L - 2;  //! = level of det
    long YSlvl = L - 2; //! = level of det
    long SXlvl = L - 1; //! =  level of (sum (YX[i]))
    long YXlvl = L ;
    long Covlvl = L ;   //! fresh
    
    //! Parameters for GWAS
    long scaleBits = 4;                //! scale factor for covariance
    long sampleDim2 = (1 << (long)ceil(log2(sampleDim)));   //! closet PoT
    long sdimBits = (long)ceil(log2(sampleDim));            //! log2(sampleDim)
    long factorDim2 = factorDim * factorDim;
    
    long nXbatching = nslots/(sampleDim2 * factorDim);     //! replicated number of a user's Xdata in a single ciphertext
    long nCovbatching = nslots/(sampleDim2 * 8);           //! replicated number of a user's (X^T * X) in a single ciphertext
    long nencsnp = (long)ceil((double)nsnp/nslots);         //! number of ciphertexts for snp encryptions
    long nterms = factorDim * (factorDim + 1)/2;
    
    long logQ = logq0 + logp * (L - 1) + logp0;
    cout << "(logN,logp,L,K,logQ,h) = ("  << logN << "," << logp << "," << L << "," << K << "," << logQ << "," << h << ")" << endl;
    cout << "(dim,nslots,nencsnp) = ("  << factorDim << "," << nslots  << "," << nencsnp  << ")" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    auto start = chrono::steady_clock::now();
    
    Context context(logN, logp, L, K, logq0, logp0, h);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    ExtScheme extscheme(secretKey, context, scheme);
    
    //! evk for multiplication
    extscheme.addDecompTwoKeys(secretKey);
    extscheme.addDecompThreeKeys(secretKey);
    
    //! evk for rotation
    extscheme.addDecompLeftRotKeys(secretKey);
    //extscheme.addDecompRightRotKeys(secretKey);
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time= " << timeElapsed << " s" << endl;
    int ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
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
    
    cipherPvals.encryptSIMDXData(encYXData, enccovData, yData, xData, factorDim, sampleDim, sampleDim2, nXbatching, nCovbatching, nterms, scaleBits, nslots, YXlvl, Covlvl) ;
    
    cipherPvals.encryptSData(encSData, encYSData, encSXData, yData, xData, sData, factorDim, sampleDim, nsnp, nencsnp, nslots, Slvl, YSlvl, SXlvl);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    //! MB
    double sizeYX = 1 * 2 * (logq0 + (YXlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizecov = 34 * 2 * (logq0 + (Covlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeS = sampleDim * nencsnp * 2 * (logq0 + (Slvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeYS = sampleDim * nencsnp * 2 * (logq0 + (YSlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeSX = sampleDim * nencsnp * factorDim * 2 * (logq0 + (Slvl - 1) * logp) * (1 << logN) / (1 << 23);
    cout << "Freshly Ciphertexts size: " << (sizeYX + sizecov + sizeS + sizeYS + sizeSX) << "(MB)" << endl;
    
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
    Ciphertext* encYX = new Ciphertext[factorDim];
    cipherPvals.aggYXData(encYX, encYXData, sdimBits, nXbatching, factorDim, nslots);
    
    
    // "+------------------------------------+"
    //! 2. Hessian
    // "+------------------------------------+"
    

    Ciphertext* encAdj = new Ciphertext[10];
    Ciphertext encDet;
    
    cipherPvals.encAdjoint(encDet, encAdj, enccovData, sdimBits, nCovbatching);

    
#if defined(__DEBUG_)
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
    

    Ciphertext encWtYX;
    cipherPvals.extSqrQuadForm(encWtYX, encYX, encAdj, factorDim);              //! L - 3 = 1
    
    Ciphertext encYnorm = extscheme.multByConstMT(encDet, (long) sampleDim);     //! L - 2 = 2
    
    scheme.modDownToAndEqual(encYnorm, encWtYX.l);
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
            scheme.modDownToAndEqual(encSX1[k], encYX[0].l);    //! L - 1
        }
        //! YSnorm = <ystar, sstar> = encYS - encWtYS
        //encYSnorm[j] = scheme.modDownTo(encYS[j], encDet.l);    //! L - 2
        encYSnorm[j] = extscheme.mult(encYS[j], encDet);
        scheme.reScaleByAndEqual(encYSnorm[j], 1);              //! L - 3
        
        cipherPvals.extQuadForm(encWtYS[j], encYX, encAdj, encSX1, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtYS[j], encYSnorm[j].l);
        scheme.subAndEqual(encYSnorm[j], encWtYS[j]);
        
        //! Snorm = <sstar, sstar> = encS - encWtS
        //encSnorm[j] = scheme.modDownTo(encS[j], encDet.l);      //! L - 2
        encSnorm[j] = extscheme.mult(encS[j], encDet);
        scheme.reScaleByAndEqual(encSnorm[j], 1);               //! L - 3
        
        cipherPvals.extSqrQuadForm(encWtSS[j], encSX1, encAdj, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtSS[j], encSnorm[j].l);    //! L - 3
        scheme.subAndEqual(encSnorm[j], encWtSS[j]);
    }
    NTL_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    //cout << "3. Norm = " << timeElapsed << " s" << endl;
    //cout << "logq: " << encYSnorm[0].l  << "," << encSnorm[0].l << endl;   // L - 3
    totalEvaltime = timeElapsed;
    cout << "Total Evaluation Timing = " << totalEvaltime << " s" << endl;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
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
    
    double* YSnorm1 = new double[nsnp];
    double* Snorm1  = new double[nsnp];
   
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        if(j != nencsnp - 1){
            long starting = j * nslots;
            for(long l = 0; l < nslots; ++l){
                YSnorm1[starting + l] = YSnorm[j][l];
                Snorm1[starting + l] = Snorm[j][l];
                zScore[starting + l] =  sqrt((pow(YSnorm[j][l],2) * deg)/ (Snorm[j][l] * (Ynorm - YSnorm[j][l])));
                pVals[starting + l] =  pnorm(zScore[starting + l]);
            }
        }
        else{
            long starting = (nencsnp - 1) * nslots;
            for(long l = 0; l < nslots1; ++l){
                YSnorm1[starting + l] = YSnorm[j][l];
                Snorm1[starting + l] = Snorm[j][l];
                zScore[starting + l] = sqrt((pow(YSnorm[j][l],2) * deg)/ (Snorm[j][l] * (Ynorm - YSnorm[j][l])));
                pVals[starting + l] =  pnorm(zScore[starting + l]);
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    printvectorToFile(YSnorm1, "LinRegResult/FastHELinReg_YSnorm.txt", nsnp);
    printvectorToFile(Snorm1, "LinRegResult/FastHELinReg_Snorm.txt", nsnp);
    //printvectorToFile(pVals, filename, nsnp);
    printPvalsToFile(pVals, snptag, filename, nsnp);
    
    delete[] encYS;
    delete[] encS;
    delete[] encSX;
    delete[] encYX;
    delete[] encAdj;
    delete[] encWtYS;
    delete[] encYSnorm;
    delete[] encWtSS;
    delete[] encSnorm;
    
    delete[] YSnorm;
    delete[] Snorm;
    delete[] YSnorm1;
    delete[] Snorm1;
}

