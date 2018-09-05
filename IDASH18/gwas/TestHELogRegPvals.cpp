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
#include "TestLogRegPvals.h"
#include "CipherLinRegPvals.h"
#include "TestHELinRegPvals.h"
#include "CipherLogRegPvals.h"
#include "TestHELogRegPvals.h"

#define ver1 1
#define ver2 0

void TestHELRPvals::testHELogReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename){
    
    struct rusage usage;
    long memoryscale = (1 << 20);   //! 2^20: linux, 2^30: mac
    
    //! Parameters for GWAS
    long YXscaleBits = 2;
    long covscaleBits = 2;      //! scale factor for covariance
    long subblocksize = 16;     //! size of subblocks for fully replicating size-n vector
    
    //! HE Parameters
    long logp = 43;                    //! all the msg are scaled by "p", logq0 - logp = (final bits of precision), logp: scaled bit
    long logq0 = 53;
    long logp0 = 61;
    
    long logN = 15;
    long nslots = (1 << (logN-1));     //! total number of plaintext slots
    long L = 23;                       //! 9 (LogReg) + 4 (Pr) + 1 (W) + 5 (Z, two-inverse) + 2 (Z^T * X) + 1 (multiplied by adj * XW2S)= 22
    long K = 1;
    long h = 170;                       //! Hamming weight of sk (logq = 1060)
    
    //! encryption level for snp data
    long Slvl = 3;
    long SXlvl = 3;    //! encW2.lvl,  fast: L - 19;

    
    //! encryption level for Xdata
    long YXlvl = L;
    long Xlvl = L - 9;      //! encBeta.lvl = 1 + (2 + log2(deg(sigmoid))) * (iter)
    long Ylvl = L - 13;     //! encPr.lvl
    long Covlvl = 6;

    //! Parameters
    long sampleDim2 = (1 << (long)ceil(log2(sampleDim)));   //! closet PoT
    long sdimBits = (long)ceil(log2(sampleDim));            //! log2(sampleDim)
    long fdimBits = (long)ceil(log2(factorDim));
    long factorDim2 = factorDim * factorDim;
    
    long nXbatching = nslots/(sampleDim2 * factorDim);     //! replicated number of a user's Xdata in a ciphertext
    long nCovbatching = nslots/(sampleDim2 * 8);           //! replicated number of a user's (X^T * X) in a ciphertext
    long bBits = (long)ceil(log2(nXbatching));             //! batching bits for XData

    long xBatchingBits = fdimBits + bBits;                 //! log2(number of slots for a single user in a ciphertext)
    long nslots1 =  (1 << xBatchingBits) * sampleDim;      //! number of packed messages in the slots

    long nencsnp = (long)ceil((double)nsnp/nslots);         //! number of ciphertexts for snp encryptions
    long nterms = factorDim * (factorDim + 1)/2;    

    //! Parameters for Logistic Regression of covariates
    long sigdeg = 3;
    long numIter = 3;
    double gammaUp = 1;
    double gammaDown = 1;
    
    long logQ = logq0 + logp * (L - 1) + logp0;
    cout << "(logN,logp,L,K,logQ,h) = ("  << logN << "," << logp << "," << L << "," << K << "," << logQ << "," << h << ")" << endl;
    cout << "(dim,nslots,nbatching,nencsnp,subblock) = ("  ;
    cout << factorDim << "," << nslots <<  "," << nXbatching << "," << nencsnp << "," << subblocksize << ")" << endl;
 
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

    //! evk for rotation (generate all the left rotation keys/ right keys for 1 and 2)
    extscheme.addDecompLeftRotKeys(secretKey);
    for(long i = 0; i < fdimBits; ++i){
        extscheme.addDecompRightRotKey(secretKey, 1 << i); //! 1, 2
    }
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time= " << timeElapsed << " s" << endl;
    int ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    

    cout << "+------------------------------------+" << endl;
    cout << "|             Encryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    CipherPvals cipherPvals(scheme, secretKey, extscheme);
    CipherLRPvals cipherLRPvals(scheme, secretKey, extscheme, cipherPvals);
    
    start = chrono::steady_clock::now();
    
    Ciphertext encYXData;   //! E( - y[1] * X[1] - , ..., - y[n] * X[n] - )
    Ciphertext encXData;
    Ciphertext encYData;
    Ciphertext* enccovData = new Ciphertext[34];

    cipherLRPvals.encXData(encYXData, encXData, encYData, enccovData, yData, xData, factorDim, sampleDim, nXbatching, nCovbatching, nterms, YXscaleBits, covscaleBits, nslots, YXlvl, Xlvl, Ylvl, Covlvl);
    
    Ciphertext** encSData = new Ciphertext*[sampleDim];     //! encSData[n][encsnp]
    Ciphertext*** encSXData = new Ciphertext**[sampleDim];  //! encSXData[n][encsnp][4]
    
    for(long i = 0; i < sampleDim; ++i){
        encSData[i] = new Ciphertext[nencsnp];
        encSXData[i] = new Ciphertext*[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSXData[i][k] = new Ciphertext[nencsnp];
        }
    }
    
    cipherLRPvals.encryptSData(encSData, encSXData, xData, sData, factorDim, sampleDim, nsnp, nencsnp, nslots, Slvl, SXlvl);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    double totalEvaltime;
    
    // "+------------------------------------+"
    //! 1. Logistic regression on covariates, 1 + (iter - 1) * 4 = 9
    // "+------------------------------------+"
    start = chrono::steady_clock::now();
    
    Ciphertext encBeta;
    uint64_t* poly = cipherLRPvals.generateNLGDAuxPoly(nslots, nslots1, factorDim);
    cipherLRPvals.encNLGDiteration(encBeta, encYXData, factorDim, sampleDim, fdimBits, sdimBits, xBatchingBits, poly, numIter);

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "1. Logistic regression = " << timeElapsed << " s" << endl;
    totalEvaltime = timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    // "+------------------------------------+"
    //! 2. Update the weights
    // "+------------------------------------+"
    start = chrono::steady_clock::now();
    
    Ciphertext encZData;  //! (L - 9) - 9 = L - 18
    Ciphertext encWData;  //! (L - 9) - 5  = L - 14
    Ciphertext encW2Data; //! (L - 9) - 6 = L - 15
    Ciphertext encZWData; //! (L - 9) - 6 = L - 15
    
    long steps = 3; //! degree for approximation of 1/p*(1-p)
    cipherLRPvals.encZWData(encZData, encWData, encW2Data, encZWData, encBeta, encXData, encYData, poly, fdimBits, steps);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "2. Update the weights = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    delete[] poly;
 #if 1
    // "+------------------------------------+"
    //! 3.1. Hessian: encAdj[4], encDet
    // "+------------------------------------+"
    
    start = chrono::steady_clock::now();
    
    Ciphertext* encAdj = new Ciphertext[10];    //! (L - 21) or (L - 20) = 3
    Ciphertext encDet;
    
    Ciphertext encWData1 = scheme.modDownBy(encWData, 3);
    cipherLRPvals.encAdjoint(encDet, encAdj, encWData1, enccovData, sdimBits, nCovbatching); //! (X^T * W * X)
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.1. (X^T * W *  X)^-1 = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    delete[] enccovData;
    
#if defined(__DEBUG_)
    double dtemp;
    cout << "1/(s^3) * Adj: [" ;
    for(long l = 0; l < nterms; ++l){
        cipherPvals.decSingleData(dtemp, encAdj[l]);
        cout << dtemp << ", " ;
    }
    cout << "] " << endl ;
    
    cipherPvals.decSingleData(dtemp, encDet);
    cout << "1/(s^3) * det: 1.16008 ?= " << dtemp << endl;
#endif

    // "+------------------------------------+"
    //! 3.2. encZXData[4], L - 21,
    // "+------------------------------------+"
    
    start= chrono::steady_clock::now();

    Ciphertext* encZX = new Ciphertext[factorDim];
    cipherLRPvals.encZXData(encZX, encXData, encZData, sdimBits, nXbatching, factorDim, nslots); // [-43.8965903, -7.693702862, -7.157768273, -18.84315905]
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.2. Z^T * X   = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;

    // "+------------------------------------+"
    //! 3.3. encSX[4][nencsnp] = enc(S^T * X),
    // "+------------------------------------+"
    
    start = chrono::steady_clock::now();
    
    Ciphertext** encSX = new Ciphertext*[factorDim];
    cipherLRPvals.encSXData(encSX, encSXData, factorDim, sampleDim, nencsnp);   //! (4*n ADD)
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.3. S^T * X  = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    // "+------------------------------------+"
    //! 3.4.0 encZWSData[nencsnp]
    // "+------------------------------------+"
    
    start = chrono::steady_clock::now();
    
    long niter = (long) ceil((double)sampleDim/subblocksize);    //! 245/16 = 16, 245/8 = 31, or 245/4 = 62
    long nstep = (long)(log2(subblocksize));
    
    long* nblock = new long[nstep];
    long* rot = new long[nstep];
    nblock[0] = subblocksize * nXbatching * factorDim;  //!  used block size in the first ciphertext
    rot[0] = (nblock[0] >> 1);
    for(long l = 1; l < nstep; ++l){
        nblock[l] = rot[l - 1];
        rot[l] = (nblock[l] >> 1);
    }
    
    uint64_t** poly0 = new uint64_t*[niter];
    uint64_t** poly1 = new uint64_t*[nstep];
    
    cipherLRPvals.generateRepAuxPoly(poly0, poly1, nslots, niter, nstep, nblock, rot);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.4.0. genpoly = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    // "+------------------------------------+"
    //! 3.4. encZWSData[nencsnp], L - 21
    // "+------------------------------------+"
   
    start = chrono::steady_clock::now();
    
    Ciphertext* encZWS = new Ciphertext[nencsnp];
    cipherLRPvals.encVecSData(encZWS, encZWData, encSData, poly0, poly1, sampleDim, nencsnp, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.4. Enc(Z^T * W * S) = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    // "+------------------------------------+"
    //! 3.5. encZWSData[nencsnp], L - 21
    // "+------------------------------------+"
    start = chrono::steady_clock::now();
    
    Ciphertext* encSWS = new Ciphertext[nencsnp];
    encWData1 = scheme.modDownBy(encWData, 1); //! (L - 14) -> (L - 15)
    cipherLRPvals.encVecSData(encSWS, encWData1, encSData, poly0, poly1, sampleDim, nencsnp, nslots, subblocksize, niter, nstep, nblock[0], rot);

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.5.  Enc(S^T * W * S)   = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    // "+------------------------------------+"
    //! 3.6. encW2SX[nencsnp][4], L - 21
    // "+------------------------------------+"
    
    start= chrono::steady_clock::now();
    
    Ciphertext** encW2SX = new Ciphertext*[nencsnp];
    cipherLRPvals.encVecMultipleSData(encW2SX, encW2Data, encSXData, poly0, poly1, factorDim, sampleDim, nencsnp, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.6. Enc(X^T * W2 * S)  = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    delete[] nblock;
    delete[] rot;
    
#if 1 // defined(__DEBUG_)
    //cout << L - encZWS[0].l << "," << L - encSWS[0].l << "," << L - encW2SX[0][0].l << endl;
    
    cout << "---------------" << endl;
    double* res = new double[nslots];
    cipherPvals.decVector(res, encZWS[0], nslots);
    cout << "Z^T * W ^ S : [" ;
    for(long l = 0; l < 10; ++l) cout << res[l] << ",";
    cout << "]" << endl;
    cout << "---------------" << endl;
    cipherPvals.decVector(res, encSWS[0], nslots);
    cout << "Z^T * W ^ S : [" ;
    for(long l = 0; l < 10; ++l) cout << res[l] << ",";
    cout << "]" << endl;
    
    cout << "---------------" << endl;
    cipherPvals.decVector(res, encSWS[0], nslots);
    cout << "S^T * W ^ S : [" ;
    for(long l = 0; l < 10; ++l){
        cout << res[l] << ",";
    }
    cout << "]" << endl;
    cout << "---------------" << endl;
    cout << "X^T * W2 ^ S : [" ;
    cipherPvals.decVector(res, encW2SX[0][0], nslots);
    for(long l = 0; l < 10; ++l){
        cout << res[l] << ",";
    }
    cout << "]" << endl;
    cout << "---------------" << endl;
#endif

    // "+------------------------------------+"
    //! 4. Norm
    // "+------------------------------------+"
    //! ZSnorm = <zstar, sstar> = (encDet * encZWS) - (encZX) * encAdj * (encXW2S)
    //! Snorm = (encDet * encSWS) - (encSX) * encAdj * (encXW2S)
    
    start= chrono::steady_clock::now();
    
    Ciphertext* encZSnorm = new Ciphertext[nencsnp];
    Ciphertext* encSnorm = new Ciphertext[nencsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        Ciphertext* encSX1 = new Ciphertext[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSX1[k] = encSX[k][j];
        }
        
        //! encZWS: (L - 21) -> encZSnorm: (L - 22)
        encZSnorm[j] = scheme.modDownTo(encDet, encZWS[j].l);
        extscheme.multAndEqual(encZSnorm[j], encZWS[j]);
        scheme.reScaleByAndEqual(encZSnorm[j], 1);
        
        //! encSWS: (L - 20) -> encZSnorm: (L - 21)
        encSnorm[j] = scheme.modDownTo(encDet, encSWS[j].l);
        extscheme.multAndEqual(encSnorm[j], encSWS[j]);
        scheme.reScaleByAndEqual(encSnorm[j], 1);
        
        Ciphertext res0; //! L - 22 = 1
        Ciphertext res1; //! L - 22 = 1

        cipherLRPvals.extQuadForm(res0, res1, encZX, encSX1, encAdj, encW2SX[j], factorDim);
        
        scheme.modDownToAndEqual(encZSnorm[j], res0.l);
        scheme.subAndEqual(encZSnorm[j], res0);
        
        scheme.modDownToAndEqual(encSnorm[j], res1.l);
        scheme.subAndEqual(encSnorm[j], res1);
    }
    NTL_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "4. Norm = " << timeElapsed << " s" << endl;
    cout << "logq: " << encZSnorm[0].l  << "," << encSnorm[0].l << endl;   // L - 3
    totalEvaltime += timeElapsed;
    cout << "Total Evaluation Timing = " << totalEvaltime << " s" << endl;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    ofstream outf(filename);
    outf.close();
    
    start = chrono::steady_clock::now();
    
    double** ZSnorm = new double*[nencsnp];
    double** Snorm = new double*[nencsnp];
    double det;
    
    cipherLRPvals.decryptResult(ZSnorm, Snorm, det, encZSnorm, encSnorm, encDet, nencsnp, nslots);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Decryption time = " << timeElapsed << " s" << endl;
    
    //! zscore = ZSnorm / sqrt(Snorm)
    zScore = new double[nsnp];   //! absolute value of beta/sigma
    pVals = new double[nsnp];
    long nslots0 = nsnp - (nencsnp - 1) * nslots;
    
    double* ZSnorm1 = new double[nsnp];
    double* Snorm1  = new double[nsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        if(j != nencsnp - 1){
            long starting = j * nslots;
            for(long l = 0; l < nslots; ++l){
                ZSnorm1[starting + l] = ZSnorm[j][l];
                Snorm1[starting + l] = Snorm[j][l];
                double beta = (ZSnorm[j][l]/Snorm[j][l]);
                double err = abs(det/Snorm[j][l]);
                zScore[starting + l] = beta / sqrt(abs(err));
                pVals[starting + l] =  pnorm(abs(zScore[starting + l]));
            }
        }
        else{
            long starting = (nencsnp - 1) * nslots;
            for(long l = 0; l < nslots0; ++l){
                double beta = (ZSnorm[j][l]/Snorm[j][l]);
                double err = abs(det/Snorm[j][l]);
                zScore[starting + l] = beta / sqrt(abs(err));
                pVals[starting + l] =  pnorm(abs(zScore[starting + l]));
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    printvectorToFile(ZSnorm1, "LogRegResult/HE_ZSnorm.txt", nsnp);
    printvectorToFile(Snorm1, "LogRegResult/HE_Snorm.txt", nsnp);
    printvectorToFile(pVals, filename, nsnp);
 
    delete[] encZX;
    delete[] encSX;
    delete[] encZWS;
    delete[] encSWS;
    delete[] encW2SX;
    
    delete[] encZSnorm;
    delete[] encSnorm;
    delete[] ZSnorm;
    delete[] Snorm;
    delete[] ZSnorm1;
    delete[] Snorm1;
#endif
}

//****************************************************************************************/
//!@ sdeg: degree of approximation polynomial of sigmoid (for testing)

void TestHELRPvals::testHELogReg_accuracy(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long numIter, long sdeg, string filename){
    
    long logN = 13;
    
    long Betalvl = 1 + 4 * (numIter - 1);
    long Prlvl;
    
    if(sdeg == 3){
        Prlvl = Betalvl + 4;
    } else{
       Prlvl = Betalvl + 5;
    }
    
    long L = Prlvl + 11;
    long K = 1;
    long h = 170;   //! Hamming weight of sk

    //------------------------------
    //! encryption level for snp data
    long Slvl = 3;
    long SXlvl = 3;    //! encW2.lvl

    //! encryption level for Xdata
    long YXlvl = L;
    long Xlvl;
    
    if(sdeg == 3){
        Xlvl = L - Betalvl;
    } else{
        Xlvl = L - Betalvl + 1;         //! need to rescale by "8" before computing IP
    }
    
    long Ylvl = L - Prlvl;
    long Covlvl = L - (Prlvl + 5);
    
    //! Parameters for GWAS
    long YXscaleBits = 2;
    long covscaleBits = 2;              //! scale factor for covariance
    long subblocksize = 16;             //! size of subblocks for fully replicating size-n vector

    //! HE Parameters
    long logp = 43;                     //! all the msg are scaled by "p", logq0 - logp = (final bits of precision), logp: scaled bit
    long logq0 = logp + 10;
    long logp0 = logq0 + 7;
    
    //! Parameters
    long nslots = (1 << (logN-1));     //! total number of plaintext slots
    long sampleDim2 = (1 << (long)ceil(log2(sampleDim)));   //! closet PoT
    long sdimBits = (long)ceil(log2(sampleDim));            //! log2(sampleDim)
    long fdimBits = (long)ceil(log2(factorDim));
    long factorDim2 = factorDim * factorDim;
    
    long nXbatching = nslots/(sampleDim2 * factorDim);     //! replicated number of a user's Xdata in a ciphertext
    long nCovbatching = nslots/(sampleDim2 * 8);           //! replicated number of a user's (X^T * X) in a ciphertext
    long bBits = (long)ceil(log2(nXbatching));             //! batching bits for XData
    
    long xBatchingBits = fdimBits + bBits;                 //! log2(number of slots for a single user in a ciphertext)
    long nslots1 =  (1 << xBatchingBits) * sampleDim;      //! number of packed messages in the slots
    
    long nencsnp = (long)ceil((double)nsnp/nslots);         //! number of ciphertexts for snp encryptions
    long nterms = factorDim * (factorDim + 1)/2;
    
    //! Parameters for Logistic Regression of covariates
    //long sigdeg = 3;
    //long numIter = 4;
    long scale = 4; //! scale for sigmoid approximation, a0 + (a1 * s) * (IP / s) + ...
    double gammaUp = 1;
    double gammaDown = 1;
    
    long logQ = logq0 + logp * (L - 1) + logp0;
    cout << "(logN,logp,L,K,logQ,h) = ("  << logN << "," << logp << "," << L << "," << K << "," << logQ << "," << h << ")" << endl;
    cout << "(dim,nslots,nbatching,nencsnp,subblock) = ("  ;
    cout << factorDim << "," << nslots <<  "," << nXbatching << "," << nencsnp << "," << subblocksize << ")" << endl;
    cout << "(niter,sdeg) = (" << numIter << "," << sdeg << ")" << endl; 
    
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
    for(long i = 0; i < fdimBits; ++i){
        extscheme.addDecompRightRotKey(secretKey, 1 << i); //! 1, 2
    }
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time= " << timeElapsed << " s" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Encryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    CipherPvals cipherPvals(scheme, secretKey, extscheme);
    CipherLRPvals cipherLRPvals(scheme, secretKey, extscheme, cipherPvals);
    
    start = chrono::steady_clock::now();
    
    Ciphertext encYXData;   //! E( - y[1] * X[1] - , ..., - y[n] * X[n] - )
    Ciphertext encXData;
    Ciphertext encYData;
    Ciphertext* enccovData = new Ciphertext[34];
    
    cipherLRPvals.encXData(encYXData, encXData, encYData, enccovData, yData, xData, factorDim, sampleDim, nXbatching, nCovbatching, nterms, YXscaleBits, covscaleBits, nslots, YXlvl, Xlvl, Ylvl, Covlvl);
    
    Ciphertext** encSData = new Ciphertext*[sampleDim];
    Ciphertext*** encSXData = new Ciphertext**[sampleDim];
    
    for(long i = 0; i < sampleDim; ++i){
        encSData[i] = new Ciphertext[nencsnp];
        encSXData[i] = new Ciphertext*[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSXData[i][k] = new Ciphertext[nencsnp];
        }
    }
    
    cipherLRPvals.encryptSData(encSData, encSXData, xData, sData, factorDim, sampleDim, nsnp, nencsnp, nslots, Slvl, SXlvl);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    double totalEvaltime;
    
    // "+------------------------------------+"
    //! 1. Logistic regression on covariates, 1 + (iter - 1) * 4 = 9
    // "+------------------------------------+"
    start = chrono::steady_clock::now();
    
    Ciphertext encBeta;
    uint64_t* poly = cipherLRPvals.generateNLGDAuxPoly(nslots, nslots1, factorDim);
    cipherLRPvals.encNLGDiteration(encBeta, encYXData, factorDim, sampleDim, fdimBits, sdimBits, xBatchingBits, poly, numIter);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "1. Logistic regression = " << timeElapsed << " s" << endl;
    totalEvaltime = timeElapsed;

    // "+------------------------------------+"
    //! 2. Update the weights
    // "+------------------------------------+"
    start = chrono::steady_clock::now();
    
    Ciphertext encZData;  //! (L - 9) - 10 = L - 19
    Ciphertext encWData;  //! (L - 9) - 5  = L - 14
    Ciphertext encW2Data; //! (L - 9) - 6 = L - 15
    Ciphertext encZWData; //! (L - 9) - 6 = L - 15
    
    long steps = 4; //! degree for approximation of 1/p*(1-p)
    cipherLRPvals.encZWData(encZData, encWData, encW2Data, encZWData, encBeta, encXData, encYData, poly, fdimBits, steps, sdeg, scale);
    
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "2. Update the weights = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
 
    // "+------------------------------------+"
    //! 3. Hessian
    // "+------------------------------------+"
    
    start = chrono::steady_clock::now();
    
    Ciphertext* encAdj = new Ciphertext[10];    //! (L - 21) or (L - 20) = 3
    Ciphertext encDet;
    
    Ciphertext encWData1 = scheme.modDownBy(encWData, 4);
    cipherLRPvals.encAdjoint(encDet, encAdj, encWData1, enccovData, sdimBits, nCovbatching); //! (X^T * W * X)
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.1. (X^T * W *  X)^-1 = " << timeElapsed << " s" << endl;
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
    cout << "1/(s^3) * det: 1.16008 ?= " << dtemp << endl;
#endif
 #if 1
    // "+------------------------------------+"
    //! encZXData[4], L - 21, [-43.8965903, -7.693702862, -7.157768273, -18.84315905],
    start= chrono::steady_clock::now();
    
    Ciphertext encZXData = scheme.modDownTo(encXData, encZData.l);
    extscheme.multAndEqual(encZXData, encZData);
    scheme.reScaleByAndEqual(encZXData, 1);
    
    Ciphertext* encZX = new Ciphertext[factorDim];
    cipherPvals.aggYXData(encZX, encZXData, sdimBits, nXbatching, factorDim, nslots);
    
#if 0
    cout << "Z: [" ;
    double* res = new double[16 * 10];
    cipherPvals.decVector(res, encZData, 16 * 10);
    for(long l = 0; l < 16 * 10; l += 16) cout << res[l] << endl;
    cout << "] " << endl ;
    
    cout << "X: [" ;
    cipherPvals.decVector(res, encXData, 16 * 10);
    for(long l = 0; l < 16 * 10; l += 16) cout << res[l + 1] << endl;
    cout << "] " << endl ;
    
    cout << "ZX: [" ;
    for(long i = 0; i < 4; ++i){
        double tmp;
        cipherPvals.decSingleData(tmp, encZX[i]);
        cout << tmp << "," ;
    }
    cout << "] " << endl ;
#endif

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.2. Z^T * X   = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    // "+------------------------------------+"
    //! encSX[nencsnp] = enc(S^T * X),
    start = chrono::steady_clock::now();
    
    Ciphertext** encSX = new Ciphertext*[factorDim];
    cipherLRPvals.encSXData(encSX, encSXData, factorDim, sampleDim, nencsnp);   //! (4*n ADD)
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.3. S^T * X  = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    // "+-----------------------------------------------------+"
    
    start = chrono::steady_clock::now();
    
    long niter = (long) ceil((double)sampleDim/subblocksize);    //! 245/16 = 16, 245/8 = 31, or 245/4 = 62
    //long nslots_block = subblocksize * nXbatching * factorDim;  //!  used block size in the first ciphertext
    long nstep = (long)(log2(subblocksize));
    
    long* nblock = new long[nstep];
    long* rot = new long[nstep];
    nblock[0] = subblocksize * nXbatching * factorDim;  //!  used block size in the first ciphertext
    rot[0] = (nblock[0] >> 1);
    for(long l = 1; l < nstep; ++l){
        nblock[l] = rot[l - 1];
        rot[l] = (nblock[l] >> 1);
    }
    
    uint64_t** poly0 = new uint64_t*[niter];
    uint64_t** poly1 = new uint64_t*[nstep];
    
    cipherLRPvals.generateRepAuxPoly(poly0, poly1, nslots, niter, nstep, nblock, rot);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.4.0. genpoly = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    
    //! 3.4. encZWData: (L - 15) -> (L - 16)
    //! encZWS[n]: (L - 21)
    start = chrono::steady_clock::now();
    
    Ciphertext* encZW = new Ciphertext[sampleDim];
    Ciphertext encZWData1 = scheme.modDownBy(encZWData, 1); // (L - 15) -> (L - 16)
    cipherLRPvals.fullReplicate(encZW, encZWData1, poly0, poly1, sampleDim, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.4.1. Enc(Z^T * W[i]) = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    start = chrono::steady_clock::now();
    
    Ciphertext* encZWS = new Ciphertext[nencsnp];
    cipherLRPvals.encZSData(encZWS, encZW, encSData, sampleDim, nencsnp);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.4.2. Enc(Z^T * W * S) = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    delete[] encZW;
    
    
    //! 3.5 encW: (L - 14) -> (L - 16)
    //! encW[n]: (L - 21)
    start = chrono::steady_clock::now();
    
    Ciphertext* encW = new Ciphertext[sampleDim];
    encWData1 = scheme.modDownBy(encWData, 2); //! (L - 14) -> (L - 16)
    cipherLRPvals.fullReplicate(encW, encWData1, poly0, poly1, sampleDim, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.5.1.  Enc(W[i])   = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    
    start = chrono::steady_clock::now();
    
    Ciphertext* encSWS = new Ciphertext[nencsnp];
    cipherLRPvals.encZSData(encSWS, encW, encSData, sampleDim, nencsnp);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.5.2.  Enc(S^T * W * S)   = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    delete[] encW;
    // "+------------------------------------+"
    //! encW2[n], (L - 15) - 4 = L - 19/  encXW2S[nencsnp][factordim], L - 20 (sub-block = 8)
    //! encW2[n], (L - 15) - 5 = L - 20/  encXW2S[nencsnp][factordim], L - 21 (sub-block = 16)
    
    start= chrono::steady_clock::now();
    
    Ciphertext* encW2 = new Ciphertext[sampleDim];
    
#if ver1
    Ciphertext encW2Data1 = scheme.modDownBy(encW2Data, 1);
    cipherLRPvals.fullReplicate(encW2, encW2Data1, poly0, poly1, sampleDim, nslots, subblocksize, niter, nstep, nblock[0], rot);
#endif
#if ver2
    cipherLRPvals.fullReplicate(encW2, encW2Data, poly0, poly1, sampleDim, nslots, subblocksize, niter, nstep, nblock[0], rot);
#endif
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.6. Enc(W2[i])  = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    start= chrono::steady_clock::now();
    
    Ciphertext** encW2SX = new Ciphertext*[nencsnp];
    cipherLRPvals.encW2SXData(encW2SX, encW2, encSXData, factorDim, sampleDim, nencsnp);
    
#if 0
    cipherLRPvals.encFastW2SXData(encW2SX, encW, encSXData, factorDim, sampleDim, nencsnp);
#endif
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.6. Enc(X^T * W2 * S)  = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    delete[] encW2;
#if 1 // defined(__DEBUG_)
    //cout << L - encZWS[0].l << "," << L - encSWS[0].l << "," << L - encW2SX[0][0].l << endl;
    
    cout << "---------------" << endl;
    double* res = new double[nslots];
    cipherPvals.decVector(res, encZWS[0], nslots);
    cout << "Z^T * W ^ S : [" ;
    for(long l = 0; l < 10; ++l){
        cout << res[l] << ",";
    }
    cout << "]" << endl;
    cout << "---------------" << endl;
    cipherPvals.decVector(res, encSWS[0], nslots);
    cout << "S^T * W ^ S : [" ;
    for(long l = 0; l < 10; ++l){
        cout << res[l] << ",";
    }
    cout << "]" << endl;
    cout << "---------------" << endl;
    cout << "X^T * W2 ^ S : [" ;
    cipherPvals.decVector(res, encW2SX[0][0], nslots);
    for(long l = 0; l < 10; ++l){
        cout << res[l] << ",";
    }
    cout << "]" << endl;
    cout << "---------------" << endl;
#endif
    
    // "+------------------------------------+"
    //! 3. Norm
    // "+------------------------------------+"
    //! ZSnorm = <zstar, sstar> = (encDet * encZWS) - (encZX) * encAdj * (encXW2S)
    //! Snorm = (encDet * encSWS) - (encSX) * encAdj * (encXW2S)
    
    start= chrono::steady_clock::now();
    
    Ciphertext* encZSnorm = new Ciphertext[nencsnp];
    Ciphertext* encSnorm = new Ciphertext[nencsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        Ciphertext* encSX1 = new Ciphertext[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSX1[k] = encSX[k][j];
        }
        
        //! encZWS: (L - 21) -> encZSnorm: (L - 22)
        encZSnorm[j] = scheme.modDownTo(encDet, encZWS[j].l);
        extscheme.multAndEqual(encZSnorm[j], encZWS[j]);
        scheme.reScaleByAndEqual(encZSnorm[j], 1);
        
        //! encSWS: (L - 20) -> encZSnorm: (L - 21)
        encSnorm[j] = scheme.modDownTo(encDet, encSWS[j].l);
        extscheme.multAndEqual(encSnorm[j], encSWS[j]);
        scheme.reScaleByAndEqual(encSnorm[j], 1);
        
        Ciphertext res0; //! L - 22 = 1
        Ciphertext res1; //! L - 22 = 1
#if ver1
        cipherLRPvals.extQuadForm(res0, res1, encZX, encSX1, encAdj, encW2SX[j], factorDim);
#endif
#if ver2
        cipherLRPvals.extQuadForm16(res0, res1, encZX, encSX1, encAdj, encW2SX[j], factorDim);
#endif
        
        scheme.modDownToAndEqual(encZSnorm[j], res0.l);
        scheme.subAndEqual(encZSnorm[j], res0);
        
        scheme.modDownToAndEqual(encSnorm[j], res1.l);
        scheme.subAndEqual(encSnorm[j], res1);
    }
    NTL_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "4. Norm = " << timeElapsed << " s" << endl;
    cout << "logq: " << encZSnorm[0].l  << "," << encSnorm[0].l << endl;   // L - 3
    totalEvaltime += timeElapsed;
    cout << "Total Evaluation Timing = " << totalEvaltime << " s" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    ofstream outf(filename);
    outf.close();
    
    start = chrono::steady_clock::now();
    
    double** ZSnorm = new double*[nencsnp];
    double** Snorm = new double*[nencsnp];
    double det;
    
    cipherLRPvals.decryptResult(ZSnorm, Snorm, det, encZSnorm, encSnorm, encDet, nencsnp, nslots);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Decryption time = " << timeElapsed << " s" << endl;
    
    //! zscore = ZSnorm / sqrt(Snorm)
    zScore = new double[nsnp];   //! absolute value of beta/sigma
    pVals = new double[nsnp];
    long nslots0 = nsnp - (nencsnp - 1) * nslots;
    
    double* ZSnorm1 = new double[nsnp];
    double* Snorm1  = new double[nsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        if(j != nencsnp - 1){
            long starting = j * nslots;
            for(long l = 0; l < nslots; ++l){
                ZSnorm1[starting + l] = ZSnorm[j][l];
                Snorm1[starting + l] = Snorm[j][l];
                double beta = (ZSnorm[j][l]/Snorm[j][l]);
                double err = abs(det/Snorm[j][l]);
                zScore[starting + l] = beta / sqrt(abs(err));
                pVals[starting + l] =  pnorm(abs(zScore[starting + l]));
            }
        }
        else{
            long starting = (nencsnp - 1) * nslots;
            for(long l = 0; l < nslots0; ++l){
                double beta = (ZSnorm[j][l]/Snorm[j][l]);
                double err = abs(det/Snorm[j][l]);
                zScore[starting + l] = beta / sqrt(abs(err));
                pVals[starting + l] =  pnorm(abs(zScore[starting + l]));
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    printvectorToFile(ZSnorm1, "LogRegResult/HE_ZSnorm.txt", nsnp);
    printvectorToFile(Snorm1, "LogRegResult/HE_Snorm.txt", nsnp);
    printvectorToFile(pVals, filename, nsnp);
    
    delete[] encZX;
    delete[] encSX;
    delete[] encZWS;
    delete[] encSWS;
    delete[] encW2SX;
    
    delete[] encZSnorm;
    delete[] encSnorm;
    delete[] ZSnorm;
    delete[] Snorm;
    delete[] ZSnorm1;
    delete[] Snorm1;
#endif
 
}


/*

void TestHELRPvals::testHELogReg_block8(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp,  string filename){
    
    //! Parameters for GWAS
    long YXscaleBits = 2;
    long covscaleBits = 2;  //! scale factor for covariance
    long subblocksize = 8;  //! size of subblocks for fully replicating size-n vector
    
    //! HE Parameters
    long logp = 43;                    //! all the msg are scaled by "p", logq0 - logp = (final bits of precision), logp: scaled bit
    long logq0 = 55;
    long logp0 = 61;
    
    long logN = 13;
    long nslots = (1<< (logN-1));      //! total number of plaintext slots
    long L = 23;                       //! 9 (LogReg) + 4 (Pr) + 1 (W) + 5 (Z, two-inverse) + 2 (Z^T * X) + 1 (multiplied by adj * XW2S)= 22
    long K = 1;
    long h = 170;                       //! Hamming weight of sk
    
    //! encryption level for snp data
    long Slvl = L - 20;
    long SXlvl = L - 19;    //! encW2.lvl,  fast: L - 18;
    
    //! encryption level for Xdata
    long YXlvl = L;
    long Xlvl = L - 9;      //! encBeta.lvl = 1 + (2 + log2(deg(sigmoid))) * (iter)
    long Ylvl = L - 13;     //! encPr.lvl
    long Covlvl = L - 17;
    
    //! Parameters
    long sampleDim2 = (1 << (long)ceil(log2(sampleDim)));   //! closet PoT
    long sdimBits = (long)ceil(log2(sampleDim));            //! log2(sampleDim)
    long fdimBits = (long)ceil(log2(factorDim));
    long factorDim2 = factorDim * factorDim;
    
    long nXbatching = nslots/(sampleDim2 * factorDim);     //! replicated number of a user's Xdata in a ciphertext
    long nCovbatching = nslots/(sampleDim2 * 8);           //! replicated number of a user's (X^T * X) in a ciphertext
    long bBits = (long)ceil(log2(nXbatching));             //! batching bits for XData
    
    long xBatchingBits = fdimBits + bBits;                 //! log2(number of slots for a single user in a ciphertext)
    long nslots1 =  (1 << xBatchingBits) * sampleDim;      //! number of packed messages in the slots
    
    long nencsnp = (long)ceil((double)nsnp/nslots);         //! number of ciphertexts for snp encryptions
    long nterms = factorDim * (factorDim + 1)/2;
    
    //! Parameters for Logistic Regression of covariates
    long sigdeg = 3;
    long numIter = 3;
    double gammaUp = 1;
    double gammaDown = 1;
    
    long logQ = logq0 + logp * (L - 1) + logp0;
    cout << "(logN,logp,L,K,logQ) = ("  << logN << "," << logp << "," << L << "," << K << "," << logQ << ")" << endl;
    cout << "(dim,nslots,nbatching,nencsnp,subblock) = ("  ;
    cout << factorDim << "," << nslots <<  "," << nXbatching << "," << nencsnp << "," << subblocksize << ")" << endl;
    
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
    extscheme.addDecompRightRotKeys(secretKey);
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time= " << timeElapsed << " s" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Encryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    CipherPvals cipherPvals(scheme, secretKey, extscheme);
    CipherLRPvals cipherLRPvals(scheme, secretKey, extscheme, cipherPvals);
    
    start = chrono::steady_clock::now();
    
    Ciphertext encYXData;   //! E( - y[1] * X[1] - , ..., - y[n] * X[n] - )
    Ciphertext encXData;
    Ciphertext encYData;
    Ciphertext* enccovData = new Ciphertext[34];
    
    
    cipherLRPvals.encXData(encYXData, encXData, encYData, enccovData, yData, xData, factorDim, sampleDim, nXbatching, nCovbatching, nterms, YXscaleBits, covscaleBits, nslots, YXlvl, Xlvl, Ylvl, Covlvl);
    
    Ciphertext** encSData = new Ciphertext*[sampleDim];
    Ciphertext*** encSXData = new Ciphertext**[sampleDim];
    
    for(long i = 0; i < sampleDim; ++i){
        encSData[i] = new Ciphertext[nencsnp];
        encSXData[i] = new Ciphertext*[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSXData[i][k] = new Ciphertext[nencsnp];
        }
    }
    
    cipherLRPvals.encryptSData(encSData, encSXData, xData, sData, factorDim, sampleDim, nsnp, nencsnp, nslots, Slvl, SXlvl);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    double totalEvaltime;
    
    // "+------------------------------------+"
    //! 1. Logistic regression on covariates, 1 + (iter - 1) * 4 = 9
    // "+------------------------------------+"
    start = chrono::steady_clock::now();
    
    Ciphertext encBeta;
    uint64_t* poly = cipherLRPvals.generateNLGDAuxPoly(nslots, nslots1, factorDim);
    cipherLRPvals.encNLGDiteration(encBeta, encYXData, factorDim, sampleDim, fdimBits, sdimBits, xBatchingBits, poly, numIter);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "1. Logistic regression = " << timeElapsed << " s" << endl;
    totalEvaltime = timeElapsed;
    
    // "+------------------------------------+"
    //! 2. Update the weights
    // "+------------------------------------+"
    start = chrono::steady_clock::now();
    
    Ciphertext encZData;  //! (L - 9) - 10 = L - 19
    Ciphertext encWData;  //! (L - 9) - 5  = L - 14
    Ciphertext encW2Data; //! (L - 9) - 6 = L - 15
    Ciphertext encZWData; //! (L - 9) - 6 = L - 15
    
    long steps = 4; //! degree for approximation of 1/p*(1-p)
    cipherLRPvals.encZWData(encZData, encWData, encW2Data, encZWData, encBeta, encXData, encYData, poly, fdimBits, steps);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "2. Update the weights = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    //cout << L - encZData.l << "," << L - encWData.l << "," << L - encW2Data.l << "," << L - encZWData.l << endl;
    
    // "+------------------------------------+"
    //! 3. Hessian
    // "+------------------------------------+"
    
    start = chrono::steady_clock::now();
    
    Ciphertext* encCov = new Ciphertext[34];    //! (L - 14) ->  (L - 17) -> (L - 18) = 4
    Ciphertext* encAdj = new Ciphertext[10];    //! (L - 20) = 3
    Ciphertext encDet;
    
    Ciphertext encWData1 = scheme.modDownBy(encWData, 3);
    
    cipherLRPvals.encWcov(encCov, encWData1, enccovData, sdimBits, nCovbatching); //! (X^T * W * X)
    cipherPvals.encSIMDAdjoint(encDet, encAdj, encCov);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.1. (X^T * W *  X)^-1 = " << timeElapsed << " s" << endl;
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
    cout << "1/(s^3) * det: 1.1623 ?= " << dtemp << endl;
#endif
    
    // "+------------------------------------+"
    //! encZXData[4], L - 21, [-43.8965903, -7.693702862, -7.157768273, -18.84315905]
    start= chrono::steady_clock::now();
    
    Ciphertext encZXData = scheme.modDownTo(encXData, encZData.l);
    extscheme.multAndEqual(encZXData, encZData);
    scheme.reScaleByAndEqual(encZXData, 1);
    
    Ciphertext* encZX = new Ciphertext[factorDim];
    cipherPvals.aggYXData(encZX, encZXData, sdimBits, nXbatching, factorDim, nslots);
    
    cout << "ZX: [" ;
    for(long i = 0; i < 4; ++i){
        double tmp;
        cipherPvals.decSingleData(tmp, encZX[i]);
        cout << tmp << "," ;
    }
    cout << "] " << endl ;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.2. Z^T * X   = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    // "+------------------------------------+"
    //! encSX[nencsnp] = enc(S^T * X)
    start = chrono::steady_clock::now();
    
    Ciphertext** encSX = new Ciphertext*[factorDim];
    cipherLRPvals.encSXData(encSX, encSXData, factorDim, sampleDim, nencsnp);   //! (4*n ADD)
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.3. S^T * X  = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    // "+-----------------------------------------------------+"
    //! encZWData: (L - 15) -> (L - 16),
    //! encZW[n]: (L - 21) or (L - 20)
    
    long niter = (long) ceil((double)sampleDim/subblocksize);    //! 245/16 = 16, 245/8 = 31, or 245/4 = 62
    //long nslots_block = subblocksize * nXbatching * factorDim;  //!  used block size in the first ciphertext
    long nstep = (long)(log2(subblocksize));
    
    long* nblock = new long[nstep];
    long* rot = new long[nstep];
    nblock[0] = subblocksize * nXbatching * factorDim;  //!  used block size in the first ciphertext
    rot[0] = (nblock[0] >> 1);
    for(long l = 1; l < nstep; ++l){
        nblock[l] = rot[l - 1];
        rot[l] = (nblock[l] >> 1);
    }
    
    
    uint64_t** poly0 = new uint64_t*[niter];
    uint64_t** poly1 = new uint64_t*[3];
    
    cipherLRPvals.generateRepAuxPoly(poly0, poly1, nslots, niter, nstep, nblock, rot);
    
    start = chrono::steady_clock::now();
    
    Ciphertext* encZW = new Ciphertext[sampleDim];
    Ciphertext encZWData1 = scheme.modDownBy(encZWData, 1); // (L - 15) -> (L - 16)
    cipherLRPvals.fullReplicate(encZW, encZWData1, poly0, poly1, sampleDim, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.4.1. Enc(Z^T * W[i]) = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    //! encZWS = (L - 22) or (L - 21)
    start = chrono::steady_clock::now();
    
    Ciphertext* encZWS = new Ciphertext[nencsnp];
    cipherLRPvals.encZSData(encZWS, encZW, encSData, sampleDim, nencsnp);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.4.2. Enc(Z^T * W * S) = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    // "+-----------------------------------------------------+"
    //! encW: (L - 14) -> (L - 16)
    //! encW[n]: (L - 21) or (L - 20)
    start = chrono::steady_clock::now();
    
    Ciphertext* encW = new Ciphertext[sampleDim];
    encWData1 = scheme.modDownBy(encWData, 2); // (L - 14) -> (L - 16)
    cipherLRPvals.fullReplicate(encW, encWData1, poly0, poly1, sampleDim, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.5.1.  Enc(W[i])   = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    
    start = chrono::steady_clock::now();
    
    Ciphertext* encSWS = new Ciphertext[nencsnp];
    cipherLRPvals.encZSData(encSWS, encW, encSData, sampleDim, nencsnp);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.5.2.  Enc(S^T * W * S)   = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    // "+------------------------------------+"
    //! encW2[n], (L - 15) - 4 = L - 19/  encXW2S[nencsnp][factordim], L - 20 (sub-block = 8)
    
    start= chrono::steady_clock::now();
    
    Ciphertext* encW2 = new Ciphertext[sampleDim];
    Ciphertext** encW2SX = new Ciphertext*[nencsnp];
    
#if 1
    cipherLRPvals.fullReplicate(encW2, encW2Data, poly0, poly1, sampleDim, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.6. Enc(W2[i])  = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    start= chrono::steady_clock::now();
    
    cipherLRPvals.encW2SXData(encW2SX, encW2, encSXData, factorDim, sampleDim, nencsnp);
#endif
    
#if 0
    cipherLRPvals.encFastW2SXData(encW2SX, encW, encSXData, factorDim, sampleDim, nencsnp);
#endif
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.6. Enc(X^T * W2 * S)  = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
#if 1 // defined(__DEBUG_)
    cout << L - encZW[0].l << "," << L - encSWS[0].l << "," << L - encW2SX[0][0].l << endl;
    
    cout << "---------------" << endl;
    cout << "ZW: [" ;
    for(long i = 0; i < 10; ++i){
        double tmp;
        cipherPvals.decSingleData(tmp, encZW[i]);
        cout << tmp << "," ;
    }
    cout << "] " << endl ;
    
    cout << "---------------" << endl;
    double* res = new double[nslots];
    cipherPvals.decVector(res, encZWS[0], nslots);
    cout << "Z^T * W ^ S : [" ;
    for(long l = 0; l < 10; ++l){
        cout << res[l] << ",";
    }
    cout << "]" << endl;
    cout << "---------------" << endl;
    
    cipherPvals.decVector(res, encSWS[0], nslots);
    cout << "S^T * W ^ S : [" ;
    for(long l = 0; l < 10; ++l){
        cout << res[l] << ",";
    }
    cout << "]" << endl;
    cout << "---------------" << endl;
    cout << "X^T * W2 ^ S : [" ;
    cipherPvals.decVector(res, encW2SX[0][0], nslots);
    for(long l = 0; l < 10; ++l){
        cout << res[l] << ",";
    }
    cout << "]" << endl;
    cout << "---------------" << endl;
#endif
    
    // "+------------------------------------+"
    //! 3. Norm
    // "+------------------------------------+"
    //! ZSnorm = <zstar, sstar> = (encDet * encZWS) - (encZX) * encAdj * (encXW2S)
    //! Snorm = (encDet * encSWS) - (encSX) * encAdj * (encXW2S)
    
    start= chrono::steady_clock::now();
    
    Ciphertext* encZSnorm = new Ciphertext[nencsnp];
    Ciphertext* encSnorm = new Ciphertext[nencsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        Ciphertext* encSX1 = new Ciphertext[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSX1[k] = encSX[k][j];
        }
        
        //! encZWS: (L - 21) -> encZSnorm: (L - 22)
        encZSnorm[j] = scheme.modDownTo(encDet, encZWS[j].l);
        extscheme.multAndEqual(encZSnorm[j], encZWS[j]);
        scheme.reScaleByAndEqual(encZSnorm[j], 1);
        
        //! encSWS: (L - 20) -> encZSnorm: (L - 21)
        encSnorm[j] = scheme.modDownTo(encDet, encSWS[j].l);
        extscheme.multAndEqual(encSnorm[j], encSWS[j]);
        scheme.reScaleByAndEqual(encSnorm[j], 1);
        
        Ciphertext res0; //! L - 22 = 1
        Ciphertext res1; //! L - 22 = 1
        
        cipherLRPvals.extQuadForm8(res0, res1, encZX, encSX1, encAdj, encW2SX[j], factorDim);
        
        scheme.modDownToAndEqual(encZSnorm[j], res0.l);
        scheme.subAndEqual(encZSnorm[j], res0);
        
        scheme.modDownToAndEqual(encSnorm[j], res1.l);
        scheme.subAndEqual(encSnorm[j], res1);
    }
    NTL_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "4. Norm = " << timeElapsed << " s" << endl;
    cout << "logq: " << encZSnorm[0].l  << "," << encSnorm[0].l << endl;   // L - 3
    totalEvaltime += timeElapsed;
    cout << "Total Evaluation Timing = " << totalEvaltime << " s" << endl;
    
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    ofstream outf(filename);
    outf.close();
    
    start = chrono::steady_clock::now();
    
    double** ZSnorm = new double*[nencsnp];
    double** Snorm = new double*[nencsnp];
    double det;
    
    cipherLRPvals.decryptResult(ZSnorm, Snorm, det, encZSnorm, encSnorm, encDet, nencsnp, nslots);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Decryption time = " << timeElapsed << " s" << endl;
    
    //! zscore = ZSnorm / sqrt(Snorm)
    zScore = new double[nsnp];   //! absolute value of beta/sigma
    pVals = new double[nsnp];
    long nslots0 = nsnp - (nencsnp - 1) * nslots;
    
    double* ZSnorm1 = new double[nsnp];
    double* Snorm1  = new double[nsnp];
    
    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        if(j != nencsnp - 1){
            long starting = j * nslots;
            for(long l = 0; l < nslots; ++l){
                ZSnorm1[starting + l] = ZSnorm[j][l];
                Snorm1[starting + l] = Snorm[j][l];
                double beta = (ZSnorm[j][l]/Snorm[j][l]);
                double err = abs(det/Snorm[j][l]);
                zScore[starting + l] = beta / sqrt(abs(err));
                pVals[starting + l] =  pnorm(abs(zScore[starting + l]));
            }
        }
        else{
            long starting = (nencsnp - 1) * nslots;
            for(long l = 0; l < nslots0; ++l){
                ZSnorm1[starting + l] = ZSnorm[j][l];
                Snorm1[starting + l] = Snorm[j][l];
                double beta = (ZSnorm[j][l]/Snorm[j][l]);
                double err = abs(det/Snorm[j][l]);
                zScore[starting + l] = beta / sqrt(abs(err));
                pVals[starting + l] =  pnorm(abs(zScore[starting + l]));
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    printvectorToFile(ZSnorm1, "LogRegResult/HE_ZSnorm.txt", nsnp);
    printvectorToFile(Snorm1, "LogRegResult/HE_Snorm.txt", nsnp);
    printvectorToFile(pVals, filename, nsnp);

}
*/

