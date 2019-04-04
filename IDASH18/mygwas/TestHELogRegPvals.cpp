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


void TestHELRPvals::new_testFastHELogReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, vector<string> snptag, string filename, long numGDIter, double gammaUp){
    

    long subblocksize;          //! size of subblocks for fully replicating size-n vector, (=2^s)
    long logN;                  //! degree of cyclotomic polynomial
    long logp;                  //! all the msg are scaled by "p", logq0 - logp = (final bits of precision), logp: scaled bit
    long logq0;
    long logp0;
    
    if(numGDIter == 1){         //! accuracy: 0.00534701 (
        subblocksize = 16;      //! 128: increasing, 64(60s), 32(60s), 16()
        logN = 15;
        logp = 43;
        logq0 = 51;
        logp0 = 60;
    }
    else if(numGDIter == 2){     //! accuracy: 0.00556404
        subblocksize = 16;      //! 73s
        logN = 15;
        logp = 43;
        logq0 = 51;
        logp0 = 60;
    }
    else if(numGDIter == 3){     //! accuracy: 0.0052
        subblocksize = 16;      //! 8(240s = 75s/45s/90s), 16 (230s = 85s/41s/70s), 32( = 90s/40s/70s), 64(250s = 110s/40s/70s)
        logN = 16;
        logp = 45;
        logq0 = 54;
        logp0 = 62;
    }
    //! (45, 53, 62) => (0.994125,0.000548787,0.00587467,0.999451), scaledbits = 1
    //! (45, 54, 62) => (0.996726,0.000658183,0.00327439,0.999342)
    
    
    long L = (2 + 4 * numGDIter) + (1 + (long)ceil(log2(subblocksize))) + 4;
    

    struct rusage usage;
    long memoryscale = (1 << 30);   //! 2^20: linux, 2^30: mac
    
    //! Parameters for Logistic Regression of covariates
    long trainsigdeg = 3;
    long testsigdeg = 3;

    //! Parameters for GWAS
    long YXscaleBits = 2;
    long covscaleBits = 1;      //! scale factor for covariance
    
    //--------------
    //! HE Parameters
    //--------------
    long nslots = (1 << (logN-1));     //! total number of plaintext slots
    long K = 1;
   
    //! encryption level for snp/Xdata
    long Betalvl = 1 + (2 + (long)ceil(log2(trainsigdeg))) * (numGDIter - 1);     // 9 when numIter = 3
    long YXlvl = L;
    long Xlvl = L - Betalvl;    //! encBeta.lvl = logreglvl
    long Ylvl = L - (Betalvl + (2 + (long)ceil(log2(testsigdeg))));     //! encY.l = encPr.l
    long Covlvl = 4 + (long)ceil(log2(factorDim - 1));  //! k = 4: Covlvl = 6
    long Covlvl1 = 4;
    long Slvl = 3;
    long SXlvl = 4;    //! encW2.lvl,  fast: L - 19;
    
    //! Parameters for encoding/encryption
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
    
    
    long logQ = logq0 + logp * (L - 1) + logp0;
    cout << "(logN,logp,logq0,logp0,L,K,logQ) = ("  << logN << "," << logp<< "," << logq0 << "," << logp0 << "," << L << "," << K << "," << logQ << ")" << endl;
    cout << "(niter,dim,nslots,nbatching,nencsnp,subblock) = ("  ;
    cout << numGDIter << "," << factorDim << "," << nslots <<  "," << nXbatching << "," << nencsnp << "," << subblocksize << ")" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    auto start = chrono::steady_clock::now();
    
    Context context(logN, logp, L, K, logq0, logp0);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    ExtScheme extscheme(secretKey, context, scheme);
    
    //! evk for multiplication
    extscheme.addDecompTwoKeys(secretKey);
    extscheme.addDecompThreeKeys(secretKey);
    
    //! evk for rotation (generate all the left rotation keys/ right keys for 1 and 2)
    extscheme.addDecompLeftRotKeys(secretKey);
    
    NTL_EXEC_RANGE(fdimBits, first, last);
    for (long i = first; i < last; ++i) {
        extscheme.addDecompRightRotKey(secretKey, 1 << i); //! 1, 2
    }
    NTL_EXEC_RANGE_END;
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time = " << timeElapsed << " s" << endl;
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
    Ciphertext* enccovData = new Ciphertext[22];
    
    cipherLRPvals.new_encXData(encYXData, encXData, encYData, enccovData, yData, xData, factorDim, sampleDim, nXbatching, nCovbatching, nterms, YXscaleBits, covscaleBits, nslots, YXlvl, Xlvl, Ylvl, Covlvl, Covlvl1);
    
    Ciphertext** encSData = new Ciphertext*[sampleDim];     //! encSData[n][encsnp]
    Ciphertext*** encSXData = new Ciphertext**[sampleDim];  //! encSXData[n][encsnp][4]
    
    for(long i = 0; i < sampleDim; ++i){
        encSData[i] = new Ciphertext[nencsnp];
        encSXData[i] = new Ciphertext*[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSXData[i][k] = new Ciphertext[nencsnp];
        }
    }
    //cipherLRPvals.encryptSData(encSData, encSXData, xData, sData, factorDim, sampleDim, nsnp, nencsnp, nslots, Slvl, SXlvl);
    cipherLRPvals.new_encSData(encSXData, xData, sData, factorDim, sampleDim, nsnp, nencsnp, nslots, SXlvl);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    //! MB
#if defined(__DEBUG_)
    double sizeX = 1 * 2 * (logq0 + (Xlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeY = 1 * 2 * (logq0 + (Ylvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeYX = 1 * 2 * (logq0 + (YXlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizecov = 18 * 2 * (logq0 + (Covlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizecov1 = 4 * 2 * (logq0 + (Covlvl1 - 1) * logp) * (1 << logN) / (1 << 23);
    
    //double sizeS = sampleDim * nencsnp * 2 * (logq0 + (Slvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeSX = sampleDim * nencsnp * factorDim * 2 * (logq0 + (SXlvl - 1) * logp) * (1 << logN) / (1 << 23);
    cout << "Freshly Ciphertexts size: " << (sizeX + sizeY + sizeYX + sizecov + sizecov1 + sizeSX) << "(MB)" << endl;
#endif
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    double totalEvaltime;
    
    cout << "+------------------------------------+" << endl;
    cout << "|          Step0: LogReg/cov         |" << endl;
    cout << "+------------------------------------+" << endl;

    start = chrono::steady_clock::now();
    
    Ciphertext encBeta;
    uint64_t* poly = cipherLRPvals.generateNLGDAuxPoly(nslots, nslots1, factorDim);
    cipherLRPvals.encNLGD(encBeta, encYXData, factorDim, sampleDim, fdimBits, sdimBits, xBatchingBits, poly, numGDIter, gammaUp);
    
    end = std::chrono::steady_clock::now();
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    totalEvaltime = timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    
    cout << "Timing = " << timeElapsed << "(s) , Mem : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

#if 1 //defined(__DEBUG_)
    double* cwData = new double[factorDim];
    cipherPvals.decVector(cwData, encBeta, factorDim);
    cout << "encBeta.l = " << encBeta.l << ", Beta = [" << cwData[0] << "," << cwData[1] << "," <<  cwData[2] << "," << cwData[3] << "]" << endl;
    cout << "--------------------------------------" << endl;
#endif
    
    cout << "+------------------------------------+" << endl;
    cout << "|     Step 1 & 2: Update weights     |" << endl;
    cout << "+------------------------------------+" << endl;
    start = chrono::steady_clock::now();
    
    Ciphertext encWZData;   //! lvl = 9
    Ciphertext encWData;    //! lvl = 8
    cipherLRPvals.encZWData(encWData, encWZData, encBeta, encXData, encYData, poly, fdimBits);
    
    end = std::chrono::steady_clock::now();
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    
    cout << "Timing = " << timeElapsed << "(s) , Mem : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "encW.l = " << encWData.l << ", encWZ.l = " << encWZData.l << endl;
    delete[] poly;
    
    

    cout << "+------------------------------------+" << endl;
    cout << "|     Step 3.0: polyGen (for rep)    |" << endl;
    cout << "+------------------------------------+" << endl;
    
    start = chrono::steady_clock::now();
    
    long niter = (long) ceil((double)sampleDim/subblocksize);    //! 245/32 = 8, 245/16 = 16, 245/8 = 31, or 245/4 = 62
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
    
    cipherLRPvals.generateRepAuxPoly(poly0, poly1, nslots, niter, nstep, nblock, rot, encWData.l);
    
    end = std::chrono::steady_clock::now();
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    cout << "Timing = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
    cout << "+------------------------------------+" << endl;
    cout << "|    Step 3 & 5: Enc(sj^T * W * z)   |" << endl;
    cout << "+------------------------------------+" << endl;

    //! Input: encSData[n][encsnp], encSXData[n][encsnp][4]
    //! encSWZ = Enc(s1^T * Wz, ...., sp^T * Wz) where (sj^T * WZ) = sum_i s[i][j] * (w[i] * z[i])
    start = chrono::steady_clock::now();
    
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < nencsnp; ++j){
            encSData[i][j] = scheme.modDownTo(encSXData[i][j][0], 3);   //! 4 -> 3
        }
    }

    Ciphertext* encSWZ = new Ciphertext[nencsnp];
    cipherLRPvals.encVecSData(encSWZ, encWZData, encSData, poly0, poly1, sampleDim, nencsnp, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
    end = std::chrono::steady_clock::now();
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout << "Timing = " << timeElapsed << "(s) , Mem : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
   
    cout << "+------------------------------------+" << endl;
    cout << "|    Step 3 & 4: Enc(sj^T * W * X)   |" << endl;
    cout << "+------------------------------------+" << endl;
    
    start = chrono::steady_clock::now();
    
    Ciphertext** encSWX = new Ciphertext*[nencsnp];
    
    cipherLRPvals.encVecSXData(encSWX, encWData, encSXData, poly0, poly1, factorDim, sampleDim, nencsnp, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
    end = std::chrono::steady_clock::now();
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout << "Timing = " << timeElapsed << "(s) , Mem : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    delete[] nblock;
    delete[] rot;
    
    
    cout << "+------------------------------------+" << endl;
    cout << "|         Step 6: X^T * W * z        |" << endl;
    cout << "+------------------------------------+" << endl;

    start= chrono::steady_clock::now();
    
    Ciphertext* encXWZ = new Ciphertext[factorDim];
    
    Ciphertext encX1Data = scheme.modDownTo(encXData, 4);
    Ciphertext encWZ1Data = scheme.modDownTo(encWZData, 4);
    
    cipherLRPvals.encZXData(encXWZ, encX1Data, encWZ1Data, sdimBits, nXbatching, factorDim, nslots);
    
    end = std::chrono::steady_clock::now();
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout << "Timing = " << timeElapsed << "(s) , Mem : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
#if 1 //defined(__DEBUG_)
    
    cout << "(lvl(X^T * W * z), lvl(sj^T * W * Z), lvl(SWX)) = " << encXWZ[0].l << "," <<  encSWZ[0].l  << "," << encSWX[0][0].l << endl;
    //cout << "---------------" << endl;
    double* res = new double[nslots];
    
    cout << "----------------------" << endl;
    cipherPvals.decVector(res, encSWZ[0], nslots);
    cout << "(Step 3&5) St * W * z : [" ;
    for(long l = 0; l < 10; ++l) cout << res[l] << ",";
    cout << "]" << endl;
    
    cout << "----------------------" << endl;
    cout << "(Step 4&5) St * W * X: [" ;
    cipherPvals.decVector(res, encSWX[0][0], nslots);
    for(long l = 0; l < 10; ++l){
        cout << res[l] << ",";
    }
    cout << "]" << endl;
    
    cout << "----------------------" << endl;
    cout << "(Step 6) Zt * W * X : [" ;
    //! plain: iter1 : -12.49899303,-2.500586912,-2.394780179,-5.459060554,
    //! iter2 : -12.48347256,-2.49475959,-2.388665939,-5.451455696
    //! iter3 : -10.95658524,-1.921525908,-1.787232471,-4.703389604
    for(long j = 0; j < 4; ++j){
        double dtemp;
        cipherPvals.decSingleData(dtemp, encXWZ[j]);
        if(j < 3) cout << dtemp << ",";
        else{
            cout << dtemp << "]" << endl;
        }
    }
    cout << "----------------------" << endl;
#endif
   
    

#if 1
    
    cout << "+------------------------------------+" << endl;
    cout << "|      Step 7: (X^T * W * X)^-1      |" << endl;
    cout << "+------------------------------------+" << endl;
    
    start = chrono::steady_clock::now();
    
    Ciphertext* encAdj = new Ciphertext[10];    //! (L - 21) or (L - 20) = 3
    Ciphertext encDet;
    Ciphertext enctmpWData = scheme.modDownTo(encWData, Covlvl);
    cipherLRPvals.new_encAdjoint(encDet, encAdj, enctmpWData, enccovData, sdimBits, nCovbatching); //! (X^T * W * X)
    
    end = std::chrono::steady_clock::now();
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout << "Timing = " << timeElapsed << "(s) , Mem : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    delete[] enccovData;
    
#if 1 // defined(__DEBUG_)
    //cout << encAdj[0].l << endl;
    double dtemp;
    cout << "1/(s^3) * Adj.lvl = " << encAdj[0].l <<   ", [" ;
    for(long l = 0; l < nterms; ++l){
        cipherPvals.decSingleData(dtemp, encAdj[l]);
        cout << dtemp << ", " ;
    }
    cout << "] " << endl ;
    
    cipherPvals.decSingleData(dtemp, encDet);
    cout << "1/(s^3) * det: 9.3 ?= " << dtemp << endl;
#endif
    
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Step 8: Norm           |" << endl;
    cout << "+------------------------------------+" << endl;
    //! ZSnorm = <zstar, sstar> = (encDet * encZWS) - (encZX) * encAdj * (encXW2S)
    //! Snorm = (encDet * encSWS) - (encSX) * encAdj * (encXW2S)
    
    start= chrono::steady_clock::now();
    
    Ciphertext* encZSnorm = new Ciphertext[nencsnp];
    Ciphertext* encSnorm = new Ciphertext[nencsnp];
    Ciphertext* encSWS = new Ciphertext[nencsnp];
    

    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        //! encZWS: 2 -> encZSnorm: 1
        encZSnorm[j] = scheme.modDownTo(encDet, encSWZ[j].l);
        extscheme.multAndEqualMT(encZSnorm[j], encSWZ[j]);
        scheme.reScaleByAndEqual(encZSnorm[j], 1);
        
        //! encSWS: 3 -> encZSnorm: 2
        encSnorm[j] = encSWX[j][0];
        scheme.modDownToAndEqual(encSnorm[j], encDet.l);
        extscheme.multAndEqualMT(encSnorm[j], encDet);
        scheme.reScaleByAndEqual(encSnorm[j], 1);
        
        Ciphertext res0; //! lvl = 1
        Ciphertext res1; //! lvl = 1
        
        cipherLRPvals.extQuadForm(res0, res1, encXWZ, encAdj, encSWX[j], factorDim);
        
        scheme.modDownToAndEqual(encZSnorm[j], res0.l);
        scheme.subAndEqual(encZSnorm[j], res0);
        
        scheme.modDownToAndEqual(encSnorm[j], res1.l);
        scheme.subAndEqual(encSnorm[j], res1);
    }
    NTL_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    
    cout << "Timing = " << timeElapsed << "(s) , Mem : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "*** Total Evaluation Timing = " << totalEvaltime << " s ***" << endl;
    //cout << "logq: " << encZSnorm[0].l  << "," << encSnorm[0].l << endl;   // L - 3
    
    
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
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout << "Timing = " << timeElapsed << "(s) , Mem : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|           Reconstruction           |" << endl;
    cout << "+------------------------------------+" << endl;
    
     start = chrono::steady_clock::now();
    
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
    
    end = std::chrono::steady_clock::now();
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout << "Timing = " << timeElapsed << "(s) , Mem : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
//    printvectorToFile(ZSnorm1, "LogRegResult/HE_ZSnorm.txt", nsnp);
//    printvectorToFile(Snorm1, "LogRegResult/HE_Snorm.txt", nsnp);
    printPvalsToFile(pVals, snptag, filename, nsnp);
    
    delete[] encXWZ;
    delete[] encSWZ;
    delete[] encSWS;
    
    delete[] encZSnorm;
    delete[] encSnorm;
    delete[] ZSnorm;
    delete[] Snorm;
    delete[] ZSnorm1;
    delete[] Snorm1;
#endif
}


//! submitted version of idash18
void TestHELRPvals::testFastHELogReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, vector<string> snptag, string filename){

    struct rusage usage;
    long memoryscale = (1 << 30);   //! 2^20: linux, 2^30: mac

    //! Parameters for GWAS
    long YXscaleBits = 2;
    long covscaleBits = 2;      //! scale factor for covariance
    long subblocksize = 16;     //! size of subblocks for fully replicating size-n vector (= 2^s)

    //! HE Parameters
    long logp = 43;                    //! all the msg are scaled by "p", logq0 - logp = (final bits of precision), logp: scaled bit
    long logq0 = 53;
    long logp0 = 61;

    long logN = 15;
    long nslots = (1 << (logN-1));     //! total number of plaintext slots
    long L = 23;                       //! 9 (LogReg) + 4 (Pr) + 1 (W) + 5 (Z, two-inverse) + 2 (Z^T * X) + 1 (multiplied by adj * XW2S)= 22
    long K = 1;

    //! encryption level for snp data
    long Slvl = 3;
    long SXlvl = 4;    //! encW2.lvl,  fast: L - 19;

    //! encryption level for Xdata
    long YXlvl = L;
    long Xlvl = L - 9;      //! encBeta.lvl = 1 + (2 + log2(deg(sigmoid))) * (iter)
    long Ylvl = L - 13;     //! encPr.lvl
    long Covlvl = 6;        //! one level: mult by W, two levels: adj -> 3
    long Covlvl1 = 4;

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
    long numGDIter = 3;
    double gammaUp = 1;
    double gammaDown = 1;

    long logQ = logq0 + logp * (L - 1) + logp0;
    cout << "(logN,logp,L,K,logQ,h) = ("  << logN << "," << logp << "," << L << "," << K << "," << logQ  << ")" << endl;
    cout << "(dim,nslots,nbatching,nencsnp,subblock) = ("  ;
    cout << factorDim << "," << nslots <<  "," << nXbatching << "," << nencsnp << "," << subblocksize << ")" << endl;

    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;

    auto start = chrono::steady_clock::now();

    Context context(logN, logp, L, K, logq0, logp0);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    ExtScheme extscheme(secretKey, context, scheme);

    //! evk for multiplication
    extscheme.addDecompTwoKeys(secretKey);
    extscheme.addDecompThreeKeys(secretKey);

    //! evk for rotation (generate all the left rotation keys/ right keys for 1 and 2)
    extscheme.addDecompLeftRotKeys(secretKey);

    NTL_EXEC_RANGE(fdimBits, first, last);
    for (long i = first; i < last; ++i) {
        extscheme.addDecompRightRotKey(secretKey, 1 << i); //! 1, 2
    }
    NTL_EXEC_RANGE_END;

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time = " << timeElapsed << " s" << endl;
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

    cipherLRPvals.encXData(encYXData, encXData, encYData, enccovData, yData, xData, factorDim, sampleDim, nXbatching, nCovbatching, nterms, YXscaleBits, covscaleBits, nslots, YXlvl, Xlvl, Ylvl, Covlvl, Covlvl1);

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
    timeElapsed = chrono::duration <double, milli> (end - start).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    //! MB
#if 0
    double sizeX = 1 * 2 * (logq0 + (Xlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeY = 1 * 2 * (logq0 + (Ylvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeYX = 1 * 2 * (logq0 + (YXlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizecov = 30 * 2 * (logq0 + (Covlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizecov1 = 4 * 2 * (logq0 + (Covlvl1 - 1) * logp) * (1 << logN) / (1 << 23);

    double sizeS = sampleDim * nencsnp * 2 * (logq0 + (Slvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeSX = sampleDim * nencsnp * factorDim * 2 * (logq0 + (SXlvl - 1) * logp) * (1 << logN) / (1 << 23);
    cout << "Freshly Ciphertexts size: " << (sizeX + sizeY + sizeYX + sizecov + sizecov1 + sizeS + sizeSX) << "(MB)" << endl;
#endif

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
    cipherLRPvals.encNLGD(encBeta, encYXData, factorDim, sampleDim, fdimBits, sdimBits, xBatchingBits, poly, numGDIter);

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "1. Logistic regression / covariates = " << timeElapsed << " s" << endl;
    totalEvaltime = timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;

    double* cwData = new double[factorDim];
    cipherPvals.decVector(cwData, encBeta, factorDim);
    cout << "[" << cwData[0] << "," << cwData[1] << "," <<  cwData[2] << "," << cwData[3] << "]" << endl;
    cout << "------------------------" << endl;

    // "+------------------------------------+"
    //! 2. Update the weights
    // "+------------------------------------+"
    start = chrono::steady_clock::now();

    Ciphertext encZWData;   //! lvl = 8
    Ciphertext encWData;    //! lvl = 7
    cipherLRPvals.encZWData(encWData, encZWData, encBeta, encXData, encYData, poly, fdimBits);

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "2. Update the weights (w, z) = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "encWData.l = " << encWData.l << ", encZWData.l = " << encZWData.l << endl;
    cout << "------------------------" << endl;

    delete[] poly;

#if 1
    // "+------------------------------------+"
    //! 3.1. Hessian: encAdj[4], encDet
    // "+------------------------------------+"

    start = chrono::steady_clock::now();

    Ciphertext* encAdj = new Ciphertext[10];    //! (L - 21) or (L - 20) = 3
    Ciphertext encDet;
    encBeta = scheme.modDownBy(encWData, 3);
    cipherLRPvals.encAdjoint(encDet, encAdj, encBeta, enccovData, sdimBits, nCovbatching); //! (X^T * W * X)

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3. (X^T * W *  X)^-1 = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "encAdj.l = " << encAdj[0].l << ", encDet.l = " << encDet.l << endl;
    cout << "------------------------" << endl;
    delete[] enccovData;

#if  1 // defined(__DEBUG_)
    cout << encAdj[0].l << endl;
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
    //! 3.2. encZWXData[4], lvl = 4 -> 2
    // "+------------------------------------+"
    // [-43.8965903, -7.693702862, -7.157768273, -18.84315905]

    start= chrono::steady_clock::now();

    Ciphertext* encZWX = new Ciphertext[factorDim];
    encBeta = scheme.modDownTo(encXData, 4);
    Ciphertext encZW1Data = scheme.modDownTo(encZWData, 4);
    cipherLRPvals.encZXData(encZWX, encBeta, encZW1Data, sdimBits, nXbatching, factorDim, nslots);

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.2. Z^T * W * X   = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    cout << "------------------------" << endl;


    // "+------------------------------------+"
    //! 3.4.0 polyGen
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

    cipherLRPvals.generateRepAuxPoly(poly0, poly1, nslots, niter, nstep, nblock, rot, encWData.l);

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.4.0. genpoly = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    cout << "------------------------" << endl;

    // "+------------------------------------+"
    //! 3.4. encZWSData[nencsnp], lvl = 2
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
    cout << "------------------------" << endl;

    // "+------------------------------------+"
    //! 3.5. encSWSData[nencsnp], lvl = 2
    //! 3.6. encSWXData[nencsnp][factorDim], lvl = 3
    // "+------------------------------------+"
    start = chrono::steady_clock::now();

    Ciphertext** encSWX = new Ciphertext*[nencsnp];

    cipherLRPvals.encVecSXData(encSWX, encWData, encSXData, poly0, poly1, factorDim, sampleDim, nencsnp, nslots, subblocksize, niter, nstep, nblock[0], rot);

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.5.  Enc(S^T * W * X)  = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "------------------------" << endl;

    delete[] nblock;
    delete[] rot;

#if 1 // defined(__DEBUG_)
    cout << "(encZWX.l, ZWS.l, SWX.l) = " << encZWX[0].l << "," <<  encZWS[0].l  << "," << encSWX[0][0].l << endl;
    //cout << "---------------" << endl;
    double* res = new double[nslots];

    cout << "Z^T * W ^ X : [" ;
    for(long j = 0; j < 4; ++j){
        double dtemp;
        cipherPvals.decSingleData(dtemp, encZWX[j]);
        cout << dtemp << ",";
    }
    cout << "]" << endl;
    cout << "---------------" << endl;
    cipherPvals.decVector(res, encZWS[0], nslots);
    cout << "Z^T * W ^ S : [" ;
    for(long l = 0; l < 10; ++l) cout << res[l] << ",";
    cout << "]" << endl;
    cout << "---------------" << endl;
    cout << "X^T * W ^ S : [" ;
    cipherPvals.decVector(res, encSWX[0][0], nslots);
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

    start = chrono::steady_clock::now();

    Ciphertext* encZSnorm = new Ciphertext[nencsnp];
    Ciphertext* encSnorm = new Ciphertext[nencsnp];
    Ciphertext* encSWS = new Ciphertext[nencsnp];


    NTL_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        //! encZWS: 2 -> encZSnorm: 1
        encZSnorm[j] = scheme.modDownTo(encDet, encZWS[j].l);
        extscheme.multAndEqualMT(encZSnorm[j], encZWS[j]);
        scheme.reScaleByAndEqual(encZSnorm[j], 1);

        //! encSWS: 3 -> encZSnorm: 2
        encSnorm[j] = encSWX[j][0];
        scheme.modDownToAndEqual(encSnorm[j], encDet.l);
        extscheme.multAndEqualMT(encSnorm[j], encDet);
        scheme.reScaleByAndEqual(encSnorm[j], 1);

        Ciphertext res0; //! lvl = 1
        Ciphertext res1; //! lvl = 1

        cipherLRPvals.extQuadForm(res0, res1, encZWX, encAdj, encSWX[j], factorDim);

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
    cout << "------------------------" << endl;

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
    //printvectorToFile(pVals, filename, nsnp);
    printPvalsToFile(pVals, snptag, filename, nsnp);

    delete[] encZWX;
    delete[] encZWS;
    delete[] encSWS;

    delete[] encZSnorm;
    delete[] encSnorm;
    delete[] ZSnorm;
    delete[] Snorm;
    delete[] ZSnorm1;
    delete[] Snorm1;
#endif
}



