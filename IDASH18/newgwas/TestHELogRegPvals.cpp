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

#include "TestLogRegPvals.h"
#include "CipherLinRegPvals.h"
#include "CipherLogRegPvals.h"
#include "TestHELogRegPvals.h"


void TestHELRPvals::new_testFastHELogReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, vector<string> snptag, string filename){
    
    struct rusage usage;
    long memoryscale = (1 << 30);   //! 2^20: linux, 2^30: mac
    
    //! Parameters for Logistic Regression of covariates
    long trainsigdeg = 3;
    long testsigdeg = 3;
    long numIter = 2;
    double gammaUp = 1;
    double gammaDown = 1;
    
    //! Parameters for GWAS
    long YXscaleBits = 2;
    long covscaleBits = 2;      //! scale factor for covariance
    long subblocksize = 16;     //! size of subblocks for fully replicating size-n vector
    
    //! HE Parameters
    long logp = 43;                    //! all the msg are scaled by "p", logq0 - logp = (final bits of precision), logp: scaled bit
    long logq0 = 51;
    long logp0 = 60;
    
    long logN = 15;
    long nslots = (1 << (logN-1));     //! total number of plaintext slots
    long L = 19;                       //! 5 (LogReg) + 4 (Pr) + 1 (W) + 5 (Z, two-inverse) + 2 (Z^T * X) + 1 (multiplied by adj * XW2S)= 22
    long K = 1;
    long h = 170;                       //! Hamming weight of sk (logq = 1103)
    
    //! encryption level for snp data
    long Slvl = 3;
    long SXlvl = 4;    //! encW2.lvl,  fast: L - 19;
    
    //! encryption level for Xdata
    long Betalvl = 1 + (2 + (long)ceil(log2(trainsigdeg))) * (numIter - 1);     // 9 when numIter = 3
    long YXlvl = L;
    long Xlvl = L - Betalvl;    //! encBeta.lvl = logreglvl
    long Ylvl = L - (Betalvl + (2 + (long)ceil(log2(testsigdeg))));     //! encPr.lvl
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
    
    
    
    long logQ = logq0 + logp * (L - 1) + logp0;
    cout << "(logN,logp,L,K,logQ,logq0,logp0) = ("  << logN << "," << logp << "," << L << "," << K << "," << logQ << "," << logq0 << "," << logp0 << ")" << endl;
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
    
    cipherLRPvals.encryptSData(encSData, encSXData, xData, sData, factorDim, sampleDim, nsnp, nencsnp, nslots, Slvl, SXlvl);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    //! MB
#if 1
    double sizeX = 1 * 2 * (logq0 + (Xlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeY = 1 * 2 * (logq0 + (Ylvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizeYX = 1 * 2 * (logq0 + (YXlvl - 1) * logp) * (1 << logN) / (1 << 23);
    double sizecov = 18 * 2 * (logq0 + (Covlvl - 1) * logp) * (1 << logN) / (1 << 23);
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
    cipherLRPvals.encNLGDiteration(encBeta, encYXData, factorDim, sampleDim, fdimBits, sdimBits, xBatchingBits, poly, numIter);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "1. Logistic regression = " << timeElapsed << " s" << endl;
    totalEvaltime = timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    double* cwData = new double[factorDim];
    cipherPvals.decVector(cwData, encBeta, factorDim);
    cout << "encBeta.l = " << encBeta.l << ", Beta = [" << cwData[0] << "," << cwData[1] << "," <<  cwData[2] << "," << cwData[3] << "]" << endl;
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
    cout << "2. Update the weights = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout << "encW.l = " << encWData.l << ", encZW.l = " << encZWData.l << endl;
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
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
    cipherLRPvals.new_encAdjoint(encDet, encAdj, encBeta, enccovData, sdimBits, nCovbatching); //! (X^T * W * X)
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3.1. (X^T * W *  X)^-1 = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    cout << "------------------------" << endl;
    delete[] enccovData;
    
#if 1// defined(__DEBUG_)
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
    
    cipherLRPvals.generateRepAuxPoly(poly0, poly1, nslots, niter, nstep, nblock, rot);
    
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
    
    cipherLRPvals.encVecSData(encSWX, encWData, encSXData, poly0, poly1, factorDim, sampleDim, nencsnp, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
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
    
    start= chrono::steady_clock::now();
    
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


void TestHELRPvals::testFastHELogReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, vector<string> snptag, string filename){
    
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
    
    long logN = 13;
    long nslots = (1 << (logN-1));     //! total number of plaintext slots
    long L = 23;                       //! 9 (LogReg) + 4 (Pr) + 1 (W) + 5 (Z, two-inverse) + 2 (Z^T * X) + 1 (multiplied by adj * XW2S)= 22
    long K = 1;
    long h = 170;                       //! Hamming weight of sk (logq = 1103)
    
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
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
    
    //! MB
#if 1
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
    cipherLRPvals.encNLGDiteration(encBeta, encYXData, factorDim, sampleDim, fdimBits, sdimBits, xBatchingBits, poly, numIter);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "1. Logistic regression = " << timeElapsed << " s" << endl;
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
    cout << "2. Update the weights = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
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
    cout << "3.1. (X^T * W *  X)^-1 = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    ret = getrusage(RUSAGE_SELF, &usage);
    cout<< "Memory Usage : " << (double) usage.ru_maxrss/(memoryscale)  << "(GB)" << endl;
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
    
    cipherLRPvals.generateRepAuxPoly(poly0, poly1, nslots, niter, nstep, nblock, rot);
    
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
    
    cipherLRPvals.encVecSData(encSWX, encWData, encSXData, poly0, poly1, factorDim, sampleDim, nencsnp, nslots, subblocksize, niter, nstep, nblock[0], rot);
    
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
    
    start= chrono::steady_clock::now();
    
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



