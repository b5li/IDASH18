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

#include <thread>

// #include "NTL/ZZX.h"
// #include <NTL/RR.h>
// #include "NTL/vec_RR.h"
// #include "NTL/mat_RR.h"
// #include <NTL/BasicThreadPool.h>

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
#include "TestHEPvals.h"

#include "threadpool.h" 

void TestHEPvals::testHELinReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp,  string filename){
    
    //! Parameters for approx-HE
    long logp = 40;
    long logN = 13;
    long nslots = (1<< (logN-1)); //! total number of plaintext slots
    long L = 5;
    long K = L;

    //! Parameters for GWAS
    long factorDim2 = factorDim * factorDim;
    long dim = (1 << (long)ceil(log2(factorDim))); //! closet PoT
    long dim2 = dim * dim;
    
    long nencsnp = (long)ceil((double)nsnp/nslots);  // number of snp encryptions
    long nterms = dim * (dim + 1)/2;    
    
    cout << "(logN,logp,L,K) = ("  << logN << "," << logp << "," << L << "," << K << ")" << endl;
    cout << "(dim,nslots,nencsnp, nterms) = ("  << dim << "," << nslots  << "," << nencsnp<< "," << nterms  << ")" << endl;
 
    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
 
    auto start= chrono::steady_clock::now();
    
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Scheme generation time= " << timeElapsed << " s" << endl;
    
    CipherPvals cipherPvals(scheme, secretKey);
 

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
    
    cipherPvals.encryptXData(encYXData, enccovData, yData, xData, factorDim, sampleDim, dim, nslots, L);
    cipherPvals.encryptSData(encSData, encYSData, encSXData, yData, xData, sData, factorDim, sampleDim, nsnp, nencsnp,nslots, L);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Encryption time (X and snp) = " << timeElapsed << " s" << endl;
    

    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;

    double totalEvaltime;
    
    //! 1. Aggregation
    //! encYS[j] = Y^T * S = sum_i encYSData[i][j]
    //! encS[j]  = sum_i SData[i][j]^2
    
    start= chrono::steady_clock::now();
    
    Ciphertext* encYS = new Ciphertext[nencsnp];
    Ciphertext* encS = new Ciphertext[nencsnp];
    
    TP_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
       encYS[j] = encYSData[0][j];
       encS[j] = encSData[0][j];
       for(long i = 1; i < sampleDim; ++i){
          scheme.addAndEqual(encYS[j], encYSData[i][j]);
          scheme.addAndEqual(encS[j], encSData[i][j]);
       }
    }
    TP_EXEC_RANGE_END;
    
    //! encYX[k] = YXvec = y^T * X[k] = [-29, -6.506140997, -6.368385672, -12.88024343]
    //! encSX[factorDim][nencsnp]
    
    Ciphertext* encYX = new Ciphertext[factorDim];
    Ciphertext** encSX = new Ciphertext*[factorDim];
    
    TP_EXEC_RANGE(factorDim, first, last);
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
    TP_EXEC_RANGE_END;
    
#if  defined(__DEBUG_)
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
    

    //! 2. Hessian
    start= chrono::steady_clock::now();
    Ciphertext encDet;
    Ciphertext* encAdj = new Ciphertext[nterms];
    cipherPvals.HesInverse(encDet, encAdj, enccovData, dim, nslots, L);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "2. (X^T * X)^-1 = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    

    //! 3. Norm
    //! encwtYX = (y^T * X) * adj(X^T * X) * (X^T * y)
    //! Ynorm = |A| <Y,Y> - wtYX = |A| * sampleDim - wtYX
    
    
    start= chrono::steady_clock::now();
    
    Ciphertext encWtYX;
    cipherPvals.SqrQuadForm(encWtYX, encYX, encAdj, factorDim);  // 12.7 * 19152/ scalefactor
    
    Ciphertext encYnorm = scheme.multByConst(encDet,  (double) sampleDim);
    scheme.reScaleByAndEqual(encYnorm, 1);
    scheme.subAndEqual(encYnorm, encWtYX);
    
    //! YSnorm = <ystar, sstar> = encYS - encWtYS
    Ciphertext* encWtYS = new Ciphertext[nencsnp];
    Ciphertext* encYSnorm = new Ciphertext[nencsnp];
    Ciphertext* encWtSS  = new Ciphertext[nencsnp];
    Ciphertext* encSnorm = new Ciphertext[nencsnp];
    
    TP_EXEC_RANGE(nencsnp, first, last);
    for(long j = first; j < last; ++j){
        Ciphertext* encSX1 = new Ciphertext[factorDim];
        for(long k = 0; k < factorDim; ++k){
            encSX1[k] = encSX[k][j];
        }

        cipherPvals.QuadForm(encWtYS[j], encYX, encAdj, encSX1, factorDim);
        encYSnorm[j] = scheme.modDownTo(encYS[j], encDet.l);
        scheme.multAndEqual(encYSnorm[j], encDet);
        scheme.reScaleByAndEqual(encYSnorm[j], 1);
        scheme.subAndEqual(encYSnorm[j], encWtYS[j]);
      
        cipherPvals.SqrQuadForm(encWtSS[j], encSX1, encAdj, factorDim); // L - 3
        encSnorm[j] = scheme.modDownTo(encS[j], encDet.l);
        scheme.multAndEqual(encSnorm[j], encDet); // L - 3
        scheme.reScaleByAndEqual(encSnorm[j], 1);
        scheme.subAndEqual(encSnorm[j], encWtSS[j]);
    }
    TP_EXEC_RANGE_END;
    
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
    
    TP_EXEC_RANGE(nencsnp, first, last);
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
    TP_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "Decryption time = " << timeElapsed << " s" << endl;
    
    //printvectorToFile(zScore, filename, nsnp);
    printvectorToFile(pVals, filename, nsnp);
}
