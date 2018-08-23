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

#ifdef USE_NTL
#include <NTL/BasicThreadPool.h>
#endif

#include "threadpool.h"

#include "NTL/ZZX.h"
#include <NTL/RR.h>
#include "NTL/vec_RR.h"
#include "NTL/mat_RR.h"
#include <NTL/BasicThreadPool.h>

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

#include "sys.h" 

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
    long scalefactor = 32.0;  //! scale factor for covariance
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


    std::cout << "allocate ciphertext" << std::endl;

    {
       MemoryUsage mem = getMemoryUsage();
       cout << "Peak memory = " << mem.vmpeak/1024 << "KB" << std::endl;
       cout << "Curr memory = " << mem.vmrss/1024  << "KB" << std::endl;
    }

    cipherPvals.encryptXData(encYXData, enccovData, yData, xData, factorDim, sampleDim, dim, nterms, scalefactor, nslots, L);
    std::cout << "encrypt X" << std::endl;

    {
       MemoryUsage mem = getMemoryUsage();
       cout << "Peak memory = " << mem.vmpeak/1024 << "KB" << std::endl;
       cout << "Curr memory = " << mem.vmrss/1024  << "KB" << std::endl;
    }
    
    cipherPvals.encryptSData(encSData, encYSData, encSXData, yData, xData, sData, factorDim, sampleDim, nsnp, nencsnp,  nslots, L);
    std::cout << "encrypt S" << std::endl;

    {
       MemoryUsage mem = getMemoryUsage();
       cout << "Peak memory = " << mem.vmpeak/1024 << "KB" << std::endl;
       cout << "Curr memory = " << mem.vmrss/1024  << "KB" << std::endl;
    }
    
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
    
    start= chrono::steady_clock::now();
    
    Ciphertext* encCov = new Ciphertext[nterms];
    Ciphertext* encYS = new Ciphertext[nencsnp];
    Ciphertext* encS = new Ciphertext[nencsnp];
    
 
    //! 10s
    TP_EXEC_RANGE(nterms, first, last);
    for(long k = first; k < last; ++k){
        encCov[k] = enccovData[0][k];
        for(long i = 1; i < sampleDim; ++i){
            scheme.addAndEqual(encCov[k], enccovData[i][k]);
        }
    }
    TP_EXEC_RANGE_END;

    for(long i = 0; i < sampleDim; ++i){
       delete [] enccovData[i];
       // for(long k = 0; k < factorDim; ++k){
       //    encSXData[i][k] = new Ciphertext[nencsnp];
       // }
    }
    delete [] enccovData;

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

    for(long i = 0; i < sampleDim; ++i){
       delete [] encSData[i];
       // for(long k = 0; k < factorDim; ++k){
       //    encSXData[i][k] = new Ciphertext[nencsnp];
       // }
    }
    delete [] encSData;

    
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

    for(long i = 0; i < sampleDim; ++i){
       delete [] encYXData[i];
       delete [] encYSData[i];
       for(long k = 0; k < factorDim; ++k){
          delete [] encSXData[i][k];
       }
       delete [] encSXData[i];
    }
    delete [] encYXData;
    delete [] encYSData;
 
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
    //cipherPvals.HesInverse(encDet, encAdj, enccovData, dim, nslots, L);
    cipherPvals.encAdjoint(encDet, encAdj, encCov);
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "2. (X^T * X)^-1 = " << timeElapsed << " s" << endl;
    totalEvaltime += timeElapsed;
    
#if 1
    // "+------------------------------------+"
    //! 3. Norm
    // "+------------------------------------+"
    //! encwtYX = (y^T * X) * adj(X^T * X) * (X^T * y)
    //! Ynorm = |A| <Y,Y> - wtYX = |A| * sampleDim - wtYX
    
    
    start= chrono::steady_clock::now();
    
    Ciphertext encWtYX;
    cipherPvals.SqrQuadForm(encWtYX, encYX, encAdj, factorDim);  // L - 3, 12.7
    
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
    
#if 1
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
        
        //! encYSnorm = encDet * (encYS[j] * scalefactor) - encWtYS
        encYSnorm[j] = scheme.multByConst(encYS[j], (double) scalefactor);
        scheme.reScaleByAndEqual(encYSnorm[j], 1);              //! L - 1
        
        scheme.modDownToAndEqual(encYSnorm[j], encDet.l);       //! L - 3
        scheme.multAndEqual(encYSnorm[j], encDet);
        scheme.reScaleByAndEqual(encYSnorm[j], 1);              //! L - 4
        
        cipherPvals.QuadForm(encWtYS[j], encYX, encAdj, encSX1, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtYS[j], encYSnorm[j].l);   //! L - 4
        scheme.subAndEqual(encYSnorm[j], encWtYS[j]);
    
        //! encSnorm
        encSnorm[j] = scheme.multByConst(encS[j], (double) scalefactor);
        scheme.reScaleByAndEqual(encSnorm[j], 1);               //! L - 1
        
        scheme.modDownToAndEqual(encSnorm[j], encDet.l);        //! L - 2
        scheme.multAndEqual(encSnorm[j], encDet);
        scheme.reScaleByAndEqual(encSnorm[j], 1);               //! L - 3
        
        cipherPvals.SqrQuadForm(encWtSS[j], encSX1, encAdj, factorDim); //! L - 3
        scheme.modDownToAndEqual(encWtSS[j], encSnorm[j].l);    //! L - 4
        scheme.subAndEqual(encSnorm[j], encWtSS[j]);

        delete [] encSX1;
    }
    TP_EXEC_RANGE_END;
    
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "3. Norm = " << timeElapsed << " s" << endl;
    cout << "logq: " << encYSnorm[0].l  << "," << encSnorm[0].l << endl;   // L - 3
    totalEvaltime += timeElapsed;
    cout << "Total Evaluation Timing = " << totalEvaltime << " s" << endl;

    delete [] encWtYS;
    delete [] encWtSS;
    delete [] encAdj;
    for(long i = 0; i < factorDim; i++) {
       delete [] encSX[i];
    }
    delete [] encSX;

    delete [] encYX;
 

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

    delete [] encYSnorm;
    delete [] encSnorm;
#endif
#endif
}
