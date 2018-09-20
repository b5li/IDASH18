#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#define USE_NTL 1
#ifdef USE_NTL
#include <NTL/BasicThreadPool.h>
#endif
#include "threadpool.h"

#include "Database.h"
#include "BasicTest.h"
#include "TestLinRegPvals.h"
#include "TestHELinRegPvals.h"

using namespace std;
using namespace NTL;

/* make testlinreg */

int main(int argc, char **argv) {
    int num_threads = (argc>3?atoi(argv[3]):4);
#ifdef USE_NTL
    SetNumThreads(num_threads);
#else
    IDASH::initThreadPool(4);
#endif
	string covariate_filename(argv[1]);     // path to covariate file
    string snp_filename(argv[2]);           // path to snp file
    
    double siglevel = 0.01;
    
    cout << endl;
    cout << "+------------------------------------+" << endl;
    cout << "|    0. Read the covariates & SNP    |" << endl;
    cout << "+------------------------------------+" << endl;
    
 
    long sampleDim = 0, ncols = 1;
    vector<string> tag;
    vector<vector<string>> covfile;
    char covfile_split_char = ',';
    DataFromFile(tag, covfile, covariate_filename, ncols, sampleDim, covfile_split_char);
    long factorDim = ncols - 1;   //! number(covariates + outcome), original Data has id numbers

    long snp_sampleDim = 0, nsnp = 0;
    vector<string> snptag;
    vector<vector<string>> snpfile;
    char snpfile_split_char = 0x20;  // SP
    DataFromFile(snptag, snpfile, snp_filename, nsnp, snp_sampleDim, snpfile_split_char);
    
    cout << "(sampleDim, factorDim, nsnp): (" << sampleDim << "," << factorDim << "," << nsnp << ")" << endl;
    
    //! Normalize the data
    double*  yData = new double[sampleDim];
    double** xData = new double*[sampleDim];   //! Xmat: covariate (Xmat[i][0]: y[i], 1<= i <<= n)
    double** sData = new double*[sampleDim];   //! sData: snp
    for(long i = 0; i < sampleDim; ++i){
        xData[i] = new double[factorDim];
        sData[i] = new double[nsnp];
    }
    ListzData(yData, xData, factorDim, sampleDim, tag, covfile);
    ListsData(sData, nsnp, sampleDim,  snpfile);

    double ** X0 = new double*[sampleDim];
    for(long i = 0; i < sampleDim; i++) {
       X0[i] = new double[factorDim-1];
       for(long j = 1; j < factorDim; j++) {
          X0[i][j-1] = xData[i][j];
       }
    }
    normalizeData(xData, xData, factorDim, sampleDim);
    normalizeData(X0, X0, factorDim-1, sampleDim);

    double ** X = new double*[sampleDim];
    for(long i = 0; i < sampleDim; i++) {
       X[i] = new double[factorDim];
       X[i][0] = 1;
       for(long j = 0; j < factorDim-1; j++) {
          X[i][j+1] = X0[i][j];
       }
    }

    //! Evaluation
    double* zScore_ct;
    double* pVals_ct;
 
    TestHEPvals::testTrivialHELinReg(zScore_ct, pVals_ct, yData, X, sData, factorDim, sampleDim, nsnp, snptag, "LinRegResult/TrivialHELinReg_Pvals.txt");

    cout << "+------------------------------------+" << endl;
    cout << "|            Quality Check           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    //! origianl semi-parallel logistic regression (plaintext)
    double* zScore_pt;
    double* pVals_pt;
    SimpleDataFromFile(pVals_pt, "LinRegResult/Plain_Pvals.txt");
    
    double TP, FP, FN, TN;
    pvalsError(TP, FP, FN, TN, pVals_ct, pVals_pt, nsnp, siglevel);
    
    cout << "Error (TP, FP, FN, TN) of pt/ct : " << TP << "," <<  FP << ","  << FN << "," << TN << endl;
    
    
    //MemoryUsage mem = getMemoryUsage();
    //cout << "Peak memory = " << mem.vmpeak/1024 << "KB" << std::endl;
    //cout << "Curr memory = " << mem.vmrss/1024  << "KB" << std::endl;

    // vec_RR zScore;
    // TestPvals::testLinReg(zScore,yOrig,xData,sData,factorDim,sampleDim,nsnp,"plain.txt");
    // double * zScore = 0;
    // TestPvals::testLinRegFast(zScore,yOrig,xData,sData,factorDim,sampleDim,nsnp,10,8,"plain.txt");

 
	return 0;
}
