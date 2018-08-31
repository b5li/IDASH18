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

#include "sys.h"

using namespace std;
using namespace NTL;

/* make testlinreg */


int main(int argc, char **argv) {
	
#ifdef USE_NTL
    SetNumThreads(4);
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
    normalizeData(xData, xData, factorDim, sampleDim);
 
    //! y ~ glm(X) <-> (2y-1) ~ glm(X)
    for(long i = 0; i < sampleDim; ++i){
        yData[i] = (2 * yData[i] - 1);
    }

    //! Evaluation
    double* zScore_ct;
    double* pVals_ct;
 
    TestHEPvals::testFastHELinReg(zScore_ct, pVals_ct, yData, xData, sData, factorDim, sampleDim, nsnp, "LinRegResult/FastHELinReg_pvals.txt");
    //TestHEPvals::testHESIMDLinReg(zScore_ct, pVals_ct, yData, xData, sData, factorDim, sampleDim, nsnp,  "Result/HELinRegpvals_MR.txt");
    
    //TestHEPvals::testHELinReg(zScore_ct, pVals_ct, yData, xData, sData, factorDim, sampleDim, nsnp,  "Result/HEpvals.txt");
    
 
 
    cout << "+------------------------------------+" << endl;
    cout << "|            Quality Check           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    //! origianl semi-parallel logistic regression (plaintext)
    double* zScore_pt;
    double* pVals_pt;
    SimpleDataFromFile(pVals_pt, "LinRegResult/pvals_parallel_logreg.txt");
    
    double TP, FP, FN, TN;
    pvalsError(TP, FP, FN, TN, pVals_ct, pVals_pt, nsnp, siglevel);
    
    cout << "Error (TP, FP, FN, TN) of pt/ct : " << TP << "," <<  FP << ","  << FN << "," << TN << endl;

    
    cout << "+------------------------------------+" << endl;
    cout << "|            Memory Check            |" << endl;
    cout << "+------------------------------------+" << endl;
    //MemoryUsage mem = getMemoryUsage();
    //cout << "Peak memory = " << mem.vmpeak/1024 << "KB" << std::endl;
    //cout << "Curr memory = " << mem.vmrss/1024  << "KB" << std::endl;
 
	return 0;
}