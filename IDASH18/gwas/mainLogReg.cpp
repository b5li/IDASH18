#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <NTL/BasicThreadPool.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Database.h"
#include "BasicTest.h"

#include "TestLinRegPvals.h"
#include "TestLogRegPvals.h"

#include "TestHELinRegPvals.h"
#include "TestHELogRegPvals.h"

using namespace std;
using namespace NTL;

/* make testlogreg */

int main(int argc, char **argv) {
	
	SetNumThreads(4);
    
	string covariate_filename(argv[1]);     // path to covariate file
    string snp_filename(argv[2]);           // path to snp file
    
    double siglevel = 0.01;
    
    cout << endl;
    cout << "+------------------------------------+" << endl;
    cout << "|    1. Read the covariates & SNP    |" << endl;
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
    
    double* zScore = new double[nsnp];
    double* pVals = new double[nsnp];
    
    long numIter = 4;   //! iteration of NLGD
    long testsigdeg = 7;

    
    TestHELRPvals::testHELogReg(zScore, pVals, yData, xData, sData, factorDim, sampleDim, nsnp, snptag, "LogRegResult/HELogReg_Pvals.txt");
   
    //TestHELRPvals::testHELogReg_accuracy(zScore, pVals, yData, xData, sData, factorDim, sampleDim, nsnp, numIter, testsigdeg, snptag, "LogRegResult/HE_Pvals_deg7.txt");
    //TestLRPvals::testLogRegGD(zScore, pVals, yData, xData, sData, factorDim, sampleDim, nsnp, 3, testsigdeg, numIter, 10,  1, "LogRegResult/Plain_Pvals_deg3.txt"); //! plain comp
    
 
    cout << "+------------------------------------+" << endl;
    cout << "|          2. Quality Check          |" << endl;
    cout << "+------------------------------------+" << endl;
    
    // origianl regression (plaintext)
    double* pVals_pt;
    SimpleDataFromFile(pVals_pt, "LogRegResult/Plain_Pvals.txt");
    
    double TP, FP, FN, TN;
    pvalsError(TP, FP, FN, TN, pVals, pVals_pt, nsnp, siglevel);
    cout << "Error (TP, FP, FN, TN) of pt/ct : " << TP << "," <<  FP << ","  << FN << "," << TN << endl;
 
	return 0;
}
