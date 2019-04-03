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
    
	string covariate_filename(argv[1]);     //! path to covariate file
    string snp_filename(argv[2]);           //! path to snp file
  
    double siglevel = 0.01;
    long pretrain_sigdeg = 3;
    long update_sigdeg = 3;
    double gammaUp = 1.0;
    double gammaDown = 1.0;
    
    
    double* pvals_semipar;
    SimpleDataFromFile(pvals_semipar, "LogRegResult/plain_pvals_semipar.txt");
    
    double* pvals_gold;
    SimpleDataFromFile(pvals_gold, "LogRegResult/plain_pvals_gold.txt");
    
    cout << endl;
    cout << "+------------------------------------+" << endl;
    cout << "|    1. Read the covariates & SNP    |" << endl;
    cout << "+------------------------------------+" << endl;
    
    long sampleDim = 0, ncols = 1;
    char covfile_split_char = ',';
    vector<string> tag;
    vector<vector<string>> covfile;
    DataFromFile(tag, covfile, covariate_filename, ncols, sampleDim, covfile_split_char);
    long factorDim = ncols - 1;   //! number(covariates + outcome), original Data has id numbers
    
    long snp_sampleDim = 0, nsnp = 0;
    char snpfile_split_char = 0x20;  // SP
    vector<string> snptag;
    vector<vector<string>> snpfile;
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
    normalizeData_Colwise(xData, xData, factorDim, sampleDim);
    
    //! y ~ glm(X) <-> (2y-1) ~ glm(X)
    for(long i = 0; i < sampleDim; ++i){
        yData[i] = (2 * yData[i] - 1);
    }
    
    
    double* zscore = new double[nsnp];
    double* pvals = new double[nsnp];

    string zfilename = "LogRegResult/plain_result/plain_zscores_GDiter";
    string pfilename = "LogRegResult/plain_result/plain_pvals_GDiter";
    string txtext = ".txt";
    
    double* TP_semipar = new double[4]; double* TP_gold = new double[4];
    double* FP_semipar = new double[4]; double* FP_gold = new double[4];
    double* FN_semipar = new double[4]; double* FN_gold = new double[4];
    double* TN_semipar = new double[4]; double* TN_gold = new double[4];
    
    for(long numGDIter = 1; numGDIter < 4; ++numGDIter){
        string zfileres = zfilename + std::to_string(numGDIter) + txtext;
        string pfileres = pfilename + std::to_string(numGDIter) + txtext;
        
        TestLRPvals::testLogRegGD(zscore, pvals, yData, xData, sData, factorDim, sampleDim, nsnp, pretrain_sigdeg, update_sigdeg, numGDIter, gammaUp, gammaDown, zfileres, pfileres);
        
        //! Quality Check
        long i = numGDIter - 1;
        pvalsError(TP_semipar[i], FP_semipar[i], FN_semipar[i], TN_semipar[i], pvals, pvals_semipar, nsnp, siglevel);
        pvalsError(TP_gold[i], FP_gold[i], FN_gold[i], TN_gold[i], pvals, pvals_gold, nsnp, siglevel);
    }
    
    cout << "+------------------------------------+" << endl;
    cout << "| 4. Quality check  (TP, FP, FN, TN) |" << endl;
    cout << "+------------------------------------+" << endl;
    
    for(long i = 0; i < 3; ++i){
        cout << "#(GD iter) = " <<  i+1 << endl;
        cout << "vs pt(semi-parallel) = (" << TP_semipar[i] << "," <<  FP_semipar[i] << ","  << FN_semipar[i] << "," << TN_semipar[i] << ")" << endl;
        cout << "vs pt(gold) = " << TP_gold[i] << "," <<  FP_gold[i] << ","  << FN_gold[i] << "," << TN_gold[i] << ")" << endl;
        cout << "=================================" << endl;
    }
    
//     long numNTIter = 5;
//    string zfilename("LogRegResult/plain_zscores_NT.txt");
//    string pfilename("LogRegResult/plain_pvals_NT.txt");
//    TestLRPvals::testLogRegNT(zscore, pvals, yData, xData, sData, factorDim, sampleDim, nsnp, numNTIter, gammaUp, zfilename, pfilename);

   

	return 0;
}
