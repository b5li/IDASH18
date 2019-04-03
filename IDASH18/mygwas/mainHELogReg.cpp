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
    long numGDIter = atoi(argv[3]);         //! iteration of NLGD
    double gammaUp = 1.0;   //! learning rate of GD (-> gammaUp/n)
    

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
    
    double* HEzscore = new double[nsnp];
    double* HEpvals = new double[nsnp];
    
    string pfilename = "LogRegResult/HE_pvals_iter";
    string txtext = ".txt";
    string pfileres = pfilename + std::to_string(numGDIter) + txtext;
    
    TestHELRPvals::new_testFastHELogReg(HEzscore, HEpvals, yData, xData, sData, factorDim, sampleDim, nsnp, snptag, pfileres, numGDIter,  gammaUp);    //! 0.00523846
    
    //TestHELRPvals::testFastHELogReg(HEzscore, HEpvals, yData, xData, sData, factorDim, sampleDim, nsnp, snptag, "LogRegResult/HE_Pvals_L23.txt");    //! 0.00523846
    

   
    cout << "+------------------------------------+" << endl;
    cout << "|          4. Quality Check          |" << endl;
    cout << "+------------------------------------+" << endl;
    
    double* pvals_semipar;
    SimpleDataFromFile(pvals_semipar, "LogRegResult/plain_pvals_semipar.txt");
    
    double* pvals_gold;
    SimpleDataFromFile(pvals_gold, "LogRegResult/plain_pvals_gold.txt");
    
    double TP, FP, FN, TN;
    double siglevel = 0.01;      //! significane level
    //pvalsError(TP, FP, FN, TN, HEpvals, pvals_semipar, nsnp, siglevel);
    pvalsError(TP, FP, FN, TN, HEpvals, pvals_semipar, nsnp, siglevel);
    cout << "vs pt(semi-parallel) = (" << TP << "," <<  FP << ","  << FN << "," << TN << ")" << endl;
    pvalsError(TP, FP, FN, TN, HEpvals, pvals_gold, nsnp, siglevel);
    cout << "vs pt(gold) = " << TP << "," <<  FP << ","  << FN << "," << TN << ")" << endl;
 
	return 0;
}
