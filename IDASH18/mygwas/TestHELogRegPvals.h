/*
 * @file       TestHELogRegPvalues.h, header file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim
 * @date       Aug. 26, 2018
 * @copyright  GNU Pub License
 */

#ifndef HEML_TESTHELRPVALUES_H_
#define HEML_TESTHELRPVALUES_H_

#include <complex>
#include <iostream>
#include <vector>
#include <string>

#include "../src/Scheme.h"
#include "../src/SecretKey.h"

using namespace std;
using namespace NTL;

//! Testing functions: logistic regression
class TestHELRPvals
{
public:

    static void testFastHELogReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, vector<string> snptag, string filename);
    
    static void new_testFastHELogReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, vector<string> snptag, string filename, long numGDIter,  const double gammaUp = 1.0);
    //static void testHELogReg_block8(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp,  string filename);

};


#endif /* TESTPVALUES_H_ */