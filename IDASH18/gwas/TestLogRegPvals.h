/*
 * @file       TestLogRegPvalues.h, header file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */

#ifndef TESTLRPVALUES_H_
#define TESTLRPVALUES_H_

#include <complex>
#include "Database.h"

#include <iostream>
#include <vector>
#include <string>

#include "NTL/RR.h"
#include "NTL/mat_RR.h"
#include "NTL/vec_RR.h"

using namespace std;
using namespace NTL;

static double sigmoid3[3] = {0.5,0.15012,-0.001593};
static double sigmoid5[4] = {0.5,0.19131,-0.0045963, 0.0000412332};   // sometimes it can be zero near 4
static double sigmoid7[5] = {0.5,0.216884,-0.00819276,0.000165861,-0.00000119581};

//! Testing functions
class TestLRPvals
{
public:
    
    static void trueNLGDiteration(double** zData, double*& wData, double*& vData, long factorDim, long sampleDim, double gamma, double eta, long kdeg);
    static void trueNewtoniteration(double** xData, double* yData, vec_RR& twData1, long factorDim, long sampleDim, double gamma);
    
    static void calculatePvals(double*& ip, double*& prob, double*& wvec, double*& zvec, double* wData, double* yData, double** zData, long kdeg, long factorDim, long sampleDim, double gamma); 
    static void calculateAUC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& AUC);
    
    
    /********************************************************************/
    //! Logistic Reg
    static void testLogRegGD(double*& zScore, double*& pvals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long trainsigdeg, long testsigdeg, long numIter, double gammaUp, double gammaDown, string filename);
    static void testLogRegNT(double*& zScore, double*& pvals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long numIter, double gammaUp, string filename);
    
    static void updateIRLS(vec_RR& zScore, double* wvec, double* zvec, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename);
    static void updateIRLSFast(double*& zScore, double* wvec, double* zvec, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename);
  
    
    
};


#endif /* TESTPVALUES_H_ */
