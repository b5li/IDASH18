/*
 * @file       TestPvalues.h, header file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */

#ifndef HEML_TESTPVALUES_H_
#define HEML_TESTPVALUES_H_

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
 
//! Testing functions
class TestPvals
{
public:
    
    static void trueNLGDiteration(double** zData, double*& wData, double*& vData, long factorDim, long sampleDim, double gamma, double eta);
    static void calculateAUC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& AUC);
    

    //! Linear Reg
    static void testLinReg(vec_RR& zScore,  double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename);
    static void testLinRegFast(double*& zScore, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long ninvstep, long ninvBits,   string filename);
    
    static void testLinRegGD(vec_RR& zScore, double** xData, double** sData, long factorDim, long sampleDim, long nsnp,long numIter, double gamma, string filename);
    
    /*
    //! Logistic Reg
    static void testLogRegGD(vec_RR& zScore, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long numIter, double gammaUp, double gammaDown, string filename);
    static void testLogRegNT(vec_RR& zScore, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long numIter, double gammaUp, string filename);
    
    static void updateIRLS(vec_RR& zScore, double* wvec, double* zvec, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename);
    static void updateIRLSFast(vec_RR& zScore, double* wvec, double* zvec, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename);
    */
    
    
};


#endif /* TESTPVALUES_H_ */
