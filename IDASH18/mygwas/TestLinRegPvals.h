/*
 * @file       TestLinRegPvalues.h, header file
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
    
    //! Linear Reg
    static void testLinReg(vec_RR& zScore,  double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename);
    static void testLinRegFast(double*& zScore, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long ninvstep, long ninvBits,   string filename);
    
    static void testLinRegGD(vec_RR& zScore, double** xData, double** sData, long factorDim, long sampleDim, long nsnp,long numIter, double gamma, string filename);
    
    
    
};


#endif /* TESTPVALUES_H_ */
