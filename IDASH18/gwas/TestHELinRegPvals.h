/*
 * @file       TestHEPvalues.h, header file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */

#ifndef HEML_TESTHEPVALUES_H_
#define HEML_TESTHEPVALUES_H_

#include <complex>

#include <iostream>
#include <vector>
#include <string>

#include "../src/Scheme.h"
#include "../src/SecretKey.h"

 
using namespace std;
using namespace NTL;


//! Testing functions
class TestHEPvals
{
public:

    //! later applied to new KS (Decomp)
    static void testFastHELinReg(double*& zScore, double*& pVals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp,  string filename);
    
    
    
};


#endif /* TESTPVALUES_H_ */
