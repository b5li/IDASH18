/*
 * @file       BasicTest.h, header file
 * @brief      defining basic functions
 *
 * @author     Miran Kim
 * @date       Aug. 7, 2018
 * @copyright  GNU Pub License
 */

#ifndef HEML_BASICTEST_H_
#define HEML_BASICTEST_H_

#include <complex>

#include <iostream>
#include <vector>
#include <string>

#include "NTL/RR.h"
#include "NTL/mat_RR.h"
#include "NTL/vec_RR.h"

using namespace std;
using namespace NTL;

void calculateAdjoint(double*& res, double& det, double* x, long scalefactor);


double trueIP(double* a, double* b, long size);

void initialWDataVDataZero(double* wData, double* vData, long factorDim);

void calculateQuadForm(RR& res, vec_RR u, Mat<RR> Data, vec_RR v);

RR calculateErr(Mat<RR> mat1, Mat<RR> mat2);

void zScoreError(double& TP, double& FP, double& FN, double& TN, double* zScore1, double* zScore2, long len, double criticalval);

void pvalsError(double& TP, double& FP, double& FN, double& TN, double* pvals1, double* pvals2, long len, double siglevel);

double pnorm(double x);

void printRvector(vec_RR& vec, const long k = 0);

void printRmatrix(Mat<RR>& mat, const long k = 0);

#endif /* BASICTEST_H_ */
