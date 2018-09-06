/*
 * @file       TestPvalues.h, header file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */

#ifndef HEML_DATABASE_H_
#define HEML_DATABASE_H_

#include <complex>

#include <iostream>
#include <vector>
#include <string>

#include "NTL/RR.h"
#include "NTL/mat_RR.h"
#include "NTL/vec_RR.h"

using namespace std;
using namespace NTL;


void SimpleDataFromFile(double*& Data,  string path);

void DataFromFile(vector<string>& tag, vector<vector<string>>& sline, string path, long& factorDim, long& sampleDim, char split_char);


void ListzData(double*& yData, double**& zData, long& factorDim, long& sampleDim, vector<string> tag, vector<vector<string>> covfile);
void ListsData(double**& sData, long& nsnp, long& sampleDim, vector<vector<string>> snpfile);

void normalizeData(double**& Data, double** zData, long factorDim, long sampleDim);


void printvector(double* vec, const long k = 0);

void printvectorToFile(double* vec, string filename, const long k = 0);

void printRvectorToFile(vec_RR& vec, string filename, const long k = 0);

void printRvectorToFile1(vec_RR& vec, string filename, const long k = 0);

void printRmatrixToFile(Mat<RR>& mat, string filename, const long k = 0);

#define Error "\033[5;31mError:\033[0m"

#endif /* TESTDATABASE_H_ */
