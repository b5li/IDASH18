/*
 * @file       BasicTest.cpp, cpp file
 * @brief      defining basic functions  
 *
 * @author     Miran Kim
 * @date       Aug. 7, 2018
 * @copyright  GNU Pub License
 */

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>

#include "NTL/ZZX.h"
#include "NTL/RR.h"
#include "NTL/vec_RR.h"
#include "NTL/mat_RR.h"
#include <NTL/BasicThreadPool.h>

#include "Database.h"
#define printout 0

using namespace std;

//!@ Input: vec_RR
//!@ Function: print the vector
//!  If k = 0, then print out all the components of an input vector

void printRvector(vec_RR& vec, long k){
    long len;
    
    if(k==0){
        len= vec.length();
    }
    else{
        len = k;
    }
    
    cout << "[" ;
    for(int i = 0; i< len; ++i){
        cout << vec[i] ;
        if(i < len -1) cout << ", ";
    }
    cout << "]" << endl;
}


//!@ Input: R-matrix
//!@ Function: print the matrix
void printRmatrix(Mat<RR>& mat, const long k){
    long rlen, clen;
    
    if(k==0){
        rlen= mat.NumRows();
        clen= mat.NumCols();
    }
    else{
        rlen = k;
        clen = k;
    }
    
    for(int i = 0; i< rlen; ++i){
        cout << "[";
        for(int j = 0; j < clen; ++j){
            cout << mat[i][j] ;
            if(j < clen-1) cout << "\t";
        }
        cout << "]\n";
    }
}


double trueIP(double* a, double* b, long size) {
    double res = 0.0;
    for(long i = 0; i < size; ++i) {
        res += a[i] * b[i];
    }
    return res;
}

void initialWDataVDataZero(double* wData, double* vData, long factorDim) {
    for (long i = 0; i < factorDim; ++i) {
        wData[i] = 0.0;
        vData[i] = 0.0;
    }
}


//! compute u^T * Data * v  (u,v: given as column vectors)
void calculateQuadForm(RR& res, vec_RR u, Mat<RR> Data, vec_RR v){
    res = to_RR("0");
    
    if((Data.NumRows()!= u.length()) || (Data.NumCols()!= v.length())){
        cout << "Error" << endl;
    }
    for(long i = 0; i < Data.NumRows(); ++i){
        for(long j = 0; j < Data.NumCols(); ++j){
            res += Data[i][j] * u[i] * v[j];
        }
    }
}

//! @ Function:
RR calculateErr(Mat<RR> mat1, Mat<RR> mat2){
    long nrows = mat1.NumRows();
    long ncols = mat1.NumCols();
    
    if((nrows != mat2.NumRows()) || (ncols != mat2.NumCols())){
        cout << "Error: cannot compare matrices" << endl;
    }
    
    RR res = to_RR("0");
    for(long i = 0; i < nrows; ++i){
        for(long j = 0; j < ncols; ++j){
            RR rtemp = abs(mat1[i][j] - mat2[i][j]);
            if(rtemp > res){
                res = rtemp;
            }
        }
    }
    
    return res;
}

void zScoreError(double& TP, double& FP, double& FN, double& TN, double* zScore1, double* zScore2, long len, double criticalval){
    TP = 0;
    FP = 0;
    FN = 0;
    TN = 0;
    long H1 = 0;   // number of reject
    
    for(long i = 0; i < len; ++i){
        //cout << zScore1[i] << "," << zScore2[i] << endl;
        if(abs(zScore1[i]) > criticalval){   // = pvals1[i] <= 0.01
            if(abs(zScore2[i]) > criticalval){
                TP++;
            }
            else{
                FN++;
            }
            H1++;
        }
        else{
            if(abs(zScore2[i]) > criticalval){
                FP++;
            }
            else{
                TN++;
            }
        }
    }
    long Type1 = (FP + TN);
    
    TP = (double)TP/H1;
    FP = (double)FP/Type1;
    FN = (double)FN/H1;
    
    long H0 = len - H1;
    TN = (double)TN/H0;
}


void pvalsError(double& TP, double& FP, double& FN, double& TN, double* pvals1, double* pvals2, long len, double siglevel){
    TP = 0; // 1441
    FP = 0;
    FN = 0;
    TN = 0;
    long H1 = 0;   // number of reject : 1528
    
    for(long i = 0; i < len; ++i){
        //cout << i << ": " << pvals1[i] << "," << pvals2[i] << endl;
        if((pvals1[i] < siglevel)||(pvals1[i] == siglevel)){
            H1++;
            
            if((pvals2[i] < siglevel)||(pvals2[i] == siglevel)){
                TP++;
            }
            else{
                FN++;
            }
        }
        else{
            if((pvals2[i] < siglevel)||(pvals2[i] == siglevel)){
                FP++;
            }
            else{
                TN++;
            }
        }
    }
    //cout << TP << "," << FP << "," << FN << "," << TN << endl;   // 1441,0,87, 9115
    long Type1 = (FP + TN);
    
    TP = (double)TP/H1;
    FP = (double)FP/Type1;
    FN = (double)FN/H1;
    
    long H0 = len - H1;
    TN = (double)TN/H0;
}

//! suppose that x >= 0
double pnorm(double x){
    return erfc(x / sqrt(2));
}


