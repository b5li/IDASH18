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

void calculateAdjoint(double*& res, double& det, double* x, long scalefactor){
    // x[0] x[1] x[2] x[3]
    //      x[4] x[5] x[6]
    //           x[7] x[8]
    //                x[10]
    
    res = new double[10];
    
    double sqrtable[4] = {3,5,6,8};   // {{0,3},{1,2},{1,3},{2,3},
    double table[20][2] = {
        {1, 5},{1, 7},{1, 8},{1, 9},
        {2, 4},{2, 5},{2, 6},{2, 9},
        {3, 5},{3, 6},{3, 8},
        {4, 7},{4, 8},{4, 9},
        {5, 6},{5, 8},{5, 9},
        {6, 7},{6, 8},{7, 9},
    };
/*
    double table1[20][2][2] = {
        {{0,1}, {1,2}},  {{0,1}, {2,2}}, {{0,1}, {2,3}},  {{0,1}, {3,3}},
        {{0,2}, {1,1}},  {{0,2}, {1,2}}, {{0,2}, {1,3}},  {{0,2}, {3,3}},
        {{0,3}, {1,2}},  {{0,3}, {1,3}}, {{0,3}, {2,3}},
        {{1,1}, {2,2}},  {{1,1}, {2,3}}, {{1,1}, {3,3}},
        {{1,2}, {1,3}},  {{1,2}, {2,3}}, {{1,2}, {3,3}},
        {{1,3}, {2,2}},  {{1,3}, {2,3}}, {{2,2}, {3,3}},
    };
*/
    
    double* sqrtemp = new double[4];
    double* temp = new double[20];  //! precomputed value
    
    for(long i = 0; i < 4; ++i){
        long j0= sqrtable[i];
        sqrtemp[i] = pow(x[j0], 2);
    }
    for(long i = 0; i < 20; ++i){
        long j0= table[i][0];
        long j1= table[i][1];
        temp[i] = x[j0] * x[j1];
    }
    
    
    double* adj = new double[30];
    
    //! 0
    adj[0] = temp[19] - sqrtemp[3];
    adj[1] = temp[18] - temp[16];
    adj[2] = sqrtemp[2];
    
    //! 1
    adj[3] = sqrtemp[3]- temp[19];
    adj[4] = - adj[1];
    adj[5] = temp[17] - temp[15];
    
    //! 2
    adj[6] = adj[4];
    adj[7] = sqrtemp[2] - temp[13];
    adj[8] = temp[12] - temp[14];
    
    //! 3
    adj[9]  = adj[5];
    adj[10] = adj[8];
    adj[11] = sqrtemp[1] - temp[11];
    
    //! 4
    adj[12] = - adj[3];
    adj[13] = temp[10] - temp[7];
    adj[14] = sqrtemp[0];
    
    //! 5
    adj[15] = adj[1];
    adj[16] = temp[3] - temp[9];
    adj[17] = temp[8] - temp[2];
    
    //! 6
    adj[18] = -adj[5];
    adj[19] = temp[6] - temp[2];
    adj[20] = temp[1] - temp[5];

    //! 7
    adj[21] = - adj[7];
    adj[22] = - adj[16];
    adj[23] = sqrtemp[0];
    
    //! 8
    adj[24] = - adj[8];
    adj[25] = - adj[19];
    adj[26] = temp[4] - temp[0];
    
    //! 9
    adj[27] = temp[11] - sqrtemp[1];
    adj[28] = -adj[20];
    adj[29] = temp[4];
    
    
    res[0] = x[4] * adj[0] +  x[5] * (adj[1] + temp[18]) - x[7] * adj[2];
    res[4] = x[0] * adj[12] + x[2] * (adj[13] + temp[10]) - x[7] * adj[14];
    res[7] = x[0] * adj[21] + x[1] * (adj[22] + temp[9]) - x[4] * adj[23];
    res[9] = x[0] * adj[27] + x[1] * (adj[28] + temp[5]) - x[2] * adj[29];
    
    res[1] = x[1] * adj[3] + x[2] * adj[4] + x[3] * adj[5];
    res[2] = x[1] * adj[6] + x[2] * adj[7] + x[3] * adj[8];
    res[3] = x[1] * adj[9] + x[2] * adj[10] + x[3] * adj[11];
    res[5] = x[0] * adj[15] + x[2] * adj[16] + x[3] * adj[17];
    res[6] = x[0] * adj[18] + x[2] * adj[19] + x[3] * adj[20];
    res[8] = x[0] * adj[24] + x[1] * adj[25] + x[3] * adj[26];
    
    
    
    det = (x[0] * res[0] + x[1] * res[1] + x[2] * res[2] + x[3] * res[3]) * scalefactor;
    
    cout << res[0] << "," << res[1]  << "," << res[2] << "," << res[3] << endl;
    cout << res[4] << "," << res[5]  << "," << res[6] << endl;
    cout << res[7] << "," << res[8] << endl;
    cout << res[9] << endl;
    
    cout << "det: " << det << endl;
}

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


