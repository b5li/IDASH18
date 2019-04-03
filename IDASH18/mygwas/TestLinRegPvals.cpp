/*
 * @file       TestLinRegPvalues.cpp, cpp file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#include "NTL/ZZX.h"
#include <NTL/RR.h>
#include "NTL/vec_RR.h"
#include "NTL/mat_RR.h"
#include <NTL/BasicThreadPool.h>

#include "BasicTest.h"
#include "TestLinRegPvals.h"


void TestPvals::testLinReg(vec_RR& zScore, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename){
    //! NX = I - X * (XT * X)^-1 * X^T
    Mat<RR> Xmat;
    Xmat.SetDims(sampleDim, factorDim);
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < factorDim; ++j){
            Xmat[i][j]  = to_RR(xData[i][j]);
            //Xmat[i][j]  = to_RR(xData[i][j] * xData[i][0]);
        }
    }
    Mat<RR> Imat;
    Imat.SetDims(sampleDim, sampleDim);
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < sampleDim; ++j){
            Imat[i][j] = to_RR("0");
        }
        Imat[i][i] = to_RR("1");
    }
    
    Mat<RR> Xtrans;
    transpose(Xtrans, Xmat);
    
    Mat<RR> Xcov;  //! k * k
    mul(Xcov, Xtrans, Xmat);
    
    Mat<RR> Xcovinv;
    inv(Xcovinv, Xcov);
    
    Mat<RR> Nx;   //! n * n
    mul(Nx, Xmat, Xcovinv);
    mul(Nx, Nx, Xtrans);

    sub(Nx, Imat, Nx);
    //printRmatrix(NX, 10);
    
    //! ystar = Nx * y
    Mat<RR> Ymat;
    Ymat.SetDims(sampleDim, 1);
    for(long i = 0; i < sampleDim; ++i){
        Ymat[i][0]  = to_RR(yData[i]);
        //Ymat[i][0]  = to_RR(xData[i][0]);
    }
    
    Mat<RR> Ystar;
    mul(Ystar, Nx, Ymat);
    
    vec_RR Yvec;
    Yvec.SetLength(sampleDim);
    for(long i = 0; i < sampleDim; ++i){
        Yvec[i] = Ystar[i][0];
    }
    
    //! sstar = Nx * s: n * p
    Mat<RR> Smat, Sstar;
    Smat.SetDims(sampleDim, nsnp);
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < nsnp; ++j){
            Smat[i][j]  = to_RR(sData[i][j]);
        }
    }
    mul(Sstar, Nx, Smat);
    
    //! compute beta, var, and alpha = |beta/sigma|
    vec_RR beta, var, RSS;
    beta.SetLength(nsnp);
    var.SetLength(nsnp);
    zScore.SetLength(nsnp);
    RSS.SetLength(nsnp);
    
    RR avgRSS = to_RR("0");
    
    RR Ynorm;
    InnerProduct(Ynorm, Yvec, Yvec);
   
    for(long j = 0; j < nsnp; ++j){
        vec_RR Svec; //! jth snp
        Svec.SetLength(sampleDim);
        for(long i = 0; i < sampleDim; ++i){
            Svec[i] = Sstar[i][j];
        }
        
        RR YS, Snorm;
        InnerProduct(YS, Yvec, Svec);
        InnerProduct(Snorm, Svec, Svec);
        
        beta[j] = YS/ Snorm;
        var[j] = (Ynorm - YS)/((sampleDim - factorDim - 1)* Snorm);   //! factordim = k + 1
        zScore[j] = (beta[j]/sqrt(var[j]));
        
        RSS[j] = Ynorm - YS;
        avgRSS += RSS[j];
    }
    
    ofstream outf(filename);
    outf.close();
    
    outf.open(filename, fstream::in | fstream::out | fstream::app);   /// open the file
    outf << "beta: " << endl;
    outf.close();
    
    printRvectorToFile(beta, filename);
    
    outf.open(filename, fstream::in | fstream::out | fstream::app);   /// open the file
    outf << "zScore: " << endl;
    outf.close();
    
    printRvectorToFile(zScore, filename);
    
    outf.open(filename, fstream::in | fstream::out | fstream::app);   /// open the file
    outf << "avg RSS: " << avgRSS/nsnp << endl;
    outf.close();
}


void TestPvals::testLinRegFast(double*& zScore, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long ninvstep, long invBits, string filename){
    
    Mat<RR> Xmat;
    Xmat.SetDims(sampleDim, factorDim);
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < factorDim; ++j){
            Xmat[i][j]  = to_RR(xData[i][j] * (yData[i]));
        }
    }
    Mat<RR> Imat;
    Imat.SetDims(sampleDim, sampleDim);
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < sampleDim; ++j){
            Imat[i][j] = to_RR("0");
        }
        Imat[i][i] = to_RR("1");
    }
    
    //! information for n users
    //! Xcov_small[i]: information of ith user
    vector<Mat<RR>> Xcov_small;
    for(long i = 0; i < sampleDim; ++i){
        Mat<RR> temp;
        temp.SetDims(factorDim, factorDim);
        for(long j = 0; j < factorDim; ++j){
            for(long l = 0; l < factorDim; ++l){
                temp[j][l] = Xmat[i][j] * Xmat[i][l];
            }
        }
        Xcov_small.push_back(temp);
    }
    
    Mat<RR> Xcov = Xcov_small[0];  //! k * k
    for(long i = 1; i < sampleDim; ++i){
        add(Xcov, Xcov, Xcov_small[i]);
    }
    //printRmatrix(Xcov);
    
    //------------------------------
    // Inverse
    //-----------------------------
    Mat<RR> Xcovinv;
    
    long invscalar = (1 << invBits);
    Mat<RR> IDmat;
    IDmat.SetDims(factorDim, factorDim);
    for(long i = 0; i < factorDim; ++i){
        for(long j = 0; j < factorDim; ++j){
            IDmat[i][j]= to_RR("0");
        }
        IDmat[i][i]= to_RR("1");
    }

    // Amat = I - X/n
    Mat<RR> Amat;
    mul(Amat, Xcov, 1./invscalar);
    sub(Amat, IDmat, Amat);
    
    Mat<RR> res;
    
    Mat<RR> invmat1;
    add(invmat1, IDmat, Amat);
    mul(res, invmat1, 1./invscalar);
    
    for(long i = 1; i< ninvstep; ++i){
        Amat = sqr(Amat);
        Mat<RR> tmp;
        add(tmp, IDmat, Amat);
        mul(invmat1, invmat1, tmp);
        
        mul(res, invmat1, 1./invscalar);
        //printRmatrix(res, nrows);
        //RR inverror1 = getError(invmat, res, nrows, ncols);
        //cout << "----------------------------------------------------" << endl;
    }
    Xcovinv = res;
    
    if(ninvstep == 0){
        inv(Xcovinv, Xcov);
    }
    
    printRmatrix(Xcovinv);
    
    //! 1. Ynorm <- <ystar, ystar> = <y, y> - (y^T * X) * Xcovinv *  (X^T * y)
    //! YXvec = y^T * X = [-29, -6.506140997, -6.368385672, -12.88024343]
    vec_RR Yvec;
    Yvec.SetLength(sampleDim);
    for(long i = 0; i < sampleDim; ++i){
        Yvec[i] = to_RR(yData[i]);
    }
    RR Ynorm = to_RR(sampleDim);
    
    vec_RR YXvec;
    YXvec.SetLength(factorDim);
    for(long j = 0; j < factorDim; ++j){
        YXvec[j] = to_RR("0");
        for(long i = 0; i < sampleDim; ++i){
            YXvec[j] +=  Xmat[i][j];

        }
    }
    //printRvector(YXvec);
    
    RR temp;
    calculateQuadForm(temp, YXvec, Xcovinv, YXvec);
    sub(Ynorm, Ynorm, temp);
    //cout << "(y^T * X) * Xcovinv *  (X^T * y): " << temp << endl;   //! 12.72937741
    
    vec_RR beta, var, RSS;
    beta.SetLength(nsnp);
    var.SetLength(nsnp);
    RSS.SetLength(nsnp);
    RR avgRSS = to_RR("0");
    zScore = new double[nsnp];
    
    ofstream outf(filename);
    outf.close();
    
    
    for(long j = 0; j < nsnp; ++j){
        vec_RR Svec; //! jth snp
        Svec.SetLength(sampleDim);
        for(long i = 0; i < sampleDim; ++i){
            Svec[i] = to_RR(sData[i][j] * (yData[i]));
        }
        
        //! SXvec = s^T * X
        vec_RR SXvec;
        SXvec.SetLength(factorDim);
        for(long k = 0; k < factorDim; ++k){
            SXvec[k] = to_RR("0");
            for(long i = 0; i < sampleDim; ++i){
                SXvec[k] +=  Svec[i] * Xmat[i][k];
            }
        }
        
        //! YS = <ystar, sstar> = <y, s> - (y^T * X) * Xcovinv *  (X^T * s)
        RR YS = to_RR("0");
        for(long i = 0; i < sampleDim; ++i){
            YS += Svec[i];
        }
        
        calculateQuadForm(temp, YXvec, Xcovinv, SXvec);
        sub(YS, YS, temp);
        
        //! Snorm =  <sstar, sstar> = <s, s> - (s^T * X) * Xcovinv *  (s^T * s)
        RR Snorm;
        InnerProduct(Snorm, Svec, Svec);   //! <s, s>
        calculateQuadForm(temp, SXvec, Xcovinv, SXvec);
        sub(Snorm, Snorm, temp);
        
        //outf.open(filename, fstream::in | fstream::out | fstream::app);
        //outf << j << ": " <<  YS << "," << Snorm << endl;
        //outf.close();
        
        beta[j] = YS/ Snorm;
        var[j] = (Ynorm - YS)/((sampleDim - factorDim - 1)* Snorm);   //! factordim = k + 1
        
        RR rtemp = (beta[j]/sqrt(var[j]));
        conv(zScore[j], rtemp);

        
        RSS[j] = Ynorm - YS;
        avgRSS += RSS[j];
    }
    
    //outf.open(filename, fstream::in | fstream::out | fstream::app);   /// open the file
    //outf << "beta: " << endl;
    //outf.close();
    

    printvectorToFile(zScore, filename, nsnp);
    
    //outf.open(filename, fstream::in | fstream::out | fstream::app);   /// open the file
    //outf << "avg RSS: " << avgRSS/nsnp << endl;
    //outf.close();
}

 

