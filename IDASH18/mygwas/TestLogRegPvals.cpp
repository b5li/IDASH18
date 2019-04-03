/*
 * @file       TestLogRegPvalues.cpp, cpp file
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

#include "Database.h"
#include "BasicTest.h"
#include "TestLogRegPvals.h"



//! ip: inner product of wData and zData
//! prob: sigmoid(ip)
//! wvec = prob * (1-prob)
//! zvec = ip + (z-p)/(prob * (1-prob))

//!@ Function: one update of NLGD
void TestLRPvals::trueNLGDiteration(double*& wData, double*& vData, double** zData, long factorDim, long sampleDim, double gamma, double eta, long kdeg) {
    double* grad = new double[factorDim]();
    
    //! miran: need to initialization
    initialVecZero(grad, factorDim);
    
    double ipmax = 0.0;
    double ipmin = 0.0;
    
    for(long j = 0; j < sampleDim; ++j) {
        double ip = trueIP(vData, zData[j], factorDim);
        
        double tmp; // = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
        switch(kdeg){
            case 3:
                tmp = (sigmoid3[0] + ip * sigmoid3[1] + pow(ip, 3) * sigmoid3[2]);
                break;
            case 5:
                tmp = (sigmoid5[0] + ip * sigmoid5[1] + pow(ip, 3) * sigmoid5[2] + pow(ip, 5) * sigmoid5[3]);
                break;
            case 7:
                tmp = (sigmoid7[0] + ip * sigmoid7[1] + pow(ip, 3) * sigmoid7[2] + pow(ip, 5) * sigmoid7[3] + pow(ip, 7) * sigmoid7[4]);
                break;
            default:
                tmp = -1. / (1. + exp(ip));
        }
        
        for(long i = 0; i < factorDim; ++i) {
            grad[i] += tmp * (double) zData[j][i];
        }
        
        if(j == 0){
            ipmax = ip;
            ipmin = ip;
        }
        
        if(ip > ipmax){
            ipmax = ip;
        }
        if(ip < ipmin){
            ipmin = ip;
        }
    }
    cout << "(IPmin, IPmax) = (" << ipmin << "," << ipmax << ")" << endl;
    
    for (long i = 0; i < factorDim; ++i) {
        double tmpw = vData[i] - gamma * grad[i];
        vData[i] = (1.0 - eta) * tmpw + eta * wData[i];
        wData[i] = tmpw;
    }
    delete[] grad;
}

//!@ Input: yData = {+1/-1} (encoded)
//!@ xcale = 1 when xData <- y * xData, 0 otherwise
void TestLRPvals::calculatePvals(double*& prob, double*& wvec, double*& zvec, double* wData, double* yData, double** xData, long kdeg, long factorDim, long sampleDim,  bool xscale){
    double* ip = new double[sampleDim];
    double ipmax = 0.0;
    double ipmin = 0.0;
    
    for(int i = 0; i < sampleDim; ++i){
        if(xscale){
            ip[i] =  xData[i][0] * trueIP(xData[i], wData, factorDim); // <z,w> = (2y-1) * <z', w>, ip(xi, w) = log (p/ (1-p))
        } else{
            ip[i] =  trueIP(xData[i], wData, factorDim); 
        }
        
        switch(kdeg){
            case 3:
                prob[i] = (sigmoid3[0] + ip[i] * sigmoid3[1] + pow(ip[i], 3) * sigmoid3[2]);
                break;
            case 5:
                prob[i] = (sigmoid5[0] + ip[i] * sigmoid5[1] + pow(ip[i], 3) * sigmoid5[2] + pow(ip[i], 5) * sigmoid5[3]);
                break;
            case 7:
                prob[i] = (sigmoid7[0] + ip[i] * sigmoid7[1] + pow(ip[i], 3) * sigmoid7[2] + pow(ip[i], 5) * sigmoid7[3] + pow(ip[i], 7) * sigmoid7[4]);
                break;
            default:
                prob[i] = 1. / (1. + exp(-ip[i]));
        }
        wvec[i] = prob[i] * (1.0 - prob[i]);
        
        double yData0;
        if(yData[i] == 1){
            yData0 = 1.0;
        } else{
            yData0 = 0.0;
        }
        zvec[i] = ip[i] + (yData0 - prob[i])/(wvec[i]);   // regulate the learning rate (alpha)
        
        //!
        if(i == 0){
            ipmax = ip[i];
            ipmin = ip[i];
        }
        
        if(ip[i] > ipmax){
            ipmax = ip[i];
        }
        if(ip[i] < ipmin){
            ipmin = ip[i];
        }
        
        //cout  << yData0 << "," << ip[i] << "," << prob[i] << "," << wvec[i] << "," << zvec[i] << endl;
    }
     cout << "(IPmin, IPmax) = (" << ipmin << "," << ipmax << ")" << endl;
}

void TestLRPvals::calculateAUC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& auc) {
    long TN = 0, FP = 0;
    
    vector<double> thetaTN;
    vector<double> thetaFP;
    
    for(int i = 0; i < sampleDim; ++i){
        if(zData[i][0] > 0){
            if(trueIP(zData[i], wData, factorDim) < 0) TN++;
            thetaTN.push_back(zData[i][0] * trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        } else{
            if(trueIP(zData[i], wData, factorDim) < 0) FP++;
            thetaFP.push_back(zData[i][0] * trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        }
    }
    
    correctness = 100.0 - (100.0 * (FP + TN) / sampleDim);
    //    cout << "Failure rate: (y = 1) " << TN << "/" << thetaTN.size() << " + (y = 0) " << FP << "/" ;
    //    cout << thetaFP.size() << " = " <<  (100.0 * (FP + TN) / sampleDim) << " %." << endl;
    cout << "Correctness: " << correctness  << " %." << endl;
    
    if(thetaFP.size() == 0 || thetaTN.size() == 0) {
        cout << "n_test_yi = 0 : cannot compute AUC" << endl;
        auc = 0.0;
    } else{
        auc = 0.0;
        for(long i = 0; i < thetaTN.size(); ++i){
            for(long j = 0; j < thetaFP.size(); ++j){
                if(thetaFP[j] <= thetaTN[i]) auc++;
            }
        }
        auc /= thetaTN.size() * thetaFP.size();
        cout << "AUC: " << auc << endl;
    }
}



//! Iterative least squares
void TestLRPvals::updateIRLSFast(double*& zScore, double*& pvals, double* wvec, double* zvec, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string zfilename, string pfilename){
    
    Mat<RR> Xmat; Xmat.SetDims(sampleDim, factorDim);
    Mat<RR> Imat; ident(Imat, sampleDim);
    
    vec_RR Wvec;
    Wvec.SetLength(sampleDim);
    for(long i = 0; i < sampleDim; ++i){
        Wvec[i] = to_RR(wvec[i]);
    }
    
    vec_RR Zvec;
    Zvec.SetLength(sampleDim);
    for(long i = 0; i < sampleDim; ++i){
        Zvec[i] = to_RR(zvec[i]);
    }
    
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < factorDim; ++j){
            Xmat[i][j]  = to_RR(xData[i][j]);
            //Xmat[i][j]  = to_RR(xData[i][j] * xData[i][0]);  // return back to the original data
        }
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
    
    //! Xcov = X^T * W * X = sum_i (w[i] * x[i]^T *  x[i])
    //! Xcov2 = X^T * W2 * X = sum_i (w[i]^2 * x[i]^T *  x[i])
    Mat<RR> Xcov, Xcov2;
    mul(Xcov, Xcov_small[0], Wvec[0]);  //! k * k
    mul(Xcov2, Xcov_small[0], Wvec[0] * Wvec[0]);
    
    Mat<RR> rtemp;
    for(long i = 1; i < sampleDim; ++i){
        mul(rtemp, Xcov_small[i], Wvec[i]);
        add(Xcov, Xcov, rtemp);
        mul(rtemp, Xcov_small[i], Wvec[i] * Wvec[i]);
        add(Xcov2, Xcov2, rtemp);
    }
    
    Mat<RR> Xcovinv;
    inv(Xcovinv, Xcov);
 
#if 1 // defined(__DEBUG_)
    cout << " (X^T * W * X): " << endl;
    printRmatrix(Xcov);
    
    cout << " (X^T * W * X)^-1: " << endl;
    printRmatrix(Xcovinv);
    
    RR det = determinant(Xcov);
    cout << "det = " << (det) << endl;
    cout << "================================" << endl;
    
    cout << " 1/2 * X^T * W * X: " << endl;
    Mat<RR> scaledXcov;
    mul(scaledXcov, Xcov, to_RR(1.0/2));
    printRmatrix(scaledXcov);
    
    cout << " (1/2 X^T * W * X)^-1: " << endl;
    Mat<RR> scaledXcovinv;
    inv(scaledXcovinv, scaledXcov);
    printRmatrix(scaledXcovinv);
    
  
    RR det1 = determinant(scaledXcov);
    cout << "4 * det = " << (to_RR("2") * det1) << endl;
    Mat<RR> scaledXcovadj;
    mul(scaledXcovadj, scaledXcovinv, det1);
    cout << "adj = " << endl;
    printRmatrix(scaledXcovadj);
#endif
    
    //! Xcovinv2 = (X^T W X) ^-1 * (X^T W^2 X) * (X^T W X) ^-1
    Mat<RR> Xcovinv2;
    mul(Xcovinv2, Xcovinv, Xcov2);
    mul(Xcovinv2, Xcovinv2, Xcovinv);
    
    //---------------------------------------------------
    vec_RR beta, var, RSS, err;
    beta.SetLength(nsnp);
    var.SetLength(nsnp);
    err.SetLength(nsnp);

    vec_RR Snorm;
    vec_RR ZSnorm;
    Snorm.SetLength(nsnp);
    ZSnorm.SetLength(nsnp);
    

    //! ZXvec = z^T * X,
    //cout << "Z^T * X: " << endl;
    vec_RR ZXvec;
    ZXvec.SetLength(factorDim);
    for(long k = 0; k < factorDim; ++k){
        ZXvec[k] = to_RR("0");
        for(long i = 0; i < sampleDim; ++i){
            ZXvec[k] +=  Zvec[i] * Xmat[i][k];
        }
        //cout << ZXvec[k] << endl;
    }
   
    //! ZWXvec = z^T * W * X
    vec_RR ZWXvec;
    ZWXvec.SetLength(factorDim);
    cout << "Zt * W * X:" ;
    for(long k = 0; k < factorDim; ++k){
        ZWXvec[k] = to_RR("0");
        for(long i = 0; i < sampleDim; ++i){
            ZWXvec[k] +=  (Wvec[i] * Zvec[i]) * Xmat[i][k];
        }
        cout << ZWXvec[k] << "," ;
    }
    cout << endl;
    
    
    //! 1. sum wi * (s*[i][j])^2 = s^T * (NXN) * s^T
    //! 2. sum wi * zi * s*[i][j]
    double** SX = new double*[nsnp];
    double** SWX = new double*[nsnp];
    double** SW2X = new double*[nsnp];
    double* SWS = new double[nsnp];
    double* ZWS = new double[nsnp];
    
    for(long j = 0; j < nsnp; ++j){
        vec_RR Svec; //! jth snp
        Svec.SetLength(sampleDim);
        for(long i = 0; i < sampleDim; ++i){
            Svec[i] = to_RR(sData[i][j]);   // s' = (2y-1) * s
        }
        
        //! SXvec = s^T * X
        vec_RR SXvec;
        SXvec.SetLength(factorDim);
        SX[j] = new double[factorDim];
        for(long k = 0; k < factorDim; ++k){
            SXvec[k] = to_RR("0");
            for(long i = 0; i < sampleDim; ++i){
                SXvec[k] +=  Svec[i] * Xmat[i][k];
            }
            conv(SX[j][k], SXvec[k]);
        }

        //! SWXvec = s^T * W * X
        vec_RR SWXvec;
        SWXvec.SetLength(factorDim);
        SWX[j] = new double[factorDim];
        for(long k = 0; k < factorDim; ++k){
            SWXvec[k] = to_RR("0");
            for(long i = 0; i < sampleDim; ++i){
                SWXvec[k] +=  (Wvec[i] * Svec[i]) * Xmat[i][k];
            }
            conv(SWX[j][k], SWXvec[k]);
            
        }
        
        
    
        //! SW2Xvec = s^T * W2 * X
        vec_RR SW2Xvec;
        SW2Xvec.SetLength(factorDim);
        SW2X[j] = new double[factorDim];
        for(long k = 0; k < factorDim; ++k){
            SW2Xvec[k] = to_RR("0");
            for(long i = 0; i < sampleDim; ++i){
                SW2Xvec[k] +=  ((Wvec[i] * Wvec[i]) * Svec[i]) * Xmat[i][k];
            }
            conv(SW2X[j][k], SW2Xvec[k]);
        }
        
        //! temp1 = S^T * W * S
        RR temp1 = to_RR("0");
        for(long i = 0; i < sampleDim; ++i){
            temp1 += Wvec[i] * Svec[i] * Svec[i];
        }
        conv(SWS[j], temp1);
    
        
        RR temp2;
        calculateQuadForm(temp2, SXvec, Xcovinv, SW2Xvec);
        
        RR temp3;
        calculateQuadForm(temp3, SWXvec, Xcovinv, SWXvec);
        
        RR temp4;
        calculateQuadForm(temp4, SXvec, Xcovinv2, SWXvec);
        
        //Snorm[j] = (temp1 - temp2);
        Snorm[j] = (temp1 - temp2 - temp3 + temp4);  // S^T * W * S
        //cout << j << ": " << temp2 << "," << temp3 <<  "," << temp4 << endl;
        
        var[j] = inv(Snorm[j]);
        
        //----------------------------
        // temp1 = Z^T * W * S
        temp1 = to_RR("0");
        for(long i = 0; i < sampleDim; ++i){
            temp1 += Wvec[i] * Zvec[i] * Svec[i];
        }
        conv(ZWS[j], temp1);
        
        calculateQuadForm(temp2, ZXvec, Xcovinv, SW2Xvec);
        calculateQuadForm(temp3, ZWXvec, Xcovinv, SWXvec);
        calculateQuadForm(temp4, ZXvec, Xcovinv2, SWXvec);
        
        //ZSnorm[j] = (temp1 - temp2);
        ZSnorm[j] = (temp1 - temp2 - temp3 + temp4);   // Z^T * W * S
        //cout << j << ": " << temp2 << "," << temp3 <<  "," << temp4 << endl;
        
        beta[j] = (ZSnorm[j]/ Snorm[j]);
        err[j] = sqrt(var[j]);
        temp1  = abs(beta[j]/err[j]);
        conv(zScore[j], temp1);
        pvals[j] = pnorm(zScore[j]);
    }
    
    printvectorToFile(zScore, zfilename, nsnp);
    printvectorToFile(pvals, pfilename, nsnp);
    
    
//    printRvectorToFile1(ZSnorm, "LogRegResult/Plain_ZSnorm.txt");
//    printRvectorToFile1(Snorm, "LogRegResult/Plain_Snorm.txt");
//    printRvectorToFile1(beta, "LogRegResult/pt_beta.txt");
//    printRvectorToFile1(err, "LogRegResult/pt_err.txt");

//    cout << "Zvec: ";
//    printRvector(Zvec, 10);
//
//    cout << "Z^T * W * S : [ " ;
//    printvector(ZWS, 10);
//
//    cout << "S^T * W * S : [ " ;
//    printvector(SWS, 10);
//
//    cout << "X^T * W2 * S : [ " ;
//    for(long i = 0; i < 10; ++i) cout << SW2X[i][0] << "," ;
//    cout << "]" << endl;
  
}

/********************************************************************/
//! Covariates: Logistic regressio using GD approach
//! Input: yData (+1/-1), xData (cov), sData (snp),
//!        pretrain_sigdeg (degree for sigmoid approximation) , numGDIter (iteration number of GD), gammaUp/gammaDown
void TestLRPvals::testLogRegGD(double*& zScore, double*& pvals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long pretrain_sigdeg, long update_sigdeg, long numGDIter,  double gammaUp, double gammaDown, string zfilename, string pfilename){
    
    cout << "+------------------------------------+" << endl;
    cout << "|       2. Pre-training / cov        |" << endl;
    cout << "+------------------------------------+" << endl;
    
    double* wData = new double[factorDim];  //! NGD
    double* vData = new double[factorDim];
    
    //! Initial step for NGD
    initialWDataVDataZero(wData, vData, factorDim);   //! Set as zero
    
    //! Iterative step for NGD
    double alpha0 = 0.01;
    double alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
    
    double** yxData = new double*[sampleDim];
    for(long i = 0; i < sampleDim; ++i){
        yxData[i] = new double[factorDim];
        for(long j = 0; j < factorDim; ++j){
            yxData[i][j]  = xData[i][j] * yData[i];
        }
    }
    
    double neweta[5] = {0.989901, 0 , -0.281783, -0.434061, -0.531076};
    
    for (long iter = 0; iter < numGDIter; ++iter) {
        double eta = (1 - alpha0) / alpha1;
        double gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDim : gammaUp / (iter - gammaDown) / sampleDim;
        trueNLGDiteration(wData, vData, yxData, factorDim, sampleDim, gamma, eta, pretrain_sigdeg);
        cout << iter << ": (gamma = " << gamma << " : w = ["  ;
        cout << wData[0] << "," << wData[1] << "," <<  wData[2] << "," << wData[3] << "], v = [" ;
        cout << vData[0] << "," << vData[1] << "," <<  vData[2] << "," << vData[3] << "]" << endl;
//        cout << "----------------------------------------------------------" << endl;
        alpha0 = alpha1;
        alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
    }
    
    double truecor, trueauc;
    calculateAUC(xData, wData, factorDim, sampleDim, truecor, trueauc);
    cout << "Logistic regression with covariates (using NLGD) / niter = " << numGDIter << ", deg(sig) = " << pretrain_sigdeg << " / " ;
    cout << ", Correctness: " << truecor << " %, AUC: " << trueauc << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|         3. Training / snps         |" << endl;
    cout << "+------------------------------------+" << endl;
    
    
//    double* wData1 = new double[factorDim + 1];  //! NGD
//    double* vData1 = new double[factorDim + 1];
//
//    for(long i = 0; i < factorDim; ++i){
//        wData1[i] = wData[i];
//        vData1[i] = vData[i];
//    }
//    wData1[factorDim] = 0.0;
//    vData1[factorDim] = 0.0;
//
//    double** yxData1 = new double*[sampleDim];
//    for(long i = 0; i < sampleDim; ++i){
//        yxData1[i] = new double[factorDim + 1];
//        for(long j = 0; j < factorDim; ++j){
//            yxData1[i][j]  = yxData[i][j];
//        }
//        yxData1[i][factorDim] = sData[i][0] * yData[i]; //! 0th column
//    }
//
//    for (long iter = 0; iter < 20; ++iter) {
//        double eta = (1 - alpha0) / alpha1;
//        double gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDim : gammaUp / (iter - gammaDown) / sampleDim;
//        trueNLGDiteration(wData1, vData1, yxData1, factorDim + 1, sampleDim, gamma, eta, pretrain_sigdeg);
//        cout << iter << ": (eta = " << eta << " : w = ["  ;
//        cout << wData1[0] << "," << wData1[1] << "," <<  wData1[2] << "," << wData1[3] << "], v = [" ;
//        cout << vData1[0] << "," << vData1[1] << "," <<  vData1[2] << "," << vData1[3] << "]" << endl;
//        cout << "beta[0] = : "<<  wData1[4] << endl;
//        cout << "----------------------------------------------------------" << endl;
//        alpha0 = alpha1;
//        alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
//    }

    
    double* prob = new double[sampleDim];   //! sigmoid(ip), ip = inner product of x and w
    double* wvec = new double[sampleDim];   //! p * (1-p)
    double* zvec = new double[sampleDim];   //! working variable for the next step (Newton)

    //! update the values
    calculatePvals(prob, wvec, zvec, wData, yData, yxData, update_sigdeg, factorDim, sampleDim);  // sigmoid -> deg5

    cout << "wData: "  ;
    printvector(wData, factorDim);

    cout << "p[i]: " ;
    printvector(prob, 10);

    cout << "wvec: " ;
    printvector(wvec, 10);

    cout << "zvec: " ;
    printvector(zvec, 10);

    //updateIRLS(zScore, wvec, zvec, xData, sData, factorDim, sampleDim, nsnp, filename);
    updateIRLSFast(zScore, pvals, wvec, zvec, yData, xData, sData, factorDim, sampleDim, nsnp, zfilename, pfilename);

}


//! Input: xData = (1|cov), yData (+1/-1)
void TestLRPvals::trueNewtoniteration(double*& wData, double** xData, double* yData, long factorDim, long sampleDim, double gamma){

    double* twData = new double[factorDim];
    for(long i = 0; i < factorDim; ++i){
        twData[i] = wData[i];
    }
    
    double* ip = new double[sampleDim];
    double* pvals = new double[sampleDim];   //! sigmoid(ip)
    double* wvec = new double[sampleDim];    //! p * (1-p)
    double* zvec = new double[sampleDim];
    
    calculatePvals(pvals, wvec, zvec, twData, yData, xData, 0, factorDim, sampleDim, false);   //! learning rate = (gamma/sampleDim)
    
    vec_RR Wvec;
    Wvec.SetLength(sampleDim);
    for(long i = 0; i < sampleDim; ++i){
        Wvec[i] = to_RR(wvec[i]);
    }
    
    vec_RR Zvec;
    Zvec.SetLength(sampleDim);
    for(long i = 0; i < sampleDim; ++i){
        Zvec[i] = to_RR(zvec[i]);
    }
    
    //! information for n users
    //! Xcov_small[i]: information of ith user
    vector<Mat<RR>> Xcov_small;
    for(long i = 0; i < sampleDim; ++i){
        Mat<RR> temp;
        temp.SetDims(factorDim, factorDim);
        for(long j = 0; j < factorDim; ++j){
            for(long l = 0; l < factorDim; ++l){
                temp[j][l] = to_RR(xData[i][j]) * to_RR(xData[i][l]);
            }
        }
        Xcov_small.push_back(temp);
    }
    
    //! Xcov = X^T * W * X
    Mat<RR> Xcov;
    mul(Xcov, Xcov_small[0], Wvec[0]);  //! k * k
    
    Mat<RR> rtemp;
    for(long i = 1; i < sampleDim; ++i){
        mul(rtemp, Xcov_small[i], Wvec[i]);
        add(Xcov, Xcov, rtemp);
    }
    
    Mat<RR> Xcovinv;
    inv(Xcovinv, Xcov);
    
    //! ZWXvec = z^T * W* X
    vec_RR ZWXvec;
    ZWXvec.SetLength(factorDim);
    for(long k = 0; k < factorDim; ++k){
        ZWXvec[k] = to_RR("0");
        for(long i = 0; i < sampleDim; ++i){
            ZWXvec[k] +=  (Wvec[i] * Zvec[i]) * to_RR(xData[i][k]);
        }
    }
    
    // ZWXvec^T = X^T * W * z
    Mat<RR> XWZmat;
    XWZmat.SetDims(factorDim, 1);
    for(long k = 0; k < factorDim; ++k){
        XWZmat[k][0] = ZWXvec[k];
    }
    
    Mat<RR> res;
    mul(res, Xcovinv, XWZmat);
    for(long k = 0; k < factorDim; ++k){
        conv(wData[k], res[k][0]);
    }
}



//! Covariates: Logistic regressio using Newton approach
void TestLRPvals::testLogRegNT(double*& zScore, double*& pvals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long numNTIter, double gammaUp, string zfilename, string pfilename){
    
    double** yxData = new double*[sampleDim];
    for(long i = 0; i < sampleDim; ++i){
        yxData[i] = new double[factorDim];
        for(long j = 0; j < factorDim; ++j){
            yxData[i][j]  = xData[i][j] * yData[i];
        }
    }
    
    cout << "+------------------------------------+" << endl;
    cout << "|          2. Pre-training           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    //! Pre-training

    double* wData = new double[factorDim];
    initialVecZero(wData, factorDim);      //! Set as zero
    
    for (long iter = 0; iter < numNTIter; ++iter) {
        trueNewtoniteration(wData, xData, yData, factorDim, sampleDim, gammaUp);
        cout << iter << ": " << endl;
        printvector(wData, factorDim);
    }
    //wData[0] = 0.6743; wData[1] = 4.044; wData[2] = 5.5659;  wData[3] = -3.7416; //!(from R)

    double truecor, trueauc;
    calculateAUC(xData, wData, factorDim, sampleDim, truecor, trueauc);
    cout << "Logistic regression with covariates (using Newton) / niter = " << numNTIter ;
    cout << ", Correctness: " << truecor << " %, AUC: " << trueauc << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|         3. Training / snps         |" << endl;
    cout << "+------------------------------------+" << endl;
    
    double* ip = new double[sampleDim];
    double* prob = new double[sampleDim];   //! sigmoid(ip)
    double* wvec = new double[sampleDim];    //! p * (1-p)
    double* zvec = new double[sampleDim];
    
    calculatePvals(prob, wvec, zvec, wData, yData, yxData, 0, factorDim, sampleDim);
    
    cout << "wData: "  ;
    printvector(wData, factorDim);
    
    cout << "p[i]: " ;
    printvector(prob, 10);
    
    cout << "wvec: " ;
    printvector(wvec, 10);
    
    cout << "zvec: " ;
    printvector(zvec, 10);
    
    //updateIRLS(zScore, wvec, zvec, xData, sData, factorDim, sampleDim, nsnp, filename);
    updateIRLSFast(zScore, pvals, wvec, zvec, yData, xData, sData, factorDim, sampleDim, nsnp, zfilename, pfilename);

}



//
////! Update beta and zScore
////! where beta[j] = (sstar[j]^T * W * sstar[j])^-1 * (sstar[j]^T * W * zstar) for 1 <= j <= nsnp
//
//void TestLRPvals::updateIRLS(vec_RR& zScore, double* wvec, double* zvec, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename){
//    Mat<RR> Wmat; Wmat.SetDims(sampleDim, sampleDim);
//    Mat<RR> Xmat; Xmat.SetDims(sampleDim, factorDim);
//    Mat<RR> Imat; Imat.SetDims(sampleDim, sampleDim);
//
//    for(long i = 0; i < sampleDim; ++i){
//        for(long j = 0; j < sampleDim; ++j){
//            Wmat[i][j] = to_RR("0");
//            Imat[i][j] = to_RR("0");
//        }
//        Wmat[i][i] = to_RR(wvec[i]);
//        Imat[i][i] = to_RR("1");
//    }
//    for(long i = 0; i < sampleDim; ++i){
//        for(long j = 0; j < factorDim; ++j){
//            Xmat[i][j]  = to_RR(xData[i][j] * xData[i][0]);
//        }
//    }
//
//    Mat<RR> Xtrans;
//    transpose(Xtrans, Xmat);
//
//    //! Xcov = X^T * W * X (: k * k)
//    Mat<RR> Xcov;
//    mul(Xcov, Xtrans, Wmat);
//    mul(Xcov, Xcov, Xmat);
//
//    Mat<RR> Xcovinv;
//    inv(Xcovinv, Xcov);
//
//    //! Nx = I - X * Xcovinv * X^T * W
//    Mat<RR> Nx;   //! n * n
//    mul(Nx, Xmat, Xcovinv);
//    mul(Nx, Nx, Xtrans);
//    mul(Nx, Nx, Wmat);
//    sub(Nx, Imat, Nx);
//
//    //printRmatrix(Nx, 10);
//
//    //! zstar = Nx * z
//    Mat<RR> Zmat;
//    Zmat.SetDims(sampleDim, 1);
//    for(long i = 0; i < sampleDim; ++i){
//        Zmat[i][0]  = to_RR(zvec[i]);
//    }
//    Mat<RR> Zstar;
//    mul(Zstar, Nx, Zmat);
//
//
//    vec_RR Zvec;
//    Zvec.SetLength(sampleDim);
//    for(long i = 0; i < sampleDim; ++i){
//        Zvec[i] = Zstar[i][0] * Wmat[i][i];    //! for simplicity, z <= z * w
//    }
//
//    //! sstar = Nx * s: n * p
//    Mat<RR> Smat, Sstar;
//    Smat.SetDims(sampleDim, nsnp);
//    for(long i = 0; i < sampleDim; ++i){
//        for(long j = 0; j < nsnp; ++j){
//            Smat[i][j]  = to_RR(sData[i][j]);
//        }
//    }
//    mul(Sstar, Nx, Smat);
//
//
//    //! compute beta, var, and alpha = |beta/sigma|
//    vec_RR beta, var;
//    beta.SetLength(nsnp);
//    var.SetLength(nsnp);
//    zScore.SetLength(nsnp);
//
//    for(long j = 0; j < nsnp; ++j){
//        vec_RR Svec; //! jth snp
//        Svec.SetLength(sampleDim);
//        for(long i = 0; i < sampleDim; ++i){
//            Svec[i] = Sstar[i][j];
//        }
//
//        vec_RR WSvec; //! jth snp
//        WSvec.SetLength(sampleDim);
//        for(long i = 0; i < sampleDim; ++i){
//            WSvec[i] = Sstar[i][j] * Wmat[i][i];
//        }
//
//        RR ZS, WS;
//        InnerProduct(ZS, Zvec, Svec);   //! sum (wi * zi) * si
//        InnerProduct(WS, WSvec, Svec);  //! sum (wi * si) * si
//
//        beta[j] = ZS/ WS;
//        var[j]  = inv(WS);   //! factordim = k + 1
//        zScore[j] = abs(beta[j]/sqrt(var[j]));
//    }
//
//    fstream outf;
//    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);
//    outf << "beta: " << endl;
//    outf.close();
//
//    printRvectorToFile(beta, filename);
//
//    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);
//    outf << "zScore: " << endl;
//    outf.close();
//    printRvectorToFile(zScore, filename);
//
//}


