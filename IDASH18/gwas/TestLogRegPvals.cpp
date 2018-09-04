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



void TestLRPvals::trueNLGDiteration(double** zData, double*& wData, double*& vData, long factorDim, long sampleDim, double gamma, double eta, long kdeg) {
    double* grad = new double[factorDim]();
    
    //! miran: need to initialization
    for(long i = 0; i < factorDim; ++i) {
        grad[i] = 0.0;
    }
    
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
    }
    
    
    for (long i = 0; i < factorDim; ++i) {
        double tmpw = vData[i] - gamma * grad[i];
        vData[i] = (1.0 - eta) * tmpw + eta * wData[i];
        wData[i] = tmpw;
    }
    delete[] grad;
}

void TestLRPvals::trueNewtoniteration(double** xData, double* yData, vec_RR& twData1, long factorDim, long sampleDim, double gamma){
    Mat<RR> Xmat;
    Xmat.SetDims(sampleDim, factorDim);
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < factorDim; ++j){
            Xmat[i][j]  = to_RR(xData[i][j] * xData[i][0]);
        }
    }
    
    double* twData = new double[factorDim];
    for(long i = 0; i < factorDim; ++i){
        conv(twData[i], twData1[i]);
    }
    
    double* ip = new double[sampleDim];
    double* pvals = new double[sampleDim];   //! sigmoid(ip)
    double* wvec = new double[sampleDim];    //! p * (1-p)
    double* zvec = new double[sampleDim];
    
    calculatePvals(ip, pvals, wvec, zvec, twData, yData, xData, 0, factorDim, sampleDim, gamma);   //! update 1/n
    
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
                temp[j][l] = Xmat[i][j] * Xmat[i][l];
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
            ZWXvec[k] +=  (Wvec[i] * Zvec[i]) * Xmat[i][k];
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
        twData1[k] = res[k][0];
    }
}

void TestLRPvals::calculatePvals(double*& ip, double*& prob, double*& wvec, double*& zvec, double* wData, double* yData, double** zData, long kdeg, long factorDim, long sampleDim, double gamma){
    
    for(int i = 0; i < sampleDim; ++i){
        ip[i] =  zData[i][0] * trueIP(zData[i], wData, factorDim); // <z,w> = (2y-1) * <z', w>, ip(xi, w) = log (p/ (1-p))
        
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
        zvec[i] = ip[i] + (gamma/sampleDim) * (yData0 - prob[i])/(wvec[i]);   // regulate the learning rate (alpha)
        if(i <= 10) cout  << yData0 << "," << ip[i] << "," << prob[i] << "," << wvec[i] << "," << zvec[i] << endl;
    }
}

void TestLRPvals::calculateAUC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& auc) {
    cout << "w:";
    for (long i = 0; i < factorDim; ++i) {
        cout << wData[i] << ",";
    }
    cout << endl;
    
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

/********************************************************************/
//! Covariates: Logistic regressio using GD approach
void TestLRPvals::testLogRegGD(double*& zScore, double*& pvals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long trainsigdeg, long testsigdeg, long numIter,  double gammaUp, double gammaDown, string filename){
    
    double* twData = new double[factorDim];  //! NGD
    double* tvData = new double[factorDim];
    initialWDataVDataZero(twData, tvData, factorDim);
    
    //! NLGD
    double alpha0, alpha1, eta, gamma;
    alpha0 = 0.01;
    alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
    
    double** yxData = new double*[sampleDim];
    for(long i = 0; i < sampleDim; ++i){
        yxData[i] = new double[factorDim];
        for(long j = 0; j < factorDim; ++j){
            yxData[i][j]  = xData[i][j] * yData[i];   // +1/-1
        }
    }
    
    for (long iter = 0; iter < numIter; ++iter) {
        eta = (1 - alpha0) / alpha1;
        gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDim : gammaUp / (iter - gammaDown) / sampleDim;
        trueNLGDiteration(yxData, twData, tvData, factorDim, sampleDim, gamma, eta, trainsigdeg);          //! true sigmoid
        cout << iter << ": " << eta << endl; 
        cout << twData[0] << "," << twData[1] << "," <<  twData[2] << "," << twData[3] << endl;
        cout << tvData[0] << "," << tvData[1] << "," <<  tvData[2] << "," << tvData[3] << endl;
        cout << "----------------------------------------------------------" << endl;
        alpha0 = alpha1;
        alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
    }

    ofstream outf(filename);
    outf.close();
    outf.open(filename, fstream::in | fstream::out | fstream::app);   /// open the file
    outf << "Logistic regression without snp (NLGD) / " ;
    outf << "niter: " << numIter << ", gammaup: " << gammaUp << endl;
    outf.close();
    
    double truecor, trueauc;
    calculateAUC(xData, twData, factorDim, sampleDim, truecor, trueauc);
    
    outf.open(filename, fstream::in | fstream::out | fstream::app);   /// open the file
    outf << "Correctness: " << truecor << " %, AUC: " << trueauc << endl;
    outf << "-------------------------------------------" << endl;
    outf.close();

 
    double* ip = new double[sampleDim];
    double* prob = new double[sampleDim];   //! sigmoid(ip)
    double* wvec = new double[sampleDim];     //! p * (1-p)
    double* zvec = new double[sampleDim];
    
    calculatePvals(ip, prob, wvec, zvec, twData, yData, yxData, testsigdeg, factorDim, sampleDim, sampleDim);  // sigmoid -> deg5
 
    updateIRLSFast(zScore, wvec, zvec, yData, xData, sData, factorDim, sampleDim, nsnp, "LogRegResult/Plain_Zscores_deg3.txt");
    
    for(long i = 0; i < nsnp; ++i){
        pvals[i] = pnorm(zScore[i]);
    }
    cout << "zScore: " ;
    printvector(zScore, 10);
    printvectorToFile(zScore, "LogRegResult/Plain_Zscores_deg3.txt", nsnp);
    
    cout << "pvals: " ;
    printvector(pvals, 10);
    printvectorToFile(pvals, filename, nsnp);
}

void TestLRPvals::updateIRLSFast(double*& zScore, double* wvec, double* zvec, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename){
    
    Mat<RR> Xmat; Xmat.SetDims(sampleDim, factorDim);
    Mat<RR> Imat; Imat.SetDims(sampleDim, sampleDim);
    
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
        for(long j = 0; j < sampleDim; ++j){
            Imat[i][j] = to_RR("0");
        }
        Imat[i][i] = to_RR("1");
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
    
    //! Xcov = X^T * W * X
    //! Xcov2 = X^T * W2 * X
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
    
    cout << " X^T * W * X: " << endl;
    printRmatrix(Xcov);

    Mat<RR> Xcovinv;
    inv(Xcovinv, Xcov);
    
    cout << " (X^T * W * X)^-1: " << endl;
    printRmatrix(Xcovinv);
    
    //! Xcovinv2 = (X^T W X) ^-1 * (X^T W^2 X) * (X^T W X) ^-1
    Mat<RR> Xcovinv2;
    mul(Xcovinv2, Xcovinv, Xcov2);
    mul(Xcovinv2, Xcovinv2, Xcovinv);
    
    vec_RR beta, var, RSS, err;
    beta.SetLength(nsnp);
    var.SetLength(nsnp);
    err.SetLength(nsnp);

    vec_RR Snorm;
    vec_RR ZSnorm;
    Snorm.SetLength(nsnp);
    ZSnorm.SetLength(nsnp);
    
    //------------------------------------
    //! ZXvec = z^T * X,
    cout << "Z^T * X: " << endl;
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
    for(long k = 0; k < factorDim; ++k){
        ZWXvec[k] = to_RR("0");
        for(long i = 0; i < sampleDim; ++i){
            ZWXvec[k] +=  (Wvec[i] * Zvec[i]) * Xmat[i][k];
        }
    }
    
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
        
        Snorm[j] = (temp1 - temp2);
        //RR Snorm = (temp1 - temp2 - temp3 + temp4);  // S^T * W * S
        //cout << j << ": " << (temp3 - temp4) << ", ";
        
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
        
        ZSnorm[j] = (temp1 - temp2);
        //beta[j] = (temp1 - temp2 - temp3 + temp4);   // Z^T * W * S
        //cout << ": " << (temp3 - temp4) << endl;
        
        beta[j] = (ZSnorm[j]/ Snorm[j]);
        err[j] = sqrt(var[j]);
        temp1  = abs(beta[j]/err[j]);
        conv(zScore[j], temp1);
    }
    
    printRvectorToFile1(ZSnorm, "LogRegResult/Plain_ZSnorm.txt");
    printRvectorToFile1(Snorm, "LogRegResult/Plain_Snorm.txt");
    
    printRvectorToFile1(beta, "LogRegResult/pt_beta.txt");
    printRvectorToFile1(err, "LogRegResult/pt_err.txt");
    //printvectorToFile(zScore, filename, nsnp);

    cout << "Zvec: [";
    for(long i = 0; i < 10; ++i) cout << Zvec[i] << "," ;
    cout << "]" << endl;
    
    cout << "Z^T * W * S : [ " ;
    for(long i = 0; i < 10; ++i) cout << ZWS[i] << "," ;
    cout << "]" << endl;
    
    cout << "S^T * W * S : [ " ;
    for(long i = 0; i < 10; ++i) cout << SWS[i] << "," ;
    cout << "]" << endl;
    
    cout << "X^T * W2 * S : [ " ;
    for(long i = 0; i < 10; ++i) cout << SW2X[i][0] << "," ;
    cout << "]" << endl;
  
}





//! Covariates: Logistic regressio using Newton approach
void TestLRPvals::testLogRegNT(double*& zScore, double*& pvals, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, long numIter, double gammaUp, string filename){
    
    vec_RR twData1; //! Newton
    twData1.SetLength(factorDim);
    for(long i = 0; i < factorDim; ++i){
        twData1[i] = to_RR("0");
    }
    
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < factorDim; ++j){
            xData[i][j]  = xData[i][j] * yData[i];
        }
    }
    
    for (long iter = 0; iter < numIter; ++iter) {
        trueNewtoniteration(xData, yData, twData1, factorDim, sampleDim, gammaUp);
        cout << iter << ": " << endl;
        printRvector(twData1);
    }
    double* twData = new double[factorDim];
    for(long k = 0; k < factorDim; ++k){
        conv(twData[k], twData1[k]);
    }
    
    //twData[0] = 0.6743;
    //twData[1] = 4.044;
    //twData[2] = 5.5659;
    //twData[3] = -3.7416;
    
    ofstream outf(filename);
    outf.close();
    
    outf.open(filename, fstream::in | fstream::out | fstream::app);   /// open the file
    outf << "Logistic regression without snp (Newton) / niter:" << numIter << endl;
    outf.close();
    
    
    double truecor, trueauc;
    calculateAUC(xData, twData, factorDim, sampleDim, truecor, trueauc);
    printvectorToFile(twData, filename, factorDim);
    
    outf.open(filename, fstream::in | fstream::out | fstream::app);   /// open the file
    outf << "Correctness: " << truecor << " %, AUC: " << trueauc << endl;
    outf << "-------------------------------------------" << endl;
    outf.close();
    
    double* ip = new double[sampleDim];
    double* prob = new double[sampleDim];   //! sigmoid(ip)
    double* wvec = new double[sampleDim];    //! p * (1-p)
    double* zvec = new double[sampleDim];
    
    calculatePvals(ip, prob, wvec, zvec, twData, yData, xData, 0, factorDim, sampleDim, sampleDim);
    
    cout << "wData: "  ;
    printvector(twData, factorDim);
    
    cout << "p[i]: " ;
    printvector(prob, 10);
    
    cout << "zvec: " ;
    printvector(zvec, 10);
    
    //updateIRLS(zScore, wvec, zvec, xData, sData, factorDim, sampleDim, nsnp, filename);
    updateIRLSFast(zScore, wvec, zvec, yData, xData, sData, factorDim, sampleDim, nsnp, filename);
    
    pvals = new double[nsnp];
    for(long i = 0; i < nsnp; ++i){
        double dtemp;
        conv(dtemp, zScore[i]);
        pvals[i] = pnorm(dtemp);
    }
}

//! Update beta and zScore
//! where beta[j] = (sstar[j]^T * W * sstar[j])^-1 * (sstar[j]^T * W * zstar) for 1 <= j <= nsnp

void TestLRPvals::updateIRLS(vec_RR& zScore, double* wvec, double* zvec, double* yData, double** xData, double** sData, long factorDim, long sampleDim, long nsnp, string filename){
    Mat<RR> Wmat; Wmat.SetDims(sampleDim, sampleDim);
    Mat<RR> Xmat; Xmat.SetDims(sampleDim, factorDim);
    Mat<RR> Imat; Imat.SetDims(sampleDim, sampleDim);
    
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < sampleDim; ++j){
            Wmat[i][j] = to_RR("0");
            Imat[i][j] = to_RR("0");
        }
        Wmat[i][i] = to_RR(wvec[i]);
        Imat[i][i] = to_RR("1");
    }
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < factorDim; ++j){
            Xmat[i][j]  = to_RR(xData[i][j] * xData[i][0]);
        }
    }
    
    Mat<RR> Xtrans;
    transpose(Xtrans, Xmat);
    
    //! Xcov = X^T * W * X (: k * k)
    Mat<RR> Xcov;
    mul(Xcov, Xtrans, Wmat);
    mul(Xcov, Xcov, Xmat);
    
    Mat<RR> Xcovinv;
    inv(Xcovinv, Xcov);
    
    //! Nx = I - X * Xcovinv * X^T * W
    Mat<RR> Nx;   //! n * n
    mul(Nx, Xmat, Xcovinv);
    mul(Nx, Nx, Xtrans);
    mul(Nx, Nx, Wmat);
    sub(Nx, Imat, Nx);
    
    //printRmatrix(Nx, 10);
    
    //! zstar = Nx * z
    Mat<RR> Zmat;
    Zmat.SetDims(sampleDim, 1);
    for(long i = 0; i < sampleDim; ++i){
        Zmat[i][0]  = to_RR(zvec[i]);
    }
    Mat<RR> Zstar;
    mul(Zstar, Nx, Zmat);
    
    
    vec_RR Zvec;
    Zvec.SetLength(sampleDim);
    for(long i = 0; i < sampleDim; ++i){
        Zvec[i] = Zstar[i][0] * Wmat[i][i];    //! for simplicity, z <= z * w
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
    vec_RR beta, var;
    beta.SetLength(nsnp);
    var.SetLength(nsnp);
    zScore.SetLength(nsnp);
    
    for(long j = 0; j < nsnp; ++j){
        vec_RR Svec; //! jth snp
        Svec.SetLength(sampleDim);
        for(long i = 0; i < sampleDim; ++i){
            Svec[i] = Sstar[i][j];
        }
        
        vec_RR WSvec; //! jth snp
        WSvec.SetLength(sampleDim);
        for(long i = 0; i < sampleDim; ++i){
            WSvec[i] = Sstar[i][j] * Wmat[i][i];
        }
        
        RR ZS, WS;
        InnerProduct(ZS, Zvec, Svec);   //! sum (wi * zi) * si
        InnerProduct(WS, WSvec, Svec);  //! sum (wi * si) * si
        
        beta[j] = ZS/ WS;
        var[j]  = inv(WS);   //! factordim = k + 1
        zScore[j] = abs(beta[j]/sqrt(var[j]));
    }
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);
    outf << "beta: " << endl;
    outf.close();
    
    printRvectorToFile(beta, filename);
    
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);
    outf << "zScore: " << endl;
    outf.close();
    printRvectorToFile(zScore, filename);
    
}


