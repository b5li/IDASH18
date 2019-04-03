
#ifndef HEML_CIPHERGD_H_
#define HEML_CIPHERGD_H_

#include <complex>
#include "../src/Scheme.h"
#include "../src/SecretKey.h"
#include "../src/ExtScheme.h"

#include "CipherLinRegPvals.h"

static double scaledsigmoid3[3] = {2,2.40192,-0.407808};  //! (2 + 2.4 * (X/4) - 0.4 * (X/4)^3) * (X/4) = sigmoid(X) * X

static double scaledsigmoid5[4] = {0.5,0.76524,-0.2941632,0.042222797};  //!  (a0 + a1 * (X/4) + a2 * (X/4)^3 + a3 * (X/4)^5) = sigmoid(X)
static double scaledsigmoid7[5] = {0.5,3.470144,-33.557545,173.917864,-320.9978};   //! scaled by "16"

using namespace std;
using namespace NTL;

class CipherLRPvals {
public:
	Scheme& scheme;
	SecretKey& secretKey;
    ExtScheme& extscheme;
    CipherPvals& cipherPvals;
    
	CipherLRPvals(Scheme& scheme, SecretKey& secretKey, ExtScheme& extscheme, CipherPvals& cipherPvals) : scheme(scheme), secretKey(secretKey), extscheme(extscheme), cipherPvals(cipherPvals) {}
    

    /********************************************************************/
    //! Encryption functions for data
    
	void encXData(Ciphertext& encYXData, Ciphertext& encXData, Ciphertext& encYData, Ciphertext*& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long nXbatching, long nCovbatching, long nterms, long YXscaleBits, long covscaleBits, long nslots, long YXlvl, long Xlvl, long Ylvl, long Covlvl, long Covlvl1);
    
    void new_encXData(Ciphertext& encYXData, Ciphertext& encXData, Ciphertext& encYData, Ciphertext*& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long nXbatching, long nCovbatching, long nterms, long YXscaleBits, long covscaleBits, long nslots, long YXlvl, long Xlvl, long Ylvl, long Covlvl, long Covlvl1);

    void encryptSData(Ciphertext**& encSData, Ciphertext***& encSXData, double** xData, double** sData,  long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots, long Slvl, long SXlvl) ;
    
    void new_encSData(Ciphertext***& encSXData, double** xData, double** sData,  long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots, long SXlvl) ;

    void decryptResult(double**& ZSnorm, double**& Snorm, double& det, Ciphertext* encZSnorm, Ciphertext* encSnorm, Ciphertext encDet, long nencsnp, long nslots);
    
	/********************************************************************/
    //! Logistic regression
    
    void encNLGD(Ciphertext& encWData, Ciphertext encYXData, long factorDim,  long sampleDim, long fdimBits, long sdimBits, long xBatchingBits, uint64_t* poly, long numIter, const double gammaUp = 1.0);
    
	uint64_t* generateNLGDAuxPoly(long nslots, long nslots1, long batch);

    void encNLGD_oneiter(Ciphertext& encWData, Ciphertext encYXData, long sdimBits, long xBatchingBits, double gamma);
    
    void encNLGDiteration0(Ciphertext*& encData, Ciphertext& encGrad, Ciphertext encYXData, long sdimBits, long xBatchingBits, double gamma, double eta);
    
    void encNLGDiteration_eta0(Ciphertext& encWData, Ciphertext& encVData, Ciphertext encYXData, Ciphertext encW0Data, Ciphertext encGrad2, uint64_t* poly, long fdimBits, long sdimBits, long xBatchingBits, double gamma);
    
    void encNLGDiteration(Ciphertext& encWData, Ciphertext& encVData, Ciphertext encYXData, Ciphertext encW0Data, Ciphertext encGrad0, Ciphertext encGrad2, uint64_t* poly, long fdimBits, long sdimBits, long xBatchingBits, double gamma, double eta);
    
    void encNLGDiterationW(Ciphertext& encWData, Ciphertext encVData, Ciphertext encYXData, Ciphertext encW0Data, Ciphertext encGrad2, uint64_t* poly, long fdimBits, long sdimBits, long xBatchingBits, double gamma);

    
	/********************************************************************/
    
    void encZWData(Ciphertext& encWData, Ciphertext& encZWData, Ciphertext encBeta, Ciphertext encXData, Ciphertext encYData, uint64_t* poly, long fdimBits, long sdeg = 3, long scale = 1);
    
    void encAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext encWData, Ciphertext* enccovData, long sdimBits, long nCovbatching);
    
    void new_encAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext encWData, Ciphertext* enccovData, long sdimBits, long nCovbatching);

    void encZXData(Ciphertext*& encZX, Ciphertext encXData, Ciphertext encZData, long sdimBits, long nbatching, long factorDim, long nslots);
    
    void encVecSData(Ciphertext*& encRes, Ciphertext encVecData, Ciphertext** encSData, uint64_t** poly0, uint64_t** poly, long sampleDim, long nencsnp,  long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot);
    
    void encVecSXData(Ciphertext**& encRes, Ciphertext encVecData, Ciphertext*** encSXData, uint64_t** poly0, uint64_t** poly, long factorDim, long sampleDim, long nencsnp, long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot);
    
   
    
    void generateRepAuxPoly(uint64_t**& poly0, uint64_t**& poly, long nslots, long niter, long nstep, long* nblock, long* rot, long polylvl);
    void fullReplicate(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long* rot, long subblocksize, long nstep);
    void sparseReplicate(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots, long nvals, long* rot, long subblocksize, long nstep);
    void fullReplicate16(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long* rot);
    void sparseReplicate16(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots, long nvals, long* rot);
    
   
    void extQuadForm(Ciphertext& res0, Ciphertext& res1, Ciphertext* encData0,  Ciphertext* encMatrix, Ciphertext* encData2, long factorDim);
    
    
    /********************************************************************/
    
    void encSXData(Ciphertext**& encSX, Ciphertext*** encSXData, long factorDim, long sampleDim, long nencsnp);
    
    void encZWData(Ciphertext& encZData, Ciphertext& encWData, Ciphertext& encW2Data, Ciphertext& encZWData, Ciphertext encBeta, Ciphertext encXData, Ciphertext encYData, uint64_t* poly, long fdimBits, long sdeg = 3, long scale = 1);
    
    void encTwoInverse(Ciphertext& encWinv, Ciphertext encPr1, Ciphertext encPr2, long steps);
    
    void encWinverse(Ciphertext&Winv, Ciphertext encPr);
    
     
    void encVecMultipleSData(Ciphertext**& encZS, Ciphertext encVData, Ciphertext*** encSData, uint64_t** poly0, uint64_t** poly, long factorDim, long sampleDim, long nencsnp, long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot);
    
    void encVecSData_new(Ciphertext*& encZS, Ciphertext encVData, Ciphertext** encSData, uint64_t** poly0, uint64_t** poly, long sampleDim, long nencsnp, long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot);
    
    void encVecMultipleSData_new(Ciphertext**& encZS, Ciphertext encVData, Ciphertext*** encSData, uint64_t** poly0, uint64_t** poly, long factorDim, long sampleDim, long nencsnp, long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot);
    
    
    void extQuadForm(Ciphertext& res0, Ciphertext& res1, Ciphertext* encData0, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim);
    
    void fullReplicate4(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long* rot);
    
    void sparseReplicate4(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots, long nvals, long* rot);
 
    
};

#endif /* CIPHERGD_H_ */
