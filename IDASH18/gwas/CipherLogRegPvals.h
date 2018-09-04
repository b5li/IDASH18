
#ifndef HEML_CIPHERGD_H_
#define HEML_CIPHERGD_H_

#include <complex>
#include "../src/Scheme.h"
#include "../src/SecretKey.h"
#include "../src/ExtScheme.h"

#include "CipherLinRegPvals.h"

static double scaledsigmoid3[3] = {2,2.40192,-0.407808};  //! (2 + 2.4 * (X/4) - 0.4 * (X/4)^3) * (X/4) = sigmoid(X) * X

static double scaledsigmoid5[4] = {0.5,0.76524,-0.2941632,0.042222797};  //!  (a0 + a1 * (X/4) + a2 * (X/4)^3 + a3 * (X/4)^5) = sigmoid(X)
static double scaledsigmoid7[5] = {0.5,0.867536,-0.5243366,0.16984166,-0.0195922};

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
    
	void encXData(Ciphertext& encYXData, Ciphertext& encXData, Ciphertext& encYData, Ciphertext*& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long nXbatching, long nCovbatching, long nterms, long YXscaleBits, long covscaleBits, long nslots, long YXlvl, long Xlvl, long Ylvl, long Covlvl);

    void encryptSData(Ciphertext**& encSData, Ciphertext***& encSXData, double** xData, double** sData,  long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots, long Slvl, long SXlvl) ;
    

    void decryptResult(double**& ZSnorm, double**& Snorm, double& det, Ciphertext* encZSnorm, Ciphertext* encSnorm, Ciphertext encDet, long nencsnp, long nslots);
    
	/********************************************************************/
    //! Logistic regression
    
    void encNLGDiteration(Ciphertext& encWData, Ciphertext encYXData, long factorDim,  long sampleDim, long fdimBits, long sdimBits, long xBatchingBits, uint64_t* poly, long numIter);
    
	uint64_t* generateNLGDAuxPoly(long nslots, long nslots1, long batch);

    void encNLGDiteration0(Ciphertext& encWData, Ciphertext& encVData, Ciphertext& encGrad, Ciphertext encYXData, long sdimBits, long xBatchingBits, double gamma, double eta);
    
    void encNLGDiteration1(Ciphertext& encWData, Ciphertext& encVData, Ciphertext encYXData, Ciphertext encW0Data, Ciphertext encGrad2, uint64_t* poly, long fdimBits, long sdimBits, long xBatchingBits, double gamma);
    
    void encNLGDiteration2(Ciphertext& encWData, Ciphertext& encVData, Ciphertext encYXData, Ciphertext encW0Data, Ciphertext encGrad0, Ciphertext encGrad2, uint64_t* poly, long fdimBits, long sdimBits, long xBatchingBits, double gamma, double eta);
    
    void encNLGDiteration_final(Ciphertext& encWData, Ciphertext& encVData, Ciphertext encYXData, Ciphertext encW0Data, Ciphertext encGrad2, uint64_t* poly, long fdimBits, long sdimBits, long xBatchingBits, double gamma);
    
    
    void encZWData(Ciphertext& encZData, Ciphertext& encWData, Ciphertext& encW2Data, Ciphertext& encZWData, Ciphertext encBeta, Ciphertext encXData, Ciphertext encYData, uint64_t* poly, long fdimBits, long steps, long sdeg = 3, long scale = 1);
    
    void encTwoInverse(Ciphertext& encWinv, Ciphertext encPr1, Ciphertext encPr2, long steps);
  
	/********************************************************************/
    
    void encWcov(Ciphertext*& encCov, Ciphertext encWData, Ciphertext* enccovData, long sdimBits, long nCovbatching);
    
    void encZXData(Ciphertext*& encZX, Ciphertext encXData, Ciphertext encZData, long sdimBits, long nbatching, long factorDim, long nslots);
    
    void encSXData(Ciphertext**& encSX, Ciphertext*** encSXData, long factorDim, long sampleDim, long nencsnp);
    
    void encVecSData(Ciphertext*& encZS, Ciphertext encVData, Ciphertext** encSData, uint64_t** poly0, uint64_t** poly, long sampleDim, long nencsnp, long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot);
    
    void encVecMultipleSData(Ciphertext**& encZS, Ciphertext encVData, Ciphertext*** encSData, uint64_t** poly0, uint64_t** poly, long factorDim, long sampleDim, long nencsnp, long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot);
    
    /********************************************************************/
    void generateRepAuxPoly(uint64_t**& poly0, uint64_t**& poly, long nslots, long niter, long nstep, long* nblock, long* rot);
    
    void fullReplicate16(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long* rot);
    
    void sparseReplicate16(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots, long nvals, long* rot);
    
    void extQuadForm(Ciphertext& res0, Ciphertext& res1, Ciphertext* encData0, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim);
    
    
    
    
    
    
    /********************************************************************/
    
    
    void encZSData(Ciphertext*& encZS, Ciphertext* encZData, Ciphertext** encSData, long sampleDim, long nencsnp);
    
    void encW2SXData(Ciphertext**& encW2SX, Ciphertext* encW2Data, Ciphertext*** encSXData, long factorDim, long sampleDim, long nencsnp);
    
    void encFastW2SXData(Ciphertext**& encW2SX, Ciphertext* encWData, Ciphertext*** encSXData, long factorDim, long sampleDim, long nencsnp);
    
    void fullReplicate(Ciphertext*& res, Ciphertext encData, uint64_t** poly0, uint64_t** poly, long sampleDim, long nslots, long subblocksize, long niter, long nstep, long nblock, long* rot);
    
    void subReplicate(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long subblocksize, long nstep, long* rot);
    

    void fullReplicate8(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots, long nslots_block);
    
    void sparseReplicate8(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots, long nvals, long* rot);
    
    void fullReplicate4(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots, long nslots_block);
    
    void sparseReplicate4(Ciphertext*& res, Ciphertext encData, uint64_t** poly, long nslots, long nvals, long* rot);
    
    void extQuadForm8(Ciphertext& res0, Ciphertext& res1, Ciphertext* encData0, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim);
    
    void extQuadForm16(Ciphertext& res0, Ciphertext& res1, Ciphertext* encData0, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim);
    
 
};

#endif /* CIPHERGD_H_ */
