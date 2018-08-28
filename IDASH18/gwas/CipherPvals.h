/*
 * @file       CipherPvalues.h, header file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */

#ifndef HEML_CIPHERPVALUES_H_
#define HEML_CIPHERPVALUES_H_

#include <complex>

#include <iostream>
#include <vector>
#include <string>

#include "NTL/RR.h"
#include "NTL/mat_RR.h"
#include "NTL/vec_RR.h"

#include "../src/Scheme.h"
#include "../src/SecretKey.h"

#include "../src/ExtScheme.h"

using namespace std;
using namespace NTL;


//! Testing functions
class CipherPvals
{
public:
    
    Scheme& scheme;
    SecretKey& secretKey;
    ExtScheme& extscheme;
    
    CipherPvals(Scheme& scheme, SecretKey& secretKey, ExtScheme& extscheme) : scheme(scheme), secretKey(secretKey), extscheme(extscheme) {}
   
    /********************************************************************/
    //! (Basic) Encryption/Decryption functions
    void encValue(Ciphertext& encData, double data, long nslots, long L);
    
    void encFullyPackedVec(Ciphertext& encData, double* data, long nslots, long L);
    
    void encSparselyPackedVec(Ciphertext& encData, double* data, long len, long nslots, long L);
    
    void decSingleData(double& Data, Ciphertext encData);
    
    void decVector(double*& Data, Ciphertext encData, long len);
    
    void computeCov(double*& covariance, double* data, long dim, long scalefactor);

    /********************************************************************/
    //! Encryption functions for data
    
    void encryptSData(Ciphertext**& encSData, Ciphertext**& encYSData, Ciphertext***& encSXData,  double* yData, double** xData, double** sData,  long factorDim, long sampleDim, long nsnp, long nencsnp, long nslots, long L) ;
    
    void encryptXData(Ciphertext**& encYXData, Ciphertext**& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long dim, long nterms, long scalefactor, long nslots, long L) ;
    
    void encryptSIMDXData(Ciphertext& encYXData, Ciphertext*& enccovData, double* yData, double** xData, long factorDim, long sampleDim, long sampleDim2, long nXbatching, long nCovbatching, long dim, long nterms, long scalefactor, long nslots, long L) ;
    
    /********************************************************************/
    
    void aggYXData(Ciphertext*& encYX, Ciphertext encYXData, long log2sampleDim, long nbatching, long factorDim, long nslots);
    
    void aggCovData(Ciphertext*& encCov, Ciphertext* enccovData,  long log2sampleDim, long nbatching);
    
    void fullReplicate4(Ciphertext*& res, Ciphertext Data, long nslots);
    
    void encSIMDAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* encData);
    
    void encSIMDAdjoint_smalldep(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* encData);
    
    void encAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* encData);
   
    void QuadForm(Ciphertext& res, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim);
    
    void SqrQuadForm(Ciphertext& res, Ciphertext* encData, Ciphertext* encMatrix, long factorDim);
    
    void extencAdjoint(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext* encData);
    
    void extQuadForm(Ciphertext& res, Ciphertext* encData1, Ciphertext* encMatrix, Ciphertext* encData2, long factorDim);
    
    void extSqrQuadForm(Ciphertext& res, Ciphertext* encData, Ciphertext* encMatrix, long factorDim);
    
    /********************************************************************/
    //! just encrypt (for testging)
    void HesInverse(Ciphertext& encDet, Ciphertext*& encAdj, Ciphertext** enccovData, long dim, long nslots, long L);
    
    
};


#endif /* CIPHERPVALUES_H_ */
