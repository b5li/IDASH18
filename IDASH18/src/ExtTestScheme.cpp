
#include "TestScheme.h"
#include "Numb.h"
#include "Context.h"
#include "SecretKey.h"
#include "Scheme.h"
#include "EvaluatorUtils.h"
#include "StringUtils.h"
#include "TimeUtils.h"
#include "SchemeAlgo.h"

#include "ExtCiphertext.h"
#include "ExtScheme.h"
#include "ExtTestScheme.h"

#include <iostream>
#include <vector>
#include <chrono>

// 
using namespace std;
using namespace chrono;

void ExtTestScheme::testBasic(long logN,  long L, long logp, long logSlots){
    cout << "!!! START TEST PROD BATCH !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    ExtScheme extscheme(secretKey, context, scheme);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = 1 << logSlots;
    
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots, bound);
    
    complex<double>* mvecNeg = new complex<double>[slots];
    complex<double>* mvecAdd = new complex<double>[slots];
    complex<double>* mvecSub = new complex<double>[slots];
    
    for(long i = 0; i < slots; i++) {
        mvecNeg[i] = - mvec1[i];
        mvecAdd[i] = mvec1[i] + mvec2[i];
        mvecSub[i] = mvec1[i] - mvec2[i];
    }
    
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    
    //! conversion
    timeutils.start("Conversion To Extended Ciphertext");
    ExtCiphertext extcipher1 =  extscheme.toExtCipher(cipher1);
    ExtCiphertext extcipher2 =  extscheme.toExtCipher(cipher2);
    timeutils.stop("Conversion To Extended Ciphertext");
    
    timeutils.start("Conversion To Normal Ciphertext");
    Ciphertext ctxt = extscheme.toNormalCipher(extcipher1);
    timeutils.stop("Conversion To Normal Ciphertext");
 
    complex<double>* dvec = scheme.decrypt(secretKey, ctxt);
    StringUtils::showcompare(mvec1, dvec, slots, "val");
    
    //! negate
    timeutils.start("Homomorphic Negate");
    ExtCiphertext extNegCipher = extscheme.negate(extcipher1);
    timeutils.stop("Homomorphic Negate");
    
    Ciphertext negCipher1 = extscheme.toNormalCipher(extNegCipher);
    complex<double>* dvecNeg1 = scheme.decrypt(secretKey, negCipher1);
    StringUtils::showcompare(mvecNeg, dvecNeg1, slots, "neg");
    
    //! negateAndEqual
    timeutils.start("Homomorphic NegateAndEqual");
    ExtCiphertext extNegCipherEqual(extcipher1);
    extscheme.negateAndEqual(extNegCipherEqual);
    timeutils.stop("Homomorphic NegateAndEqual");
    
    Ciphertext negCipher2 = extscheme.toNormalCipher(extNegCipherEqual);
    complex<double>* dvecNeg2 = scheme.decrypt(secretKey, negCipher2);
    StringUtils::showcompare(mvecNeg, dvecNeg2, slots, "neg");
    
    //! addition
    timeutils.start("Homomorphic Addition");
    ExtCiphertext extAddCipher = extscheme.add(extcipher1, extcipher2);
    timeutils.stop("Homomorphic Addition");
    
    Ciphertext addCipher1 = extscheme.toNormalCipher(extAddCipher);
    complex<double>* dvecAdd1 = scheme.decrypt(secretKey, addCipher1);
    StringUtils::showcompare(mvecAdd, dvecAdd1, slots, "add");
    
    timeutils.start("Homomorphic AdditionAndEqual");
    ExtCiphertext extAddCipherEqual = extcipher1;
    extscheme.addAndEqual(extAddCipherEqual, extcipher2);
    timeutils.stop("Homomorphic AdditionAndEqual");
    
    Ciphertext addCipher2 = extscheme.toNormalCipher(extAddCipherEqual);
    complex<double>* dvecAdd2 = scheme.decrypt(secretKey, addCipher2);
    StringUtils::showcompare(mvecAdd, dvecAdd2, slots, "add");
    
    
    //! subtraction
    timeutils.start("Homomorphic Subtraction");
    ExtCiphertext extSubCipher = extscheme.sub(extcipher1, extcipher2);
    timeutils.stop("Homomorphic Subtraction");
    
    Ciphertext subCipher1 = extscheme.toNormalCipher(extSubCipher);
    complex<double>* dvecSub1 = scheme.decrypt(secretKey, subCipher1);
    StringUtils::showcompare(mvecSub, dvecSub1, slots, "sub");
    
    timeutils.start("Homomorphic SubtractionAndEqual");
    ExtCiphertext extSubCipherEqual = extcipher1;
    extscheme.subAndEqual(extSubCipherEqual, extcipher2);
    timeutils.stop("Homomorphic SubtractionAndEqual");
    
    Ciphertext subCipher2 = extscheme.toNormalCipher(extSubCipherEqual);
    complex<double>* dvecSub2 = scheme.decrypt(secretKey, subCipher2);
    StringUtils::showcompare(mvecSub, dvecSub2, slots, "sub");

}


void ExtTestScheme::testThreeProd(long logN,  long L, long logp, long logSlots){
    cout << "!!! START TEST BASIC !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L + 1;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    ExtScheme extscheme(secretKey, context, scheme);
    
    extscheme.addThreeMultKey(secretKey);   // Ps^3 -> s
    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec3 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvecMult = new complex<double>[slots];
    complex<double>* mvecSqr = new complex<double>[slots];

    for(long i = 0; i < slots; i++) {
        mvecMult[i] = mvec1[i] * mvec2[i] * mvec3[i];
        mvecSqr[i] = mvec1[i] * mvec1[i];
    }
    
    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    Ciphertext cipher3 = scheme.encrypt(mvec3, slots, L);
    timeutils.stop("Encrypt two batch");
    
#if 0
    //! Naive (c1*c2, *c3)
    timeutils.start("Homomorphic Multiplication & Rescaling (Naive)");
    Ciphertext naivemultCipher = scheme.mult(cipher1, cipher2);
    scheme.multAndEqual(naivemultCipher, cipher3);
    
    auto start= chrono::steady_clock::now();
    scheme.reScaleByAndEqual(naivemultCipher, 2);
    timeutils.stop("Homomorphic Multiplication & Rescaling");
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count();
    cout << "RS time= " << timeElapsed << " ms" << endl;
#endif
    //! try1
    timeutils.start("Homomorphic Multiplication & Rescaling (Three in once)");
    ExtCiphertext extmultCipher = extscheme.rawmult3(cipher1, cipher2, cipher3);
    timeutils.stop("Homomorphic Multiplication & Rescaling");
    
    Ciphertext multCipher = extscheme.ModRaiseKeySwitch(extmultCipher);
    scheme.reScaleByAndEqual(multCipher, 2);
    
    //! try2
    timeutils.start("Homomorphic Multiplication & Rescaling (Three in once)");
    ExtCiphertext extmultCipher0 = extscheme.rawmult3_(cipher1, cipher2, cipher3);
    timeutils.stop("Homomorphic Multiplication & Rescaling");
    
    Ciphertext multCipher0 = extscheme.ModRaiseKeySwitch(extmultCipher0);
    scheme.reScaleByAndEqual(multCipher0, 2);
    
    //! try3
    timeutils.start("Homomorphic Multiplication & Rescaling (two and one)");
    ExtCiphertext extmultCipher1 = extscheme.rawmult(cipher1, cipher2);
    ExtCiphertext extmultCipher2 = extscheme.rawmult(extmultCipher1, cipher3);
    timeutils.stop("Homomorphic Multiplication & Rescaling");
    
    
    timeutils.start("Homomorphic Key Switching (s^3)");
    Ciphertext multCipher1 = extscheme.ModRaiseKeySwitch(extmultCipher2);
    timeutils.stop("Homomorphic Key Switching");
    
    scheme.reScaleByAndEqual(multCipher1, 2);
    

    timeutils.start("Homomorphic Square");
    ExtCiphertext extsqrCiphext = extscheme.rawsquare(cipher1);
    timeutils.stop("Homomorphic Square");
    
    
    Ciphertext sqrCipher = extscheme.ModRaiseKeySwitch(extsqrCiphext);
    scheme.reScaleByAndEqual(sqrCipher, 1);
    
    //complex<double>* dvecMultNaive = scheme.decrypt(secretKey, naivemultCipher);
    complex<double>* dvecMult = scheme.decrypt(secretKey, multCipher);
    complex<double>* dvecMult0 = scheme.decrypt(secretKey, multCipher0);
    complex<double>* dvecMult1 = scheme.decrypt(secretKey, multCipher1);
    complex<double>* dvecSqr = scheme.decrypt(secretKey, sqrCipher);

    //StringUtils::showcompare(mvecMult, dvecMultNaive, slots, "naivemult");
    StringUtils::showcompare(mvecMult, dvecMult, slots, "raw0");
    StringUtils::showcompare(mvecMult, dvecMult0, slots, "raw1");
    StringUtils::showcompare(mvecMult, dvecMult1, slots, "two/one");
    StringUtils::showcompare(mvecSqr,  dvecSqr, slots, "sqr");
}



void ExtTestScheme::testDecompKS(long logN,  long L, long logp, long logSlots){
    cout << "!!! START TEST KEY SWITCHING !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = 1;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    
    ExtScheme extscheme(secretKey, context, scheme);
    
    extscheme.addDecompTwoKeys(secretKey);
    extscheme.addDecompThreeKeys(secretKey);
    
    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec3 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvecMult1 = new complex<double>[slots];
    complex<double>* mvecMult2 = new complex<double>[slots];
    
    for(long i = 0; i < slots; i++) {
        mvecMult1[i] = mvec1[i] * mvec2[i] ;
        mvecMult2[i] = mvec1[i] * mvec2[i] * mvec3[i] ;
    }
    
    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    Ciphertext cipher3 = scheme.encrypt(mvec3, slots, L);
    timeutils.stop("Encrypt two batch");
    
 
    timeutils.start("Homomorphic Multiplication ");
    ExtCiphertext extmultCipher1 = extscheme.rawmult(cipher1, cipher2);
    ExtCiphertext extmultCipher2 = extscheme.rawmult(extmultCipher1, cipher3);
    timeutils.stop("Homomorphic Multiplication");
    
    timeutils.start("Homomorphic Key Switching (s^2)");
    Ciphertext multCipher1 = extscheme.DecompKeySwitch(extmultCipher1);
    timeutils.stop("Homomorphic Key Switching");
    
    scheme.reScaleByAndEqual(multCipher1, 1);
    
    complex<double>* dvecMult1 = scheme.decrypt(secretKey, multCipher1);
    StringUtils::showcompare(mvecMult1, dvecMult1, slots, "raw0");
    
    timeutils.start("Homomorphic Multiplication ");
    Ciphertext multCipher2 = extscheme.mult(cipher1, cipher2);
    timeutils.stop("Homomorphic Multiplication");
    scheme.reScaleByAndEqual(multCipher2, 1);
    
    dvecMult1 = scheme.decrypt(secretKey, multCipher2);
    StringUtils::showcompare(mvecMult1, dvecMult1, slots, "raw0");
    
    timeutils.start("Homomorphic Multiplication (AndEqual) ");
    extscheme.multAndEqual(cipher1, cipher2);
    timeutils.stop("Homomorphic Multiplication");
    scheme.reScaleByAndEqual(cipher1, 1);
    
    dvecMult1 = scheme.decrypt(secretKey, cipher1);
    StringUtils::showcompare(mvecMult1, dvecMult1, slots, "raw0");
    
    timeutils.start("Homomorphic Key Switching (s^2 & s^3)");
    Ciphertext multCipher3 = extscheme.DecompKeySwitch(extmultCipher2);
    timeutils.stop("Homomorphic Key Switching");
    
    scheme.reScaleByAndEqual(multCipher3, 2);
    
    complex<double>* dvecMult2 = scheme.decrypt(secretKey, multCipher3);
    StringUtils::showcompare(mvecMult2, dvecMult2, slots, "raw0");
}

//! test for consecutive multiplications using decompoision KS
void ExtTestScheme::testProdDecompKS(long logN,  long L, long logp, long logSlots){
    cout << "!!! START TEST KEY SWITCHING !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L + 1;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    
    ExtScheme extscheme(secretKey, context, scheme);
    
    extscheme.addDecompTwoKeys(secretKey);
    extscheme.addDecompThreeKeys(secretKey);
    
    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvecMult = new complex<double>[slots];
  
    for(long i = 0; i < slots; i++) {
        mvecMult[i] = mvec1[i];
    }
    
    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    timeutils.stop("Encrypt two batch");

    Ciphertext multCipher = cipher1;
    
    for(long i = 0; i < L - 2; ++i){
        Ciphertext tmp = scheme.modDownTo(cipher2, multCipher.l);
        multCipher = extscheme.mult(multCipher, tmp);
        scheme.reScaleByAndEqual(multCipher, 1);
        
        for(long j = 0; j < slots; j++) {
            mvecMult[j] *= mvec2[j];
        }
        cout << i + 1 << "-level mult" << endl;
        complex<double>* dvecMult = scheme.decrypt(secretKey, multCipher);
        StringUtils::showcompare(mvecMult, dvecMult, slots, "mul");
    }
}

void ExtTestScheme::testDecompRotate(long logN,  long L, long logRotSlots, long RotSlots, long logp, long logSlots){
    
    cout << "!!! START TEST Rotation !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = 1;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    ExtScheme extscheme(secretKey, context, scheme);
    
    extscheme.addDecompTwoKeys(secretKey);
    extscheme.addDecompThreeKeys(secretKey);
    
    extscheme.addDecompLeftRotKeys(secretKey);
    extscheme.addDecompRightRotKeys(secretKey);
    
    long slots = (1 << logSlots);

    complex<double>* mvec = EvaluatorUtils::randomComplexArray(slots);
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots);
    complex<double>* mvec3 = EvaluatorUtils::randomComplexArray(slots);
    complex<double>* mvec4 = EvaluatorUtils::randomComplexArray(slots);
    
    complex<double>* mvec5 = EvaluatorUtils::randomComplexArray(slots);
    cout << "[" ;
    for(long i = 0; i < slots; ++i){
        cout << mvec[i] << "," ;
        mvec1[i] = mvec[i];
        mvec2[i] = mvec[i];
        mvec3[i] = mvec[i];
        mvec4[i] = mvec[i];
        
        mvec5[i] = mvec[i] * mvec[i];
        
    }
    cout << "]" << endl;
    
    timeutils.start("Encrypt two batch");
    Ciphertext cipher = scheme.encrypt(mvec, slots, L);
    timeutils.stop("Encrypt two batch");
    
#if 0
    Ciphertext multCipher = extscheme.mult(cipher, cipher);
    scheme.reScaleByAndEqual(multCipher, 1);
    
    long RotSlots2 = (1 << logRotSlots);
    Ciphertext multrotCipher = extscheme.leftRotateFast(multCipher, RotSlots2);
    complex<double>* dvec = scheme.decrypt(secretKey, multrotCipher);
    EvaluatorUtils::leftRotateAndEqual(mvec5, slots, RotSlots2);
    StringUtils::showcompare(mvec5, dvec, slots, "mul/rot");
#endif
#if 1
    /**********************************************************/
    timeutils.start("Homomorphic Left Rotation Fast");
    long RotSlots1 = (1 << logRotSlots);
    Ciphertext rotcipher1 = extscheme.leftRotateFast(cipher, RotSlots1);
    timeutils.stop("Homomorphic Left Rotation Fast");
    
    complex<double>* dvec1 = scheme.decrypt(secretKey, rotcipher1);
    
    EvaluatorUtils::leftRotateAndEqual(mvec1, slots, RotSlots1);
    StringUtils::showcompare(mvec1, dvec1, slots, "rot");

 
    timeutils.start("Homomorphic Left Rotation (PoT)");
    rotcipher1 = extscheme.leftRotateByPo2(cipher, logRotSlots);
    timeutils.stop("Homomorphic Left Rotation (PoT)");
    
    dvec1 = scheme.decrypt(secretKey, rotcipher1);
    StringUtils::showcompare(mvec1, dvec1, slots, "rot");
    
    

    timeutils.start("Homomorphic Left Rotation ");
    Ciphertext rotcipher2 = extscheme.leftRotate(cipher, RotSlots);
    timeutils.stop("Homomorphic Left Rotation");
    
    complex<double>* dvec2 = scheme.decrypt(secretKey, rotcipher2);
    
    EvaluatorUtils::leftRotateAndEqual(mvec2, slots, RotSlots);
    StringUtils::showcompare(mvec2, dvec2, slots, "rot");
 
    /**********************************************************/
    timeutils.start("Homomorphic Right Rotation Fast");
    RotSlots1 = (1 << logRotSlots);
    rotcipher1 = extscheme.rightRotateFast(cipher, RotSlots1);
    timeutils.stop("Homomorphic Left Rotation Fast");
    

    dvec1 = scheme.decrypt(secretKey, rotcipher1);
    
    EvaluatorUtils::rightRotateAndEqual(mvec3, slots, RotSlots1);
    StringUtils::showcompare(mvec3, dvec1, slots, "rot");
    
    
    timeutils.start("Homomorphic Right Rotation (PoT)");
    rotcipher1 = extscheme.rightRotateByPo2(cipher, logRotSlots);
    timeutils.stop("Homomorphic Right Rotation (PoT)");
    
    dvec1 = scheme.decrypt(secretKey, rotcipher1);
    
    StringUtils::showcompare(mvec3, dvec1, slots, "rot");
    
    timeutils.start("Homomorphic Right Rotation ");
    rotcipher2 = extscheme.rightRotate(cipher, RotSlots);
    timeutils.stop("Homomorphic Right Rotation");
    
    dvec2 = scheme.decrypt(secretKey, rotcipher2);
    
    EvaluatorUtils::rightRotateAndEqual(mvec4, slots, RotSlots);
    StringUtils::showcompare(mvec4, dvec2, slots, "rot");
#endif
}
 


