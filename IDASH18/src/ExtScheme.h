
#ifndef HEAANNTT_EXTSCHEME_H_
#define HEAANNTT_EXTSCHEME_H_

#include <map>
#include <chrono>

#include "Common.h"
#include "Ciphertext.h"
#include "Context.h"
#include "ExtCiphertext.h"
#include "Plaintext.h"
#include "SecretKey.h"
#include "Key.h"
#include "Numb.h"
#include "Scheme.h"

using namespace std;


//! miran

static long THREEMULTIPLICATION = 3;
//static long DECOMPTWOKey = 4;
//static long DECOMPTHREEKeay = 5;

class ExtScheme {
public:

    Context& context;
    
    Scheme& scheme;

	map<long, Key> keyMap; ///< contain THREEMultiplication

    map<long, Key> decompTwoKeyMap;
    map<long, Key> decompThreeKeyMap;
    
    map<pair<long, long>, Key> decompLeftRotKeyMap;
    //map<long, long, Key> decompRightRotKeyMap;
    
	ExtScheme(SecretKey& secretKey, Context& context, Scheme& scheme);

    /********************************************************************/
    //! add keys for KS (Ps^3 -> s)
    void addThreeMultKey(SecretKey& secretKey);
    
    //! add key for Word-decomp
    void addDecompTwoKeys(SecretKey& secretKey); ///< contain decomposition keys, if generated
    void addDecompThreeKeys(SecretKey& secretKey); ///< contain decomposition keys, if generated
    
    void addDecompLeftRotKey(SecretKey& secretKey, long rot);
    void addDecompLeftRotKeys(SecretKey& secretKey);
    void addDecompRightRotKeys(SecretKey& secretKey);
    
    Ciphertext decompKeySwitch(ExtCiphertext& cipher);
    
    Ciphertext leftRotateFast(Ciphertext& cipher, long rotSlots);
    Ciphertext rightRotateFast(Ciphertext& cipher, long rotSlots);
    
    /********************************************************************/
    //! conversion to a normal ciphertext (b, a)
    ExtCiphertext toExtCipher(Ciphertext& cipher);
    Ciphertext toNormalCipher(ExtCiphertext& extcipher);
    
    // Homomorphic Negation
    ExtCiphertext negate(ExtCiphertext& cipher);
    void negateAndEqual(ExtCiphertext& cipher);
    
    // Homomorphic Addition
    ExtCiphertext add(ExtCiphertext& cipher1, ExtCiphertext& cipher2);
    void addAndEqual(ExtCiphertext& cipher1, ExtCiphertext& cipher2);
    
    // Homomorphic subtraction
    ExtCiphertext sub(ExtCiphertext& cipher1, ExtCiphertext& cipher2);
    void subAndEqual(ExtCiphertext& cipher1, ExtCiphertext& cipher2);
    void sub2AndEqual(ExtCiphertext& cipher1, ExtCiphertext& cipher2);
    
    // Raw multiplication
    ExtCiphertext rawmult(Ciphertext& cipher1, Ciphertext& cipher2);
    ExtCiphertext rawmult(ExtCiphertext& cipher1, Ciphertext& cipher2);
    ExtCiphertext rawsquare(Ciphertext& cipher);
    
    ExtCiphertext rawmult3(Ciphertext& cipher1, Ciphertext& cipher2, Ciphertext& cipher3);
    ExtCiphertext rawmult3_(Ciphertext& cipher1, Ciphertext& cipher2, Ciphertext& cipher3);
    
    //! Homomorphic Key-switching using "scale-up"
    Ciphertext keySwitch(ExtCiphertext& cipher);
    
    /********************************************************************/
    void reScaleByAndEqual(ExtCiphertext& cipher, long dl);
    void reScaleToAndEqual(ExtCiphertext& cipher, long l);
    
    ExtCiphertext modDownBy(ExtCiphertext& cipher, long dl);
    void modDownByAndEqual(ExtCiphertext& cipher, long dl);
    ExtCiphertext modDownTo(ExtCiphertext& cipher, long dl);
    void modDownToAndEqual(ExtCiphertext& cipher, long dl);
 
};

#endif
