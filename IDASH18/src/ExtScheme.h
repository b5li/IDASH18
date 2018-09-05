/*
 * @file       ExtScheme.h, header file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim, Yongsoo Song
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */


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


static long THREEMULTIPLICATION = 0;


class ExtScheme {
public:

    Context& context;
    
    Scheme& scheme;

	map<long, Key> keyMap; ///< contain THREEMultiplication
    map<long, Key> decompTwoKeyMap; ///< contain DecompTwoMultiplication
    map<long, Key> decompThreeKeyMap; ///< contain DecompThreeMultiplication
    
    map<pair<long, long>, Key> decompLeftRotKeyMap;
    
	ExtScheme(SecretKey& secretKey, Context& context, Scheme& scheme);

    uint64_t* p0InvModqi;
    
    /********************************************************************/
    void addThreeMultKey(SecretKey& secretKey); ///< contain scale up keys, if generated
    
    void addDecompTwoKeys(SecretKey& secretKey); ///< contain decomposition keys, if generated
    void addDecompThreeKeys(SecretKey& secretKey); ///< contain decomposition keys, if generated
    
    void addDecompLeftRotKey(SecretKey& secretKey, long rot);
    void addDecompLeftRotKeys(SecretKey& secretKey); ///< contain all the decomposition keys of rotations, if generated
    
    void addDecompRightRotKey(SecretKey& secretKey, long rot);
    void addDecompRightRotKeys(SecretKey& secretKey);
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
    
    /********************************************************************/
    //! Homomorphic Key-switching using "modulus raising "
    Ciphertext ModRaiseKeySwitch(ExtCiphertext& cipher);
    
    //! Homomorphic Key-swithcing using "decomposition"
    Ciphertext DecompKeySwitch(ExtCiphertext& cipher);
    
    void rnsDecomp(uint64_t*& res, uint64_t* a, long i, long l);
    void mulDecompKey(uint64_t*& axtmp, uint64_t*& bxtmp, uint64_t* axi, Key keys, long l);
    void modDownByp0(uint64_t*& a, long l);
    
    //! Rotation using "Decomposition Key Switching:
    Ciphertext leftRotateFast(Ciphertext& cipher, long rotSlots);
    Ciphertext leftRotateByPo2(Ciphertext& cipher, long logRotSlots);
    Ciphertext leftRotate(Ciphertext& cipher, long rotSlots);
    
    Ciphertext rightRotateFast(Ciphertext& cipher, long rotSlots);
    Ciphertext rightRotateByPo2(Ciphertext& cipher, long logRotSlots);
    Ciphertext rightRotate(Ciphertext& cipher, long logRotSlots);
    
    /********************************************************************/
    void reScaleByAndEqual(ExtCiphertext& cipher, long dl);
    void reScaleToAndEqual(ExtCiphertext& cipher, long l);
    
    ExtCiphertext modDownBy(ExtCiphertext& cipher, long dl);
    void modDownByAndEqual(ExtCiphertext& cipher, long dl);
    ExtCiphertext modDownTo(ExtCiphertext& cipher, long dl);
    void modDownToAndEqual(ExtCiphertext& cipher, long dl);
 
    /********************************************************************/
    //! Raw multiplication + decomposition KS
    Ciphertext square(Ciphertext& cipher);
    void squareAndEqual(Ciphertext& cipher);
    
    Ciphertext mult(Ciphertext& cipher1, Ciphertext& cipher2);
    void multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);
};

#endif
