

#ifndef HEAANNTT_EXTCIPHERTEXT_H_
#define HEAANNTT_EXTCIPHERTEXT_H_

#include "Common.h"

class ExtCiphertext {

public:
    
	long N; ///< Dimension of Ring

	long slots; ///< The length of plaintext vector

	long l; ///< The level of this ciphertext

    long deg; //! The degree of ciphertext (1,s) : degree 1, (1,s,s2): degree 2
    
    uint64_t* bx;  // 1
    
    uint64_t** ax;
    
	// Default constructor
	ExtCiphertext();

	// Constructor
    ExtCiphertext(uint64_t** ax1, uint64_t* bx, long N, long slots, long l, long deg);
	//ExtCiphertext(uint64_t* ax3, uint64_t* ax2, uint64_t* ax1, uint64_t* bx, long N, long slots, long l);

	// Copy constructor
	ExtCiphertext(const ExtCiphertext& cipher);
	ExtCiphertext& operator=(const ExtCiphertext &o);

    // Move constructor/assignment operator
    ExtCiphertext(ExtCiphertext&& cipher);
	ExtCiphertext& operator=(ExtCiphertext&& o);

    // Destructor
    ~ExtCiphertext();
	
};

#endif
