

#include "ExtCiphertext.h"


ExtCiphertext::ExtCiphertext() : ax(nullptr), bx(nullptr),  N(0), slots(0), l(0), deg(1) {
    /*
    long deg1 = 1;
    ax = new uint64_t*[deg1];
    for(long i = 0; i < deg1; ++i){
        ax[i] = nullptr;
    }*/
    
}

ExtCiphertext::ExtCiphertext(uint64_t** ax1, uint64_t* bx, long N, long slots, long l, long deg) : bx(bx), N(N), slots(slots), l(l), deg(deg){
    ax = new uint64_t*[deg];
    for(long i = 0; i < deg; ++i){
        ax[i] = ax1[i];
    }
}

ExtCiphertext::ExtCiphertext(const ExtCiphertext& cipher) : N(cipher.N), slots(cipher.slots), l(cipher.l), deg(cipher.deg) {
	bx = new uint64_t[N * l];

    ax = new uint64_t*[cipher.deg];
    for(long j = 0; j < cipher.deg; ++j){
        ax[j] = new uint64_t[N * l];
    }
    
    for (long i = 0; i < N * l; ++i) {
        bx[i] = cipher.bx[i];
        for(long j = 0; j < cipher.deg; ++j){
            ax[j][i] = cipher.ax[j][i];
        }
    }
}

ExtCiphertext& ExtCiphertext::operator=(const ExtCiphertext& o) {
  
    if(this == &o){
        return *this; // handling of self assignment, thanks for your advice, arul.
    }
    delete[] bx;
    delete[] ax;
    N = o.N;
    l = o.l;
    slots = o.slots;
    deg = o.deg;
    
    bx = new uint64_t[N * l];
    
    ax = new uint64_t*[deg];
    for(long j = 0; j < o.deg; ++j){
        ax[j] = new uint64_t[N * l];
    }
    
    for (long i = 0; i < N * l; ++i) {
        bx[i] = o.bx[i];
        for(long j = 0; j < o.deg; ++j){
            ax[j][i] = o.ax[j][i];
        }
    }
    
    return *this;
}

// Move constructor
ExtCiphertext::ExtCiphertext(ExtCiphertext&& cipher) : N(cipher.N), slots(cipher.slots), l(cipher.l), deg(cipher.deg) {
   ax = cipher.ax;              // obtain the pointers from other
   bx = cipher.bx;
   cipher.ax = nullptr;         // move is done so we clear other
   cipher.bx = nullptr;
}

// Move copy assignment operator
ExtCiphertext& ExtCiphertext::operator=(ExtCiphertext&& o) {
    if(this == &o){
        return *this; // handling of self assignment, thanks for your advice, arul.
    }

    delete[] bx;
    delete[] ax;

    ax = o.ax;                  // obtain the pointers from other
    bx = o.bx;
    N = o.N;
    l = o.l;
    slots = o.slots;
    deg = o.deg;
    
    o.ax = nullptr;             // clear the other
    o.bx = nullptr;
    
    return *this;
}

ExtCiphertext::~ExtCiphertext() {
   if(ax) {
      for(long j = 0; j < deg; j++) {
         if(ax[j]) {
            delete [] ax[j];
         }
      }
      delete [] ax;
   }
   if(bx) {
      delete [] bx;
   }
}

 
