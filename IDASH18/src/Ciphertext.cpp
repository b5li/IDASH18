/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Ciphertext.h"

Ciphertext::Ciphertext() : ax(nullptr), bx(nullptr), N(0), slots(0), l(0) {}

Ciphertext::Ciphertext(uint64_t* _ax, uint64_t* _bx, long _N, long _slots, long _l) : ax(_ax), bx(_bx), N(_N), slots(_slots), l(_l) { }

Ciphertext::Ciphertext(const Ciphertext& cipher) : N(cipher.N), slots(cipher.slots), l(cipher.l) {
	ax = new uint64_t[N * l];
	bx = new uint64_t[N * l];
	for (long i = 0; i < N * l; ++i) {
		ax[i] = cipher.ax[i];
		bx[i] = cipher.bx[i];
	}
}

Ciphertext& Ciphertext::operator=(const Ciphertext& o) {
	if(this == &o) return *this; // handling of self assignment, thanks for your advice, arul.
	delete[] ax;
	delete[] bx;
	N = o.N;
	l = o.l;
	slots = o.slots;
	ax = new uint64_t[N * l];
	bx = new uint64_t[N * l];
	for (long i = 0; i < N * l; ++i) {
		ax[i] = o.ax[i];
		bx[i] = o.bx[i];
	}
	return *this;
}


Ciphertext::Ciphertext(Ciphertext&& cipher) : N(cipher.N), slots(cipher.slots), l(cipher.l) {
  ax = cipher.ax;               // obtain the array pointer from other
  bx = cipher.bx;
  cipher.ax = nullptr;          // move is done so we clear the other
  cipher.bx = nullptr;
}

Ciphertext& Ciphertext::operator=(Ciphertext&& o) {
  if(this == &o) return *this;
  delete[] ax;
  delete[] bx;

  ax = o.ax;                    // obtain the pointers and other fields from other
  bx = o.bx;
  N = o.N;
  l = o.l;
  slots = o.slots;

  o.ax = nullptr;               // clear the other
  o.bx = nullptr;

  return *this;
}

  

Ciphertext::~Ciphertext() {
   if(ax) {
      delete [] ax;
   }
   if(bx) {
      delete [] bx;
   }
}
