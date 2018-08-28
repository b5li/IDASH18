

#ifndef HEAANNTT_EXTTESTSCHEME_H_
#define HEAANNTT_EXTTESTSCHEME_H_

#include <iostream>

using namespace std;

class ExtTestScheme {
public:
	
    //! miran
    static void testBasic(long logN,  long L, long logp, long logSlots);

    static void testThreeProd(long logN,  long L, long logp, long logSlots);
    
    //! KS
    static void testDecompKS(long logN,  long L, long logp, long logSlots);
    
    static void testDecompRotate(long logN,  long L, long logRotSlots, long logp, long logSlots);
    
};

#endif /* TESTSCHEME_H_ */
