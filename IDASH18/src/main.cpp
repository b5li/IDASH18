

// g++ -std=c++11 -O2 -I/usr/local/include -pthread main.cpp ../lib/libHEAAN.a -o foo -L/usr/local/lib -lntl -lgmp -lm

#include "ExtTestScheme.h"
#include "TestScheme.h"
#include "Numb.h"

int main() {
    
    
//    ExtTestScheme::testBasic(14, 6, 50, 2);  //! extended ciphertext

//    ExtTestScheme::testThreeProd(15, 6, 35, 1);
 
//    ExtTestScheme::testDecompKS(15, 3, 40, 1);
    
//    ExtTestScheme::testProdDecompKS(13, 5, 50, 0);
    
   ExtTestScheme::testDecompRotate(14, 4, 1, 3, 40, 3);   // rot by (1 << 1), rot by 3
    

	return 0;
}
