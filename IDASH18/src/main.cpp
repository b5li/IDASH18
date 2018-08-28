/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

// g++ -std=c++11 -O2 -I/usr/local/include -pthread main.cpp ../lib/libHEAAN.a -o foo -L/usr/local/lib -lntl -lgmp -lm

#include "ExtTestScheme.h"
#include "TestScheme.h"
#include "Numb.h"

int main() {
    
    
//    ExtTestScheme::testBasic(14, 6, 50, 2);  //! extended ciphertext

//   ExtTestScheme::testThreeProd(15, 6, 50, 1);
 
//    ExtTestScheme::testDecompKS(14, 6, 50, 1);
    
    ExtTestScheme::testDecompRotate(14, 6, 1, 50, 3);
    
    //-----------------------------------------
//	TestScheme::testEncodeSingle(14, 1, 55);

//	TestScheme::testEncodeBatch(15, 6, 55, 3);

//	TestScheme::testBasic(15, 11, 55, 3);

//	TestScheme::testConjugateBatch(15, 6, 55, 1);

//	TestScheme::testRotateByPo2Batch(16, 26, 40, 1, 4, false);

//	TestScheme::testRotateBatch(15, 6, 55, 3, 4, true);

//	TestScheme::testimultBatch(16, 16, 55, 2);

//	TestScheme::testPowerOf2Batch(16, 15, 50, 2, 3);

//	TestScheme::testInverseBatch(14, 5, 55, 4, 3);

//	TestScheme::testExponentBatch(14, 5, 55, 7, 3);

//	TestScheme::testSigmoidBatch(16, 15, 55, 3, 3);

//	TestScheme::testSlotsSum(16, 15, 40, 3);

//	TestScheme::testMeanVariance(14, 3, 55, 13);

//	TestScheme::testHEML("data/uis.txt", 0, 5);

	return 0;
}
