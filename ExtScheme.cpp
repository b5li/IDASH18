

#include "ExtScheme.h"

//ExtScheme::ExtScheme(Context& context) : context(context) {}

ExtScheme::ExtScheme(SecretKey& secretKey, Context& context, Scheme& scheme) : context(context), scheme(scheme) {
    //! (p0^-1 mod ql)
    p0InvModqi = new uint64_t[context.L]();
    for(long i = 0; i < context.L; ++i){
        p0InvModqi[i] = invMod(context.pVec[0], context.qVec[i]);
    }
}

/**************************************************************************************/

//!@ Generate the evaluation key = Enc(P * sk^3)
void ExtScheme::addThreeMultKey(SecretKey& secretKey) {
    uint64_t* ex = new uint64_t[(context.L + context.K) << context.logN]();
    uint64_t* ax = new uint64_t[(context.L + context.K) << context.logN]();
    uint64_t* bx = new uint64_t[(context.L + context.K) << context.logN]();
    uint64_t* sx3 = new uint64_t[(context.L + context.K) << context.logN]();
    
    context.mul(sx3, secretKey.sx, secretKey.sx, context.L);
    context.mulAndEqual(sx3, secretKey.sx, context.L);
    
    context.evalAndEqual(sx3, context.L);
    
    context.sampleGauss(ex, context.L, context.K);
    context.NTTAndEqual(ex, context.L, context.K);
    
    context.addAndEqual(ex, sx3, context.L);
    
    context.sampleUniform(ax, context.L, context.K);
    context.mul(bx, ax, secretKey.sx, context.L, context.K);
    context.sub2AndEqual(ex, bx, context.L, context.K);
    
    delete[] ex;
    delete[] sx3;
    
    keyMap.insert(pair<long, Key>(THREEMULTIPLICATION, Key(ax, bx)));
}

//!@ Generate the evaluation key = Enc(pVec[0] * RNS(Ei) * s^2)
void ExtScheme::addDecompTwoKeys(SecretKey& secretKey){
    uint64_t** ax = new uint64_t*[context.L];
    uint64_t** bx = new uint64_t*[context.L];
    uint64_t* p0sx2 = new uint64_t[context.L << context.logN];
    
    //! Compute pVec[0] * s^2
    context.mul(p0sx2, secretKey.sx, secretKey.sx, context.L);
    context.mulConstAndEqual(p0sx2, context.pVec[0], context.L);
   
    //! SwitchKeyGen Enc(pVec[0] * RNS(Ei) * s^2)
    uint64_t* ex = new uint64_t[(context.L + 1) << context.logN];
    
    for(long i = 0; i < context.L; ++i){
        ax[i] = new uint64_t[(context.L + 1) << context.logN];
        bx[i] = new uint64_t[(context.L + 1) << context.logN];
        
        context.sampleGauss(ex, context.L, 1);
        context.NTTAndEqual(ex, context.L, 1);
        
        //! Add pVec[0] * s^2 to i-th poly
        context.qiAddAndEqual(ex + (i << context.logN), p0sx2 + (i << context.logN), i);
        
        context.sampleUniform(ax[i], context.L, 1);
        context.mul(bx[i], ax[i], secretKey.sx, context.L, 1);
        context.sub2AndEqual(ex, bx[i], context.L, 1);
        
        decompTwoKeyMap.insert(pair<long, Key>(i, Key(ax[i], bx[i])));
    }
    delete[] ex;
    delete[] p0sx2;
}

//!@ Generate the evaluation key = Enc(pVec[0] * RNS(Ei) * s^3)
void ExtScheme::addDecompThreeKeys(SecretKey& secretKey){
    uint64_t** ax = new uint64_t*[context.L];
    uint64_t** bx = new uint64_t*[context.L];
 
    uint64_t* p0sx3 = new uint64_t[context.L << context.logN];
    
    // Compute pVec[0] * s^3
    context.mul(p0sx3, secretKey.sx, secretKey.sx, context.L);
    context.mulAndEqual(p0sx3, secretKey.sx, context.L);
    context.mulConstAndEqual(p0sx3, context.pVec[0], context.L);
 
    // SwitchKeyGen Enc(pVec[0] * RNS(Ei) * s^3)
    uint64_t* ex = new uint64_t[(context.L + 1) << context.logN];
 
    for(long i = 0; i < context.L; ++i){
        ax[i] = new uint64_t[(context.L + 1) << context.logN];
        bx[i] = new uint64_t[(context.L + 1) << context.logN];
 
        context.sampleGauss(ex, context.L, 1);
        context.NTTAndEqual(ex, context.L, 1);
        
        // Add pVec[0] * s^2 to i-th poly
        context.qiAddAndEqual(ex + (i << context.logN), p0sx3 + (i << context.logN), i);
 
        context.sampleUniform(ax[i], context.L, 1);
        context.mul(bx[i], ax[i], secretKey.sx, context.L, 1);
        context.sub2AndEqual(ex, bx[i], context.L, 1);
 
        decompThreeKeyMap.insert(pair<long, Key>(i, Key(ax[i], bx[i])));
    }
    delete[] ex;
    delete[] p0sx3;
}

//!@ Generate the evaluation key = Enc(pVec[0] * RNS(Ei) * s^3)
void ExtScheme::addDecompLeftRotKey(SecretKey& secretKey, long rot){
    uint64_t** ax = new uint64_t*[context.L];
    uint64_t** bx = new uint64_t*[context.L];
    
    uint64_t* p0sxrot = new uint64_t[context.L << context.logN];
    
    //! Compute pVec[0] * sx(x^rot)
    context.leftRot(p0sxrot, secretKey.sx, context.L, rot);
    context.mulConstAndEqual(p0sxrot, context.pVec[0], context.L);
    
    //! SwitchKeyGen Enc(pVec[0] * RNS(0,..1,..,0) * s(x^rot))
    uint64_t* ex = new uint64_t[(context.L + 1) << context.logN];
    
    for(long i = 0; i < context.L; ++i){
        ax[i] = new uint64_t[(context.L + 1) << context.logN];
        bx[i] = new uint64_t[(context.L + 1) << context.logN];
        
        context.sampleGauss(ex, context.L, 1);
        context.NTTAndEqual(ex, context.L, 1);
        
        //! Add pVec[0] * s^2 to i-th poly
        context.qiAddAndEqual(ex + (i << context.logN), p0sxrot + (i << context.logN), i);
        
        context.sampleUniform(ax[i], context.L, 1);
        context.mul(bx[i], ax[i], secretKey.sx, context.L, 1);
        context.sub2AndEqual(ex, bx[i], context.L, 1);
        
        decompLeftRotKeyMap.insert(pair< pair<long, long>, Key> (make_pair(rot, i), Key(ax[i], bx[i])));
    }
    delete[] ex;
    delete[] p0sxrot;
}

//!@ Generate all the left rotation key (PoT rotations)
void ExtScheme::addDecompLeftRotKeys(SecretKey& secretKey){
    for (long i = 0; i < context.logNh; ++i) {
        long idx = 1 << i;
        if(decompLeftRotKeyMap.find(make_pair(idx,0)) == decompLeftRotKeyMap.end()) {
            addDecompLeftRotKey(secretKey, idx);
        }
    }
}


//!@ Generate all the right rotation key
void ExtScheme::addDecompRightRotKeys(SecretKey& secretKey){
    for (long i = 0; i < context.logNh; ++i) {
        long idx = context.Nh - (1 << i);
        if(decompLeftRotKeyMap.find(make_pair(idx,0)) == decompLeftRotKeyMap.end()) {
            addDecompLeftRotKey(secretKey, idx);
        }
    }
}

//!@ Generate the right rotation key by "rot" positions
void ExtScheme::addDecompRightRotKey(SecretKey& secretKey, long rot){
    long idx = context.Nh - rot;
    if(decompLeftRotKeyMap.find(make_pair(idx,0)) == decompLeftRotKeyMap.end()) {
        addDecompLeftRotKey(secretKey, idx);
    }
    //cout << "rot:" << rot << "," << "idx: " << idx << endl; 
}

/**************************************************************************************/
//!@ convert to cipher to extended
//!@  bx + ax * s  -> bx + ax * s  + 0 * s^2 + ...
ExtCiphertext ExtScheme::toExtCipher(Ciphertext& cipher){
    long deg = 1;
    uint64_t** ax = new uint64_t*[deg];
    ax[0] = cipher.ax;
    return ExtCiphertext(ax, cipher.bx, context.N, cipher.slots, cipher.l, deg);
}


//!@ suppose (b, a[0], 0,,,,0)
Ciphertext ExtScheme::toNormalCipher(ExtCiphertext& extcipher){
    return Ciphertext(extcipher.ax[0], extcipher.bx, context.N, extcipher.slots, extcipher.l);
}

/**************************************************************************************/
//!@ Basic functions on extended ciphertexts

ExtCiphertext ExtScheme::negate(ExtCiphertext& cipher) {
    uint64_t* bxres = new uint64_t[cipher.l << context.logN];
    uint64_t** axres = new uint64_t*[cipher.deg];

    context.negate(bxres, cipher.bx, cipher.l);
    for(long j = 0; j < cipher.deg; ++j){
        axres[j] = new uint64_t[cipher.l << context.logN];
        context.negate(axres[j], cipher.ax[j], cipher.l);
    }

    return ExtCiphertext(axres, bxres, context.N, cipher.slots, cipher.l, cipher.deg);
}

void ExtScheme::negateAndEqual(ExtCiphertext& cipher) {
    long shift = 0;
    for (long i = 0; i < cipher.l; ++i) {
        for(long j = 0; j < cipher.deg; ++j){
            context.qiNegateAndEqual(cipher.ax[j] + shift, i);
        }
        context.qiNegateAndEqual(cipher.bx + shift, i);
        shift += context.N;
    }
}

ExtCiphertext ExtScheme::add(ExtCiphertext& cipher1, ExtCiphertext& cipher2) {
    if((cipher1.l != cipher2.l)||(cipher1.deg != cipher2.deg)) {
        throw invalid_argument("Ciphertexts are not comparable");
    }
    uint64_t* bxres = new uint64_t[cipher1.l << context.logN];
    uint64_t** axres = new uint64_t*[cipher1.deg];
    
    context.add(bxres, cipher1.bx, cipher2.bx, cipher1.l);
    for(long j = 0; j < cipher1.deg; ++j){
        axres[j] = new uint64_t[cipher1.l << context.logN];
        context.add(axres[j], cipher1.ax[j], cipher2.ax[j], cipher1.l);
    }

    return ExtCiphertext(axres, bxres, context.N, cipher1.slots, cipher1.l, cipher1.deg);
}

void ExtScheme::addAndEqual(ExtCiphertext& cipher1, ExtCiphertext& cipher2) {
    if((cipher1.l != cipher2.l)||(cipher1.deg != cipher2.deg)) {
        throw invalid_argument("Ciphertexts are not comparable");
    }
    
    context.addAndEqual(cipher1.bx, cipher2.bx, cipher1.l);
    for(long j = 0; j < cipher1.deg; ++j){
        context.addAndEqual(cipher1.ax[j], cipher2.ax[j], cipher1.l);
    }
}

ExtCiphertext ExtScheme::sub(ExtCiphertext& cipher1, ExtCiphertext& cipher2) {
    if((cipher1.l != cipher2.l)||(cipher1.deg != cipher2.deg)) {
        throw invalid_argument("Ciphertexts are not comparable");
    }
    
    uint64_t* bxres = new uint64_t[cipher1.l << context.logN];
    uint64_t** axres = new uint64_t*[cipher1.deg];
    
    context.sub(bxres, cipher1.bx, cipher2.bx, cipher1.l);
    for(long j = 0; j < cipher1.deg; ++j){
        axres[j] = new uint64_t[cipher1.l << context.logN];
        context.sub(axres[j], cipher1.ax[j], cipher2.ax[j], cipher1.l);
    }
    return ExtCiphertext(axres, bxres, context.N, cipher1.slots, cipher1.l, cipher1.deg);
}

void ExtScheme::subAndEqual(ExtCiphertext& cipher1, ExtCiphertext& cipher2) {
    if((cipher1.l != cipher2.l)||(cipher1.deg != cipher2.deg)) {
        throw invalid_argument("Ciphertexts are not comparable");
    }
    
    context.subAndEqual(cipher1.bx, cipher2.bx, cipher1.l);
    for(long j = 0; j < cipher1.deg; ++j){
        context.subAndEqual(cipher1.ax[j], cipher2.ax[j], cipher1.l);
    }
}

void ExtScheme::sub2AndEqual(ExtCiphertext& cipher1, ExtCiphertext& cipher2) {
    if((cipher1.l != cipher2.l)||(cipher1.deg != cipher2.deg)) {
        throw invalid_argument("Ciphertexts are not comparable");
    }
    
    context.sub2AndEqual(cipher1.bx, cipher2.bx, cipher1.l);
    for(long j = 0; j < cipher1.deg; ++j){
        context.sub2AndEqual(cipher1.ax[j], cipher2.ax[j], cipher1.l);
    }
}

void ExtScheme::reScaleByAndEqual(ExtCiphertext& cipher, long dl) {
    for (long i = 0; i < dl; ++i) {
        context.reScaleAndEqual(cipher.bx, cipher.l);
        for(long j = 0; j < cipher.deg; ++j){
            context.reScaleAndEqual(cipher.ax[j], cipher.l);
        }
        cipher.l -= 1;
    }
}

void ExtScheme::reScaleToAndEqual(ExtCiphertext& cipher, long l) {
    long dl = cipher.l - l;
    reScaleByAndEqual(cipher, dl);
}


ExtCiphertext ExtScheme::modDownBy(ExtCiphertext& cipher, long dl) {
    uint64_t** ax =  new uint64_t*[cipher.deg];
    
    uint64_t* bx = context.modDown(cipher.bx, cipher.l, dl);
    for(long j = 0; j < cipher.deg; ++j){
        ax[j] = context.modDown(cipher.ax[j], cipher.l, dl);
    }
    return ExtCiphertext(ax, bx, context.N, cipher.slots, cipher.l - dl, cipher.deg);
}

void ExtScheme::modDownByAndEqual(ExtCiphertext& cipher, long dl) {
    context.modDownAndEqual(cipher.bx, cipher.l, dl);
    for(long j = 0; j < cipher.deg; ++j){
        context.modDownAndEqual(cipher.ax[j], cipher.l, dl);
    }
    cipher.l -= dl;
}

ExtCiphertext ExtScheme::modDownTo(ExtCiphertext& cipher, long l) {
    long dl = cipher.l - l;
    return modDownBy(cipher, dl);
}

void ExtScheme::modDownToAndEqual(ExtCiphertext& cipher, long l) {
    long dl = cipher.l - l;
    modDownByAndEqual(cipher, dl);
}

/**************************************************************************************/
//!@ Multiplication

//!@ expand (b1 + a1 * s) * (b2 + a2 * s)
//!@ 3M + 4A
ExtCiphertext ExtScheme::rawmult(Ciphertext& cipher1, Ciphertext& cipher2){

    uint64_t* axbx1 = new uint64_t[cipher1.l << context.logN]();
    uint64_t* axbx2 = new uint64_t[cipher1.l << context.logN]();
    
    uint64_t* axax = new uint64_t[cipher1.l << context.logN]();
    uint64_t* bxbx = new uint64_t[cipher1.l << context.logN]();
    
    context.mul(bxbx, cipher1.bx, cipher2.bx, cipher1.l);   //! b1 * b2
    context.mul(axax, cipher1.ax, cipher2.ax, cipher1.l);   //! a1 * a2
    
    context.add(axbx1, cipher1.ax, cipher1.bx, cipher1.l);  //! (a1 + b1)
    context.add(axbx2, cipher2.ax, cipher2.bx, cipher2.l);  //! (a2 + b2)
    context.mulAndEqual(axbx1, axbx2, cipher1.l);           //! (a1 + b1)*(a2 + b2)
    
    context.subAndEqual(axbx1, bxbx, cipher1.l);
    context.subAndEqual(axbx1, axax, cipher1.l);
    
    long deg = 2;
    uint64_t** extax = new uint64_t*[deg];
    for(long j = 0; j < deg; ++j){
        extax[j] = new uint64_t[context.N * cipher1.l];
    }
    extax[0] = axbx1;
    extax[1] = axax;

    return ExtCiphertext(extax, bxbx, context.N, cipher1.slots, cipher1.l, deg);
}


ExtCiphertext ExtScheme::rawsquare(Ciphertext& cipher){
    
    uint64_t* axbx1 = new uint64_t[cipher.l << context.logN]();

    uint64_t* axax = new uint64_t[cipher.l << context.logN]();
    uint64_t* bxbx = new uint64_t[cipher.l << context.logN]();
    
    context.square(bxbx, cipher.bx, cipher.l);   //! b1 * b2
    context.square(axax, cipher.ax, cipher.l);   //! a1 * a2
    
    context.add(axbx1, cipher.ax, cipher.bx, cipher.l);  //! (a1 + b1)
    context.squareAndEqual(axbx1, cipher.l);             //! (a1 + b1)*(a2 + b2) = ()^2
    
    context.subAndEqual(axbx1, bxbx, cipher.l);
    context.subAndEqual(axbx1, axax, cipher.l);
    
    long deg = 2;
    uint64_t** extax = new uint64_t*[deg];
    for(long j = 0; j < deg; ++j){
        extax[j] = new uint64_t[context.N * cipher.l];
    }
    extax[0] = axbx1;
    extax[1] = axax;

    return ExtCiphertext(extax, bxbx, context.N, cipher.slots, cipher.l, deg);
}

//!@ expand (b1 + a1[0] * s + a1[1] * s^2)  * (b2 + a2 * s)
//!@ 5M + 5A
ExtCiphertext ExtScheme::rawmult(ExtCiphertext& cipher1, Ciphertext& cipher2){
    
    uint64_t* axbx1 = new uint64_t[cipher1.l << context.logN]();
    uint64_t* axbx2 = new uint64_t[cipher1.l << context.logN]();
    
    uint64_t* axax = new uint64_t[cipher1.l << context.logN]();
    uint64_t* bxbx = new uint64_t[cipher1.l << context.logN]();
    
    uint64_t* ax2 = new uint64_t[cipher1.l << context.logN]();
    uint64_t* ax3 = new uint64_t[cipher1.l << context.logN]();
    
    context.mul(bxbx, cipher1.bx, cipher2.bx, cipher1.l);       //! -> bx
    context.mul(axax, cipher1.ax[0], cipher2.ax, cipher1.l);    //! -> (ax[1] + (*)) s2
    
    context.add(axbx1, cipher1.ax[0], cipher1.bx, cipher1.l);
    context.add(axbx2, cipher2.ax, cipher2.bx, cipher2.l);
    context.mulAndEqual(axbx1, axbx2, cipher1.l);
    context.subAndEqual(axbx1, bxbx, cipher1.l);
    context.subAndEqual(axbx1, axax, cipher1.l);              //! -> (ax[0]) * s
    
    context.mul(ax2, cipher1.ax[1], cipher2.bx, cipher1.l);   //! (ax[1] + a1[1] * b2) s2
    context.addAndEqual(ax2, axax, cipher1.l);
    
    context.mul(ax3, cipher1.ax[1], cipher2.ax, cipher1.l);   //! (a1[1] * a2) s3
    
    long deg = 3;
    uint64_t** extax = new uint64_t*[deg];
    for(long j = 0; j < deg; ++j){
        extax[j] = new uint64_t[context.N * cipher1.l];
    }
    extax[0] = axbx1;
    extax[1] = ax2;
    extax[2] = ax3;
    
    return ExtCiphertext(extax, bxbx, context.N, cipher1.slots, cipher1.l, deg);
}


//!@ expand (b1 + a1 * s) * (b2 + a2 * s)
//!@ 9M, 6A
ExtCiphertext ExtScheme::rawmult3(Ciphertext& cipher1, Ciphertext& cipher2, Ciphertext& cipher3) {
    uint64_t* bx = new uint64_t[cipher1.l << context.logN]();
    uint64_t* ax1 = new uint64_t[cipher1.l << context.logN]();  //! * s
    uint64_t* ax2 = new uint64_t[cipher1.l << context.logN]();  //! * s2
    uint64_t* ax3 = new uint64_t[cipher1.l << context.logN]();  //! * s3
    
    uint64_t* ax12 = new uint64_t[cipher1.l << context.logN]();
    uint64_t* bx12 = new uint64_t[cipher1.l << context.logN]();
    
    uint64_t* axbx1 = new uint64_t[cipher1.l << context.logN]();
    uint64_t* axbx2 = new uint64_t[cipher1.l << context.logN]();
    
    //! bx = (b1 * b2) * b3
    context.mul(bx12, cipher1.bx, cipher2.bx, cipher2.l);
    context.mul(bx, bx12, cipher3.bx, cipher1.l);
    
    //! ax3 = (a1 * a2) * a3
    context.mul(ax12, cipher1.ax, cipher2.ax, cipher2.l);
    context.mul(ax3, ax12, cipher3.ax, cipher1.l);
    
    //! abbx1 = a1 * b2 + a2 * b1 = (a1 + b1) * (a2 + b2)
    context.add(axbx1, cipher1.ax, cipher1.bx, cipher1.l);
    context.add(axbx2, cipher2.ax, cipher2.bx, cipher2.l);
    context.mulAndEqual(axbx1, axbx2, cipher1.l);
    context.subAndEqual(axbx1, ax12, cipher1.l);
    context.subAndEqual(axbx1, bx12, cipher1.l);
    
    //! ax1 = (axbx1) * b3 + (b1 * b2) * a3
    context.mul(ax1, axbx1, cipher3.bx, cipher1.l);
    context.mulAndEqual(bx12, cipher3.ax, cipher1.l);
    context.addAndEqual(ax1, bx12, cipher1.l);
    
    //! ax2 = (axbx1) * a3 + (a1 * a2) * b3
    context.mul(ax2, axbx1, cipher3.ax, cipher1.l);
    context.mulAndEqual(ax12, cipher3.bx, cipher1.l);
    context.addAndEqual(ax2, ax12, cipher1.l);
    
    long deg = 3;
    uint64_t** extax = new uint64_t*[deg];
    for(long j = 0; j < deg; ++j){
        extax[j] = new uint64_t[context.N * cipher1.l];
    }
    extax[0] = ax1;
    extax[1] = ax2;
    extax[2] = ax3;
    
    return ExtCiphertext(extax, bx, context.N, cipher1.slots, cipher1.l, deg);
}


//!@ expand (b1 + a1 * s) * (b2 + a2 * s):
//!@ 8M, 8A
ExtCiphertext ExtScheme::rawmult3_(Ciphertext& cipher1, Ciphertext& cipher2, Ciphertext& cipher3) {
    uint64_t* bx = new uint64_t[cipher1.l << context.logN]();
    uint64_t* ax1 = new uint64_t[cipher1.l << context.logN]();  //! * s
    uint64_t* ax2 = new uint64_t[cipher1.l << context.logN]();  //! * s2
    uint64_t* ax3 = new uint64_t[cipher1.l << context.logN]();  //! * s3
    
    uint64_t* ax12 = new uint64_t[cipher1.l << context.logN]();
    uint64_t* bx12 = new uint64_t[cipher1.l << context.logN]();
    
    uint64_t* axbx1 = new uint64_t[cipher1.l << context.logN]();
    uint64_t* axbx2 = new uint64_t[cipher1.l << context.logN]();
    
    //! bx = (b1 * b2) * b3)
    context.mul(bx12, cipher1.bx, cipher2.bx, cipher2.l);
    context.mul(bx, bx12, cipher3.bx, cipher1.l);
    
    //! ax3 = (a1 * a2) * a3
    context.mul(ax12, cipher1.ax, cipher2.ax, cipher2.l);
    context.mul(ax3, ax12, cipher3.ax, cipher1.l);
   
    //! abbx1 = a1 * b2 + a2 * b1
    context.add(axbx1, cipher1.ax, cipher1.bx, cipher1.l);
    context.add(axbx2, cipher2.ax, cipher2.bx, cipher2.l);
    context.mulAndEqual(axbx1, axbx2, cipher1.l);
    context.sub(ax1, axbx1, ax12, cipher1.l);
    context.sub(axbx1, ax1, bx12, cipher1.l);
    
    //! note: (c0, c1, c2) = (b1 * b2, axbx1, a1 * a2) = (bx12, axbx1, ax12)
    //! ax2 = c1 * a3
    context.mul(ax2, axbx1, cipher3.ax, cipher1.l);
    
    //! ax1 = c0 * a3  + c1 * b3
    //!      = (c0 + c1) * (a3 + b3) - (c0 * b3) - (c1 * a3)
    context.add(axbx2, cipher3.ax, cipher3.bx, cipher1.l); //! a3 + b3
    context.mulAndEqual(ax1, axbx2, cipher1.l);
    context.subAndEqual(ax1, bx, cipher1.l);
    context.subAndEqual(ax1, ax2, cipher1.l);
    
    //! ax2 = c1 * a3 + c2 (=a1 * a2) * b3
    context.mul(axbx2, ax12, cipher3.bx, cipher1.l);
    context.addAndEqual(ax2, axbx2, cipher1.l); 
 
    long deg = 3;
    uint64_t** extax = new uint64_t*[deg];
    for(long j = 0; j < deg; ++j){
        extax[j] = new uint64_t[context.N * cipher1.l];
    }
    extax[0] = ax1;
    extax[1] = ax2;
    extax[2] = ax3;
    
    return ExtCiphertext(extax, bx, context.N, cipher1.slots, cipher1.l, deg);
}


/**************************************************************************************/
//! Key Switching Operations
//! (bx, ax[0]) + (bxmult2, axmult2) + (bxmult3, axmult3)

Ciphertext ExtScheme::ModRaiseKeySwitch(ExtCiphertext& cipher) {
    uint64_t* axmult2 = new uint64_t[(cipher.l + context.K) << context.logN]();
    uint64_t* bxmult2 = new uint64_t[(cipher.l + context.K) << context.logN]();
    uint64_t* axmult3 = new uint64_t[(cipher.l + context.K) << context.logN]();
    uint64_t* bxmult3 = new uint64_t[(cipher.l + context.K) << context.logN]();
    
    uint64_t* axres = new uint64_t[cipher.l << context.logN];
    uint64_t* bxres = new uint64_t[cipher.l << context.logN];
    
    //! s^2 -> s
    if(cipher.deg > 1){
        context.raiseAndEqual(cipher.ax[1], cipher.l);
        Key keys2 = scheme.keyMap.at(MULTIPLICATION);  //! generated from scheme.cpp
        
        context.mulKey(axmult2, cipher.ax[1], keys2.ax, cipher.l);
        context.mulKey(bxmult2, cipher.ax[1], keys2.bx, cipher.l);
        
        context.backAndEqual(axmult2, cipher.l);
        context.backAndEqual(bxmult2, cipher.l);
        
        context.add(axres, cipher.ax[0], axmult2, cipher.l);
        context.add(bxres, cipher.bx, bxmult2, cipher.l);
    }
    
    //! s^3 -> s
    if(cipher.deg > 2){
        context.raiseAndEqual(cipher.ax[2], cipher.l);
        Key keys3 = keyMap.at(THREEMULTIPLICATION);
        
        context.mulKey(axmult3, cipher.ax[2], keys3.ax, cipher.l);
        context.mulKey(bxmult3, cipher.ax[2], keys3.bx, cipher.l);
        
        context.backAndEqual(axmult3, cipher.l);
        context.backAndEqual(bxmult3, cipher.l);
        
        context.addAndEqual(axres, axmult3, cipher.l);
        context.addAndEqual(bxres, bxmult3, cipher.l);
    }

    delete[] axmult2;
    delete[] axmult3;
    delete[] bxmult2;
    delete[] bxmult3;
    
    return Ciphertext(axres, bxres, context.N, cipher.slots, cipher.l);
    
}

//!@ Functions for DecompKS 
void ExtScheme::rnsDecomp(uint64_t*& res, uint64_t* a, long i, long l){
    uint64_t* ai = new uint64_t[context.N];
    copy(a + (i << context.logN), a + ((i + 1) << context.logN), ai);
    context.qiINTTAndEqual(ai, i);
    
    for(long j = 0; j < l; ++j){
        uint64_t* resj = res + (j << context.logN);
        context.qiNTT(resj, ai, j);
    }
    context.piNTT(res + (l << context.logN), ai, 0);
    delete[] ai;
}

void ExtScheme::mulDecompKey(uint64_t*& axtmp, uint64_t*& bxtmp, uint64_t* axi, Key keys, long l){
    context.mul(axtmp, axi, keys.ax, l);
    context.mul(bxtmp, axi, keys.bx, l);
    long lN = (l << context.logN);
    context.piMul(axtmp + lN, axi + lN, keys.ax + (context.L << context.logN), 0);
    context.piMul(bxtmp + lN, axi + lN, keys.bx + (context.L << context.logN), 0);
}

void ExtScheme::modDownByp0(uint64_t*& a, long l){
    uint64_t* tmp = new uint64_t[context.N];
    uint64_t* al = a + (l << context.logN);
    context.piINTTAndEqual(al, 0);
    
    for(long i = 0; i < l; ++i){
        uint64_t* ai = a + (i << context.logN);
        context.qiNTT(tmp, al, i);
        context.qiSubAndEqual(ai, tmp, i);
        context.qiMulConstAndEqual(ai, p0InvModqi[i], i);
    }
    delete[] tmp;
}

Ciphertext ExtScheme::DecompKeySwitch(ExtCiphertext& cipher) {
    uint64_t* axi = new uint64_t[(cipher.l + 1) << context.logN];
    uint64_t* axtmp = new uint64_t[(cipher.l + 1) << context.logN];
    uint64_t* bxtmp = new uint64_t[(cipher.l + 1) << context.logN];
    uint64_t* axmult2 = new uint64_t[(cipher.l + 1) << context.logN]();
    uint64_t* bxmult2 = new uint64_t[(cipher.l + 1) << context.logN]();
    uint64_t* axmult3 = new uint64_t[(cipher.l + 1) << context.logN]();
    uint64_t* bxmult3 = new uint64_t[(cipher.l + 1) << context.logN]();
    
    uint64_t* axres = new uint64_t[cipher.l << context.logN];
    uint64_t* bxres = new uint64_t[cipher.l << context.logN];
    
    //! s^2 -> s
    if(cipher.deg > 1){
        
        for(long i = 0; i < cipher.l; ++i){
            Key keys = decompTwoKeyMap.at(i);
            
            //! Generate RNS(axi) mod p_0, q_l
            rnsDecomp(axi, cipher.ax[1], i, cipher.l);
            
            //! multiply by the evaluation key
            mulDecompKey(axtmp, bxtmp, axi, keys, cipher.l);
            
            //! aggregate
            context.addAndEqual(axmult2, axtmp, cipher.l, 1);
            context.addAndEqual(bxmult2, bxtmp, cipher.l, 1);
        }
        
        //! ModDown by pVec[0]
        modDownByp0(axmult2, cipher.l);
        modDownByp0(bxmult2, cipher.l);
        
        //! aggregate
        context.add(axres, axmult2, cipher.ax[0], cipher.l);
        context.add(bxres, bxmult2, cipher.bx, cipher.l);
    }
    
    //! s^3 -> s
    if(cipher.deg > 2){
        for(long i = 0; i < cipher.l; ++i){
            Key keys = decompThreeKeyMap.at(i);
    
            //! Generate RNS(axi) mod p_0, q_l
            rnsDecomp(axi, cipher.ax[2], i, cipher.l);
            
            //! multiply by the evaluation key
            mulDecompKey(axtmp, bxtmp, axi, keys, cipher.l);
            
            //! aggregate
            context.addAndEqual(axmult3, axtmp, cipher.l, 1);
            context.addAndEqual(bxmult3, bxtmp, cipher.l, 1);
        }
        
        //! ModDown by pVec[0]
        modDownByp0(axmult3, cipher.l);
        modDownByp0(bxmult3, cipher.l);
        
        //! aggregate
        context.addAndEqual(axres, axmult3, cipher.l);
        context.addAndEqual(bxres, bxmult3, cipher.l);
    }

    
    delete[] axi;
    delete[] axmult2;
    delete[] bxmult2;
    delete[] axmult3;
    delete[] bxmult3;
    delete[] axtmp;
    delete[] bxtmp;
    
    return Ciphertext(axres, bxres, context.N, cipher.slots, cipher.l);
}

/**************************************************************************************/
//! Rotation by rotSlots (rotSlots: PoT)
Ciphertext ExtScheme::leftRotateFast(Ciphertext& cipher, long rotSlots){
    
    uint64_t* axres = new uint64_t[cipher.l << context.logN]();
    uint64_t* bxres = new uint64_t[cipher.l << context.logN]();
    
    uint64_t* axi = new uint64_t[(cipher.l + 1) << context.logN];
    uint64_t* axtmp = new uint64_t[(cipher.l + 1) << context.logN];
    uint64_t* bxtmp = new uint64_t[(cipher.l + 1) << context.logN];
    uint64_t* axmult = new uint64_t[(cipher.l + 1) << context.logN]();
    uint64_t* bxmult = new uint64_t[(cipher.l + 1) << context.logN]();
    
    //! (b(X), a(X)) -> (b(X^k), a(X^k))
    context.leftRot(bxres, cipher.bx, cipher.l, rotSlots);
    context.leftRot(axres, cipher.ax, cipher.l, rotSlots);
    
    //! Key Switching of a(X^k)
    
    for(long i = 0; i < cipher.l; ++i){
        Key keys = decompLeftRotKeyMap.at(make_pair(rotSlots, i));
        
        //! Generate RNS(axi) mod p_0, q_l
        rnsDecomp(axi, axres, i, cipher.l);
        
        //! multiply by the evaluation key
        mulDecompKey(axtmp, bxtmp, axi, keys, cipher.l);
        
        //! aggregate
        context.addAndEqual(axmult, axtmp, cipher.l, 1);
        context.addAndEqual(bxmult, bxtmp, cipher.l, 1);
    }
    
    // ModDown by pVec[0]
    modDownByp0(axmult, cipher.l);
    modDownByp0(bxmult, cipher.l);
    
    context.addAndEqual(bxres, bxmult, cipher.l);
    
    return Ciphertext(axmult, bxres, context.N, cipher.slots, cipher.l);
}


//! Rotation by (1 << logRotSlots)
Ciphertext ExtScheme::leftRotateByPo2(Ciphertext& cipher, long logRotSlots){
    long rotSlots = (1 << logRotSlots);
    return leftRotateFast(cipher, rotSlots);
}


//! Rotation by rotSlots
Ciphertext ExtScheme::leftRotate(Ciphertext& cipher, long rotSlots){
    Ciphertext res = cipher;
    long remrotSlots = (rotSlots % cipher.slots);
    long logrotSlots = log2((double)remrotSlots) + 1;
    for (long i = 0; i < logrotSlots; ++i) {
        if(remrotSlots & 1 << i)  {
            res = leftRotateFast(res, (1 << i));
        }
    }
    return res; 
}

//! Right rotation
Ciphertext ExtScheme::rightRotateFast(Ciphertext& cipher, long rotSlots) {
    long rotSlots1 = context.Nh - rotSlots;
    return leftRotateFast(cipher, rotSlots1);
}

Ciphertext ExtScheme::rightRotateByPo2(Ciphertext& cipher, long logRotSlots){
    long rotSlots = context.Nh -  (1 << logRotSlots);
    return leftRotateFast(cipher, rotSlots);
}

Ciphertext ExtScheme::rightRotate(Ciphertext& cipher, long rotSlots){
    Ciphertext res = cipher;
    long remrotSlots = (rotSlots % cipher.slots);
    long logrotSlots = log2((double)remrotSlots) + 1;
    for (long i = 0; i < logrotSlots; ++i) {
        if(remrotSlots & 1 << i)  {
            long rotSlots1 = context.Nh -  (1 << i);
            res = leftRotateFast(res, rotSlots1);
        }
    }
    return res;
}

/**************************************************************************************/
//! raw multiplication + decomposition KS
Ciphertext ExtScheme::mult(Ciphertext& cipher1, Ciphertext& cipher2){
    ExtCiphertext extcipher = rawmult(cipher1, cipher2);
    return DecompKeySwitch(extcipher);
}

void ExtScheme::multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2){
    ExtCiphertext extcipher = rawmult(cipher1, cipher2);
    cipher1 = DecompKeySwitch(extcipher);
}

Ciphertext ExtScheme::square(Ciphertext& cipher){
    ExtCiphertext extcipher = rawsquare(cipher);
    return DecompKeySwitch(extcipher);
}

void ExtScheme::squareAndEqual(Ciphertext& cipher){
    ExtCiphertext extcipher = rawsquare(cipher);
    cipher = DecompKeySwitch(extcipher);
}
