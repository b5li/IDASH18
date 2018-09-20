#ifndef MATRIX_H_
#define MATRIX_H_
#include <assert.h>

// n-by-m matrix
class Matrix {
 public:
   size_t n;
   size_t m;
   double ** data;               // [1...m,m+1...2m],...,[(n-1)m+1...nm]

   explicit Matrix(size_t _n, size_t _m) : n(_n), m(_m), data(0) {
      data = new double*[n];

      for(size_t i = 0; i < n; i++) {
         data[i] = new double[m];
      }
   }

   explicit Matrix(size_t _n, size_t _m, double ** _dat) : Matrix(_n,_m) {
      for(size_t i = 0; i < _n; i++) {
         for(size_t j = 0; j < m; j++) {
            data[i][j] = _dat[i][j];
         }
      }
   }

   explicit Matrix(size_t _n, double * _dat) : Matrix(_n,1) { // column vector
      for(size_t i = 0; i < _n; i++) {
         data[i][0] = _dat[i];
      }
   }

   Matrix(Matrix const& o) : n(o.n), m(o.m), data(0) {
      data = new double*[n];
      for(size_t i = 0; i < n; i++) {
         data[i] = new double[m];
         for(size_t j = 0; j < m; j++) {
            data[i][j] = o.data[i][j];
         }
      }
   }

   inline Matrix& operator=(Matrix const& o) {
      if(&o==this) {
         return *this;
      }
      for(size_t i = 0; i < n; i++) {
         delete [] data[i];
      }
      delete [] data;

      data = new double*[n];

      for(size_t i = 0; i < n; i++) {
         data[i] = new double[m];
         for(size_t j = 0; j < m; j++) {
            data[i][j] = o.data[i][j];
         }
      }
      return *this;
   }


   ~Matrix() {
      if(data) {
         for(size_t i = 0; i < n; i++) {
            delete [] data[i];
         }
         delete [] data;
      }
   }

   inline size_t indexOf(size_t i, size_t j) const {
      assert(0<=i && i<n);
      assert(0<=j && j<m);
      return i*m + j;
   }

   // 0 <= i < n, 0 <= j < m
   inline double at(size_t i, size_t j) const {
      return data[i][j];
   }

   inline void set(size_t i, size_t j, double x) {
      data[i][j] = x;
   }
   inline void plusAndEqual(size_t i, size_t j, double x) {
      data[i][j] += x;
   }

   inline Matrix adj() const {
      assert(n==4 && m==4);
      Matrix M(n,m);
      for(size_t i = 0; i < n; i++) {
         for(size_t j = 0; j < m; j++) {
            M.data[i][j] = coFactor(i,j);
         }
      }
      return M;
   }

   inline Matrix adjT() const {
      assert(n==4 && m==4);
      Matrix M(m,n);
      for(size_t i = 0; i < n; i++) {
         for(size_t j = 0; j < m; j++) {
            M.data[j][i] = coFactor(i,j);
         }
      }
      return M;
   }

   inline Matrix coMatrix(size_t i, size_t j) const {
      assert(0<=i && i<n);
      assert(0<=j && j<m);

      Matrix M(n-1,m-1);
      size_t curr = 0;

      size_t p = 0;             // index over M
      size_t k = 0;             // index over comatrix
      while(k<n) {
         if(k==i) { k++; continue; }
         size_t q = 0;
         size_t l = 0;
         while(l<m) {
            if(l==j) { l++; continue; }
            M.data[p][q] = data[k][l];
            l++;
            q++;
         }
         k++;
         p++;
      }
      return M;
   }

   inline double coFactor(size_t i, size_t j) const {
      Matrix Mij = coMatrix(i,j);
      return ((i+j)%2 == 0) ? Mij.det() : -Mij.det(); // (-1)^(i+j) * |Mj|
   }

   inline double det() const {
      assert(n==m);
      if(n==3) {
         return data[0][0] * data[1][1] * data[2][2]
            +   data[1][0] * data[2][1] * data[0][2]
            +   data[2][0] * data[0][1] * data[1][2]
            -   data[2][0] * data[1][1] * data[0][2]
            -   data[0][1] * data[1][0] * data[2][2]
            -   data[0][0] * data[2][1] * data[1][2];
      } else if(n==4) {
         double res = 0;
         for(size_t j = 0; j < n; j++) { // expand along the 1th column
            res += data[0][j] * coFactor(0,j);
         }
         return res;
      }
      assert(!"unknown case");
      return 0;
   }

   inline Matrix inv() const {
      assert(n==m);
      Matrix M(n,n);
      Matrix const C = adjT();
      double const dinv = 1.0/det();
      for(size_t i = 0; i < n; i++) {
         for(size_t j = 0; j < n; j++) {
            M.data[i][j] = C.data[i][j] * dinv;
         }
      }
      return M;
   }

   static inline Matrix mul(Matrix const& A, Matrix const& B) {
      assert(A.m == B.n);
      Matrix M(A.n, B.m);
      for(size_t i = 0; i < M.n; i++) {
         for(size_t j = 0; j < M.m; j++) {
            double x = 0;
            for(size_t k = 0; k < A.m; k++) { x += A.data[i][k] * B.data[k][j]; }

            M.data[i][j] = x;
         }
      }
      return M;
   }

   static inline Matrix add(Matrix const& A, Matrix const& B) {
      assert(A.n == B.n);
      assert(A.m == B.m);
      Matrix M(A.n, A.m);
      for(size_t i = 0; i < A.n; i++) {
         for(size_t j = 0; j < A.m; j++) {
            M.data[i][j] = A.data[i][j] + B.data[i][j];
         }
      }
      return M;
   }

   static inline Matrix sub(Matrix const& A, Matrix const& B) {
      assert(A.n == B.n && A.m == B.m);

      Matrix M(A.n, A.m);
      for(size_t i = 0; i < A.n; i++) {
         for(size_t j = 0; j < A.m; j++) {
            M.data[i][j] = A.data[i][j] - B.data[i][j];
         }
      }
      return M;
   }

   void setData(double * a) {
      size_t curr = 0;
      for(size_t i = 0; i < n; i++) {
         for(size_t j = 0; j < m; j++) {
            data[i][j] = a[curr++];
         }
      }
   }
};

#endif



