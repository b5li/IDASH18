/*
 * @file       TestPvalues.cpp, cpp file
 * @brief      defining functions for computing pvalues of model parameters
 *
 * @author     Miran Kim
 * @date       July. 10, 2018
 * @copyright  GNU Pub License
 */

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#include "NTL/ZZX.h"
#include "NTL/RR.h"
#include "NTL/vec_RR.h"
#include "NTL/mat_RR.h"
#include <NTL/BasicThreadPool.h>

#include "Database.h"


#define printout 0

using namespace std;

//! @Input:   path (file name)
//! @fn:      Update tag (name of covariates), sline (data list), factorDim (number of covariates), sampleDim(number of samples)

void DataFromFile(vector<string>& tag, vector<vector<string>>& sline, string path, long& factorDim, long& sampleDim, char split_char) {
 
    //factorDim = 1;              // dimension of x
    //sampleDim = 0;              // number of samples
    
    //char split_char =  ',';     // = '\t' , splitting character between each string
    //char split_char1 = 0x22;
    
    ifstream openFile(path.data());
    if(openFile.is_open()) {
        string line, temp;
        getline(openFile, line);
        
        //! read the first line and make the tag list
        for(long i = 0; i < line.length(); ++i){
            if(line[i] == split_char){
                factorDim++;
            }
        }
        istringstream split(line);
        for (string each; getline(split, each, split_char); tag.push_back(each));
        
        //! final tag: .... (13)
        //tag[factorDim-1] = tag[factorDim-1].substr(0, tag[factorDim-1].size()-1);

        //! read the sample data
        while(getline(openFile, line)){
            vector<string> vecsline;
            
            istringstream split(line);
            vector<string> tokens;
            for (string each; getline(split, each, split_char); tokens.push_back(each));
            
            sline.push_back(tokens);
            sampleDim++;
        }
    } else {
        cout << "Error: cannot read file" << endl;
    }
}


void SimpleDataFromFile(double*& Data, string path) {
    vector<string> sline;
    long nrows = 0;
    
    ifstream openFile(path.data());
 
    if(openFile.is_open()) {
        string line;
        while(getline(openFile, line)){
            sline.push_back(line);
            nrows++;
        }
    } else {
        cout << "Error: cannot read file" << endl;
    }
    
    
    Data = new double[nrows];
    for(long i = 0; i < nrows; ++i){
        Data[i] = atof(sline[i].c_str());
    }
 
}


//! @Input:   tag, covfile/ snptag, snpfile
//! @fn:      Store the data in zData[sampleDim][factorDim]

void ListzData(double*& yData, double**& zData, long& factorDim, long& sampleDim, vector<string> tag, vector<vector<string>> covfile){
    long i, j , k;
    long idtag, ytag;
    
    //! find the idtag and ytag
    i = 0;
    while(1){
        if(tag[i].compare("human_id") == 0){
            idtag = i;
            break;
        }
        i++;
    }
    
    i = 0;
    while(1){
        if(tag[i].compare("condition(y)") == 0){
            ytag = i;
            break;
        }
        i++;
    }
    
    //! Read the covariates
    double** covData = new double*[sampleDim];
    long numcov = factorDim - 1;
    for(i = 0; i < sampleDim; ++i){
        covData[i] = new double[numcov];
    }
    
    j = 0;
    k = 0;
    
    while(1){
        if((j == idtag)||(j == ytag)){
            j++;
        }
        else{
            //! jth covariate
            vector<long> missing_index;
            double mean = 0.0;
            
            for(i = 0; i < sampleDim; ++i){
                if(covfile[i][j].compare("NA") != 0){
                    covData[i][k] = atof(covfile[i][j].c_str());
                    mean += covData[i][k];
                }
                else{
                    missing_index.push_back(i);
                }
            }
            mean /= (sampleDim - missing_index.size());
            
            for(long l = 0; l < missing_index.size(); ++l){
                long index = missing_index[l];
                covData[index][k] = mean;
            }
            j++;
            k++;
        }
        
        if(j == tag.size()){
            break;
        }
    }
  
#if 0
    for(i = 0; i < sampleDim; ++i){
        cout << i << ": " ;
        for(j = 0; j < numcov; ++j){
            cout << covData[i][j] << "\t" ;
        }
        cout << endl;
    }
#endif
    for(i = 0; i < sampleDim; ++i){
        yData[i] = atof(covfile[i][ytag].c_str());
        //zData[i][0] = 2 * atof(covfile[i][ytag].c_str()) - 1;
        
        zData[i][0] = 1.0;
        for(j = 1; j < factorDim; ++j){
            zData[i][j] = covData[i][j - 1];
            //zData[i][j] = zData[i][0] * covData[i][j - 1];
        }
        //cout << endl;
    }
 
}
 

void ListsData(double**& sData, long& nsnp, long& sampleDim, vector<vector<string>> snpfile){
    for(long i = 0; i < sampleDim; ++i){
        for(long j = 0; j < nsnp; ++j){
            if(snpfile[i][j].compare("NA") != 0){
                sData[i][j] =  atof(snpfile[i][j].c_str());
            }
            else{
                cout << "Error: cannot read the snp" << endl;
            }
        }
    }
}

//! @Input: zData
//! @ft: Data (normalized)
void normalizeData(double**& Data, double** zData, long factorDim, long sampleDim){
    Data = new double*[sampleDim];
    for(long i = 0; i < sampleDim; ++i){
        Data[i] = new double[factorDim];
    }
    for (long j = 0; j < factorDim; ++j) {
        double m0 = abs(zData[0][j]);
        double m1 = abs(zData[0][j]);
        for (long i = 1; i < sampleDim; ++i) {
            m0 = max(m0, abs(zData[i][j]));
            m1 = min(m1, abs(zData[i][j]));
        }
        if(m0 < 1e-10) continue;
    
        if(m0 != m1){
            for (long i = 0; i < sampleDim; ++i) {
                Data[i][j] = (zData[i][j] - m1)/(m0 - m1);
            }
        } else{
            for (long i = 0; i < sampleDim; ++i) {
                Data[i][j] = zData[i][j];
            }
        }
    }
}

void printvector(double* vec, long k){

    cout << "[" ;
    for(int i = 0; i< k; ++i){
        cout << vec[i] ;
        if(i < k -1) cout << ", ";
        if((i + 1) % 10 == 0) cout << endl;
    }
    cout << "]" << endl;
}

//!@ Input: vec_RR
//!@ Function: print the vector
//!  If k = 0, then print out all the components of an input vector
#if 0
void printRvector(vec_RR& vec, long k){
    long len;
    
    if(k==0){
        len= vec.length();
    }
    else{
        len = k;
    }
    
    cout << "[" ;
    for(int i = 0; i< len; ++i){
        cout << vec[i] ;
        if(i < len -1) cout << ", ";
        if((i + 1) % 10 == 0) cout << endl;
    }
    cout << "]" << endl;
}


//!@ Input: R-matrix
//!@ Function: print the matrix

void printRmatrix(Mat<RR>& mat, const long k){
    long rlen, clen;
    
    if(k==0){
        rlen= mat.NumRows();
        clen= mat.NumCols();
    }
    else{
        rlen = k;
        clen = k;
    }
    
    for(int i = 0; i< rlen; ++i){
        cout << "[";
        for(int j = 0; j < clen; ++j){
            cout << mat[i][j] ;
            if(j < clen-1) cout << "\t";
        }
        cout << "]\n";
    }
}
#endif

void printvectorToFile(double* vec, string filename, long k){
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   /// open the file
    
    for(int i = 0; i< k; ++i){
        outf << vec[i]  << endl;
    }

    outf.close();
}


void printRvectorToFile(vec_RR& vec, string filename, long k){
    long len;
    
    if(k == 0){
        len = vec.length();
    }
    else{
        len = k;
    }
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   /// open the file
    
    outf << "[" ;
    for(int i = 0; i< len; ++i){
        outf << vec[i] ;
        if(i < len -1) outf << ", ";
        if((i + 1) % 10 == 0) outf << endl;
    }
    outf << "]" << endl;
    outf.close();
}

void printRvectorToFile1(vec_RR& vec, string filename, long k){
    long len;
    
    if(k == 0){
        len = vec.length();
    }
    else{
        len = k;
    }
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   /// open the file
    
    for(int i = 0; i< len; ++i){
        outf << vec[i]  << endl;
    }
    outf.close();
}

void printRmatrix(Mat<RR>& mat, string filename, const long k){
    long rlen, clen;
    
    if(k==0){
        rlen= mat.NumRows();
        clen= mat.NumCols();
    }
    else{
        rlen = k;
        clen = k;
    }
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   /// open the file
    
    for(int i = 0; i< rlen; ++i){
        outf << "[";
        for(int j = 0; j < clen; ++j){
            outf << mat[i][j] ;
            if(j < clen-1) outf << "\t";
        }
        outf << "]\n";
    }
    outf.close();
}


