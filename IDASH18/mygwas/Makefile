###############################################################
#
#  Choose C++ compiler, compiler flags
#
###############################################################

PREFIX = $(HOME)
CC = g++
AR = ar
CFLAGS= -g -O2 -std=c++11 -pthread

# path to the library files
LDLIBS  = -L/usr/local/lib -lntl -lgmp -lm  

# link in a library       
LDFLAGS = -I/usr/local/include      
LIBS    = -lntl -lm


###############################################################
#
#  Files
#
###############################################################

HEADER = Database.h  BasicTest.h TestLinRegPvals.h CipherLinRegPvals.h TestHELinRegPvals.h threadpool.h CipherLogRegPvals.h TestLogRegPvals.h TestHELogRegPvals.h

SRC = $(HEADER:.h=.cpp)   

OBJ= $(SRC:.cpp=.o)


###############################################################
#
#  Library
#
###############################################################

all: libHEGWAS.a

clean:
	rm *.o libHEGWAS.a  || echo nothing to clean
	
new:
	make clean
	make all

###############################################################
#
#  Test
#
###############################################################

testlinreg:
	g++ -std=c++11 -O2 -I/usr/local/include -pthread mainLinReg.cpp libHEGWAS.a -o foo -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	./foo "idashdata/covariates.txt" "idashdata/snpMat.txt"

testlogreg:
	g++ -std=c++11 -O2 -I/usr/local/include -pthread mainLogReg.cpp libHEGWAS.a -o foo -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	./foo "idashdata/covariates.txt" "idashdata/snpMat.txt" 

testHElogreg_iter1:
	g++ -std=c++11 -O2 -I/usr/local/include -pthread mainHELogReg.cpp libHEGWAS.a -o foo -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	./foo "idashdata/covariates.txt" "idashdata/snpMat.txt" 1

testHElogreg_iter2:
	g++ -std=c++11 -O2 -I/usr/local/include -pthread mainHELogReg.cpp libHEGWAS.a -o foo -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	./foo "idashdata/covariates.txt" "idashdata/snpMat.txt" 2

testHElogreg_iter3:
	g++ -std=c++11 -O2 -I/usr/local/include -pthread mainHELogReg.cpp libHEGWAS.a -o foo -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	./foo "idashdata/covariates.txt" "idashdata/snpMat.txt" 3

obj: $(OBJ)

%.o: %.cpp $(HEADER)
	 $(CC) $(CFLAGS) -c $(LDFLAGS)  $<  
  
libHEGWAS.a: $(OBJ)
	$(AR) -q libHEGWAS.a  $(OBJ)



