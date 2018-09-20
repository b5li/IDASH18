###############################################################
#
#  Choose C++ compiler, compiler flags
#
###############################################################

PREFIX = $(HOME)
CC = g++
AR = ar
CFLAGS= -O3 -std=c++11 -pthread $(TPFLAGS)

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

HEADER = Database.h  BasicTest.h TestLinRegPvals.h CipherLinRegPvals.h TestHELinRegPvals.h threadpool.h CipherLogRegPvals.h TestLogRegPvals.h TestHELogRegPvals.h sys.h

SRC = $(HEADER:.h=.cpp)   

OBJ= $(SRC:.cpp=.o)


###############################################################
#
#  Library
#
###############################################################

all: libHEGWAS.a

clean:
	rm -f *.o libHEGWAS.a testlinreg testlogreg testtrivial || echo nothing to clean
	
new:
	make clean
	make all

###############################################################
#
#  Test
#
###############################################################

testlinreg: mainLinReg.cpp libHEGWAS.a
	g++ -std=c++11 -O3 -I/usr/local/include $(TPFLAGS) -pthread mainLinReg.cpp  libHEGWAS.a -o testlinreg -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	/usr/bin/time ./testlinreg "idashdata/covariates.txt" "idashdata/snpMat.txt"

testlogreg: mainLogReg.cpp libHEGWAS.a
	g++ -std=c++11 -O3 -I/usr/local/include $(TPFLAGS) -pthread mainLogReg.cpp  libHEGWAS.a -o testlogreg -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	/usr/bin/time ./testlogreg "idashdata/covariates.txt" "idashdata/snpMat.txt"

testtrivial: mainTrivialLinReg.cpp libHEGWAS.a
	g++ -std=c++11 -O3 -I/usr/local/include $(TPFLAGS) -pthread mainTrivialLinReg.cpp  libHEGWAS.a -o testtrivial -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	/usr/bin/time ./testtrivial "idashdata/covariates.txt" "idashdata/snpMat.txt"

obj: $(OBJ)


%.o: %.cpp $(HEADER)
	 $(CC) $(CFLAGS) -c $(LDFLAGS)  $<  
  
libHEGWAS.a: $(OBJ)
	$(AR) -q libHEGWAS.a  $(OBJ)



