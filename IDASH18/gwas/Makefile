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
	g++ -std=c++11 -O3 -I/usr/local/include $(TPFLAGS) -pthread mainLinReg.cpp  libHEGWAS.a -o foo -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	 /usr/bin/time ./foo "idashdata/covariates.txt" "idashdata/snpMat.txt"

testlogreg:
	g++ -std=c++11 -O3 -I/usr/local/include $(TPFLAGS) -pthread mainLogReg.cpp  libHEGWAS.a -o foo -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm
	 /usr/bin/time ./foo "idashdata/covariates.txt" "idashdata/snpMat.txt"

obj: $(OBJ)


%.o: %.cpp $(HEADER)
	 $(CC) $(CFLAGS) -c $(LDFLAGS)  $<  
  
libHEGWAS.a: $(OBJ)
	$(AR) -q libHEGWAS.a  $(OBJ)



