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

HEADER = Database.h  BasicTest.h TestPvals.h CipherPvals.h TestHEPvals.h

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

test:
	g++ -std=c++11 -O2 -I/usr/local/include -pthread main.cpp ../lib/libHEAAN.a libHEGWAS.a -o foo -L/usr/local/lib -lntl -lgmp -lm


obj: $(OBJ)


%.o: %.cpp $(HEADER)
	 $(CC) $(CFLAGS) -c $(LDFLAGS)  $<  
  
libHEGWAS.a: $(OBJ)
	$(AR) -q libHEGWAS.a  $(OBJ)



