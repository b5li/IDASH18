###############################################################
#
#  Choose C++ compiler, compiler flags
#
###############################################################

PREFIX = $(HOME)
CC = g++
AR = ar
CFLAGS= -O0 -g -std=c++11 -pthread $(TPFLAGS)

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

HEADER = Database.h  BasicTest.h TestPvals.h CipherPvals.h TestHEPvals.h threadpool.h

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
	g++ -std=c++11 -O0 -g -I/usr/local/include $(TPFLAGS) -pthread main.cpp sys.cpp libHEGWAS.a -o foo -L../lib -L/usr/local/lib -lHEAAN -lntl -lgmp -lm


obj: $(OBJ)


%.o: %.cpp $(HEADER)
	 $(CC) $(CFLAGS) -c $(LDFLAGS)  $<  
  
libHEGWAS.a: $(OBJ)
	$(AR) -q libHEGWAS.a  $(OBJ)



