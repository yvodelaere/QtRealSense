DEBUG=-g
#PROFILE= -pg
#OPT=-O4
CFLAGS=`pkg-config --cflags opencv` $(DEBUG) $(OPT) $(PROFILE)
CPPFLAGS=$(CFLAGS)
LIBS=`pkg-config --libs opencv` -lm -ljpeg -lpthread -lrt

CPP=g++
CC=gcc

all: lda reconstruct

.c.o: 
	$(CC) $(CFLAGS) -c $<

.cc.o: 
	$(CPP) $(CPPFLAGS) -c $<

standardLDA: standardLDA.o asciimatrixio.o beematrixio.o matrixio.o
	g++ -o $@ standardLDA.o asciimatrixio.o beematrixio.o matrixio.o $(LIBS)

lda: main.o lda.o asciimatrixio.o beematrixio.o matrixio.o
	g++ $(CPPFLAGS) -o $@ main.o lda.o asciimatrixio.o beematrixio.o matrixio.o $(LIBS)
	
tstlda: tstlda.o lda.o asciimatrixio.o beematrixio.o matrixio.o
	g++ $(CPPFLAGS) -o $@ tstlda.o lda.o asciimatrixio.o beematrixio.o matrixio.o $(LIBS)
	
reconstruct: reconstruct.o lda.o asciimatrixio.o beematrixio.o matrixio.o
	g++ -o $@ reconstruct.o lda.o asciimatrixio.o beematrixio.o matrixio.o $(LIBS)
	
install: lda
	cp lda ../bin

clean:
	rm -f *.o

