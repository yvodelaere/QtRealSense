REGISTRATION=../registration_v3
LDA=../lda

DEBUG=-g
#PROFILE= -pg
#OPT=-O4
#CFLAGS=-I$(LDA) -I$(REGISTRATION) -I/usr/local/include/opencv $(DEBUG) $(OPT) $(PROFILE)
CFLAGS=-I$(LDA) -I$(REGISTRATION) `pkg-config --cflags opencv` -I/usr/include/x86_64-linux-gnu $(DEBUG) $(OPT) $(PROFILE)
#CFLAGS=-I$(LDA) -I$(REGISTRATION) -I/opt/local/include -I/opt/local/include/opencv $(DEBUG) $(OPT) $(PROFILE)
CPPFLAGS=$(CFLAGS)
#LIBS=-L/usr/local/lib -lcv -lcxcore -lhighgui -lcvaux -lm
#LIBS=-L/opt/local/lib -lopencv_core -lopencv_highgui -lopencv_contrib -lopencv_imgproc -lm
LIBS=`pkg-config --libs opencv` -lm -ljpeg -lpthread -lrt
CPP=g++
CC=gcc

REGOBJ=$(REGISTRATION)/orderedpointset.o\
	$(REGISTRATION)/unorderedpointset.o\
	$(REGISTRATION)/symmetry.o\
	$(REGISTRATION)/opt.o\
	$(REGISTRATION)/point3d.o\
	$(REGISTRATION)/rangeimage.o\
	$(REGISTRATION)/nose.o\
	$(REGISTRATION)/filter.o\
	$(REGISTRATION)/holes.o\
	$(REGISTRATION)/register.o\
	$(REGISTRATION)/roi.o\
	$(REGISTRATION)/coordlist.o\
	$(REGISTRATION)/spline.o\

LDAOBJ=$(LDA)/lda.o\
	$(LDA)/asciimatrixio.o\
	$(LDA)/beematrixio.o\
	$(LDA)/matrixio.o\

all: faceut3d faceut3dbi

.c.o: 
	$(CC) $(CFLAGS) -c $<

.cc.o: 
	$(CPP) $(CPPFLAGS) -c $<

faceut3dbi: faceut3d.o main.o builtinclassifier.o
	$(MAKE) -C $(REGISTRATION) register
	$(MAKE) -C $(LDA) lda
	g++ -o $@ main.o faceut3d.o builtinclassifier.o $(REGOBJ) $(LDAOBJ) $(LIBS)

faceut3d: faceut3d.o main.o nobuiltinclassifier.o
	$(MAKE) -C $(REGISTRATION) register
	$(MAKE) -C $(LDA) lda
	g++ -o $@ main.o faceut3d.o nobuiltinclassifier.o $(REGOBJ) $(LDAOBJ) $(LIBS)

lda:
	PWD=`pwd`
	cd ../lda
	$(MAKE)
	cd $(PWD)

install: faceut3d faceut3dbi
	cp faceut3d faceut3dbi ../bin

clean:
	rm -f *.o

