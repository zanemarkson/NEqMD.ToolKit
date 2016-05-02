# GNU Makefile
CC=gcc -g
CFLAGS=-Wall -O2 -lm -lgfortran 

EXE=g_neqalign_vdwradii
O_FILE=neqalign.o
SOURCE=neqalign.c



HEADERS=mutil.h gutil.h futil.h
UTILDIR=../util.lib
UTILFLAG=-lzmutil
I_PATH=../util.lib/source/

all	:	$(EXE)

$(EXE)	:	$(I_PATH)mutil.h $(I_PATH)futil.h $(I_PATH)gutil.h
	$(CC) $(CFLAGS) -o ../$(EXE) $(SOURCE) -L$(UTILDIR) $(UTILFLAG) -I$(I_PATH) -lgfortran




.PHONY	:	clean
clean	:
	-rm $(O_FILE)


