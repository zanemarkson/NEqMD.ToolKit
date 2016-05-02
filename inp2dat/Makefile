# GNU Makefile
CC=gcc -g
CFLAGS=-Wall -O2 -lm -lgfortran 

EXE=g_inp2dat_d
O_FILE=inp2dat.o
SOURCE=inp2dat.c



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


