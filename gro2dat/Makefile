# GNU Makefile
CC=gcc -g
CFLAGS=-Wall -O2 -lm -lgfortran 

EXE=g_gro2dat_d
O_FILE=gro2dat.o
SOURCE=gro2dat.c



HEADERS=futil.h
UTILDIR=../util.lib
UTILFLAG=-lzmutil
I_PATH=../util.lib/source/

all	:	$(EXE)

$(EXE)	:	$(I_PATH)futil.h
	$(CC) $(CFLAGS) -o ../$(EXE) $(SOURCE) -L$(UTILDIR) $(UTILFLAG) -I$(I_PATH) -lgfortran




.PHONY	:	clean
clean	:
	-rm $(O_FILE)


