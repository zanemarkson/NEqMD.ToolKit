#!/bin/bash

rm ../libzmutil.a

gfortran  -g -O3 -fdefault-double-8 -fdefault-real-8 -c mutilexe.f gutilexe.f90 utility.f90

gcc -g -O2 -c *.c -lm


ar -cvq ../libzmutil.a *.o



rm *.o

