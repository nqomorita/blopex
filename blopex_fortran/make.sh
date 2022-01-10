#!/bin/bash

rm serial_fortran_driver
rm *.o
rm *.mod

gcc -Wall -g -O2 \
-I../blopex_abstract/include \
-I../blopex_serial_double/multivector \
-I../blopex_serial_double/pcg_multi \
-o blopex_fortran_wrapper.o -c blopex_fortran_wrapper.c

gfortran -O2 \
-o blopex_fortran.o -c blopex_fortran.f90
#-O2 -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \

gfortran -O2 \
-o serial_fortran_driver.o -c serial_fortran_driver.f90
#-O2 -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \

gfortran -O2 -o serial_fortran_driver blopex_fortran.o serial_fortran_driver.o blopex_fortran_wrapper.o \
../blopex_serial_double/multivector/multi_vector.o \
../blopex_serial_double/pcg_multi/pcg_multi.o \
../blopex_serial_double/matmultivec/matmultivec.o \
-L../blopex_abstract/lib -lBLOPEX -llapack
