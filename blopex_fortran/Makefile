
CC = gcc
FC = gfortran

CFLAGS = -O2 -Wall -g
FFLAGS = -O2

LOBPCG_ROOT_DIR   =  ../blopex_abstract
LOBPCG_SERIAL_DIR =  ../blopex_serial_double

INCLUDE = \
-I$(LOBPCG_ROOT_DIR)/include \
-I$(LOBPCG_SERIAL_DIR)/multivector \
-I$(LOBPCG_SERIAL_DIR)/matmultivec \
-I$(LOBPCG_SERIAL_DIR)/pcg_multi

LIB_BLOPEX = -L$(LOBPCG_ROOT_DIR)/lib -lBLOPEX
LIB_LAPACK = -llapack
LIBS = $(LIB_BLOPEX) $(LIB_LAPACK)

OBJS = blopex_fortran_wrapper.o blopex_fortran.o serial_fortran_driver.o

all: serial_fortran_driver

serial_fortran_driver: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

clean: FORCE
	rm -f *.o *.mod

FORCE:

#gcc -Wall -g \
#-I../blopex_abstract/include \
#-I../blopex_serial_double/multivector \
#-I../blopex_serial_double/pcg_multi \
#-o main_c.o -c main_c.c
#gfortran main.f90 main_c.o ../blopex_serial_double/multivector/multi_vector.o ../blopex_serial_double/pcg_multi/pcg_multi.o ../blopex_serial_double/matmultivec/matmultivec.o -L../blopex_abstract/lib -lBLOPEX -llapack
