OUT = murman_cole.exe
OBJ = generategrid.f90 param.f90 boundaryConditions.f90 solvers.f90 main.f90
DBG = -Wall -Wextra -Wconversion -pedantic -fbounds-check -ffpe-trap=zero,overflow,underflow -fbacktrace -g
PRC = -fdefault-real-8 -fdefault-integer-8
#PRC = -fdefault-integer-8
compile:
	gfortran $(OBJ) $(PRC) -O3 -o $(OUT)
dbg:
	gfortran $(PRC) $(DBG) $(OBJ) -o $(OUT)
clean:
	rm -fv *.out *.mod *.MOD $(OUT)
