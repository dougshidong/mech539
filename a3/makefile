OUT = murman_cole.exe
OBJ = globals.f90 grid.f90 boundaryConditions.f90 solvers.f90 output.f90 main.f90
DBG = -Wall -Wextra -Wconversion -pedantic -fbounds-check -ffpe-trap=zero,overflow,underflow -fbacktrace -pg
#PRC = -fdefault-real-8 -fdefault-integer-8
PRC = 
#PRC = -fdefault-integer-8
compile:
	gfortran $(PRC) $(OBJ) -O3 -o $(OUT)
dbg:
	gfortran $(PRC) $(DBG) $(OBJ) -o $(OUT)
clean:
	rm -fv *.out *.mod *.MOD $(OUT)
