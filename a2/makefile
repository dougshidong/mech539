OUT = a2.exe
OBJ = prec.f90 a2.f90
DBG = -g
compile:
	gfortran $(OBJ) -O3 -o $(OUT)
dbg:
	gfortran -Wall -Wextra -Wconversion -pedantic -fbounds-check -ffpe-trap=zero,overflow,underflow -g -fbacktrace $(OBJ) -o $(OUT)
clean:
	rm -fv *.out *.mod *.MOD $(OUT)
