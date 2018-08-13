# makefile for 3D CFD Repository
#
# Written by Matt Blomquist
# Last Update: 2018-07-18 (YYYY-MM-DD)
#
# This file compiles and links the cfd-3d repository source code into the
# executable file main3d.out
#
main3d.out : convergence3d.o geometry3d.o initialize3d.o solver3d.o main3d.o output3d_results.o pressure3d.o simpler3d.o temperature3d.o velocity3d.o
	ifort -o build/main3d.out -mkl -O3 build/convergence3d.o build/solver3d.o build/geometry3d.o build/initialize3d.o build/main3d.o build/output3d_results.o build/pressure3d.o build/simpler3d.o build/temperature3d.o build/velocity3d.o

convergence3d.o : src/convergence3d.f90
	ifort -o build/convergence3d.o -O3 -c src/convergence3d.f90

geometry3d.o : src/geometry3d.f90
	ifort -o build/geometry3d.o -O3 -c src/geometry3d.f90

initialize3d.o : src/initialize3d.f90
	ifort -o build/initialize3d.o -O3 -c src/initialize3d.f90

main3d.o : src/main3d.f90
	ifort -o build/main3d.o -O3 -c src/main3d.f90

output3d_results.o : src/output3d_results.f90
	ifort -o build/output3d_results.o -O3 -c src/output3d_results.f90

pressure3d.o : src/pressure3d.f90
	ifort -o build/pressure3d.o -O3 -c src/pressure3d.f90

simpler3d.o : src/simpler3d.f90
	ifort -o build/simpler3d.o -O3 -c src/simpler3d.f90

solver3d.o : src/solver3d.f90
	ifort -o build/solver3d.o -O3 -c src/solver3d.f90

temperature3d.o : src/temperature3d.f90
	ifort -o build/temperature3d.o -O3 -c src/temperature3d.f90

velocity3d.o : src/velocity3d.f90
	ifort -o build/velocity3d.o -O3 -c src/velocity3d.f90

clean :
	rm build/convergence3d.o build/geometry3d.o build/solver3d.o build/initialize3d.o build/main3d.o build/output3d_results.o build/pressure3d.o build/simpler3d.o build/temperature3d.o build/velocity3d.o
