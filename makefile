# makefile for 3D CFD Repository
#
# Written by Matt Blomquist
# Last Update: 2018-07-18 (YYYY-MM-DD)
#
# This file compiles and links the cfd-3d repository source code into the
# executable file main3d.out
#
main3d.out : convergence3d.o geometry3d.o initalize3d.o main3d.o output3d_results.o pressure3d.o pressure3d_correct.o pressure3d_solve.o pseudo3d_solve.o simpler3d.o solver1d_tdma.o solver1d_tdma.o solver3d_bcicgstab.o solver3d_bcicgstab2.o solver3d_gmres.o solver3d_paradiso.o solver3d_tdma.o temperature3d.o temperature3d_boundary.o temperature3d_solve.o temperature3d_source.o velocity3d.o velocity3d_boundary.o velocity3d_correct.o velocity3d_solve.o velocity3d_source.o
	ifort -o build/main3d.out -mkl build/convergence3d.o build/geometry3d.o build/initalize3d.o build/main3d.o build/output3d_results.o build/pressure3d.o build/pressure3d_correct.o build/pressure3d_solve.o build/pseudo3d_solve.o build/simpler3d.o build/solver1d_tdma.o build/solver3d_bcicgstab.o build/solver3d_bcicgstab2.o build/solver3d_gmres.o build/solver3d_paradiso.o build/solver3d_tdma.o build/temperature3d.o build/temperature3d_boundary.o build/temperature3d_solve.o build/temperature3d_source.o build/velocity3d.o build/velocity3d_boundary.o build/velocity3d_correct.o build/velocity3d_solve.o build/velocity3d_source.o

convergence2d.o : src/convergence3d.f90
	ifort -o build/convergence3d.o -c src/convergence3d.f90

geometry3d.o : src/geometry3d.f90
	ifort -o build/geometry3d.o -c src/geometry3d.f90

initalize3d.o : src/initalize3d.f90
	ifort -o build/initalize3d.o -c src/initalize3d.f90

main3d.o : src/main3d.f90
	ifort -o build/main3d.o -c src/main3d.f90

output3d_results.o : src/output3d_results.f90
	ifort -o build/output3d_results.o -c src/output3d_results.f90

pressure3d.o : src/pressure3d.f90
	ifort -o build/pressure3d.o -c src/pressure3d.f90

pressure3d_correct.o : src/pressure3d_correct.f90
	ifort -o build/pressure3d_correct.o -c src/pressure3d_correct.f90

pressure3d_solve.o : src/pressure3d_solve.f90
	ifort -o build/pressure3d_solve.o -c src/pressure3d_solve.f90

pseudo3d_solve.o : src/pseudo3d_solve.f90
	ifort -o build/pseudo3d_solve.o -c src/pseudo3d_solve.f90

simpler3d.o : src/simpler3d.f90
	ifort -o build/simpler3d.o -c src/simpler3d.f90

solver1d_tdma.o : src/solver1d_tdma.f90
	ifort -o build/solver1d_tdma.o -c src/solver1d_tdma.f90

solver3d_bcicgstab.o : src/solver3d_bcicgstab.f90
	ifort -o build/solver3d_bcicgstab.o -c src/solver3d_bcicgstab.f90

solver3d_bcicgstab2.o : src/solver3d_bcicgstab2.f90
	ifort -o build/solver3d_bcicgstab2.o -c src/solver3d_bcicgstab2.f90

solver3d_gmres.o : src/solver3d_gmres.f90
	ifort -o build/solver3d_gmres.o -c src/solver3d_gmres.f90

solver3d_paradiso.o : src/solver3d_paradiso.f90
	ifort -o build/solver3d_paradiso.o -c src/solver3d_paradiso.f90

solver3d_tdma.o : src/solver3d_tdma.f90
	ifort -o build/solver3d_tdma.o -c src/solver3d_tdma.f90

temperature3d.o : src/temperature3d.f90
	ifort -o build/temperature3d.o -c src/temperature3d.f90

temperature3d_boundary.o : src/temperature3d_boundary.f90
	ifort -o build/temperature3d_boundary.o -c src/temperature3d_boundary.f90

temperature3d_solve.o : src/temperature3d_solve.f90
	ifort -o build/temperature3d_solve.o -c src/temperature3d_solve.f90

temperature3d_source.o : src/temperature3d_source.f90
	ifort -o build/temperature3d_source.o -c src/temperature3d_source.f90

velocity3d.o : src/velocity3d.f90
	ifort -o build/velocity3d.o -c src/velocity3d.f90

velocity3d_boundary.o : src/velocity3d_boundary.f90
	ifort -o build/velocity3d_boundary.o -c src/velocity3d_boundary.f90

velocity3d_correct.o : src/velocity3d_correct.f90
	ifort -o build/velocity3d_correct.o -c src/velocity3d_correct.f90

velocity3d_solve.o : src/velocity3d_solve.f90
	ifort -o build/velocity3d_solve.o -c src/velocity3d_solve.f90

velocity3d_source.o : src/velocity3d_source.f90
	ifort -o build/velocity3d_source.o -c src/velocity3d_source.f90

clean :
	rm build/convergence3d.o build/geometry3d.o build/initalize3d.o build/main3d.o build/output3d_results.o build/pressure3d.o build/pressure3d_correct.o build/pressure3d_solve.o build/pseudo3d_solve.o build/simpler3d.o build/solver1d_tdma.o build/solver3d_bcicgstab.o build/solver3d_bcicgstab2.o build/solver3d_gmres.o build/solver3d_paradiso.o build/solver3d_tdma.o build/temperature3d.o build/temperature3d_boundary.o build/temperature3d_solve.o build/temperature3d_source.o build/velocity3d.o build/velocity3d_boundary.o build/velocity3d_correct.o build/velocity3d_solve.o build/velocity3d_source.o
