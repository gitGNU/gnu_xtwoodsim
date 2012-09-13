# FC = the compiler to use
FC=elmerf90
#FC=gfortran


# Compiler options
#FFLAGS=-O0 -g  # No optimization and gdb support
#FFLAGS=-O -g -Wall -Iusr/share/elmersolver/include/
#FFLAGS=-O -g -Wall  -I/build/buildd/elmerfem-5.5.0.svn.4458.dfsg/debian/tmp/usr/include -fPIC -m64 -fPIC -fPIC -I. -Ibinio -I/usr/include/mpi -m64 -fPIC -I/usr/share/elmersolver/include -shared 
#FFLAGS=-O -g -Wall -fPIC -m64 -I/usr/share/elmersolver/include -shared 
FFLAGS= -Wall 


# List libraries used by the program here
LIBS=


# Suffix-rules:  Begin by throwing away all old suffix-
# rules, and then create new ones for compiling 
# *.f90-files.
.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FFLAGS) $<

## Include the dependency-list created by makedepf90 below
#include .depend

FOBJ=BauteilClass.o ElementClass.o Eriksson2006DGL_Class.o Eriksson2006Solve.o GeneralSolverPackage.o LokaleGekoppelteMatrizenClass.o Materialeigenschaften.o Materialfunktionen.o MyStdModules.o 

Eriksson2006Solve: $(FOBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(FOBJ) $(LIBS)

BauteilClass.o : BauteilClass.f90 ElementClass.o LokaleGekoppelteMatrizenClass.o 
ElementClass.o : ElementClass.f90 MyStdModules.o Eriksson2006DGL_Class.o GeneralSolverPackage.o GeneralSolverPackage.o GeneralSolverPackage.o GeneralSolverPackage.o GeneralSolverPackage.o Materialeigenschaften.o Materialeigenschaften.o LokaleGekoppelteMatrizenClass.o 
Eriksson2006DGL_Class.o : Eriksson2006DGL_Class.f90 MyStdModules.o Materialeigenschaften.o Materialeigenschaften.o GeneralSolverPackage.o LokaleGekoppelteMatrizenClass.o 
Eriksson2006Solve.o : Eriksson2006Solve.f90 BauteilClass.o 
GeneralSolverPackage.o : GeneralSolverPackage.f90 MyStdModules.o 
LokaleGekoppelteMatrizenClass.o : LokaleGekoppelteMatrizenClass.f90 GeneralSolverPackage.o 
Materialeigenschaften.o : Materialeigenschaften.f90 
Materialfunktionen.o : Materialfunktionen.f90 
MyStdModules.o : MyStdModules.f90 
# target 'clean' for deleting object- *.mod- and other 
# unwanted files
clean:
	rm -f *.o *.mod core


# Create a dependency list using makedepf90.  All files 
# that needs to be compiled to build the program, 
# i.e all source files except include files, should 
# be given on the command line to makedepf90.  
#
# The argument to the '-o' option will be the name of the
# resulting program when running 'make', in this case 
# 'foobar'
depend .depend:
	/home/ihd/Software/makedepf90/bin/makedepf90 -o Eriksson2006Solve *.f90 > .depend

#dependMaterial .dependMaterial:
#	/home/ihd/Software/makedepf90/bin/makedepf90 -o Materialfunktionen *.f90 > .depend
