FC	    = gfortran
MPIFLG  = -O3 -fcheck=all

# compile
all: clean wannier cleanx

wannier: commvar 
	$(FC) -o wannier.x *f90 -fopenmp
commvar:
	$(FC) -c variable.f90

clean:
	rm -rf *o *mod wannier.x

cleanx:
	rm -rf *o *mod  
