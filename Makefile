FC = ifort
OPT = -O3
DLIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95_lp64
MKLINCLUDE = ${MKLROOT}/include

MKLPATH = ${MKLROOT}/lib/intel64
MKLINCLUDELAPACK = ${MKLROOT}/include/intel64/lp64

all : dcd2vcvmatrix.exe pca_imkl.exe dihxy2vcvmatrix.exe

dcd2vcvmatrix.exe :
	$(FC) $(OPT) -o $@ dcd2vcvmatrix.F90

dihxy2vcvmatrix.exe :
	$(FC) $(OPT) -o $@ dihxy2vcvmatrix.F90

pca_imkl.exe :
	$(FC) $(OPT) -o $@ pca_imkl.f90 -L$(MKLPATH) -I$(MKLINCLUDE) -I$(MKLINCLUDELAPACK) $(DLIB)

clean :
	rm -rf dcd2vcvmatrix.exe pca_imkl.exe dihxy2vcvmatrix.exe

#ifort pca_imkl.f90 -o pca_imkl -O3 -L$MKLROOT/lib/em64t -I$MKLROOT/include -I$MKLROOT/include/em64t/lp64 -lmkl_lapack95_lp64 -Wl,--start-group $MKLROOT/lib/em64t/libmkl_intel_lp64.a $MKLROOT/lib/em64t/libmkl_intel_thread.a $MKLROOT/lib/em64t/libmkl_core.a -Wl,--end-group -liomp5 -lpthread
