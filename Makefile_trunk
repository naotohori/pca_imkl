FC = ifort
OPT = -O3
DLIB = -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95
MKLINCLUDE = ${MKLROOT}/include

# @local
MKLPATH = ${MKLROOT}/lib/ia32
MKLINCLUDELAPACK = ${MKLROOT}/include/32

# @rei
#MKLPATH = ${MKLROOT}/lib/em64t
#MKLINCLUDELAPACK = ${MKLROOT}/include/em64t/lp64

all : dcd2vcvmatrix.exe pca_imkl.exe

dcd2vcvmatrix.exe :
	$(FC) $(OPT) -o $@ dcd2vcvmatrix.F90

pca_imkl.exe :
	$(FC) $(OPT) -o $@ pca_imkl.f90 -L$(MKLPATH) -I$(MKLINCLUDE) -I$(MKLINCLUDELAPACK) $(DLIB)

clean :
	rm -rf dcd2vcvmatrix pca_imkl

#ifort pca_imkl.f90 -o pca_imkl -O3 -L$MKLROOT/lib/em64t -I$MKLROOT/include -I$MKLROOT/include/em64t/lp64 -lmkl_lapack95_lp64 -Wl,--start-group $MKLROOT/lib/em64t/libmkl_intel_lp64.a $MKLROOT/lib/em64t/libmkl_intel_thread.a $MKLROOT/lib/em64t/libmkl_core.a -Wl,--end-group -liomp5 -lpthread
