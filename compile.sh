ifort dcd2vcvmatrix.F90 -O3 -o dcd2vcvmatrix

ifort pca_imkl.f90 -o pca_imkl -O3 -L$MKLROOT/lib/em64t -I$MKLROOT/include -I$MKLROOT/include/em64t/lp64 -lmkl_lapack95_lp64 -Wl,--start-group $MKLROOT/lib/em64t/libmkl_intel_lp64.a $MKLROOT/lib/em64t/libmkl_intel_thread.a $MKLROOT/lib/em64t/libmkl_core.a -Wl,--end-group -liomp5 -lpthread
