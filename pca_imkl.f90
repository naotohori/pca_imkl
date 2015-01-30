program pca_imkl

   USE mkl95_LAPACK, ONLY: SYEVX

   implicit none
   integer, parameter :: PREC=8
   integer, parameter :: luninp = 10, lunout = 11
   integer :: iopen_status
   integer :: iarg, iargc

   integer :: i, j, m, n
   integer :: n_eigen
   integer :: info
   real(PREC), allocatable :: a(:,:), w(:), z(:,:)
   character(100) :: cfile_input, cfile_output, cn_eigen

   iarg = iargc()
   if (iarg /= 3) then
      write(*,*) 'Usage: % PROGRAM [input matrix file] [# of eigen values] [output file]'
      stop
   endif

   call getarg(1, cfile_input)
   call getarg(2, cn_eigen)
   read(cn_eigen,'(i)') n_eigen
   write(*,*) 'n_eigen=',n_eigen
   call getarg(3, cfile_output)

   open(luninp, file=cfile_input, status='OLD', action='READ', iostat=iopen_status)
   if (iopen_status > 0) then
      write(*,*) 'cannot open the file in input '//cfile_input
      stop
   endif
   open(lunout, file=cfile_output, status='NEW', action='WRITE', iostat=iopen_status)
   if (iopen_status > 0) then
      write(*,*) 'cannot open the file in output '//cfile_output
      stop
   endif

   read(luninp,*) n
   write(*,*) 'n=',n

   allocate( a(n,n))
   allocate( w(n) )
   allocate( z(n,n_eigen) )
   a(:,:) = 0.0e0_PREC
   w(:) = 0.0e0_PREC
   z(:,:) = 0.0e0_PREC

   do i = 1, n
      do j = i, n
         read(luninp,*) a(j,i)
         !write(*,*) i,j,a(i,j)
      enddo
   enddo
   close(luninp)

   !do i = 1, n
   !   do j = 1, i-1
   !      a(j,i) = a(i,j)
   !   enddo
   !enddo
   m = 0

   write(*,*) 'start SYEVX'
   call SYEVX( a, w, uplo='L', il=n-n_eigen+1, iu=n, m=m, z=z, info=info)

   write(lunout, *) '#eigenvalue'
   do i = m, 1, -1
      write(lunout,'(F22.10)') w(i)
   enddo

   do i = m, 1, -1
      write(lunout, *) '#eigenvector ', i
      do j = 1, n
         write(lunout, '(F22.10)') z(j,i)
      enddo
   enddo
   close(lunout)
  
   deallocate(a, w, z)

   if (info /= 0) then
      write(*,*) 'Error: info='
      write(*,*) info
      stop
   endif

endprogram pca_imkl
