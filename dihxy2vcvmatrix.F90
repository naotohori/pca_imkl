program dihxy2vcv

  implicit none

  ! --------------------------------------------------------------------
  integer, parameter :: MXNUM = 22000
  integer :: ihead = 1
  integer :: n
  integer :: i, j, idih, jdih
  integer :: i1, i2, j1, j2
  integer :: iarg, iargc, infile, ioutfile, ifltfile
  integer :: iopen_status, input_status
  integer :: imodel, ndih, num
  integer :: ndih_pca, ipca
  integer :: idih2pca(MXNUM)
  logical :: flg_filter
  real(4) :: x(MXNUM), y(MXNUM)
  !real(8) :: x(MXNUM), y(MXNUM)
  real(8), allocatable :: covariance(:,:), variance(:)
  real(8) :: r_num
  real(8) :: val
  real(8) :: x_i, y_i, x_j, y_j
  real(8) :: ix, iy, jx, jy
  character(80) :: cinfile, cfltfile, coutfile
   
  infile = 11
  ioutfile = 13

  iarg = iargc()

  if (iarg == 2) then
     call getarg(1, cinfile)
     call getarg(2, coutfile)
     flg_filter = .False.
  else if (iarg == 3) then
     call getarg(1, cinfile)
     call getarg(2, cfltfile)
     call getarg(3, coutfile)
     flg_filter = .True.
  else
     write (*, *) 'Usage: % PROGRAM [INPUT_FILE] [OUTPUT_FILE]'
     write (*, *) 'Usage: % PROGRAM [INPUT_FILE] [FILTER FILE] [OUTPUT_FILE]'
     stop
  endif

  ! --------------------------------------------------------------------
  ! open output file
  ! --------------------------------------------------------------------
  open(ioutfile, file = coutfile, status = 'UNKNOWN', &
       action = 'WRITE', iostat = iopen_status)
  if(iopen_status > 0) then
     write (*, *) 'Error: cannot open the output file'
     stop
  end if

  ! --------------------------------------------------------------------
  ! open and read filter file
  ! --------------------------------------------------------------------
  if (flg_filter) then

     open(ifltfile, file = cfltfile, status = 'OLD', action = 'READ', iostat = iopen_status)
     if(iopen_status > 0) then
        write (*, *) 'Error: cannot open the filter file'
        stop
     end if

     i = 0
     ipca = 0
     idih2pca(:) = 0
     do 
        read (ifltfile, *, iostat = input_status) num
        if(input_status /= 0) then
           if(input_status < 0) then
              close(ifltfile)
              exit
           else
              write (*, *) 'Error: input error of filter file'
              stop
           end if
        end if
        
        i = i + 1
        if (i > MXNUM) then
           write(*,*) 'Error: number of lines exceeds MXNUM'
           stop
        endif

        if (num == 0) then
           idih2pca(i) = 0
        else
           ipca = ipca + 1
           idih2pca(i) = ipca
        endif

     enddo

     ndih_pca = ipca

  else
     do i = 1, MXNUM
        idih2pca(i) = i
     enddo
  endif

  ! --------------------------------------------------------------------
  ! open input dihxy file
  ! --------------------------------------------------------------------
  open(infile, file = cinfile, status = 'OLD', &
       action = 'READ', iostat = iopen_status, &
#ifdef UNFORMATTED
       form='unformatted', access='transparent') 
#else
       form='binary')
#endif
  if(iopen_status > 0) then
     write (*, *) 'Error: cannot open the input file'
     stop
  end if

  imodel = 1
  read(infile) ndih
  write (*, *) 'ndih = ', ndih


  if(ndih * 2 > MXNUM) then
     write (*, *) 'Error: should be increase MXNUM'
     stop
  end if


  !  End of the dcd header part
  if (flg_filter) then
     if (ndih < ndih_pca) then
        write (*, *) 'Error: ndih in dihxy < ndih_pca in the filter file'
        stop
     endif
  else
     ndih_pca = ndih
  endif

  
  allocate(covariance(ndih_pca*2, ndih_pca*2))  ! x (=cos(dih)) and y (=sin(dih))
  allocate(variance(ndih_pca*2))
  covariance(:, :) = 0.0e0_8
  variance(:) = 0.0e0_8

  do   
     read (infile, iostat = input_status) (x(n), y(n), n=1,ndih)
     if(input_status /= 0) then
        write (*, *) 'Error: input dihxy error'
        stop
     end if

     do idih = 1, ndih
        if (idih2pca(idih) == 0) then
           cycle
        endif

        i1 = 2 * (idih2pca(idih) - 1) + 1
        i2 = 2 * (idih2pca(idih) - 1) + 2
        ix = x(idih)
        iy = y(idih)

        do jdih = 1, ndih
           if (idih2pca(jdih) == 0) then
              cycle
           endif

           j1 = 2 * (idih2pca(jdih) - 1) + 1
           j2 = 2 * (idih2pca(jdih) - 1) + 2
           jx = x(jdih)
           jy = y(jdih)

           covariance(j1,i1) = covariance(j1,i1) + ix * jx
           covariance(j2,i1) = covariance(j2,i1) + ix * jy
           covariance(j1,i2) = covariance(j1,i2) + iy * jx
           covariance(j2,i2) = covariance(j2,i2) + iy * jy
        enddo

        variance(i1) = variance(i1) + ix
        variance(i2) = variance(i2) + iy
     enddo
     
     imodel = imodel + 1
     !write (*, *) 'frame number = ', imodel

     read (infile, iostat=input_status) num
     if(input_status /= 0) then
        if(input_status < 0) then
           
           close(infile)
           write (*, *) 'dihxy finished '
           exit
        else
           write (*, *) 'Error: input dihxy error'
           stop
        end if
     end if
     if (num /= ndih) then
        write(*,*) 'Error: num /= ndih'
        stop
     endif
  end do

  r_num = real(imodel - 1, kind=8)
  write(ioutfile, '(i)') ndih_pca*2
  do i = 1, ndih_pca*2
     do j = i, ndih_pca*2
        val = covariance(j,i) / r_num - variance(i) * variance(j) / (r_num ** 2) 
        write(ioutfile, '(f22.10)') val
     enddo
  enddo
  close(ioutfile)
end program dihxy2vcv
