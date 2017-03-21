program dcd2vcvmatrix

  implicit none

  ! --------------------------------------------------------------------
  integer, parameter :: MXNUM = 22000
  integer :: ihead = 1
  integer :: n
  integer :: i, j, jatom
  integer :: i1, i2, i3, j1, j2, j3
  integer :: iarg, iargc, infile, ioutfile, ifltfile
  integer :: iopen_status, input_status
  integer :: iunit, iunit2, nunit, ntitle, nblock_size
  integer :: imodel, iatom, natom, iat, num, inumber
  integer :: natom_pca, ipca
  integer :: nset, istrt, nsavc, nstep, nver
!  integer :: lunit2mp(MXNUM)
  integer :: idcd2pca(MXNUM)
  logical :: flg_filter
  real(4) :: delta
  real(4) :: x(MXNUM), y(MXNUM), z(MXNUM)
  real(8), allocatable :: covariance(:,:), variance(:)
  real(8) :: tempk
  real(8) :: r_num
  real(8) :: value
  real(8) :: x_i, y_i, z_i, x_j, y_j, z_j
  real(8) :: ix, iy, iz, jx, jy, jz
  character(80) :: cinfile, cfltfile, coutfile, title
  character(4) :: hdr
  character(26) :: chain
  character(1)  :: chainid
   
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
     idcd2pca(:) = 0
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
           idcd2pca(i) = 0
        else
           ipca = ipca + 1
           idcd2pca(i) = ipca
        endif

     enddo

     natom_pca = ipca

  else
     do i = 1, MXNUM
        idcd2pca(i) = i
     enddo
  endif

  ! --------------------------------------------------------------------
  ! open input DCD file
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
  
  ! --------------------------------------------------------------------
  ! read header
  ! --------------------------------------------------------------------
  ! ... block-size
  read (infile) nblock_size
  write (*, *) 'block size = ', nblock_size

  ! ... 'CORD' for coordinate, 'VELD' for velocity
  read (infile) hdr
  write (*, *) 'coodinate (CORD) or velocity (VELD) = ', hdr

  ! ... the number of frames
  nset = 1
  read (infile) nset
  write (*, *) 'the number of the frames =', nset

  ! ... starting step number
  read (infile) istrt

  ! ... step interval
  read (infile) nsavc

  ! ... the number of steps 
  read (infile) nstep

  ! ... the number of unit
  read (infile) nunit

  read (infile) num
  read (infile) num
  read (infile) num

  ! ... the number of free atoms, where it is set to be 0.
  read (infile) num

  ! ... time-step
  read (infile) delta

  ! ... unit-cell information
  read (infile) num

  ! ... read int x eight times
  read (infile) num
  read (infile) num
  read (infile) num
  read (infile) num
  read (infile) num
  read (infile) num
  read (infile) num
  read (infile) num

  ! version if CHARMm
  read (infile) nver

  ! block-size
  read (infile) nblock_size

  ! block-size
  read (infile) nblock_size

  ! the line number of title lines
  read (infile) ntitle

  ! title text
  read (infile) title
  read (infile) title

  ! read temperature
  read (infile) title
  read (title, *) tempk

  ! read lunit2mp
  do iunit = 1, nunit
     read (infile) title
     !read (title, '(i6)') lunit2mp(iunit)
  end do

  ! block-size
  read (infile) nblock_size


  ! block-size
  read (infile) nblock_size

  ! the number of atoms
  read (infile) natom

  ! block-size
  read (infile) nblock_size
  

  if(natom > MXNUM) then
     write (*, *) 'Error: should be increase MXNUM'
     stop
  end if

  !  End of the dcd header part
  if (flg_filter) then
     if (natom < natom_pca) then
        write (*, *) 'Error: natom in DCD < natom_pca in the filter file'
        stop
     endif
  else
     natom_pca = natom
  endif

  
  allocate(covariance(natom_pca*3, natom_pca*3))
  allocate(variance(natom_pca*3))
  covariance(:, :) = 0.0e0_8
  variance(:) = 0.0e0_8

  do   
     read (infile, iostat = input_status) num
     if(input_status /= 0) then
        if(input_status < 0) then
           
           close(infile)
           write (*, *) 'dcd finished '
           exit
        else
           write (*, *) 'Error: input DCD error'
           stop
        end if
     end if

     iat = num / 4     
     if (iat /= natom) then
        write(*,*) 'Error: iat /= natom'
        stop
     endif
     write (*, *) 'frame number = ', imodel
   
     read (infile) (x(n), n = 1, iat)
     read (infile) num
     read (infile) num
     read (infile) (y(n), n = 1, iat)
     read (infile) num
     read (infile) num
     read (infile) (z(n), n = 1, iat)
     read (infile) num
     

     do iatom = 1, natom
        if (idcd2pca(iatom) == 0) then
           cycle
        endif

        i1 = 3 * (idcd2pca(iatom) - 1) + 1
        i2 = 3 * (idcd2pca(iatom) - 1) + 2
        i3 = 3 * (idcd2pca(iatom) - 1) + 3
        ix = x(iatom)
        iy = y(iatom)
        iz = z(iatom)

        do jatom = 1, natom
           if (idcd2pca(jatom) == 0) then
              cycle
           endif

           j1 = 3 * (idcd2pca(jatom) - 1) + 1
           j2 = 3 * (idcd2pca(jatom) - 1) + 2
           j3 = 3 * (idcd2pca(jatom) - 1) + 3
           jx = x(jatom)
           jy = y(jatom)
           jz = z(jatom)

           covariance(j1,i1) = covariance(j1,i1) + ix * jx
           covariance(j2,i1) = covariance(j2,i1) + ix * jy
           covariance(j3,i1) = covariance(j3,i1) + ix * jz
           covariance(j1,i2) = covariance(j1,i2) + iy * jx
           covariance(j2,i2) = covariance(j2,i2) + iy * jy
           covariance(j3,i2) = covariance(j3,i2) + iy * jz
           covariance(j1,i3) = covariance(j1,i3) + iz * jx
           covariance(j2,i3) = covariance(j2,i3) + iz * jy
           covariance(j3,i3) = covariance(j3,i3) + iz * jz
        enddo

        variance(i1) = variance(i1) + ix
        variance(i2) = variance(i2) + iy
        variance(i3) = variance(i3) + iz
     enddo
     
     imodel = imodel + 1

  end do

  r_num = real(imodel - 1, kind=8)
  write(ioutfile, '(i)') natom_pca*3
  do i = 1, natom_pca*3
     do j = i, natom_pca*3
        value = covariance(j,i) / r_num - variance(i) * variance(j) / (r_num ** 2) 
        write(ioutfile, '(f22.10)') value
        !write(ioutfile, '(f22.10)') covariance(j,i) / r_num - variance(i) * variance(j) / (r_num ** 2)
     enddo
  enddo
  close(ioutfile)
end program dcd2vcvmatrix
