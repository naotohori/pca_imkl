program dcd2vcvmatrix

  implicit none

  ! --------------------------------------------------------------------
  integer, parameter :: MXNUM = 22000
  integer :: ihead = 1
  integer :: n
  integer :: i, j, jatom
  integer :: iarg, iargc, infile, ioutfile
  integer :: iopen_status, input_status
  integer :: iunit, iunit2, nunit, ntitle, nblock_size
  integer :: imodel, iatom, iat, num, inum, inumber
  integer :: nset, istrt, nsavc, nstep, nver
!  integer :: lunit2mp(MXNUM)
  real(4) :: delta
  real(4) :: x(MXNUM), y(MXNUM), z(MXNUM)
  real(8), allocatable :: covariance(:,:), variance(:)
  real(8) :: tempk
  real(8) :: r_num
  real(8) :: value
  real(8) :: x_i, x_j
  character(80) :: cinfile, coutfile, title
  character(4) :: hdr
  character(26) :: chain
  character(1)  :: chainid
   
  ! --------------------------------------------------------------------
  ! open files
  ! --------------------------------------------------------------------
  infile = 11
  ioutfile = 13

  iarg = iargc()
  if (iarg /= 2) then
     write (*, *) 'Usage: % PROGRAM [INPUT_FILE] [OUTPUT_FILE]'
     stop
  end if
  call getarg(1, cinfile)
  call getarg(2, coutfile)

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
  open(ioutfile, file = coutfile, status = 'UNKNOWN', &
       action = 'WRITE', iostat = iopen_status)
  if(iopen_status > 0) then
     write (*, *) 'Error: cannot open the output file'
     stop
  end if

  ! --------------------------------------------------------------------
  imodel = 1
  inum = 1
  iat = 0
  
  if (ihead == 1) then  
     do  
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
        read (infile) iat

        ! block-size
        read (infile) nblock_size
        

        if(iat > MXNUM) then
           write (*, *) 'Error: should be increase MXNUM'
           stop
        end if

        !  End of the dcd header part
        ihead = 2
        exit
     end do
  end if
  
  allocate(covariance(iat*3, iat*3))
  allocate(variance(iat*3))
  covariance(1:iat*3, 1:iat*3) = 0.0e0_8
  variance(1:iat*3) = 0.0e0_8

  do   
     read (infile, iostat = input_status) num
     if(input_status /= 0) then
        if(input_status < 0) then
           
           close(infile)
           write (*, *) 'dcd finished '
           exit
        else
           write (*, *) 'Error: input error in pdb2crd'
           stop
        end if
     end if

     iat = num / 4     
     write (*, *) 'frame number = ', imodel
   
     read (infile) (x(n), n = 1, iat)
     read (infile) num
     read (infile) num
     read (infile) (y(n), n = 1, iat)
     read (infile) num
     read (infile) num
     read (infile) (z(n), n = 1, iat)
     read (infile) num
     

     do i = 1, iat*3
        iatom = int((i-1) / 3) + 1
        selectcase (mod(i,3))
           case(1)
              x_i = x(iatom)
           case(2) 
              x_i = y(iatom)
           case(0)
              x_i = z(iatom)
        endselect
        do j = i, iat*3
           jatom = int((j-1) / 3) + 1
           selectcase (mod(j,3))
              case(1)
                 x_j = x(jatom)
              case(2)
                 x_j = y(jatom)
              case(0)
                 x_j = z(jatom)
           endselect

           covariance(j,i) = covariance(j,i) + x_i * x_j
        enddo
        variance(i) = variance(i) +  x_i
     enddo
     
     imodel = imodel + 1

  end do

  r_num = real(imodel - 1, kind=8)
  write(ioutfile, '(i)') iat*3
  do i = 1, iat*3
     do j = i, iat*3
        value = covariance(j,i) / r_num - variance(i) * variance(j) / (r_num ** 2) 
        write(ioutfile, '(f22.10)') value
        !write(ioutfile, '(f22.10)') covariance(j,i) / r_num - variance(i) * variance(j) / (r_num ** 2)
     enddo
  enddo
  close(ioutfile)
end program dcd2vcvmatrix
