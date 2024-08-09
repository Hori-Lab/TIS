!setp_dtrna15
!> @brief
!>       

subroutine setp_pmf()

  use const_index
  use const_maxsize
  use var_io,  only: iopen_lunnum
  use var_setp,    only : inpmf
#ifdef MPI_PAR
  use mpiconst
#endif
  implicit none
  
  ! -------------------------------------------------------------
  integer :: i, itype
  real(PREC) :: r, v, r_pre
  real(PREC) :: Rbin, Rmin
  integer :: Nbin
  integer :: input_status
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  ! Currently only for Mg-P
  itype = PMFTYPE%MG_P

  iopen_lunnum = iopen_lunnum + 1
  open(iopen_lunnum, file=inpmf%path(itype), status='old', action='read', iostat=input_status)
  if (input_status > 0) then
     error_message = 'Error: cannot open the file: ' // trim(inpmf%path(itype))
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  ! First line
  read (iopen_lunnum, *, iostat = input_status) r, v
  if(input_status /= 0) then
     error_message = 'Error: cannot read the first line of PMF file: ' // inpmf%path(itype)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  Rmin = r

  ! Second line
  read (iopen_lunnum, *, iostat = input_status) r, v
  if(input_status /= 0) then
     error_message = 'Error: cannot read the second line of PMF file: ' // inpmf%path(itype)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  Rbin = r - Rmin

  Nbin = 2
  r_pre = r
  do
     read (iopen_lunnum, *, iostat = input_status) r, v
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: cannot read a line of PMF file: ' // inpmf%path(itype)
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     Nbin = Nbin + 1

     ! Check the bin width
     if (abs(r - r_pre - Rbin) > 0.0001) then
        error_message = 'Error: bin width is not uniform in PMF file: '// inpmf%path(itype)
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     r_pre = r
  enddo

  if (Nbin > MXPMFBIN) then
     error_message = 'Error: Nbin > MXPMFBIN in setp_pmf'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  inpmf%Nbin(itype) = Nbin
  inpmf%Rmin(itype) = Rmin
  inpmf%Rbin(itype) = Rbin
  inpmf%Rbininv(itype) = 1.0_PREC / Rbin
  inpmf%Rmax(itype) = r

  rewind(iopen_lunnum)

  i = 0
  do
     read (iopen_lunnum, *, iostat = input_status) r, v
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: cannot read a line of PMF file: ' // inpmf%path(itype)
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     i = i + 1
     inpmf%pmf(i, itype) = v
  enddo

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inpmf, inpmf%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif
endsubroutine setp_pmf
