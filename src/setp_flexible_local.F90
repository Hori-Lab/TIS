! setp_flexible_local
!> @brief This subroutine is to read and set the parameters for flexible local 
!>        potential

! *****************************************************************************
subroutine setp_flexible_local()

  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile
  use var_setp, only : inflp
  use var_struct, only: nres, nfba, nfdih, &
       fba_para_x, fba_para_y, fba_para_y2, fdih_para

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  !----------------------------------------------------------------------------
  ! Local variables
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat, icol
  integer :: ier
  integer :: i_ang1, i_ang2, i_dih
  integer :: iba, idih, ip

  character(4) :: kfind
  character(5) :: header
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  real(PREC), allocatable, save :: fba_mult_para_y1(:,:)  !(10, nres)
  real(PREC), allocatable, save :: fba_mult_para_y2(:,:)  !(10, nres)
  integer,    allocatable, save :: fba_mult_resID1(:)     !(nres)
  integer,    allocatable, save :: fba_mult_resID2(:)     !(nres)
  real(PREC), allocatable, save :: fdih_mult_para(:,:) !(7, nres)
  integer,    allocatable, save :: fdih_mult_resID(:)  !(nres)

  !----------------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  
#ifdef MPI_PAR
  if (myrank == 0) then
#endif
     
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'flexible_local  ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)

     if (kfind /= 'FIND') then
        error_message = 'Error: cannot find "flexible_local" field in setp_flexible_local'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     error_message = 'failed in memory allocation at allocate_nativestruct, PROGRAM STOP'

     allocate( fba_mult_resID1(nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fba_mult_resID1(:) = 0
     allocate( fba_mult_resID2(nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fba_mult_resID2(:) = 0
     allocate( fdih_mult_resID(nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fdih_mult_resID(:) = 0
     allocate( fba_mult_para_y1(10, nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fba_mult_para_y1(:,:) = 0
     allocate( fba_mult_para_y2(10, nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fba_mult_para_y2(:,:) = 0
     allocate( fdih_mult_para(7, nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fdih_mult_para(:,:) = 0

     i_ang1 = 1
     i_ang2 = 1
     i_dih  = 1
     
     do iline = 1, nlines
        ctmp00 = cwkinp(iline)
        if (ctmp00(1:5) == 'ANGY1') then
           read (ctmp00, *) header, fba_mult_resID1(i_ang1), &
                fba_mult_para_y1(1,  i_ang1), &
                fba_mult_para_y1(2,  i_ang1), &
                fba_mult_para_y1(3,  i_ang1), &
                fba_mult_para_y1(4,  i_ang1), &
                fba_mult_para_y1(5,  i_ang1), &
                fba_mult_para_y1(6,  i_ang1), &
                fba_mult_para_y1(7,  i_ang1), &
                fba_mult_para_y1(8,  i_ang1), &
                fba_mult_para_y1(9,  i_ang1), &
                fba_mult_para_y1(10, i_ang1)
           i_ang1 = i_ang1 + 1
        else if (ctmp00(1:5) == 'ANGY2') then
           read (ctmp00, *) header, fba_mult_resID2(i_ang2), &
                fba_mult_para_y2(1,  i_ang2), &
                fba_mult_para_y2(2,  i_ang2), &
                fba_mult_para_y2(3,  i_ang2), &
                fba_mult_para_y2(4,  i_ang2), &
                fba_mult_para_y2(5,  i_ang2), &
                fba_mult_para_y2(6,  i_ang2), &
                fba_mult_para_y2(7,  i_ang2), &
                fba_mult_para_y2(8,  i_ang2), &
                fba_mult_para_y2(9,  i_ang2), &
                fba_mult_para_y2(10, i_ang2)
           i_ang2 = i_ang2 + 1
        else if (ctmp00(1:4) == 'DIH') then
           read (ctmp00, *) header, fdih_mult_resID(i_dih), &
                fdih_mult_para(1,  i_dih), &
                fdih_mult_para(2,  i_dih), &
                fdih_mult_para(3,  i_dih), &
                fdih_mult_para(4,  i_dih), &
                fdih_mult_para(5,  i_dih), &
                fdih_mult_para(6,  i_dih), &
                fdih_mult_para(7,  i_dih)
           i_dih = i_dih + 1
        else 
           call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
           
           do iequa = 1, nequat
              cvalue = 'k_ang'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inflp%k_ang, cvalue)
              
              cvalue = 'k_dih'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inflp%k_dih, cvalue)
           end do
        end if
     end do

     do iba = 1, nres
        if (fba_mult_resID1(iba) /= 0) then
           do ip = 1, 10
              fba_para_y(ip, fba_mult_resID1(iba)) = fba_mult_para_y1(ip, iba)
           end do
        end if
     end do
     
     do iba = 1, nres
        if (fba_mult_resID2(iba) /= 0) then
           do ip = 1, 10
              fba_para_y2(ip, fba_mult_resID2(iba)) = fba_mult_para_y2(ip, iba)
           end do
        end if
     end do

     do idih = 1, nres
        if (fdih_mult_resID(idih) /= 0) then
           do ip = 1, 7
              fdih_para(ip, fdih_mult_resID(idih)) = fdih_mult_para(ip, idih)
           end do
        end if
     end do
     

#ifdef MPI_PAR
  end if

!  call MPI_Bcast (inflp%k_ang, 1, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
!  call MPI_Bcast (inflp%k_dih, 1, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast (inflp, inflp%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_para_x, 10*nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_para_y, 10*nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_para_y2, 10*nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fdih_para, 7*nfdih, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
   
#endif

end subroutine setp_flexible_local
