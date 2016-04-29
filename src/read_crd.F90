! read_crd
!> @brief Read coordinates in CARD format

!subroutine read_crd(n_coor_velo, istep, velo_mp)
subroutine read_crd(i_record_file, velo_mp)

  use const_maxsize
  use const_index
  use const_physical
  use var_inp, only : infile
  use var_struct, only : nmp_real, nunit_real, lunit2mp, xyz_mp_rep
  use var_replica, only : n_replica_mpi, irep2grep

  implicit none

  ! --------------------------------------------------------------------
  ! arguments
  !integer, intent(in) :: n_coor_velo
  !integer(L_INT), intent(out) :: istep
!  real(PREC), intent(out), optional :: velo_mp(:,:,:)
  integer, intent(in) :: i_record_file
  real(PREC), intent(out) :: velo_mp(:,:,:)


  ! --------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------
  integer :: imp, iunit
  integer :: irep, grep, ifile(MXREPLICA)
  integer :: iline
  integer :: n_a, n_b, n_c, n_d
  real(PREC) :: r_a
  real(PREC) :: xyz(SDIM)
  character(1) :: ctype
  character(8) :: cha8a, cha8b, cha8c
  character(140) :: cline
  character(CARRAY_MSG_ERROR) :: error_message
  logical :: flg_velo

!  if (present(velo_mp)) then
  if (i_record_file == RECORD_FILE%VELO) then
     ! Reading for velocities
     flg_velo = .true.
  else
     ! Reading for coordinates
     flg_velo = .false.
  endif

  ! --------------------------------------------------------------------

  do irep = 1, n_replica_mpi

     grep = irep2grep(irep)

     if(flg_velo) then
        ifile(grep) = infile%velo(grep)
     else
        error_message = 'Error: Reading from CRD is not implemented.'
        call util_error(ERROR%STOP_ALL, error_message)
        !ifile(grep) = infile%crd(grep)
     end if
  end do

  ! ---------------------------------------------------------------------
  do irep = 1, n_replica_mpi

     iline = 0
     grep = irep2grep(irep)

     ! ----------------------------------------------------------------
     ! Read the header part in the CRD file
     ! ----------------------------------------------------------------
     ! First line
     read (ifile(grep), '(a8)') cha8a
     iline = iline + 1
     if(flg_velo) then
        if (cha8a /= '*  VELOS') then
           error_message = 'Error: Header is not "*  VELOS" in .velo file.'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
     else
        if (cha8a /= '* COORDS') then
           error_message = 'Error: Header is not "* COORDS" in .crd file.'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
     end if

     ! Skip 2 lines (*)
     read (ifile(grep), *) cline
     read (ifile(grep), *) cline
     iline = iline + 2
     
     ! Read nmp_real and check it
     read (ifile(grep), '(i10,2x,3x)') n_a
     iline = iline + 1
     if (n_a /= nmp_real) then
        error_message = 'Error: nmp_real is not consistent in .crd file.'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
   
     ! ----------------------------------------------------------------
     ! Read the data part in the CRD file
     ! ----------------------------------------------------------------
     do iunit = 1, nunit_real
        do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)

           !##### write_crd.F90 #####
           ! write (ifile(grep), &
           !      "(2i10, 2(2xa8), 3f20.10, (2xa8), (2xi3.3), a1, i4.4, f20.10)") &
           !      imp, iresnum, cha8a, cha8b,  &
           !      xyz_mp_rep(1, imp, irep), &
           !      xyz_mp_rep(2, imp, irep), &
           !      xyz_mp_rep(3, imp, irep), &
           !      cha8c, imodel, ctype, iunit, 0.0
           read (ifile(grep), &
                "(2i10, 2(2xa8), 3f20.10, (2xa8), (2xi3.3), a1, i4.4, f20.10)") &
                n_a, n_b, cha8a, cha8b, xyz(1), xyz(2), xyz(3),     &
                cha8c, n_c, ctype, n_d, r_a
           iline = iline + 1

           if (n_a /= imp) then
              write(error_message,*) 'Error: imp is invalid in .crd file. line:', iline 
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           if (flg_velo) then
              velo_mp(:, imp, irep) = xyz(:)
           else
              xyz_mp_rep(:, imp, irep) = xyz(:)
           end if

        end do
     end do
  end do

end subroutine read_crd
