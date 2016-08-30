!read_pdbatom
!> @brief Read pdb file

subroutine read_pdbatom(lun, & ![i ] target I/O unit
     nunit_atom,  & ![i ] initial and last unit
     lunit2atom,  & ![ o] correspondence list (unit -> atom)
     pdb_atom)      ![ o] string in pdbfile

  use const_maxsize
  use const_index
  use const_physical
  use var_io,   only : pdbatom
  implicit none

  ! ---------------------------------------------------------------------
  integer,       intent(in)    :: lun
  integer,       intent(in)    :: nunit_atom(2)
  integer,       intent(out)   :: lunit2atom(2, MXUNIT)
  type(pdbatom), intent(out)   :: pdb_atom(:)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: iatom, iunit
  integer :: input_status
  logical :: flg_reading         ! flag for
  character(80) :: char80
  character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
  write(*,*) '#### start read_pdbatom'
#endif

  ! ---------------------------------------------------------------------
  flg_reading   = .false.
  iunit = nunit_atom(1)
  iatom = 0

  ! ---------------------------------------------------------------------
  rewind(lun)

  do
     read (lun, '(a80)', iostat = input_status) char80
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: input error in read_pdbatom'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     if(flg_reading .and. &
          (char80(1:3) == 'TER' .or. char80(1:3) == 'END' .or. &
           char80(1:5) == 'MODEL' .or. char80(1:2) == '>>')) then
        lunit2atom(2, iunit) = iatom
        iunit = iunit + 1
        flg_reading = .false.
     end if

     if(char80(1:4) == 'ATOM') then

        ! if using H atom, this if statement should be coment out
!        if(c_name(1:1) == 'H' .or. c_name(2:2) == 'H') then
!        if(char80(13:13) == 'H' .or. char80(14:14) == 'H') then
!           cycle
!        end if

!        if(c_altloc /= " " .and. c_altloc /= "A") then
!        if(char80(17:17) /= " " .and. char80(17:17) /= "A") then
!           cycle
!        end if

        iatom = iatom + 1
        if (iatom > MXPDBATOM) then
           write(error_message,*) 'Error: iatom > MXPDBATOM, in read_pdbatom. iatom=', iatom, ' MXPDBATOM=',MXPDBATOM
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        read (char80, '(a6, i5, 1x, a4, a1, a3, 1x, a1, i4, a1, 3x, 3f8.3, 2f6.2, 10x, 2a2)', &
             iostat = input_status) &
             pdb_atom(iatom)%c_recname, pdb_atom(iatom)%i_serial, &
             pdb_atom(iatom)%c_name, pdb_atom(iatom)%c_altloc, &
             pdb_atom(iatom)%c_resname, pdb_atom(iatom)%c_chainid, &
             pdb_atom(iatom)%i_resseq, pdb_atom(iatom)%c_icode, &
             pdb_atom(iatom)%x, pdb_atom(iatom)%y, pdb_atom(iatom)%z, &
             pdb_atom(iatom)%occupancy, pdb_atom(iatom)%tempfactor, &
             pdb_atom(iatom)%c_element, pdb_atom(iatom)%c_charge
        flg_reading = .true.
     end if
  end do
  lunit2atom(2, iunit) = iatom
  if(iunit > nunit_atom(1)) then
     if(lunit2atom(2, iunit) == lunit2atom(2, iunit - 1)) then
        iunit = iunit - 1
     end if
  end if

  if(iunit /= nunit_atom(2)) then
     write(error_message, *) 'Error: the number of unit is wrong in read_pdbatom nunit(pdb file) =', iunit-nunit_atom(1)+1, ' nunit(input) = ', nunit_atom(2)-nunit_atom(1)+1
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  lunit2atom(1, nunit_atom(1)) = 1
  do iunit = nunit_atom(1) + 1, nunit_atom(2)
     lunit2atom(1, iunit) = lunit2atom(2, iunit - 1) + 1
  end do

!  write (*, *) nunit_atom(1), nunit_atom(2), lunit2atom(1, nunit_atom(1)), lunit2atom(2, nunit_atom(2))

#ifdef _DEBUG
  write(*,*) '#### end read_pdbatom'
#endif

end subroutine read_pdbatom
