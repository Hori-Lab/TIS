! read_pdb_lipid
!> @brief Read pdb file of coarse-graind lipid molecules

! ***********************************************************************
subroutine read_pdb_lipid(lun, nunit, lunit2mp, nmp, nres, ires_mp, &
     xyz_mp, cmp2seq, cmp2atom)

  use const_maxsize
  use const_index

  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lun
  integer, intent(inout) :: nmp, nunit, nres
  integer, intent(out) :: lunit2mp(2, MXUNIT)
  integer, intent(out) :: ires_mp(MXMP)
  real(PREC), intent(out) :: xyz_mp(3, MXMP)
  character(3), intent(out) :: cmp2seq(MXMP)
  character(4), intent(out) :: cmp2atom(MXMP)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: input_status, line
  integer :: imunitemp, imp, iatom, iunit, ires
  integer :: iok
  real(PREC) :: x, y, z, tempfactor, occupancy
  character(72) :: char72
  character(1) :: multistruct
  character(4) :: nameofatom
  character(6) :: nameid
  character(3) :: nameofmp
  character(2) :: chainid
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
  imp = nmp
  iunit = nunit
  ires = nres
  line = 0 
  iok = 0

  ! ---------------------------------------------------------------------
  rewind(lun)

  do
     read (lun, '(a72)', iostat = input_status) char72
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: input error in read_xyzfile_lipid'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     line = line + 1 
     if(iok == 1 .and. char72(1:5) == 'MODEL') then
        lunit2mp(2, iunit) = imp
        iunit = iunit + 1
        ires = ires + 1
        iok = 0
     end if

     if(char72(1:4) == 'ATOM') then
        read (char72, '(a6, i5, 1x, a4, a1, a3, a2, i4, 4x, 3f8.3, 2f6.2)', &
             iostat = input_status) &
             nameid, iatom, nameofatom, multistruct, &
             nameofmp, chainid, &
             imunitemp, x, y, z, occupancy, tempfactor
        if(input_status > 0) then
           error_message = 'Error: cannot read pdb file in read_xyzfile_lipid'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        
        if(nameofmp /= 'COR' .and. nameofmp /= 'INT' .and. nameofmp /= 'TAI') then
           error_message = 'Error: non lipid in read_xyzfile_lipid'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! phosphate
        if(nameofatom == ' C  ') then
           imp = imp + 1
           ires = ires + 1
           ires_mp(imp) = ires
           cmp2seq(imp) = nameofmp  
           cmp2atom(imp) = ' C  '
           xyz_mp(1, imp) = x
           xyz_mp(2, imp) = y
           xyz_mp(3, imp) = z  
        end if
        
        ! sugar
        if(nameofatom == ' O  ') then
           imp = imp + 1
           ires_mp(imp) = ires
           cmp2seq(imp) = nameofmp
           cmp2atom(imp) = ' O  '
           xyz_mp(1, imp) = x
           xyz_mp(2, imp) = y
           xyz_mp(3, imp) = z  
        end if
        
        ! base
        if(nameofatom == ' N  ') then
           imp = imp + 1
           ires_mp(imp) = ires
           cmp2seq(imp) = nameofmp
           cmp2atom(imp) = ' N  '
           xyz_mp(1, imp) = x
           xyz_mp(2, imp) = y
           xyz_mp(3, imp) = z  
        end if
     end if
  end do

  lunit2mp(2, iunit) = imp
  if(iunit > 1) then
     if(lunit2mp(2, iunit) == lunit2mp(2, iunit - 1)) then
        iunit = iunit - 1
     end if
  end if

  nunit = iunit
  nmp = imp
  nres = ires

end subroutine read_pdb_lipid
