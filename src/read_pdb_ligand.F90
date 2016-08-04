! read_pdb_ligand
!> @brief This subroutine is to read the pdb-structure for explicit ligand.

! ***********************************************************************
subroutine read_pdb_ligand(lun, nunit, lunit2mp, nmp, nres, ires_mp, &
     xyz_mp, cmp2seq, cmp2atom, iatomnum, xyz)

  use const_maxsize
  use const_index
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lun
  integer, intent(inout) :: nmp, nunit, nres
  integer, intent(out) :: lunit2mp(2, MXUNIT), iatomnum(MXMP)
  integer, intent(out) :: ires_mp(MXMP)
  real(PREC), intent(out) :: xyz(3, MXATOM_MP, MXMP), xyz_mp(3, MXMP)
  character(3), intent(out) :: cmp2seq(MXMP)
  character(4), intent(out) :: cmp2atom(MXMP)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: input_status, line
  integer :: imunitemp, imp, oldires, oldnmp, iatom, iunit, ires
  integer :: iok
  real(PREC) :: x, y, z, tempfactor, occupancy
  character(72) :: char72
  character(1) :: multistruct
  character(4) :: nameofatom
  character(6) :: nameid
  character(3) :: nameofmp
  character(2) :: chainid, dummy
  character(CARRAY_MSG_ERROR) :: error_message
#ifdef _DEBUG
  write(*,*) '#### start read_pdb_ligand'
#endif

  ! ---------------------------------------------------------------------
  imp = nmp
  oldnmp = imp
  iunit = nunit
  ires = nres
  line = 0 
  iok = 0
  oldires = 0

  ! ---------------------------------------------------------------------
  rewind(lun)

  do
     read (lun, '(a72)', iostat = input_status) char72
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: input error in read_pdb_ligand'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     line = line + 1 
     if(iok == 1 .and. char72(1:5) == 'MODEL') then
        lunit2mp(2, iunit) = imp
        oldnmp = imp
        iunit = iunit + 1
        iok = 0
     end if

     if(char72(1:4) == 'ATOM' .or. char72(1:6) == 'HETATM') then
        read (char72, '(a6, i5, 1x, a4, a1, a3, a2, i4, 4x, 3f8.3, 2f6.2)', &
             iostat = input_status) &
             nameid, iatom, nameofatom, multistruct, &
             nameofmp, chainid, &
             imunitemp, x, y, z, occupancy, tempfactor
        if(input_status > 0) then
           error_message = 'Error: cannot read pdb file in read_pdb_ligand'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        
        if(iok == 0) then
           dummy = chainid
           ires = ires + 1
           iok = 1
           oldires = imunitemp
        else if(iok == 1 .and. chainid /= dummy) then
           lunit2mp(2, iunit) = imp
           oldnmp = imp
           iunit = iunit + 1
           ires = ires + 1
           dummy = chainid
           oldires = imunitemp
        end if
        
        if(multistruct /= " " .and. multistruct /= "A") then
           cycle
        end if

        ! no check of atom
        imp = imp + 1
        if(imunitemp > oldires) then
           oldires = imunitemp
           ires = ires + 1
        end if
        ires_mp(imp) = ires
        iatomnum(imp) = 0
        cmp2seq(imp) = nameofmp  
        cmp2atom(imp) = nameofatom
        xyz_mp(1, imp) = x
        xyz_mp(2, imp) = y
        xyz_mp(3, imp) = z  
        
        if(imp <= oldnmp) cycle

        ! if using H atom, this if statement should be coment out
        if(nameofatom(1:1) == 'H' .or. nameofatom(2:2) == 'H') cycle

        iatomnum(imp) = iatomnum(imp) + 1
        xyz(1, iatomnum(imp), imp) = x
        xyz(2, iatomnum(imp), imp) = y
        xyz(3, iatomnum(imp), imp) = z
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
#ifdef _DEBUG
  write(*,*) '#### end read_pdb_ligand'
#endif

end subroutine read_pdb_ligand

