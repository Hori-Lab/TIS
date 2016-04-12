! read_pdb_dna
!> @brief Read the PDB file for DNA

subroutine read_pdb_dna(lun, nunit, lunit2mp, nmp, nres, ires_mp, &
     xyz_mp, iontype_mp, cmp2seq, cmp2atom, imp2type, iatomnum, xyz)

  use const_maxsize
  use const_index
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lun
  integer, intent(inout) :: nmp, nunit, nres
  integer, intent(inout) :: lunit2mp(2, MXUNIT), iatomnum(MXMP)
  integer, intent(inout) :: ires_mp(MXMP), iontype_mp(MXMP)
  real(PREC), intent(inout) :: xyz(3, MXATOM_MP, MXMP), xyz_mp(3, MXMP)
  character(3), intent(inout) :: cmp2seq(MXMP)
  character(4), intent(inout) :: cmp2atom(MXMP)
  integer, intent(inout)      :: imp2type(MXMP)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: input_status, line
  integer :: imunitemp, imp, oldnmp, iatom, iunit, ires
  integer :: iok
  integer :: ipupi = 1
  real(PREC) :: x, y, z, tempfactor, occupancy
  character(72) :: char72
  character(1) :: multistruct
  character(4) :: nameofatom
  character(6) :: nameid
  character(3) :: nameofmp
  character(2) :: chainid, dummy
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
  imp = nmp
  iatom = 0
  oldnmp = imp
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
        error_message = 'Error: input error in read_xyzfile_dna'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     
     line = line + 1 
     if(iok == 1 .and. char72(1:5) == 'MODEL') then
        lunit2mp(2, iunit) = imp
        oldnmp = imp
        iunit = iunit + 1
        ires = ires + 1
        iok = 0
     end if
        
     !else if (iok == 1 .and. char72(1:3) == 'TER') then

     !   lunit2mp(2, iunit) = imp
     !   oldnmp = imp
     !   iunit = iunit + 1
     !   write(*, *) 'read_pdb_dna: iunit: ', iunit
     !   ires = ires + 1
        
     if(char72(1:4) == 'ATOM') then
        read (char72, '(a6, i5, 1x, a4, a1, a3, a2, i4, 4x, 3f8.3, 2f6.2)', &
             iostat = input_status) &
             nameid, iatom, nameofatom, multistruct, &
             nameofmp, chainid, &
             imunitemp, x, y, z, occupancy, tempfactor
        if(input_status > 0) then
           error_message = 'Error: cannot read pdb file in read_xyzfile_dna'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        
        if(nameofmp == ' DA' .or. nameofmp == ' DG') then
           ipupi = 1
        else if(nameofmp == ' DT' .or. nameofmp == ' DC') then
           ipupi = 2
        else
           error_message = 'Error: non dna chain in read_xyzfile_dna'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        if(iok == 0) then
           dummy = chainid
           iok = 1
        else if(iok == 1 .and. chainid /= dummy) then
           lunit2mp(2, iunit) = imp
           oldnmp = imp
           iunit = iunit + 1
           ires = ires + 1
           dummy = chainid
        end if

        !if(iok == 0) then
        !   iok = 1
        !end if
        
        if(multistruct /= " " .and. multistruct /= "A") then
           cycle
        end if

        ! phosphate
        !if (nameofatom == ' P  ' .or. \
        !    nameofatom == ' OP1' .or. \
        !    nameofatom == ' OP2') then
            
        !    imp = imp + 1
        !    iatom = iatom + 1
            
        !    imp_atom(iatom) = imp
        !    write(*, *) 'read_pdb_dna: char72: ', char72
            
         !end if

        
        if(nameofatom == ' P  ') then
           imp = imp + 1
           ires = ires + 1
           ires_mp(imp) = ires
           iatomnum(imp) = 0
           iontype_mp(imp) = IONTYPE%P
           cmp2seq(imp) = nameofmp  
           cmp2atom(imp) = ' O  '
           imp2type(imp) = MPTYPE%DNA_PHOS
           xyz_mp(1, imp) = x
           xyz_mp(2, imp) = y
           xyz_mp(3, imp) = z  
        end if
        
        ! sugar
        if(nameofatom == " O5'") then
           imp = imp + 1
           ires_mp(imp) = ires
           iatomnum(imp) = 0
           cmp2seq(imp) = nameofmp
           cmp2atom(imp) = ' S  '
           imp2type(imp) = MPTYPE%DNA_SUGAR
           xyz_mp(1, imp) = x
           xyz_mp(2, imp) = y
           xyz_mp(3, imp) = z  
        end if 
        
        ! base
        if((ipupi == 1 .and. nameofatom == ' N9 ') .or. &
           (ipupi == 2 .and. nameofatom == ' N1 ')) then
           imp = imp + 1
           ires_mp(imp) = ires
           iatomnum(imp) = 0
           cmp2seq(imp) = nameofmp
           cmp2atom(imp) = ' N  '
           imp2type(imp) = MPTYPE%DNA_BASE
        end if

        if((ipupi == 1 .and. nameofatom == ' N1 ') .or. &
           (ipupi == 2 .and. nameofatom == ' N3 ')) then
           xyz_mp(1, imp) = x
           xyz_mp(2, imp) = y
           xyz_mp(3, imp) = z  
        end if

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

end subroutine read_pdb_dna

