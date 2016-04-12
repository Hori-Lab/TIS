! read_pdbatom_dna
!> @brief Read the pdb_atom for DNA

subroutine read_pdbatom_dna(pdb_atom, lunit2atom, nunit_atom, &
     lunit2mp, nmp, nres, ires_mp, xyz_mp, iontype_mp, &
     cmp2seq, cmp2atom, imp2type, iatomnum, xyz)

  use const_maxsize
  use const_index
  use const_physical
  use var_setp,   only : pdbatom
  implicit none

  ! ---------------------------------------------------------------------
  type(pdbatom), intent(in)    :: pdb_atom(:)
  integer,       intent(in)    :: lunit2atom(2, MXUNIT)
  integer,       intent(in)    :: nunit_atom(2)
  integer,       intent(inout) :: lunit2mp(2, MXUNIT)
  integer,       intent(inout) :: nmp, nres
  integer,       intent(inout) :: ires_mp(MXMP)
  real(PREC),    intent(inout) :: xyz_mp(3, MXMP)
  integer,       intent(inout) :: iontype_mp(MXMP)
  character(3),  intent(inout) :: cmp2seq(MXMP)
  character(4),  intent(inout) :: cmp2atom(MXMP)
  integer,       intent(inout) :: imp2type(MXMP)
  integer,       intent(inout) :: iatomnum(MXMP)
  real(PREC),    intent(inout) :: xyz(3, MXATOM_MP, MXMP)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: iatom, imp, iunit, ibase
  integer :: ires, ires_ini, ires_last
  integer :: ipos, ipos_ini, ipos_last
  integer :: idx, ndx(3)
  integer :: ibtype, ibtype_res(MXMP)
  integer :: lpos2atom(4, MXMP), ipos2atom(MXMP)
  logical :: flg_reading         ! flag for 
  character(CARRAY_MSG_ERROR) :: error_message

  integer      :: i_resseq_save ! to store previous i_resseq
  character(1) :: c_icode_save  ! to store previous c_icode_save

  ! ---------------------------------------------------------------------
  imp = nmp
  ires = nres

  ! ---------------------------------------------------------------------
  do iunit = nunit_atom(1), nunit_atom(2)
     ires_ini = ires + 1
     ibtype_res(:) = 0
     lpos2atom(:,:) = 0
     ipos2atom(:) = 0
     flg_reading = .false.

     do iatom = lunit2atom(1, iunit), lunit2atom(2, iunit)
        
        if(pdb_atom(iatom)%c_altloc /= " " .and. &
             pdb_atom(iatom)%c_altloc /= "A") cycle

        ! if using H atom, this statement should be coment out
        if(pdb_atom(iatom)%c_name(1:1) == "H" .or. &
             pdb_atom(iatom)%c_name(2:2) == "H") cycle

        ! check new residue or not
        if(flg_reading .and. &
           pdb_atom(iatom)%i_resseq == i_resseq_save .and. &
           pdb_atom(iatom)%c_icode == c_icode_save) then

        else
           ires = ires + 1
           call sub_check_basetype(pdb_atom(iatom)%c_resname, ibtype, ndx)
           ibtype_res(ires) = ibtype
           if(ires == ires_ini) then
              lpos2atom(1, ires) = 1
           else
              lpos2atom(1, ires) = lpos2atom(4, ires - 1) + 1
           end if
           lpos2atom(2, ires) = lpos2atom(1, ires) + ndx(1)
           lpos2atom(3, ires) = lpos2atom(1, ires) + ndx(2)
           lpos2atom(4, ires) = lpos2atom(1, ires) + ndx(3) - 1

           i_resseq_save = pdb_atom(iatom)%i_resseq
           c_icode_save = pdb_atom(iatom)%c_icode
           flg_reading = .true.
        end if
        
        call sub_select_idx(ibtype, pdb_atom(iatom)%c_name, idx)
        ipos = idx + lpos2atom(1, ires) - 1
        if(ipos2atom(ipos) == 0) then
           ipos2atom(ipos) = iatom
        else
           write (error_message, *) 'Error: invalid idx in reading read_pdbatom_dna', iatom, ires, idx, ipos, ipos2atom(ipos)
           call util_error(ERROR%STOP_ALL, error_message)
        end if
       
     end do
     ires_last = ires

     ! check ini pos
     ipos_ini = 0
     do ipos = lpos2atom(1, ires_ini), lpos2atom(4, ires_ini)
        if(ipos2atom(ipos) /= 0) then
           ipos_ini = ipos
           exit
        end if
     end do
     ! start from phosphate or sugar (not base)
     if(ipos_ini /= lpos2atom(1, ires_ini) .and. &
          ipos_ini /= lpos2atom(2, ires_ini)) then
        write (error_message, *) 'Error: invalid ipos_ini in reading read_pdbatom_dna', &
             ires_ini, ipos_ini, lpos2atom(1, ires_ini), lpos2atom(2, ires_ini)
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     ! check last pos
     ipos_last = 0
     do ipos = lpos2atom(1, ires_last), lpos2atom(4, ires_last)
        if(ipos2atom(ipos) /= 0) then
           ipos_last = ipos
        end if
     end do
     ! end to phosphate, sugar, or base
     if(ipos_last /= lpos2atom(2, ires_last) - 1 .and. &
          ipos_last /= lpos2atom(3, ires_last) - 1 .and. &
          ipos_last /= lpos2atom(4, ires_last)) then
        write (error_message, *) 'Error: invalid ipos_last in reading read_pdbatom_dna', &
             ires_last, ipos_last, lpos2atom(2, ires_last) - 1, &
             lpos2atom(3, ires_last) - 1, lpos2atom(4, ires_last)
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     ! check ipos2atom
     do ipos = ipos_ini, ipos_last
        if(ipos2atom(ipos) == 0) then
           write (error_message, *) 'Error: invalid ipos2atom in reading read_pdbatom_dna', ipos, ipos2atom(ipos)
           call util_error(ERROR%STOP_ALL, error_message)
        end if
     end do

     ! store information such as coorinates, residue name, and atom name
     lunit2mp(1, iunit) = imp + 1
     do ires = ires_ini, ires_last
        do ipos = lpos2atom(1, ires), lpos2atom(4, ires)
           if(ipos < ipos_ini .or. ipos > ipos_last) cycle

           iatom = ipos2atom(ipos)
           if(ipos == lpos2atom(1, ires) .or. &
                ipos == lpos2atom(2, ires) .or. &
                ipos == lpos2atom(3, ires)) then
              imp = imp + 1
              if (imp > MXMP) then
                 write(error_message,*) 'Error: imp > MXMP, in read_pdbatom_pro imp=',imp,' MXMP=',MXMP
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
              
              ires_mp(imp) = ires
              if(ibtype_res(ires) == CHEMICALTYPE%DA) then
                 cmp2seq(imp) = ' DA'
              else if(ibtype_res(ires) == CHEMICALTYPE%DG) then
                 cmp2seq(imp) = ' DG'
              else if(ibtype_res(ires) == CHEMICALTYPE%DT) then
                 cmp2seq(imp) = ' DT'
              else if(ibtype_res(ires) == CHEMICALTYPE%DC) then
                 cmp2seq(imp) = ' DC'
              end if
              iatomnum(imp) = 0
              xyz(1:3, :, imp) = INVALID_VALUE
           end if

           if(ipos == lpos2atom(1, ires)) then
              cmp2atom(imp) = " O  "
              imp2type(imp) = MPTYPE%DNA_PHOS
              iontype_mp(imp) = IONTYPE%P
              xyz_mp(1, imp) = pdb_atom(iatom)%x
              xyz_mp(2, imp) = pdb_atom(iatom)%y
              xyz_mp(3, imp) = pdb_atom(iatom)%z  
           end if
        
           ! sugar
           if(ipos == lpos2atom(2, ires)) then
              cmp2atom(imp) = " S  "
              imp2type(imp) = MPTYPE%DNA_SUGAR
              xyz_mp(1, imp) = pdb_atom(iatom)%x
              xyz_mp(2, imp) = pdb_atom(iatom)%y
              xyz_mp(3, imp) = pdb_atom(iatom)%z  
           end if
        
           ! base
           ibase = 0
           if(ibtype_res(ires) == CHEMICALTYPE%DA .or. &
                ibtype_res(ires) == CHEMICALTYPE%DG) then
              if(pdb_atom(iatom)%c_name == " N1 " .or. pdb_atom(iatom)%c_name == "N1  ") then
                 ibase = 1
              end if
           else if(ibtype_res(ires) == CHEMICALTYPE%DT .or. &
                ibtype_res(ires) == CHEMICALTYPE%DC) then
              if(pdb_atom(iatom)%c_name == " N3 " .or. pdb_atom(iatom)%c_name == "N3  ") then
                 ibase = 1
              end if
           else
              write(error_message,*) 'Error: invalid base type in read_pdbatom_pro', ires, ibtype_res(ires)
              call util_error(ERROR%STOP_ALL, error_message)
           end if
           if(ibase == 1) then
              cmp2atom(imp) = " N  "
              imp2type(imp) = MPTYPE%DNA_BASE
              xyz_mp(1, imp) = pdb_atom(iatom)%x
              xyz_mp(2, imp) = pdb_atom(iatom)%y
              xyz_mp(3, imp) = pdb_atom(iatom)%z
           end if
           
           iatomnum(imp) = iatomnum(imp) + 1
           xyz(1, iatomnum(imp), imp) = pdb_atom(iatom)%x
           xyz(2, iatomnum(imp), imp) = pdb_atom(iatom)%y
           xyz(3, iatomnum(imp), imp) = pdb_atom(iatom)%z
        end do
     end do
     lunit2mp(2, iunit) = imp
  end do

  nres = ires
  nmp = imp

!###########################################################################
contains

  ! check residue type
  subroutine sub_check_basetype(c_resname, ibasetype, ndx)

    implicit none

    character(3), intent(in) :: c_resname
    integer, intent(out)     :: ibasetype
    integer, intent(out)     :: ndx(3)

    ndx(1) = 3
    ndx(2) = 11

     ! acceptable residues
    if ( c_resname(1:1) == 'A' .OR. c_resname(2:2) == 'A' .OR. c_resname(3:3) == 'A' ) then
       ibasetype = CHEMICALTYPE%DA
       ndx(3) = 21
    else if ( c_resname(1:1) == 'G' .OR. c_resname(2:2) == 'G' .OR. c_resname(3:3) == 'G' ) then
       ibasetype = CHEMICALTYPE%DG
       ndx(3) = 22
    else if ( c_resname(1:1) == 'T' .OR. c_resname(2:2) == 'T' .OR. c_resname(3:3) == 'T' ) then
       ibasetype = CHEMICALTYPE%DT
       ndx(3) = 20
    else if ( c_resname(1:1) == 'C' .OR. c_resname(2:2) == 'C' .OR. c_resname(3:3) == 'C' ) then
       ibasetype = CHEMICALTYPE%DC
       ndx(3) = 19
    else
       write(error_message,*) 'Error: non DNA residue in read_pdbatom_dna ,', c_resname
       call util_error(ERROR%STOP_ALL, error_message)
    end if

  end subroutine sub_check_basetype


  ! return idx
  subroutine sub_select_idx(ibasetype, c_name, idx)

    implicit none

    integer, intent(in)      :: ibasetype
    character(4), intent(in) :: c_name
    integer, intent(out)     :: idx

    idx = 0
    ! phosphate
    if(c_name == ' P  ' .or. c_name == 'P   ') then
       idx = 1
    else if (c_name == ' OP1' .or. c_name =='OP1 ' .or. &
         c_name == ' O1P' .or. c_name == 'O1P ') then
       idx = 2
    else if (c_name == ' OP2' .or. c_name == 'OP2 ' .or. &
         c_name ==  ' O2P' .or. c_name == 'O2P ') then
       idx = 3
       
    ! sugar
    else if (c_name == " O5'" .or. c_name == "O5' ") then
       idx = 4
    else if (c_name == " C5'" .or. c_name == "C5' ") then
       idx = 5
    else if (c_name == " C4'" .or. c_name == "C4' ") then
       idx = 6
    else if (c_name == " O4'" .or. c_name == "O4' " .or. &
         c_name == " O1'" .or. c_name == "O1' ") then
       idx = 7
    else if (c_name == " C3'" .or. c_name == "C3' ") then
       idx = 8
    else if (c_name == " O3'" .or. c_name == "O3' ") then
       idx = 9
    else if (c_name == " C2'" .or. c_name == "C2' ") then
       idx = 10
    else if (c_name == " C1'" .or. c_name == "C1' ") then
       idx = 11
       
    ! base
    else if(ibasetype == CHEMICALTYPE%DA) then
       if (c_name == ' N9 ' .or. c_name == 'N9  ') then
          idx = 12
       else if (c_name == ' C8 ' .or. c_name == 'C8  ') then
          idx = 13
       else if (c_name == ' N7 ' .or. c_name == 'N7  ') then
          idx = 14
       else if (c_name == ' C5 ' .or. c_name == 'C5  ') then
          idx = 15
       else if (c_name == ' C6 ' .or. c_name == 'C6  ') then
          idx = 16
       else if (c_name == ' N6 ' .or. c_name == 'N6  ') then
          idx = 17
       else if (c_name == ' N1 ' .or. c_name == 'N1  ') then
          idx = 18
       else if (c_name == ' C2 ' .or. c_name == 'C2  ') then
          idx = 19
       else if (c_name == ' N3 ' .or. c_name == 'N3  ') then
          idx = 20
       else if (c_name == ' C4 ' .or. c_name == 'C4  ') then
          idx = 21
       end if

    else if(ibasetype == CHEMICALTYPE%DG) then
       if (c_name == ' N9 ' .or. c_name == 'N9  ') then
          idx = 12
       else if (c_name == ' C8 ' .or. c_name == 'C8  ') then
          idx = 13
       else if (c_name == ' N7 ' .or. c_name == 'N7  ') then
          idx = 14
       else if (c_name == ' C5 ' .or. c_name == 'C5  ') then
          idx = 15
       else if (c_name == ' C6 ' .or. c_name == 'C6  ') then
          idx = 16
       else if (c_name == ' O6 ' .or. c_name == 'O6  ') then
          idx = 17
       else if (c_name == ' N1 ' .or. c_name == 'N1  ') then
          idx = 18
       else if (c_name == ' C2 ' .or. c_name == 'C2  ') then
          idx = 19
       else if (c_name == ' N2 ' .or. c_name == 'N2  ') then
          idx = 20
       else if (c_name == ' N3 ' .or. c_name == 'N3  ') then
          idx = 21
       else if (c_name == ' C4 ' .or. c_name == 'C4  ') then
          idx = 22
       end if

    else if(ibasetype == CHEMICALTYPE%DT) then
       if (c_name == ' N1 ' .or. c_name == 'N1  ') then
          idx = 12
       else if (c_name == ' C2 ' .or. c_name == 'C2  ') then
          idx = 13
       else if (c_name == ' O2 ' .or. c_name == 'O2  ') then
          idx = 14
       else if (c_name == ' N3 ' .or. c_name == 'N3  ') then
          idx = 15
       else if (c_name == ' C4 ' .or. c_name == 'C4  ') then
          idx = 16
       else if (c_name == ' O4 ' .or. c_name == 'O4  ') then
          idx = 17
       else if (c_name == ' C5 ' .or. c_name == 'C5  ') then
          idx = 18
       else if (c_name == ' C7 ' .or. c_name == 'C7  ') then
          idx = 19
       else if (c_name == ' C6 ' .or. c_name == 'C6  ') then
          idx = 20
       end if

    else if(ibasetype == CHEMICALTYPE%DC) then
       if (c_name == ' N1 ' .or. c_name == 'N1  ') then
          idx = 12
       else if (c_name == ' C2 ' .or. c_name == 'C2  ') then
          idx = 13
       else if (c_name == ' O2 ' .or. c_name == 'O2  ') then
          idx = 14
       else if (c_name == ' N3 ' .or. c_name == 'N3  ') then
          idx = 15
       else if (c_name == ' C4 ' .or. c_name == 'C4  ') then
          idx = 16
       else if (c_name == ' N4 ' .or. c_name == 'N4  ') then
          idx = 17
       else if (c_name == ' C5 ' .or. c_name == 'C5  ') then
          idx = 18
       else if (c_name == ' C6 ' .or. c_name == 'C6  ') then
          idx = 19
       end if
    end if

    if(idx == 0) then
       write(error_message,*) 'Error: undefined atom in read_pdbatom_dna', ibasetype, c_name, idx
       call util_error(ERROR%STOP_ALL, error_message)
    end if

  end subroutine sub_select_idx

end subroutine read_pdbatom_dna

