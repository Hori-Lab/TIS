! read_pdb_dna2
!> @brief Read the PDB file for DNA2

subroutine read_pdb_dna2(lun,      & ![i ] target I/O unit
                         nunit,    & ![io] number of unit
                         nmp,      & ![io] number of mp (mass point)
                         nres,     & ![io] the number of residue
                         lunit2mp, & ![ o] correspondence list (unit -> mp)
                         ires_mp,  & ![ o] residue index       (mp   -> res )
                         cmp2seq,  & ![ o] correspondence list (mp   -> seq )
                         cmp2atom, & ![ o] correspondence list (mp   -> atom)
                         imp2type, & ![ o] correspondence list (mp   -> type)
                         imp2base, & ![ o] correspondence list (mp   -> base)
                         iatomnum, & ![ o] the number of atom  (mp   -> #   )
                         xyz       ) ![ o] coordinate of all atom

  use const_physical
  use const_maxsize
  use const_index
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in)       :: lun
  integer, intent(inout)    :: nunit, nmp, nres
  integer, intent(out)      :: lunit2mp(2, MXUNIT)
  integer, intent(out)      :: ires_mp(MXMP)
  character(3), intent(out) :: cmp2seq(MXMP)
  character(4), intent(out) :: cmp2atom(MXMP)
  integer, intent(out)      :: imp2type(MXMP)
  integer, intent(out)      :: imp2base(MXMP)
  integer, intent(out)      :: iatomnum(MXMP)
  real(PREC), intent(out)   :: xyz(SPACE_DIM, MXATOM_MP, MXMP)

  ! ---------------------------------------------------------------------
  ! local variables
  logical :: flg_reading
  logical :: flg_phosphate, flg_sugar, flg_base
  integer :: input_status
  integer :: imp, iunit, ires
  integer :: idx_pos
  character(72) :: char72
  character(CARRAY_MSG_ERROR) :: error_message
  type residue
     integer :: PHOSPHATE
     integer :: SUGAR
     integer :: BASE
     integer :: PHOSPHATE_NEXT
  endtype residue
  type(residue), parameter :: RES = residue(1,2,3,4)
  integer    :: type_res

  ! 1:Phospahte, 2:Sugar, 3:Base
  real(PREC) :: tmp_xyz(SPACE_DIM, MXATOM_MP, 3)
  integer    :: tmp_iatomnum(3)

  ! PDB ATOM Record Format (See http://www.wwpdb.org/docs.html )
  character(4) :: c_name       ! 13-16  Atom name.
  character(1) :: c_altloc     ! 17     Alternate location indicator.
  character(3) :: c_resname    ! 18-20  Residue name.
  character(1) :: c_chainid    ! 22     Chain identifier.
  integer      :: i_resseq     ! 23-26  Residue sequence number.
  character(1) :: c_icode      ! 27     Code for insertion of residues.
  real(PREC)   :: coord(3)     ! 31-54  Orthogonal coordinates in Angstroms.
  real(PREC)   :: occupancy    ! 55-60  Occupancy. (Not used)
  real(PREC)   :: tempfactor   ! 61-66  Temperature factor. (Not used)

  ! To store previous
  character(3) :: c_resname_save
  character(1) :: c_chainid_save
  integer      :: i_resseq_save
  character(1) :: c_icode_save

  ! ---------------------------------------------------------------------
  imp = nmp
  iunit = nunit - 1
  ires = nres

  c_icode_save   = '-'
  c_chainid_save = '-'
  i_resseq_save = 0
  flg_reading   = .false.
  flg_phosphate = .false.
  flg_sugar     = .false.
  flg_base      = .false.
  tmp_iatomnum(:) = 0
  tmp_xyz(1:SPACE_DIM,1:MXATOM_MP,1:3) = INVALID_VALUE

  rewind(lun)

  do
     read (lun, '(a72)', iostat = input_status) char72
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: input error in read_pdb_dna2'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     ! chain terminates
     if (flg_reading .and. (char72(1:3) == 'TER' .or. char72(1:3) == 'END') ) then
        call chain_terminate()
        flg_reading   = .false.
        flg_phosphate = .false.
        flg_sugar     = .false.
        flg_base      = .false.
        c_icode_save   = '-'
        c_chainid_save = '-'
        i_resseq_save = 0
        tmp_iatomnum(:) = 0
        tmp_xyz(1:SPACE_DIM,1:MXATOM_MP,1:3) = INVALID_VALUE
        cycle
     end if

     if(char72(1:4) == 'ATOM') then
        read (char72, '(6x,5x,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2)',iostat=input_status) &
             c_name, c_altloc, c_resname, c_chainid, i_resseq, c_icode, &
             coord(1),coord(2),coord(3), occupancy, tempfactor

        if (c_resname == ' DA' .or. c_resname == 'DA ') then
           c_resname = 'DA '
        else if (c_resname == ' DT' .or. c_resname == 'DT ') then
           c_resname = 'DT '
        else if (c_resname == ' DG' .or. c_resname == 'DG ') then
           c_resname = 'DG '
        else if (c_resname == ' DC' .or. c_resname == 'DC ') then
           c_resname = 'DC '
        end if
        !write(*, *) 'read_pdb_dna2: c_resname: ', c_resname

        !=== check ======================
        if(input_status > 0) then
           error_message = 'Error: cannot read pdb file in read_pdb_dna2'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        call sub_check_residue()

        if(c_altloc /= " " .AND. c_altloc /= "A") then
           error_message = 'Error: c_altloc in read_pdb_dna2'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        ! Ignoring hydrogen atom
        if ( c_name(1:1) == 'H'   .OR. c_name(2:2) == 'H'     .OR. &
             c_name(1:3) == " 1H" .OR. c_name(1:3) == " 2H" ) then
           cycle
        endif
        !================================

        ! New chain
        if (.not. flg_reading) then
           c_resname_save = c_resname
           c_chainid_save = c_chainid
           i_resseq_save = i_resseq
           c_icode_save = c_icode
           flg_reading = .true.
           imp = imp + 1
           call check_maxmp()
           ires = ires + 1
           iunit = iunit + 1

        ! otherwise must stay current chain
        else if (c_chainid /= c_chainid_save) then
           error_message = 'Error: c_chainid is not receptible in read_pdb_dna2'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! New residue
        if (i_resseq /= i_resseq_save .or. c_icode /= c_icode_save) then

           if (i_resseq /= i_resseq_save + 1) then
              write(error_message,*) &
              'Warning: chain is broken in read_pdb_dna2. unit=',iunit,&
              ' residueID=',i_resseq
              call util_error(ERROR%WARN_ALL, error_message)
           endif

           i_resseq_save = i_resseq
           c_icode_save  = c_icode

           ! for previous residue
           if (flg_phosphate) then
              cmp2seq(imp)    = c_resname_save
              cmp2atom(imp)   = 'DP  '
              imp2type(imp)   = MPTYPE%DNA2_PHOS
              imp2base(imp)   = BASETYPE%P
              ires_mp(imp)    = ires

              iatomnum(imp)   = tmp_iatomnum(RES%PHOSPHATE)
              xyz(1:SPACE_DIM,1:MXATOM_MP,imp)   = tmp_xyz(1:SPACE_DIM,1:MXATOM_MP,RES%PHOSPHATE)
              xyz(:,:,imp+1) = tmp_xyz(:,:,RES%SUGAR)
              xyz(:,:,imp+2) = tmp_xyz(:,:,RES%BASE)

              imp = imp + 1
              call check_maxmp()
           endif

           if (flg_sugar) then
              if (.not.flg_base) then
                 error_message = 'Error: base is needed. in read_pdb_dna2'
                 call util_error(ERROR%STOP_ALL, error_message)
              endif
              cmp2seq(imp)    = c_resname_save
              cmp2seq(imp+1)  = c_resname_save
              cmp2atom(imp)   = 'DS  '
              cmp2atom(imp+1) = 'DB  '
              imp2type(imp)   = MPTYPE%DNA2_SUGAR
              imp2base(imp)   = BASETYPE%S
              imp2type(imp+1) = MPTYPE%DNA2_BASE
              if (c_resname_save == ' DA' .or. c_resname_save == 'DA ') then
                 imp2base(imp+1) = BASETYPE%A
              else if (c_resname_save == ' DT' .or. c_resname_save == 'DT ') then
                 imp2base(imp+1) = BASETYPE%T
              else if (c_resname_save == ' DG' .or. c_resname_save == 'DG ') then
                 imp2base(imp+1) = BASETYPE%G
              else if (c_resname_save == ' DC' .or. c_resname_save == 'DC ') then
                 imp2base(imp+1) = BASETYPE%C
              end if
              ires_mp(imp)    = ires
              ires_mp(imp+1)  = ires

              iatomnum(imp)   = tmp_iatomnum(RES%SUGAR)
              iatomnum(imp+1) = tmp_iatomnum(RES%BASE)
              xyz(:,:,imp)   = tmp_xyz(:,:,RES%SUGAR)
              xyz(:,:,imp+1) = tmp_xyz(:,:,RES%BASE)

              imp = imp + 2
              call check_maxmp()
           endif

           if ((.not. flg_phosphate) .AND. (.not. flg_sugar)) then
              error_message = 'Error: failed in read_pdb_dna2 (1)'
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           c_resname_save = c_resname
           ires = ires + 1

           tmp_iatomnum(:) = 0
           tmp_xyz(1:SPACE_DIM,1:MXATOM_MP,1:3) = INVALID_VALUE
           flg_phosphate = .false.
           flg_sugar     = .false.
           flg_base      = .false.

        ! Not new residue
        else
           if (c_resname /= c_resname_save) then
              write(error_message,*) &
              'Error: failed in read_pdb_dna2. c_resname /= c_resname_save;', &
              ' iunit=',iunit, ' residueID=',i_resseq, &
              ' c_resname=',c_resname,&
              ' c_resname_save=',c_resname_save
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        endif  ! (i_resseq /= i_resseq_save .or. c_icode /= c_icode_save)

        call sub_select_idx()

        selectcase (type_res)
        ! phosphate
        case (RES%PHOSPHATE)
           flg_phosphate = .true.
           tmp_iatomnum(RES%PHOSPHATE) = tmp_iatomnum(RES%PHOSPHATE) + 1
           tmp_xyz(:, idx_pos, RES%PHOSPHATE) = coord(:)

        ! sugar
        case (RES%SUGAR)
           flg_sugar = .true.
           tmp_iatomnum(RES%SUGAR) = tmp_iatomnum(RES%SUGAR) + 1
           tmp_xyz(:, idx_pos, RES%SUGAR) = coord(:)
        ! base
        case (RES%BASE)
           flg_base = .true.
           tmp_iatomnum(RES%BASE) = tmp_iatomnum(RES%BASE) + 1
           tmp_xyz(:, idx_pos, RES%BASE) = coord(:)

        case default
           error_message = 'Error: failed in read_pdb_dna2 (2)'
           call util_error(ERROR%STOP_ALL, error_message)
        endselect

     end if ! (char72(1:4) == 'ATOM')
  end do ! do

  if (flg_reading) then
     call chain_terminate()
  endif

  nunit = iunit
  nres = ires



!###########################################################################
contains

subroutine chain_terminate()
     implicit none

     if (flg_sugar) then
        if (.NOT. flg_phosphate) then
           error_message = 'Error: phosphate is needed. in read_pdb_dna2'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        if (.NOT. flg_base) then
           error_message = 'Error: base is needed. in read_pdb_dna2'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        cmp2seq(imp)    = c_resname_save
        cmp2seq(imp+1)  = c_resname_save
        cmp2seq(imp+2)  = c_resname_save
        cmp2atom(imp)   = 'DP  '
        cmp2atom(imp+1) = 'DS  '
        cmp2atom(imp+2) = 'DB  '
        imp2type(imp)   = MPTYPE%DNA2_PHOS
        imp2base(imp)   = BASETYPE%P
        imp2type(imp+1) = MPTYPE%DNA2_SUGAR
        imp2base(imp+1) = BASETYPE%S
        imp2type(imp+2) = MPTYPE%DNA2_BASE
        if (c_resname_save == ' DA' .or. c_resname_save == 'DA ') then
           imp2base(imp+2) = BASETYPE%A
        else if (c_resname_save == ' DT' .or. c_resname_save == 'DT ') then
           imp2base(imp+2) = BASETYPE%T
        else if (c_resname_save == ' DG' .or. c_resname_save == 'DG ') then
           imp2base(imp+2) = BASETYPE%G
        else if (c_resname_save == ' DC' .or. c_resname_save == 'DC ') then
           imp2base(imp+2) = BASETYPE%C
        end if

        ires_mp(imp)    = ires
        ires_mp(imp+1)  = ires
        ires_mp(imp+2)  = ires

        iatomnum(imp)   = tmp_iatomnum(RES%PHOSPHATE)
        iatomnum(imp+1) = tmp_iatomnum(RES%SUGAR)
        iatomnum(imp+2) = tmp_iatomnum(RES%BASE)
        xyz(:,:,imp)    = tmp_xyz(:,:,RES%PHOSPHATE)
        xyz(:,:,imp+1)  = tmp_xyz(:,:,RES%SUGAR)
        xyz(:,:,imp+2)  = tmp_xyz(:,:,RES%BASE)

        imp = imp + 2
        nmp = imp

     else if (flg_phosphate) then
        cmp2seq(imp)    = c_resname_save
        cmp2atom(imp)   = 'DP  '
        imp2type(imp)   = MPTYPE%DNA2_PHOS
        imp2base(imp)   = BASETYPE%P
        ires_mp(imp)    = ires
        iatomnum(imp)   = tmp_iatomnum(RES%PHOSPHATE)
        xyz(:,:,imp)    = tmp_xyz(:,:,RES%PHOSPHATE)

        nmp = imp

     else if (flg_base) then
        !! Warning
        write(error_message,*) 'Warning: base is isolated. imp =',imp
        call util_error(ERROR%WARN_ALL, error_message)

        cmp2seq(imp)    = c_resname_save
        cmp2atom(imp)   = 'DB  '
        if (c_resname_save == ' DA' .or. c_resname_save == 'DA ') then
           imp2base(imp) = BASETYPE%A
        else if (c_resname_save == ' DT' .or. c_resname_save == 'DT ') then
           imp2base(imp) = BASETYPE%T
        else if (c_resname_save == ' DG' .or. c_resname_save == 'DG ') then
           imp2base(imp) = BASETYPE%G
        else if (c_resname_save == ' DC' .or. c_resname_save == 'DC ') then
           imp2base(imp) = BASETYPE%C
        end if

        imp2type(imp)   = MPTYPE%DNA2_BASE
        ires_mp(imp)    = ires

        iatomnum(imp) = tmp_iatomnum(RES%BASE)
        xyz(:,:,imp)  = tmp_xyz(:,:,RES%BASE)

        nmp = imp

     else
        error_message = 'Error: failed in read_pdb_dna2 (3)'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     lunit2mp(2, iunit) = nmp


  endsubroutine chain_terminate

  subroutine sub_select_idx()
     implicit none

     selectcase (c_name)
     ! phosphate
     case (' P  ')
        flg_phosphate = .true.
        type_res = RES%PHOSPHATE
        idx_pos = 1
     case ('P   ')
        flg_phosphate = .true.
        type_res = RES%PHOSPHATE
        idx_pos = 1
     case (' OP1')
        type_res = RES%PHOSPHATE
        idx_pos = 2
     case ('OP1 ')
        type_res = RES%PHOSPHATE
        idx_pos = 2
     case (' O1P')
        type_res = RES%PHOSPHATE
        idx_pos = 2
     case ('O1P ')
        type_res = RES%PHOSPHATE
        idx_pos = 2
     case (' OP2')
        type_res = RES%PHOSPHATE
        idx_pos = 3
     case ('OP2 ')
        type_res = RES%PHOSPHATE
        idx_pos = 3
     case (' O2P')
        type_res = RES%PHOSPHATE
        idx_pos = 3
     case ('O2P ')
        type_res = RES%PHOSPHATE
        idx_pos = 3

     ! sugar
     case (" S  ")
        type_res = RES%SUGAR
        idx_pos = 1
     case ("S   ")
        type_res = RES%SUGAR
        idx_pos = 1
     case (" O5'")
        type_res = RES%SUGAR
        idx_pos = 1
     case ("O5' ")
        type_res = RES%SUGAR
        idx_pos = 1
     case (" C5'")
        type_res = RES%SUGAR
        idx_pos = 2
     case ("C5' ")
        type_res = RES%SUGAR
        idx_pos = 2
     case (" C4'")
        type_res = RES%SUGAR
        idx_pos = 3
     case ("C4' ")
        type_res = RES%SUGAR
        idx_pos = 3
     case (" O1'")
        type_res = RES%SUGAR
        idx_pos = 4
     case ("O1' ")
        type_res = RES%SUGAR
        idx_pos = 4
     case (" O4'")
        type_res = RES%SUGAR
        idx_pos = 4
     case ("O4' ")
        type_res = RES%SUGAR
        idx_pos = 4
     case (" C3'")
        type_res = RES%SUGAR
        idx_pos = 5
     case ("C3' ")
        type_res = RES%SUGAR
        idx_pos = 5
     case (" O3'")
        type_res = RES%SUGAR
        idx_pos = 6
     case ("O3' ")
        type_res = RES%SUGAR
        idx_pos = 6
     case (" C2'")
        type_res = RES%SUGAR
        idx_pos = 7
     case ("C2' ")
        type_res = RES%SUGAR
        idx_pos = 7
     case (" C1'")
        type_res = RES%SUGAR
        idx_pos = 8
     case ("C1' ")
        type_res = RES%SUGAR
        idx_pos = 8

     ! base
     case (" B  ")
        type_res = RES%BASE
        idx_pos = 1
     case ("B   ")
        type_res = RES%BASE
        idx_pos = 1
     case (" N1 ")
        type_res = RES%BASE
        idx_pos = 1
     case ("N1  ")
        type_res = RES%BASE
        idx_pos = 1
     case (" N2 ")
        type_res = RES%BASE
        idx_pos = 2
     case ("N2  ")
        type_res = RES%BASE
        idx_pos = 2
     case (" N3 ")
        type_res = RES%BASE
        idx_pos = 3
     case ("N3  ")
        type_res = RES%BASE
        idx_pos = 3
     case (" N4 ")
        type_res = RES%BASE
        idx_pos = 4
     case ("N4  ")
        type_res = RES%BASE
        idx_pos = 4
     case (" N6 ")
        type_res = RES%BASE
        idx_pos = 5
     case ("N6  ")
        type_res = RES%BASE
        idx_pos = 5
     case (" N7 ")
        type_res = RES%BASE
        idx_pos = 6
     case ("N7  ")
        type_res = RES%BASE
        idx_pos = 6
     case (" N9 ")
        type_res = RES%BASE
        idx_pos = 7
     case ("N9  ")
        type_res = RES%BASE
        idx_pos = 7
     case (" C2 ")
        type_res = RES%BASE
        idx_pos = 8
     case ("C2  ")
        type_res = RES%BASE
        idx_pos = 8
     case (" C4 ")
        type_res = RES%BASE
        idx_pos = 9
     case ("C4  ")
        type_res = RES%BASE
        idx_pos = 9
     case (" C5 ")
        type_res = RES%BASE
        idx_pos = 10
     case ("C5  ")
        type_res = RES%BASE
        idx_pos = 10
     case (" C6 ")
        type_res = RES%BASE
        idx_pos = 11
     case ("C6  ")
        type_res = RES%BASE
        idx_pos = 11
     case (" C7 ")
        type_res = RES%BASE
        idx_pos = 12
     case ("C7  ")
        type_res = RES%BASE
        idx_pos = 12
     case (" C8 ")
        type_res = RES%BASE
        idx_pos = 13
     case ("C8  ")
        type_res = RES%BASE
        idx_pos = 13
     case (" O2 ")
        type_res = RES%BASE
        idx_pos = 14
     case ("O2  ")
        type_res = RES%BASE
        idx_pos = 14
     case (" O4 ")
        type_res = RES%BASE
        idx_pos = 15
     case ("O4  ")
        type_res = RES%BASE
        idx_pos = 15
     case (" O6 ")
        type_res = RES%BASE
        idx_pos = 16
     case ("O6  ")
        type_res = RES%BASE
        idx_pos = 16
     case default
        write(error_message,*) 'Error: undefined atom in read_pdb_dna2;',c_name
        call util_error(ERROR%STOP_ALL, error_message)
     endselect

  endsubroutine sub_select_idx

  subroutine sub_check_residue()
     ! acceptable residues
    if ( c_resname == ' DT' .OR. &
         c_resname == ' DA' .OR. &
         c_resname == ' DC' .OR. &
         c_resname == ' DG' .OR. &
         c_resname == 'DT ' .OR. &
         c_resname == 'DA ' .OR. &
         c_resname == 'DC ' .OR. &
         c_resname == 'DG ') then
       continue
    else
       write(error_message,*) 'Error: non DNA residue in read_pdb_dna2 ,', c_resname
       call util_error(ERROR%STOP_ALL, error_message)
    endif
  endsubroutine sub_check_residue

  subroutine check_maxmp()
    if (imp > MXMP) then
       write(error_message,*) 'Error: imp > MXMP, in read_pdb_dna2. imp=',imp,' MXMP=',MXMP
       call util_error(ERROR%STOP_ALL, error_message)
    endif
  endsubroutine check_maxmp

endsubroutine read_pdb_dna2
