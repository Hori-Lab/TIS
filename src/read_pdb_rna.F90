! read_pdb_rna
!> @brief Read the PDB file for RNA

subroutine read_pdb_rna(lun,      & ![i ] target I/O unit
                        nunit,    & ![io] number of unit
                        nmp,      & ![io] number of mp (mass point)
                        nres,     & ![io] the number of residue
                        lunit2mp, & ![ o] correspondence list (unit -> mp)
                        ires_mp,  & ![ o] residue index       (mp   -> res )
                        cmp2seq,  & ![ o] correspondence list (mp   -> seq )
                        cmp2atom, & ![ o] correspondence list (mp   -> atom)
                        imp2type, & ![ o] correspondence list (mp   -> type)
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
  integer, intent(out)      :: iatomnum(MXMP)
  real(PREC), intent(out)   :: xyz(SDIM, MXATOM_MP, MXMP)

  ! ---------------------------------------------------------------------
  ! local variables
  logical :: flg_reading
  logical :: flg_phosphate, flg_sugar, flg_base
  logical :: flg_tmp_O3, flg_next_O3
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
  real(PREC) :: tmp_xyz(SDIM, MXATOM_MP, 3)
  real(PREC) :: tmp_xyz_O3(SDIM)
  real(PREC) :: next_xyz_O3(SDIM)
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
  flg_next_O3   = .false.
  flg_tmp_O3    = .false.
  tmp_iatomnum(:) = 0
  tmp_xyz(1:SDIM,1:MXATOM_MP,1:3) = INVALID_VALUE
  tmp_xyz_O3(1:SDIM) = INVALID_VALUE
  next_xyz_O3(1:SDIM) = INVALID_VALUE

  rewind(lun)

  do
     read (lun, '(a72)', iostat = input_status) char72
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: input error in read_pdb_rna'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     ! chaine terminates
     if (flg_reading .and. (char72(1:3) == 'TER' .or. char72(1:3) == 'END') ) then
        call chain_terminate()
        flg_reading   = .false.
        flg_phosphate = .false.
        flg_sugar     = .false.
        flg_base      = .false.
        flg_next_O3   = .false.
        flg_tmp_O3    = .false.
        c_icode_save   = '-'
        c_chainid_save = '-'
        i_resseq_save = 0
        tmp_iatomnum(:) = 0
        tmp_xyz(1:SDIM,1:MXATOM_MP,1:3) = INVALID_VALUE
        tmp_xyz_O3(1:SDIM) = INVALID_VALUE
        next_xyz_O3(1:SDIM) = INVALID_VALUE
        cycle
     end if

     if(char72(1:4) == 'ATOM') then
        read (char72, '(6x,5x,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2)',iostat=input_status) &
             c_name, c_altloc, c_resname, c_chainid, i_resseq, c_icode, &
             coord(1),coord(2),coord(3), occupancy, tempfactor

        !=== check ======================
        if(input_status > 0) then
           error_message = 'Error: cannot read pdb file in read_pdb_rna'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        call sub_check_residue()

        if(c_altloc /= " " .AND. c_altloc /= "A") then
           error_message = 'Error: c_altloc in read_pdb_rna'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        ! hydrogen atom
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
           error_message = 'Error: c_chainid is not receptible in read_pdb_rna'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        ! New residue
        if (i_resseq /= i_resseq_save .or. c_icode /= c_icode_save) then

           if (i_resseq /= i_resseq_save + 1) then
              write(error_message,*) &
              'Warning: chain is broken in read_pdb_rna. unit=',iunit,&
              ' residueID=',i_resseq
              call util_error(ERROR%WARN_ALL, error_message)
           endif

           i_resseq_save = i_resseq
           c_icode_save  = c_icode

           ! for previous residue
           if (flg_phosphate) then 
              cmp2seq(imp)    = c_resname_save
              cmp2atom(imp)   = ' P  '
              imp2type(imp)   = MPTYPE%RNA_PHOS
              ires_mp(imp)    = ires

              iatomnum(imp)   = tmp_iatomnum(RES%PHOSPHATE)
              xyz(1:SDIM,1:MXATOM_MP,imp)   = tmp_xyz(1:SDIM,1:MXATOM_MP,RES%PHOSPHATE)
              if (flg_next_O3) then
                 xyz(1:SDIM,1,imp)   = next_xyz_O3(1:SDIM)
                 iatomnum(imp) = iatomnum(imp) + 1
              endif
              xyz(:,:,imp+1) = tmp_xyz(:,:,RES%SUGAR)
              xyz(:,:,imp+2) = tmp_xyz(:,:,RES%BASE)

              imp = imp + 1
              call check_maxmp()
           endif

           if (flg_sugar) then
              if (.not.flg_base) then
                 error_message = 'Error: base is needed. in read_pdb_rna'
                 call util_error(ERROR%STOP_ALL, error_message)
              endif
              cmp2seq(imp)    = c_resname_save
              cmp2seq(imp+1)  = c_resname_save
              cmp2atom(imp)   = ' S  '
              cmp2atom(imp+1) = ' B  '
              imp2type(imp)   = MPTYPE%RNA_SUGAR
              imp2type(imp+1) = MPTYPE%RNA_BASE
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
              error_message = 'Error: failed in read_pdb_rna (1)'
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           flg_next_O3 = .false.
           if (flg_tmp_O3) then
              next_xyz_O3(1:SDIM) = tmp_xyz_O3(1:SDIM)
              flg_next_O3 = .true.
           endif

           c_resname_save = c_resname
           ires = ires + 1 

           tmp_iatomnum(:) = 0
           tmp_xyz(1:SDIM,1:MXATOM_MP,1:3) = INVALID_VALUE
           tmp_xyz_O3(1:SDIM) = INVALID_VALUE
           flg_phosphate = .false.
           flg_sugar     = .false.
           flg_base      = .false.
           flg_tmp_O3    = .false.

        ! Not new residue
        else
           if (c_resname /= c_resname_save) then
              write(error_message,*) &
              'Error: failed in read_pdb_rna. c_resname /= c_resname_save;', &
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

        ! O3' --> phosphate (next residue)
        case (RES%PHOSPHATE_NEXT)
           flg_tmp_O3 = .true.
           tmp_xyz_O3(:) = coord(:)

        case default
           error_message = 'Error: failed in read_pdb_rna (2)'
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
           error_message = 'Error: phosphate is needed. in read_pdb_rna'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        if (.NOT. flg_base) then
           error_message = 'Error: base is needed. in read_pdb_rna'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        cmp2seq(imp)    = c_resname_save
        cmp2seq(imp+1)  = c_resname_save
        cmp2seq(imp+2)  = c_resname_save
        cmp2atom(imp)   = ' P  '
        cmp2atom(imp+1) = ' S  '
        cmp2atom(imp+2) = ' B  '
        imp2type(imp)   = MPTYPE%RNA_PHOS
        imp2type(imp+1) = MPTYPE%RNA_SUGAR
        imp2type(imp+2) = MPTYPE%RNA_BASE
        ires_mp(imp)    = ires
        ires_mp(imp+1)  = ires
        ires_mp(imp+2)  = ires
   
        iatomnum(imp)   = tmp_iatomnum(RES%PHOSPHATE)
        iatomnum(imp+1) = tmp_iatomnum(RES%SUGAR)
        iatomnum(imp+2) = tmp_iatomnum(RES%BASE)
        xyz(:,:,imp)    = tmp_xyz(:,:,RES%PHOSPHATE)
        xyz(:,:,imp+1)  = tmp_xyz(:,:,RES%SUGAR)
        xyz(:,:,imp+2)  = tmp_xyz(:,:,RES%BASE)
        if (flg_next_O3) then
           xyz(:,1,imp)    = next_xyz_O3(:)
           iatomnum(imp) = iatomnum(imp) + 1
        endif
   
        imp = imp + 2
        nmp = imp
   
     else if (flg_phosphate) then
        cmp2seq(imp)    = c_resname_save
        cmp2atom(imp)   = ' P  '
        imp2type(imp)   = MPTYPE%RNA_PHOS
        ires_mp(imp)    = ires
   
        iatomnum(imp)   = tmp_iatomnum(RES%PHOSPHATE)
        xyz(:,:,imp)    = tmp_xyz(:,:,RES%PHOSPHATE)
        if (flg_next_O3) then
           xyz(:,1,imp)    = next_xyz_O3(:)
           iatomnum(imp) = iatomnum(imp) + 1
        endif
   
        nmp = imp
   
     else if (flg_base) then
        !! Warning
        write(error_message,*) 'Warning: base is isolated. imp =',imp
        call util_error(ERROR%WARN_ALL, error_message)
   
        cmp2seq(imp)    = c_resname_save
        cmp2atom(imp)   = ' B  '
        imp2type(imp)   = MPTYPE%RNA_BASE
        ires_mp(imp)    = ires
   
        iatomnum(imp) = tmp_iatomnum(RES%BASE)
        xyz(:,:,imp)  = tmp_xyz(:,:,RES%BASE)
        if (flg_next_O3) then
           xyz(:,1,imp)    = next_xyz_O3(:)
           iatomnum(imp) = iatomnum(imp) + 1
        endif
   
        nmp = imp
   
     else
        error_message = 'Error: failed in read_pdb_rna (3)'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     lunit2mp(2, iunit) = nmp

  endsubroutine chain_terminate
  
  subroutine sub_select_idx()
     implicit none

     selectcase (c_name)
     ! phosphate
     case (' OP3')
        type_res = RES%PHOSPHATE
        idx_pos = 1
     case (' O3P')
        type_res = RES%PHOSPHATE
        idx_pos = 1
     case (' P  ')
        flg_phosphate = .true.
        type_res = RES%PHOSPHATE
        idx_pos = 2
     case (' OP1') 
        type_res = RES%PHOSPHATE
        idx_pos = 3
     case (' O1P') 
        type_res = RES%PHOSPHATE
        idx_pos = 3
     case (' OP2')
        type_res = RES%PHOSPHATE
        idx_pos = 4
     case (' O2P')
        type_res = RES%PHOSPHATE
        idx_pos = 4
     case (" O5'")
        type_res = RES%PHOSPHATE
        idx_pos = 5

     ! phosphate (next residue)
     case (" O3'")
        type_res = RES%PHOSPHATE_NEXT
        idx_pos = 1

     ! sugar
     case (" C5'")
        type_res = RES%SUGAR
        idx_pos = 1
     case (" C4'")
        type_res = RES%SUGAR
        idx_pos = 2
     case (" O4'")
        type_res = RES%SUGAR
        idx_pos = 3
     case (" C3'")
        type_res = RES%SUGAR
        idx_pos = 4
     case (" C2'")
        type_res = RES%SUGAR
        idx_pos = 5
     case (" O2'")
        type_res = RES%SUGAR
        idx_pos = 6
     case (" C1'")
        type_res = RES%SUGAR
        idx_pos = 7
     case ("O2'C")  ! atom name is " CM2" in PDB, one atom of residue "OMU"
        type_res = RES%SUGAR
        idx_pos = 8
     case (" C1X")
        type_res = RES%SUGAR
        idx_pos = 9
     case (" O1X")
        type_res = RES%SUGAR
        idx_pos = 10 
     case (" C2X")
        type_res = RES%SUGAR
        idx_pos = 11 
     case (" O2X")
        type_res = RES%SUGAR
        idx_pos = 12 
     case (" C3X")
        type_res = RES%SUGAR
        idx_pos = 13 
     case (" O3X")
        type_res = RES%SUGAR
        idx_pos = 14 
     case (" C4X")
        type_res = RES%SUGAR
        idx_pos = 15 
     case (" C5X")
        type_res = RES%SUGAR
        idx_pos = 16 
     case (" O5X")
        type_res = RES%SUGAR
        idx_pos = 17 
     case (" PX ")
        type_res = RES%SUGAR
        idx_pos = 18 
     case ("O1PX")
        type_res = RES%SUGAR
        idx_pos = 19 
     case ("O2PX")
        type_res = RES%SUGAR
        idx_pos = 20 
     case ("O3PX")
        type_res = RES%SUGAR
        idx_pos = 21

     ! base
     case (" N1 ")
        type_res = RES%BASE
        idx_pos = 1
     case (" CM1")
        type_res = RES%BASE
        idx_pos = 2
     case (" C2 ")
        type_res = RES%BASE
        idx_pos = 3
     case (" N2 ")
        type_res = RES%BASE
        idx_pos = 4 
     case (" O2 ")
        type_res = RES%BASE
        idx_pos = 5 
     case (" CM2")
        type_res = RES%BASE
        idx_pos = 6 
     case (" N3 ")
        type_res = RES%BASE
        idx_pos = 7 
     case (" C3 ")
        type_res = RES%BASE
        idx_pos = 8 
     case (" C4 ")
        type_res = RES%BASE
        idx_pos = 9 
     case (" N4 ")
        type_res = RES%BASE
        idx_pos = 10 
     case (" O4 ")
        type_res = RES%BASE
        idx_pos = 11 
     case (" C5 ")
        type_res = RES%BASE
        idx_pos = 12 
     case (" CM5")
        type_res = RES%BASE
        idx_pos = 13
     case (" C5M")
        type_res = RES%BASE
        idx_pos = 14
     case (" C6 ")
        type_res = RES%BASE
        idx_pos = 15 
     case (" N6 ")
        type_res = RES%BASE
        idx_pos = 16
     case (" O6 ")
        type_res = RES%BASE
        idx_pos = 17
     case (" N7 ")
        type_res = RES%BASE
        idx_pos = 18 
     case (" CM7")
        type_res = RES%BASE
        idx_pos = 19 
     case (" C8 ")
        type_res = RES%BASE
        idx_pos = 20
     case (" N9 ")
        type_res = RES%BASE
        idx_pos = 21
     case (" C10")
        type_res = RES%BASE
        idx_pos = 22
     case (" C11")
        type_res = RES%BASE
        idx_pos = 23
     case (" C12")
        type_res = RES%BASE
        idx_pos = 24
     case (" C13")
        type_res = RES%BASE
        idx_pos = 25
     case (" C14")
        type_res = RES%BASE
        idx_pos = 26
     case (" C15")
        type_res = RES%BASE
        idx_pos = 27
     case (" C16")
        type_res = RES%BASE
        idx_pos = 28
     case (" O17")
        type_res = RES%BASE
        idx_pos = 29
     case (" O18")
        type_res = RES%BASE
        idx_pos = 30 
     case (" C19")
        type_res = RES%BASE
        idx_pos = 31
     case (" N20")
        type_res = RES%BASE
        idx_pos = 32
     case (" C21")
        type_res = RES%BASE
        idx_pos = 33
     case (" O22")
        type_res = RES%BASE
        idx_pos = 34
     case (" O23")
        type_res = RES%BASE
        idx_pos = 35
     case (" C24")
        type_res = RES%BASE
        idx_pos = 36
     case (" C3U")
        type_res = RES%BASE
        idx_pos = 37
     case (" CM4")
        type_res = RES%BASE
        idx_pos = 38
     case (" C9 ")
        type_res = RES%BASE
        idx_pos = 39
     case (" C1 ")
        type_res = RES%BASE
        idx_pos = 40
     case ("BR  ")
        type_res = RES%BASE
        idx_pos = 41
     case (" CM6")
        type_res = RES%BASE
        idx_pos = 42
     case (" O10")
        type_res = RES%BASE
        idx_pos = 43
     case (" N11")
        type_res = RES%BASE
        idx_pos = 44
     case (" O14")
        type_res = RES%BASE
        idx_pos = 45
     case (" ODA")
        type_res = RES%BASE
        idx_pos = 46
     case (" ODB")
        type_res = RES%BASE
        idx_pos = 47
     case (" O13")
        type_res = RES%BASE
        idx_pos = 48
     case (" C7 ")
        type_res = RES%BASE
        idx_pos = 49
     case (" F5 ")
        type_res = RES%BASE
        idx_pos = 50
     case (" S10")
        type_res = RES%BASE
        idx_pos = 51
     case default
        write(error_message,*) 'Error: undefined atom in read_pdb_rna;',c_name
        call util_error(ERROR%STOP_ALL, error_message)
     endselect 

  endsubroutine sub_select_idx

  subroutine sub_check_residue()
     ! acceptable residues
     if ( c_resname == '  A' .OR. &
          c_resname == ' RA' .OR. &
          c_resname == 'A  ' .OR. &
          c_resname == 'RA ' .OR. &
          c_resname == 'RA5' .OR. &
          c_resname == 'RA3' .OR. &
          c_resname == '  G' .OR. &
          c_resname == ' RG' .OR. &
          c_resname == 'G  ' .OR. &
          c_resname == 'RG ' .OR. &
          c_resname == 'RG5' .OR. &
          c_resname == 'RG3' .OR. &
          c_resname == 'OMG' .OR. &
          c_resname == '2MG' .OR. &
          c_resname == 'M2G' .OR. &
          c_resname == 'YYG' .OR. &
          c_resname == '1MG' .OR. &
          c_resname == '7MG' .OR. &
          c_resname == '  U' .OR. &
          c_resname == ' RU' .OR. &  ! AMBER
          c_resname == 'U  ' .OR. &
          c_resname == 'RU ' .OR. &  ! AMBER
          c_resname == 'RU5' .OR. &  ! AMBER
          c_resname == 'RU3' .OR. &  ! AMBER
          c_resname == 'OMU' .OR. &
          c_resname == 'PSU' .OR. &
          c_resname == 'H2U' .OR. &
          c_resname == '  C' .OR. &
          c_resname == ' RC' .OR. &  ! AMBER
          c_resname == 'C  ' .OR. &
          c_resname == 'RC ' .OR. &  ! AMBER
          c_resname == 'RC5' .OR. &  ! AMBER
          c_resname == 'RC3' .OR. &  ! AMBER
          c_resname == 'OMC' .OR. &
          c_resname == '5MC' .OR. &
          c_resname == '5MU' .OR. &
          c_resname == '1MA' .OR. &
          c_resname == 'UR3' .OR. &
          c_resname == 'MA6' .OR. &
          c_resname == '4OC' .OR. &
          c_resname == '  I' .OR. &
          c_resname == 'PYY' .OR. &
          c_resname == 'CB2' .OR. &
          c_resname == 'AET' .OR. &
          c_resname == 'QUO' .OR. &
          c_resname == 'FHU' .OR. &
          c_resname == 'T6A' .OR. &
          c_resname == 'RIA' .OR. &
          c_resname == 'MIA' .OR. &
          c_resname == 'PRF' ) then

        continue
     else
        write(error_message,*) 'Error: non rna base in read_pdb_rna ,',c_resname
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  endsubroutine sub_check_residue

  subroutine check_maxmp()
     if (imp > MXMP) then
        write(error_message,*) 'Error: imp > MXMP, in read_pdb_rna. imp=',imp,' MXMP=',MXMP
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  endsubroutine check_maxmp

endsubroutine read_pdb_rna
