!read_pdb_pro
!> @brief Reads structure information about proteins from PDB files.

subroutine read_pdb_pro(lun,      & ![i ] target I/O unit
                        nunit,    & ![io] the number of unit
                        nmp,      & ![io] the number of mp(mass point)
                        nres,     & ![io] the number of residue
                        lunit2mp, & ![ o] correspondence list (unit -> mp  )
                        ires_mp,  & ![ o] residue index       (mp   -> res )
                        cmp2seq,  & ![ o] correspondence list (mp   -> seq )
                        imp2type, & ![ o] correspondence list (mp   -> type)
                        iatomnum, & ![ o] the number of atom           (mp   -> #   )
                        xyz,      & ![ o] coordinate of all atom 
                        cname_ha)   ![ o] name of heavy atoms  ! aicg

  use const_maxsize
  use const_index
  use const_physical
  implicit none

  ! ---------------------------------------------------------------------
  integer,      intent(in)    :: lun
  integer,      intent(inout) :: nmp, nunit, nres
  integer,      intent(out)   :: lunit2mp(2, MXUNIT)
  integer,      intent(out)   :: ires_mp(MXMP)
  character(3), intent(out)   :: cmp2seq(MXMP)
  integer,      intent(out)   :: imp2type(MXMP)
  integer,      intent(out)   :: iatomnum(MXMP)
  real(PREC),   intent(out)   :: xyz(3, MXATOM_MP, MXMP)
  character(4), intent(out)   :: cname_ha(MXATOM_MP, MXMP)  ! aicg

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: iatom
  integer :: input_status
  integer :: imp, iunit, ires, ii, ii2, nha
  logical :: flg_reading         ! flag for 
  character(72) :: char72
  character(CARRAY_MSG_ERROR) :: error_message

  ! PDB ATOM Record Format (See http://www.wwpdb.org/docs.html )
  character(4) :: c_name       ! 13-16  Atom name.
  character(1) :: c_altloc     ! 17     Alternate location indicator.
  character(3) :: c_resname    ! 18-20  Residue name.
  character(1) :: c_chainid    ! 22     Chain identifier.
  integer      :: i_resseq     ! 23-26  Residue sequence number.
  character(1) :: c_icode      ! 27     Code for insertion of residues.
  real(PREC)   :: x, y, z      ! 31-54  Orthogonal coordinates in Angstroms.
  real(PREC)   :: occupancy    ! 55-60  Occupancy. (Not used)
  real(PREC)   :: tempfactor   ! 61-66  Temperature factor. (Not used)

  integer      :: i_resseq_save ! to store previous i_resseq
  character(1) :: c_icode_save  ! to store previous c_icode_save

  integer :: ifunc_seq2hanum
#ifdef _DEBUG
  write(*,*) '#### start read_pdb_pro'
#endif

  ! ---------------------------------------------------------------------
  flg_reading   = .false.
  i_resseq_save = 0
  c_icode_save  = ' '
  imp = nmp
  iunit = nunit
  ires = nres

  ! ---------------------------------------------------------------------
  rewind(lun)

  do
     read (lun, '(a72)', iostat = input_status) char72
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: input error in read_pdb_pro'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     if(flg_reading .and. (char72(1:3) == 'TER' .or.       &
                           char72(1:3) == 'END' )) then
        lunit2mp(2, iunit) = imp
        iunit = iunit + 1
        flg_reading = .false.
     end if

     if(char72(1:4) == 'ATOM') then
        read (char72, '(6x, 5x, 1x, a4, a1, a3, 1x, a1, i4, a1, 3x, 3f8.3, 2f6.2)', &
             iostat = input_status) &
             c_name, c_altloc, c_resname, &
             c_chainid, i_resseq, c_icode, &
             x, y, z, occupancy, tempfactor

        if(input_status > 0) then
           error_message = 'Error: cannot read pdb file in read_pdb_pro'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        
        if(c_altloc /= " " .and. c_altloc /= "A") then
           cycle
        end if

        ! if using H atom, this if statement should be coment out
        if(c_name(1:1) == 'H' .or. c_name(2:2) == 'H') then
           cycle
        end if
        
        ! New residue
        if((.not. flg_reading) .or. &
           i_resseq/=i_resseq_save .or. c_icode /= c_icode_save) then

           imp = imp + 1
           if (imp > MXMP) then
              write(error_message,*) 'Error: imp > MXMP, in read_pdb_pro. imp=',imp,' MXMP=',MXMP
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           ires = ires + 1
           ires_mp(imp) = ires
           iatomnum(imp) = 5
           cmp2seq(imp) = c_resname
           imp2type(imp) = MPTYPE%PRO

           i_resseq_save = i_resseq
           c_icode_save = c_icode
           flg_reading = .true.
           xyz(1:3, :, imp) = INVALID_VALUE
        end if

        if(c_name == ' N  ') then
           iatom = 1
     
        else if(c_name == ' CA ') then
           iatom = 2

        else if(c_name == ' C  ') then
           iatom = 3

        else if(c_name == ' O  ') then
           iatom = 4

        else if(c_name == ' CB ') then
           iatom = 5

        else
           iatomnum(imp) = iatomnum(imp) + 1
           if (iatomnum(imp) > MXATOM_MP) then
              write(error_message,*) 'Error: too many atom in ', c_resname,' of unit ',iunit
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           iatom = iatomnum(imp)
        end if

        xyz(1, iatom, imp) = x
        xyz(2, iatom, imp) = y
        xyz(3, iatom, imp) = z

        cname_ha(iatom, imp) = c_name  ! aicg
     end if
  end do

  lunit2mp(2, iunit) = imp
  if(iunit > 1) then
     if(lunit2mp(2, iunit) == lunit2mp(2, iunit - 1)) then
        iunit = iunit - 1
     end if
  end if

  ii_loop: do ii = nmp + 1, imp
     c_resname = cmp2seq(ii)

     if (xyz(1,1,ii) > INVALID_JUDGE) then  ! N
        write (*, '(a, i6, a, a3, a)') &
             'Waring: N atom does not exist in the', ii,'th residue(', c_resname,')'
        call util_error(ERROR%WARN_ALL, error_message)
     endif
     if (xyz(1,2,ii) > INVALID_JUDGE) then  ! CA
        write (*, '(a, i6, a, a3, a)') &
             'Waring: CA atom does not exist in the', ii,'th residue(', c_resname,')'
        call util_error(ERROR%WARN_ALL, error_message)
     endif
     if (xyz(1,3,ii) > INVALID_JUDGE) then  ! C
        write (*, '(a, i6, a, a3, a)') &
             'Waring: C atom does not exist in the', ii,'th residue(', c_resname,')'
        call util_error(ERROR%WARN_ALL, error_message)
     endif
     if (xyz(1,4,ii) > INVALID_JUDGE) then  ! O
        write (*, '(a, i6, a, a3, a)') &
             'Waring: O atom does not exist in the', ii,'th residue(', c_resname,')'
        call util_error(ERROR%WARN_ALL, error_message)
     endif

     !!!! GLY
     if (c_resname == 'GLY') then

        if (iatomnum(ii) > 5) then

           if (iatomnum(ii) == 6) then
              do ii2 = 1, iunit
                 if(ii == lunit2mp(2, ii2)) then ! OXT
                    cycle ii_loop
                 end if
              end do
           end if

           write (error_message, '(2(a, i6), a, a3, a)') &
                'Error: too many heavy atoms', iatomnum(ii), &
                ' in the ', ii, 'th residue(', c_resname, ') '
           call util_error(ERROR%STOP_ALL, error_message)
        endif

     !!!! Other residues
     else
        if (xyz(1,5,ii) > INVALID_JUDGE) then  ! CB
           write (*, '(a, i6, a, a3, a)') &
                'Waring: CB atom does not exist in the', ii,'th residue(', c_resname,')'
           call util_error(ERROR%WARN_ALL, error_message)
        endif
     
        nha = ifunc_seq2hanum(c_resname)
        if (nha < 0) then
           write(error_message,*) &
           'Error:',c_resname,'is undefined at residue',ii-nmp,&
           ' of unit',iunit,', in read_pdb_pro.'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        if(iatomnum(ii) > nha) then
           if(iatomnum(ii) == nha + 1) then
              do ii2 = 1, iunit
                 if(ii == lunit2mp(2, ii2)) then ! OXT
                    cycle ii_loop
                 end if
              end do
           end if

           write (error_message, '(2(a, i6), (a, a3), (a, i6), a)') &
                'Error: too many heavy atoms', iatomnum(ii), &
                ' in the ', ii, 'th residue(', c_resname, ') ', nha, &
                ' in read_pdb_pro'
           call util_error(ERROR%STOP_ALL, error_message)

        else if(iatomnum(ii) < nha) then
           write (*, '(2(a, i6), (a, a3), (a, i6), a)') &
                'Waring: less heavy atoms', iatomnum(ii), &
                ' in the ', ii, 'th residue(', c_resname, ') ', nha, &
                ' in read_pdb_pro'
           call util_error(ERROR%WARN_ALL, error_message)
        end if
     endif
  end do ii_loop

  nunit = iunit
  nmp   = imp
  nres  = ires

#ifdef _DEBUG
  write(*,*) '#### end read_pdb_pro'
#endif

end subroutine read_pdb_pro
