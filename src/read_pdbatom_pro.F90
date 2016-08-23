!read_pdbatom_pro
!> @brief Reads structure information about proteins from pdb_atom.

subroutine read_pdbatom_pro(pdb_atom,  & ![i ] pdb atom information
                        lunit2atom,    & ![i ] correspondence list (unit -> atom)
                        nunit_atom, & ![i ] initial and last unit
                        nmp,      & ![io] the number of mp (mass point)
                        nres,     & ![io] the number of residue
                        lunit2mp, & ![ o] correspondence list (unit -> mp  )
                        ires_mp,  & ![ o] residue index       (mp   -> res )
                        cmp2seq,  & ![ o] correspondence list (mp   -> seq )
                        imp2type, & ![ o] correspondence list (mp   -> type)
                        iatomnum, & ![ o] the number of atom  (mp   -> #   )
                        xyz,      & ![ o] coordinate of all atom 
                        cname_ha)   ![ o] name of heavy atoms  ! aicg

  use const_maxsize
  use const_index
  use const_physical
  use var_io, only : pdbatom
  implicit none

  ! ---------------------------------------------------------------------
  type(pdbatom), intent(in)    :: pdb_atom(:)
  integer,       intent(in)    :: lunit2atom(2, MXUNIT)
  integer,       intent(in)    :: nunit_atom(2)
  integer,       intent(inout) :: nmp, nres
  integer,       intent(out)   :: lunit2mp(2, MXUNIT)
  integer,       intent(out)   :: ires_mp(MXMP)
  character(3),  intent(out)   :: cmp2seq(MXMP)
  integer,       intent(out)   :: imp2type(MXMP)
  integer,       intent(out)   :: iatomnum(MXMP)
  real(PREC),    intent(out)   :: xyz(3, MXATOM_MP, MXMP)
  character(4),  intent(out)   :: cname_ha(MXATOM_MP, MXMP)  ! aicg

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: iatom, iatomn
  integer :: imp, iunit, ires, ii, ii2, nha
  logical :: flg_reading         ! flag for 
  character(CARRAY_MSG_ERROR) :: error_message

  integer      :: i_resseq_save ! to store previous i_resseq
  character(1) :: c_icode_save  ! to store previous c_icode_save
  character(3) :: c_resnamel

  integer :: ifunc_seq2hanum
#ifdef _DEBUG
  write(*,*) '#### start read_pdbatom_pro'
#endif

  ! ---------------------------------------------------------------------
  c_icode_save = ' '
  imp = nmp
  ires = nres

  ! ---------------------------------------------------------------------
  do iunit = nunit_atom(1), nunit_atom(2)
     lunit2mp(1, iunit) = imp + 1
     flg_reading   = .false.

     do iatom = lunit2atom(1, iunit), lunit2atom(2, iunit)

        if(pdb_atom(iatom)%c_altloc /= " " .and. &
             pdb_atom(iatom)%c_altloc /= "A") cycle

        ! if using H atom, if this statement should be coment out
        if(pdb_atom(iatom)%c_name(1:1) == 'H' .or. &
             pdb_atom(iatom)%c_name(2:2) == 'H') cycle
        
        ! New residue
        if((.not. flg_reading) .or. &
           pdb_atom(iatom)%i_resseq /= i_resseq_save .or. &
           pdb_atom(iatom)%c_icode /= c_icode_save) then

           imp = imp + 1
           if (imp > MXMP) then
              write(error_message,*) 'Error: imp > MXMP, in read_pdbatom_pro. imp=',imp,' MXMP=',MXMP
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           ires = ires + 1
           ires_mp(imp) = ires
           iatomnum(imp) = 5
           cmp2seq(imp) = pdb_atom(iatom)%c_resname
           imp2type(imp) = MPTYPE%PRO

           i_resseq_save = pdb_atom(iatom)%i_resseq
           c_icode_save = pdb_atom(iatom)%c_icode
           flg_reading = .true.
           xyz(1:3, :, imp) = INVALID_VALUE

        end if

        if(pdb_atom(iatom)%c_name == ' N  ') then
           iatomn = 1     
        else if(pdb_atom(iatom)%c_name == ' CA ') then
           iatomn = 2
        else if(pdb_atom(iatom)%c_name == ' C  ') then
           iatomn = 3
        else if(pdb_atom(iatom)%c_name == ' O  ') then
           iatomn = 4
        else if(pdb_atom(iatom)%c_name == ' CB ') then
           iatomn = 5
        else
           iatomnum(imp) = iatomnum(imp) + 1
           if (iatomnum(imp) > MXATOM_MP) then
              write(error_message,*) 'Error: too many atom in ', pdb_atom(iatom)%c_resname,' of unit ', iunit
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           iatomn = iatomnum(imp)
        end if

        xyz(1, iatomn, imp) = pdb_atom(iatom)%x
        xyz(2, iatomn, imp) = pdb_atom(iatom)%y
        xyz(3, iatomn, imp) = pdb_atom(iatom)%z

        cname_ha(iatomn, imp) = pdb_atom(iatom)%c_name  ! aicg
     end do

     lunit2mp(2, iunit) = imp
  end do


  ! ---------------------------------------------------------------------
  ii_loop: do ii = nmp + 1, imp
     c_resnamel = cmp2seq(ii)

     if (xyz(1,1,ii) > INVALID_JUDGE) then  ! N
        write (*, '(a, i6, a, a3, a)') &
             'Waring: N atom does not exist in the', ii,'th residue(', c_resnamel,')'
        call util_error(ERROR%WARN_ALL, error_message)
     endif
     if (xyz(1,2,ii) > INVALID_JUDGE) then  ! CA
        write (*, '(a, i6, a, a3, a)') &
             'Waring: CA atom does not exist in the', ii,'th residue(', c_resnamel,')'
        call util_error(ERROR%WARN_ALL, error_message)
     endif
     if (xyz(1,3,ii) > INVALID_JUDGE) then  ! C
        write (*, '(a, i6, a, a3, a)') &
             'Waring: C atom does not exist in the', ii,'th residue(', c_resnamel,')'
        call util_error(ERROR%WARN_ALL, error_message)
     endif
     if (xyz(1,4,ii) > INVALID_JUDGE) then  ! O
        write (*, '(a, i6, a, a3, a)') &
             'Waring: O atom does not exist in the', ii,'th residue(', c_resnamel,')'
        call util_error(ERROR%WARN_ALL, error_message)
     endif

     !!!! GLY
     if (c_resnamel == 'GLY') then

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
                ' in the ', ii, 'th residue(', c_resnamel, ') '
           call util_error(ERROR%STOP_ALL, error_message)
        endif

     !!!! Other residues
     else
        if (xyz(1,5,ii) > INVALID_JUDGE) then  ! CB
           write (*, '(a, i6, a, a3, a)') &
                'Waring: CB atom does not exist in the', ii,'th residue(', c_resnamel,')'
           call util_error(ERROR%WARN_ALL, error_message)
        endif
     
        nha = ifunc_seq2hanum(c_resnamel)
        if (nha < 0) then
           write(error_message,*) &
           'Error:',c_resnamel,'is undefined at residue',ii-nmp,&
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
                ' in the ', ii, 'th residue(', c_resnamel, ') ', nha, &
                ' in read_pdb_pro'
           call util_error(ERROR%STOP_ALL, error_message)

        else if(iatomnum(ii) < nha) then
           write (*, '(2(a, i6), (a, a3), (a, i6), a)') &
                'Waring: less heavy atoms', iatomnum(ii), &
                ' in the ', ii, 'th residue(', c_resnamel, ') ', nha, &
                ' in read_pdb_pro'
           call util_error(ERROR%WARN_ALL, error_message)
        end if
     endif
  end do ii_loop

  nmp   = imp
  nres  = ires

#ifdef _DEBUG
  write(*,*) '#### end read_pdbatom_pro'
#endif

end subroutine read_pdbatom_pro
