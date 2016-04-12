!setp_native_dih 
!> @brief Constructs "native" information about dihedral angles.  &
!>        Related coefficients are also calculated.

subroutine setp_native_dih(xyz_mp_init, xyz_dna)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inpro, inmisc, indna, inrna, inligand, insimu
  use var_struct, only : nunit_all, lunit2mp, iclass_unit, &
                         cmp2seq, ndih, coef_dih, factor_dih, &
                         coef_dih_gauss, wid_dih_gauss, &
                         idih2mp, dih_nat, dih_sin_nat, dih_cos_nat, &
                         imp2type, idih2type, iunit2dih, ires_mp, &
                         nrna_st, irna_st2mp, irna_st2unit, &
                         coef_rna_st, coef_rna_st_a, coef_rna_st_fD, &
                         factor_rna_st, rna_st_nat, rna_st_nat2, &
                         rna_base_type, fdih_para, &
                         imp2unit, &
                         aicg14_nat, factor_aicg14 !aicg2
!#ifdef MPI_PAR
!  use mpiconst
!#endif
  
  implicit none

  ! ----------------------------------------------------------------------
  real(PREC), intent(in) :: xyz_mp_init(SPACE_DIM, MXMP)
  real(PREC), intent(in) :: xyz_dna(SPACE_DIM, MXMP)
  ! intent(out) :: ndih, factor_dih, idih2mp, dih_nat

  ! ----------------------------------------------------------------------
  ! function
!  integer :: ifunc_seq2id
!  real(PREC) :: rfunc_propensity

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: idih, iunit, imp, impmod, imp1, imp2, imp3, imp4, lmp
  integer :: irna_st
  integer :: impmod_P, impmod_S, impmod_B
  integer :: imp_start
  integer :: isumdih(10)
  real(PREC) :: sumdih(10)
!  integer :: id_mp, itype
!  integer :: idel, ini, las
!  integer :: ier
!  integer :: id1, id2
!  integer :: ip
!  real(PREC) :: ave, sum
!  real(PREC) :: pre(4), pre_dih(MXDIH) 
!  character(3) :: char3
  character(CARRAY_MSG_ERROR) :: error_message

  ! ----------------------------------------------------------------------
  ! calc native dih angle
  idih = 0
  irna_st = 0
  isumdih(1:10) = 0
  sumdih(1:10) = 0.0
  do iunit = 1, nunit_all
     
!     if (mod(inmisc%itype_nlocal(iunit, iunit), INTERACT%ENM) == 0) cycle
     if (inmisc%flag_nlocal_unit(iunit, iunit, INTERACT%ENM)) cycle

     lmp = lunit2mp(1, iunit)
     iunit2dih(1, iunit) = idih + 1

     if(iclass_unit(iunit) == CLASS%PRO) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_GO) .or. &
             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG1) .or. &
             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .or. &
             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS) .or. &
             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_FLP)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit) - 3
              imp1 = imp
              imp2 = imp + 1
              imp3 = imp + 2
              imp4 = imp + 3
              call nat_dih(idih, imp1, imp2, imp3, imp4, xyz_mp_init)        
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%LIP) then

     else if(iclass_unit(iunit) == CLASS%DNA) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_BDNA)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit) - 3
!              impmod = mod(imp - lmp, 3)

!              if(impmod == 0) then
              if(imp2type(imp) == MPTYPE%DNA_SUGAR) then
                 ! b-S(3')-P-(5')S A:-22.69 T:-33.42 G:-22.30 C:-32.72
                 imp1 = imp + 1
                 imp2 = imp
                 imp3 = imp + 2
                 imp4 = imp + 3
                 call nat_dih_dna(idih, imp1, imp2, imp3, imp4, xyz_dna)
                 if(cmp2seq(imp1) == ' DA') then
                    sumdih(3) = sumdih(3) + dih_nat(idih)
                    isumdih(3) = isumdih(3) + 1
                 else if(cmp2seq(imp1) == ' DT') then
                    sumdih(4) = sumdih(4) + dih_nat(idih)
                    isumdih(4) = isumdih(4) + 1
                 else if(cmp2seq(imp1) == ' DG') then 
                    sumdih(5) = sumdih(5) + dih_nat(idih)
                    isumdih(5) = isumdih(5) + 1
                 else if(cmp2seq(imp1) == ' DC') then
                    sumdih(6) = sumdih(6) + dih_nat(idih)
                    isumdih(6) = isumdih(6) + 1
                 else
                    error_message = 'Error: not DNA sequence in setp_native_dih'
                    call util_error(ERROR%STOP_ALL, error_message)
                 end if
                 
                 if(imp - 1 < lunit2mp(1, iunit)) cycle
                 ! P-(5')S(3')-P-(5')S -154.80
                 imp1 = imp - 1
                 imp2 = imp
                 imp3 = imp + 2
                 imp4 = imp + 3
                 call nat_dih_dna(idih, imp1, imp2, imp3, imp4, xyz_dna)
                 if(dih_nat(idih) < 0.0) then
                    sumdih(1) = sumdih(1) + dih_nat(idih)
                 else
                    sumdih(1) = sumdih(1) + dih_nat(idih) - 2*F_PI
                 end if
                 isumdih(1) = isumdih(1) + 1
                 
!              else if(impmod == 1) then
              else if(imp2type(imp) == MPTYPE%DNA_BASE) then
                 if(imp - 1 < lunit2mp(1, iunit)) cycle
                 ! S(3')-P-(5')S-b A:50.69 T:54.69 G:50.66 C:54.50
                 imp1 = imp - 1
                 imp2 = imp + 1
                 imp3 = imp + 2
                 imp4 = imp + 3
                 call nat_dih_dna(idih, imp1, imp2, imp3, imp4, xyz_dna)
                 if(cmp2seq(imp4) == ' DA') then
                    sumdih(7) = sumdih(7) + dih_nat(idih)
                    isumdih(7) = isumdih(7) + 1
                 else if(cmp2seq(imp4) == ' DT') then
                    sumdih(8) = sumdih(8) + dih_nat(idih)
                    isumdih(8) = isumdih(8) + 1
                 else if(cmp2seq(imp4) == ' DG') then
                    sumdih(9) = sumdih(9) + dih_nat(idih)
                    isumdih(9) = isumdih(9) + 1
                 else if(cmp2seq(imp4) == ' DC') then
                    sumdih(10) = sumdih(10) + dih_nat(idih)
                    isumdih(10) = isumdih(10) + 1
                 else
                    error_message = 'Error: not DNA sequence in setp_native_dih'
                    call util_error(ERROR%STOP_ALL, error_message)
                 end if
                 
              else
                 if(imp - 2 < lunit2mp(1, iunit)) cycle
                 ! S(3')-P-(5')S(3')-P -179.17
                 imp1 = imp - 2
                 imp2 = imp
                 imp3 = imp + 1
                 imp4 = imp + 3
                 call nat_dih_dna(idih, imp1, imp2, imp3, imp4, xyz_dna)
                 if(dih_nat(idih) < 0.0) then
                    sumdih(2) = sumdih(2) + dih_nat(idih)
                 else
                    sumdih(2) = sumdih(2) + dih_nat(idih) - 2*F_PI
                 end if
                 isumdih(2) = isumdih(2) + 1
              end if
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%DNA2) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_DNA2) .OR. &
           inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_DNA2C)) then
           select case (imp2type(lmp))
           case (MPTYPE%DNA2_PHOS) ! chain starts with P
              impmod_P = 0
              impmod_S = 1
              impmod_B = 2
              imp_start = lmp + 4
           case (MPTYPE%DNA2_SUGAR) ! chain starts with S
              impmod_S = 0
              impmod_B = 1
              impmod_P = 2
              imp_start = lmp + 3
           case (MPTYPE%DNA2_BASE) ! chain must NOT start with B
              error_message = 'Error: DNA sequence is broken in setp_native_dih'
              call util_error(ERROR%STOP_ALL, error_message)
           case default
              error_message = 'Error: not DNA sequence in setp_native_dih'
              call util_error(ERROR%STOP_ALL, error_message)
           end select
           
           do imp = imp_start, lunit2mp(2, iunit)
              impmod = mod(imp - lmp, 3)
              
              if(impmod == impmod_S) then
                 if (imp-4 >= lmp) then
                    ! P-S-P-S
                    imp1 = imp - 4
                    imp2 = imp - 3
                    imp3 = imp - 1
                    imp4 = imp 
                    call nat_dih_dna2(idih, imp1, imp2, imp3, imp4, DIHTYPE%DNA2_PSPS, xyz_mp_init)
                 endif
                 
              else if (impmod == impmod_P) then
                 ! S-P-S-P 
                 imp1 = imp - 5
                 imp2 = imp - 3
                 imp3 = imp - 2
                 imp4 = imp 
                 call nat_dih_dna2(idih, imp1, imp2, imp3, imp4, DIHTYPE%DNA2_SPSP, xyz_mp_init)
                 
              else if (impmod == impmod_B) then
                 cycle
              else 
                 error_message = 'Error: logical defect in setp_native_dih'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end do
        end if

        
     else if(iclass_unit(iunit) == CLASS%RNA) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_GO)) then
           select case (imp2type(lmp))
           case (MPTYPE%RNA_PHOS) ! chain starts with P
              impmod_P = 0
              impmod_S = 1
              impmod_B = 2
              imp_start = lmp + 4
           case (MPTYPE%RNA_SUGAR) ! chain starts with S
              impmod_S = 0
              impmod_B = 1
              impmod_P = 2
              imp_start = lmp + 3
           case (MPTYPE%RNA_BASE) ! chain must NOT start with B
              error_message = 'Error: RNA sequence is broken in setp_native_dih'
              call util_error(ERROR%STOP_ALL, error_message)
           case default
              error_message = 'Error: not RNA sequence in setp_native_dih'
              call util_error(ERROR%STOP_ALL, error_message)
           end select
           
           do imp = imp_start, lunit2mp(2, iunit)
              impmod = mod(imp - lmp, 3)
              
              if(impmod == impmod_S) then
                 if (imp-4 >= lmp) then
                    ! P-(5')S(3')-P-(5')S
                    imp1 = imp - 4
                    imp2 = imp - 3
                    imp3 = imp - 1
                    imp4 = imp 
                    call nat_dih_rna(idih, imp1, imp2, imp3, imp4, DIHTYPE%RNA_PSPS, xyz_mp_init)
                 endif
                 
                 ! b-S(3')-P-(5')S
                 imp1 = imp - 2
                 imp2 = imp - 3
                 imp3 = imp - 1
                 imp4 = imp 
                 if (rna_base_type(imp1) == 'R') then
                    call nat_dih_rna(idih, imp1, imp2, imp3, imp4, DIHTYPE%RNA_RSPS, xyz_mp_init)
                 else if (rna_base_type(imp1) == 'Y') then
                    call nat_dih_rna(idih, imp1, imp2, imp3, imp4, DIHTYPE%RNA_YSPS, xyz_mp_init)
                 else 
                    error_message = 'Error: logical defect in setp_native_bangle, invalid BSPS'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif
                 
              else if(impmod == impmod_B) then
                 if (imp-4 >= lmp) then
                    ! B-S-S-B (stack)
                    imp1 = imp - 3
                    imp2 = imp - 4
                    imp3 = imp - 1
                    imp4 = imp
                    if ((rna_base_type(imp1) == 'R' .AND. rna_base_type(imp4) == 'R') .OR. &
                        (rna_base_type(imp1) == 'R' .AND. rna_base_type(imp4) == 'Y') .OR. &
                        (rna_base_type(imp1) == 'Y' .AND. rna_base_type(imp4) == 'R') .OR. &
                        (rna_base_type(imp1) == 'Y' .AND. rna_base_type(imp4) == 'Y') )then
                       call nat_stack_rna(imp1, imp2, imp3, imp4, xyz_mp_init)
                    else
                       error_message = 'Error: logical defect in setp_native_bangle, invalid BSSB'
                       call util_error(ERROR%STOP_ALL, error_message)
                    endif
                 endif
                 
                 ! S(3')-P-(5')S-b
                 imp1 = imp - 4
                 imp2 = imp - 2
                 imp3 = imp - 1
                 imp4 = imp 
                 if (rna_base_type(imp4) == 'R') then
                    call nat_dih_rna(idih, imp1, imp2, imp3, imp4, DIHTYPE%RNA_SPSR, xyz_mp_init)
                 else if (rna_base_type(imp4) == 'Y') then
                    call nat_dih_rna(idih, imp1, imp2, imp3, imp4, DIHTYPE%RNA_SPSY, xyz_mp_init)
                 else
                    error_message = 'Error: logical defect in setp_native_bangle, invalid SPSB'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif
                 
              else if (impmod == impmod_P) then
                 ! S(3')-P-(5')S(3')-P 
                 imp1 = imp - 5
                 imp2 = imp - 3
                 imp3 = imp - 2
                 imp4 = imp 
                 call nat_dih_rna(idih, imp1, imp2, imp3, imp4, DIHTYPE%RNA_SPSP, xyz_mp_init)
                 
              else 
                 error_message = 'Error: logical defect in setp_native_bangle'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end do
        end if
           
     else if(iclass_unit(iunit) == CLASS%LIG) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_RIGID_LIG)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit) - 3
              imp1 = imp
              imp2 = imp + 1
              imp3 = imp + 2
              imp4 = imp + 3
              if(ires_mp(imp1) == ires_mp(imp4)) then
                 call nat_dih_ligand(idih, imp1, imp2, imp3, imp4, xyz_mp_init)
              end if
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%ION) then

     else 
        write(error_message,*) 'Error: unit',iunit,&
        'has undefined class in setp_native_bangle'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     iunit2dih(2, iunit) = idih

  end do
  ndih = idih
  nrna_st = irna_st
!  if(iclass_unit(1) == CLASS%DNA) then
!     sumdih(1:10) = 180.0/F_PI*sumdih(1:10)/isumdih(1:10)
!     write (*, *) "P-(5')S(3')-P-(5')S(-154.80):", sumdih(1)
!     write (*, *) "S(3')-P-(5')S(3')-S(-179.17):", sumdih(2)
!     write (*, *) "Ab-S(3')-P-(5')(-22.60):", sumdih(3)
!     write (*, *) "Tb-S(3')-P-(5')(-33.42):", sumdih(4)
!     write (*, *) "Gb-S(3')-P-(5')(-22.30):", sumdih(5)
!     write (*, *) "Cb-S(3')-P-(5')(-32.72):", sumdih(6)
!     write (*, *) "S(3')-P-(5')S-Ab(50.69):", sumdih(7)
!     write (*, *) "S(3')-P-(5')S-Tb(54.69):", sumdih(8)
!     write (*, *) "S(3')-P-(5')S-Gb(50.66):", sumdih(9)
!     write (*, *) "S(3')-P-(5')S-Cb(54.50):", sumdih(10)
!  end if

  ! ----------------------------------------------------------------------
  ! bond angle coefficient 
!  if(imodel(2) == 2) then       
!     sum = 0.0e0_PREC
!     do idih = 1, ndih
!        if(iclass_mp(idih2mp(1, idih)) /= CLASS%PRO) cycle
!        do i = 1, 4
!           imp1 = idih2mp(i, idih)       
!           char3 = cmp2seq(imp1)
!           id_mp = ifunc_seq2id(char3)        
!           itype = istype_mp(imp1)         
!           pre(i) = rfunc_propensity(id_mp, itype)
!           
!           pre_dih(idih) = 1.0e0_PREC / (-(log10(pre(1)) + log10(pre(2)) + &
!                log10(pre(3)) + log10(pre(4))))              
!               
!           pre_dih(idih) = (pre(1) * 0.5e0_PREC) * pre(2) * pre(3) * (pre(4) * 0.5e0_PREC)
!           sum = sum + pre_dih(idih)
!        end do
!     end do
!     ave = sum / real(ndih, PREC)
!     do idih = 1, ndih
!        if(iclass_mp(idih2mp(1, idih)) /= CLASS%PRO) cycle
!        factor_dih(idih) = factor_dih(idih) * pre_dih(idih) * (1.0e0_PREC / ave)
!     end do
!  end if

  !-----------------------------------------------------------------------
  ! AICG
!  if (inmisc%force_flag(INTERACT%AICG2) .OR. inmisc%force_flag_local(LINTERACT%L_AICG2)) then
  if (inmisc%force_flag_local(LINTERACT%L_AICG2) .or. &
      inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
     do idih = 1, ndih
        iunit = imp2unit(idih2mp(1, idih))
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .or. &
           inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then
           imp1 = idih2mp(1, idih)
           imp2 = idih2mp(2, idih)
           imp3 = idih2mp(3, idih)
           imp4 = idih2mp(4, idih)
           aicg14_nat(idih) = sqrt((xyz_mp_init(1, imp1) - xyz_mp_init(1, imp4))**2  &
                                 + (xyz_mp_init(2, imp1) - xyz_mp_init(2, imp4))**2  &
                                 + (xyz_mp_init(3, imp1) - xyz_mp_init(3, imp4))**2)
           factor_aicg14(idih) = factor_dih(idih)
           factor_dih(idih) = 0.0e0_PREC
           coef_dih(1, idih) = 0.0e0_PREC
           coef_dih(2, idih) = 0.0e0_PREC
        end if
     end do
   end if


  ! ----------------------------------------------------------------------
  ! del dihedral angle interaction
!  if(inmisc%ndel_ba_dih > 0) then
!     do idel = 1, inmisc%ndel_ba_dih
!        ini = inmisc%idel_ba_dih(1, idel)
!        las = inmisc%idel_ba_dih(2, idel)
!
!        do idih = 1, ndih
!           imp1 = idih2mp(1, idih)
!           imp2 = idih2mp(2, idih)
!           imp3 = idih2mp(3, idih)
!           imp4 = idih2mp(4, idih) 
!           if(imp1 <= las .and. imp4 >= ini) then
!              factor_dih(idih) = 0.0e0_PREC
!              coef_dih(1, idih) = 0.0e0_PREC
!              coef_dih(2, idih) = 0.0e0_PREC
              !nfdih = nfdih + 1
              !ifdih2mp(1, nfdih) = imp1
              !ifdih2mp(2, nfdih) = imp2
              !ifdih2mp(3, nfdih) = imp3
              !ifdih2mp(4, nfdih) = imp4
!           end if
!        end do
!
!     end do
!  end if

  !------------------------------------------------------------------
  ! AICG
!  if (inmisc%force_flag(INTERACT%AICG2) .OR. inmisc%force_flag_local(LINTERACT%L_AICG2)) then
!     do idih = 1, ndih
!        iunit = imp2unit(idih2mp(1, idih))
!        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2))then
!           factor_aicg14(idih) = factor_dih(idih)
!           factor_dih(idih) = 0.0e0_PREC
!           coef_dih(1, idih) = 0.0e0_PREC
!           coef_dih(2, idih) = 0.0e0_PREC
!        end if
!     end do
!  end if

  !if (inmisc%i_add_int == 1) then
     ! allocate flexible dihedral potential parameter array
  !allocate( fdih_para(7, nfdih), stat=ier)
   !  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
    ! fdih_para(:,:) = 0
     
     ! set parameter of flexible local potential
     !do idih = 1, nfdih
!  id1 = ifunc_seq2id( cmp2seq( ifdih2mp(2, idih) ) )
!        id2 = ifunc_seq2id( cmp2seq( ifdih2mp(3, idih) ) )
!        do ip = 1, 7
!           fdih_para(ip, idih) = inflp%dih_para(id1, id2, ip)
!        end do
!     end do
     
     !inflp%coeff = BOLTZC * insimu%tempk

!#ifdef MPI_PAR
!     call MPI_Bcast(inflp, inflp%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
!     call MPI_Bcast(fdih_para, 7*nfdih, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
!#endif
     
!  end if
contains

  subroutine nat_dih(idih, imp1, imp2, imp3, imp4, xyz_dih)

    ! -------------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3, imp4
    integer, intent(inout) :: idih
    real(PREC), intent(in) :: xyz_dih(3, MXMP)

    ! -------------------------------------------------------------------
    ! local variables
    real(PREC) :: dih_angle, si_dih, co_dih

    ! -------------------------------------------------------------------
    idih = idih + 1
    idih2type(idih) = DIHTYPE%PRO
    idih2mp(1, idih) = imp1
    idih2mp(2, idih) = imp2
    idih2mp(3, idih) = imp3
    idih2mp(4, idih) = imp4
    call util_dihangle(imp1, imp2, imp3, imp4, dih_angle, co_dih, si_dih, xyz_dih)
    dih_nat(idih) = dih_angle
    dih_sin_nat(idih) = si_dih
    dih_cos_nat(idih) = co_dih
    factor_dih(idih) = inmisc%factor_local_unit(iunit, iunit)
    coef_dih(1, idih) = factor_dih(idih) * inpro%cdih_1

    if (inmisc%i_triple_angle_term == 0) then
       coef_dih(2, idih) = 0.0e0_PREC
    else if (inmisc%i_triple_angle_term == 1 .or. inmisc%i_triple_angle_term == 2) then
       coef_dih(2, idih) = factor_dih(idih) * inpro%cdih_3
    else 
       error_message = 'Error: invalid value for i_triple_angle_term'
       call util_error(ERROR%STOP_ALL, error_message)
    endif
  end subroutine nat_dih

  subroutine nat_dih_dna(idih, imp1, imp2, imp3, imp4, xyz_dih)

    ! -------------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3, imp4
    integer, intent(inout) :: idih
    real(PREC), intent(in) :: xyz_dih(3, MXMP)

    ! -------------------------------------------------------------------
    ! local variables
    real(PREC) :: dih_angle, si_dih, co_dih

    ! -------------------------------------------------------------------
    idih = idih + 1
    idih2mp(1, idih) = imp1
    idih2mp(2, idih) = imp2
    idih2mp(3, idih) = imp3
    idih2mp(4, idih) = imp4
    call util_dihangle(imp1, imp2, imp3, imp4, dih_angle, co_dih, si_dih, xyz_dih)
    dih_nat(idih) = dih_angle
    dih_sin_nat(idih) = si_dih
    dih_cos_nat(idih) = co_dih
    factor_dih(idih) = inmisc%factor_local_unit(iunit, iunit)
    coef_dih(1, idih) = factor_dih(idih) * indna%cdih_1_dna

    if (inmisc%i_triple_angle_term == 0) then
       coef_dih(2, idih) = 0.0e0_PREC
    else if (inmisc%i_triple_angle_term == 1 .or. inmisc%i_triple_angle_term == 2) then
       coef_dih(2, idih) = factor_dih(idih) * indna%cdih_3_dna
    else
       error_message = 'Error: invalid value for i_triple_angle_term'
       call util_error(ERROR%STOP_ALL, error_message)
    endif
  end subroutine nat_dih_dna

  subroutine nat_dih_dna2(idih, imp1, imp2, imp3, imp4, i_type_dih, xyz_dih)

     use var_setp,   only : indna2
    ! -------------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3, imp4
    integer, intent(inout) :: idih
    integer, intent(in) :: i_type_dih
    real(PREC), intent(in) :: xyz_dih(SPACE_DIM, MXMP)

    ! -------------------------------------------------------------------
    ! local variables
    !real(PREC) :: dih_angle, si_dih, co_dih

    ! -------------------------------------------------------------------
    
    idih = idih + 1
    idih2type(idih) = i_type_dih
    idih2mp(1, idih) = imp1
    idih2mp(2, idih) = imp2
    idih2mp(3, idih) = imp3
    idih2mp(4, idih) = imp4
    !dih_nat(idih) = dih_angle
    !dih_sin_nat(idih) = si_dih
    !dih_cos_nat(idih) = co_dih
    factor_dih(idih) = inmisc%factor_local_unit(iunit, iunit)

    select case (i_type_dih)
    case (DIHTYPE%DNA2_PSPS)
       !coef_dih(1, idih) = factor_dih(idih) * indna2%cdih_PSPS
       !coef_dih(2, idih) = indna2%sdih_PSPS
       coef_dih_gauss(idih) = factor_dih(idih) * indna2%cdih_PSPS
       wid_dih_gauss(idih)  = indna2%sdih_PSPS
       
       dih_nat(idih) = indna2%ndih_PSPS * F_PI / 180.0e0_PREC
    case (DIHTYPE%DNA2_SPSP)
       !coef_dih(1, idih) = factor_dih(idih) * indna2%cdih_SPSP
       !coef_dih(2, idih) = indna2%sdih_SPSP
       coef_dih_gauss(idih) = factor_dih(idih) * indna2%cdih_SPSP
       wid_dih_gauss(idih)  = indna2%sdih_SPSP

       dih_nat(idih) = indna2%ndih_SPSP * F_PI / 180.0e0_PREC
    end select
    
  end subroutine nat_dih_dna2

  subroutine nat_dih_rna(idih, imp1, imp2, imp3, imp4, i_type_dih, xyz_dih)

     use var_setp,   only : inrna
    ! -------------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3, imp4
    integer, intent(inout) :: idih
    integer, intent(in) :: i_type_dih
    real(PREC), intent(in) :: xyz_dih(SPACE_DIM, MXMP)

    ! -------------------------------------------------------------------
    ! local variables
    real(PREC) :: dih_angle, si_dih, co_dih

    ! -------------------------------------------------------------------
     
    call util_dihangle(imp1, imp2, imp3, imp4, dih_angle, co_dih, si_dih, xyz_dih)

    idih = idih + 1
    idih2type(idih) = i_type_dih
    idih2mp(1, idih) = imp1
    idih2mp(2, idih) = imp2
    idih2mp(3, idih) = imp3
    idih2mp(4, idih) = imp4
    dih_nat(idih) = dih_angle
    dih_sin_nat(idih) = si_dih
    dih_cos_nat(idih) = co_dih
    factor_dih(idih) = inmisc%factor_local_unit(iunit, iunit)

    select case (i_type_dih)
    case (DIHTYPE%RNA_RSPS)
       coef_dih(1, idih) = factor_dih(idih) * inrna%cdih_1_RSPS
       coef_dih(2, idih) = factor_dih(idih) * inrna%cdih_3_RSPS
    case (DIHTYPE%RNA_YSPS)
       coef_dih(1, idih) = factor_dih(idih) * inrna%cdih_1_YSPS
       coef_dih(2, idih) = factor_dih(idih) * inrna%cdih_3_YSPS
    case (DIHTYPE%RNA_PSPS)
       coef_dih(1, idih) = factor_dih(idih) * inrna%cdih_1_PSPS
       coef_dih(2, idih) = factor_dih(idih) * inrna%cdih_3_PSPS
    case (DIHTYPE%RNA_SPSR)
       coef_dih(1, idih) = factor_dih(idih) * inrna%cdih_1_SPSR
       coef_dih(2, idih) = factor_dih(idih) * inrna%cdih_3_SPSR
    case (DIHTYPE%RNA_SPSY)
       coef_dih(1, idih) = factor_dih(idih) * inrna%cdih_1_SPSY
       coef_dih(2, idih) = factor_dih(idih) * inrna%cdih_3_SPSY
    case (DIHTYPE%RNA_SPSP)
       coef_dih(1, idih) = factor_dih(idih) * inrna%cdih_1_SPSP
       coef_dih(2, idih) = factor_dih(idih) * inrna%cdih_3_SPSP
    end select 

    if (inmisc%i_triple_angle_term == 0) then
       coef_dih(2, idih) = 0.0e0_PREC
    else if (inmisc%i_triple_angle_term /= 1) then
       error_message = 'Error: invalid value for i_triple_angle_term'
       call util_error(ERROR%STOP_ALL, error_message)
    endif
  end subroutine nat_dih_rna

  subroutine nat_stack_rna(imp1, imp2, imp3, imp4, xyz_dih)
     use var_setp,   only : inrna
     integer, intent(in) :: imp1, imp2, imp3, imp4
     real(PREC), intent(in) :: xyz_dih(SPACE_DIM, MXMP)
 
     real(PREC) :: dih_angle, si_dih, co_dih
     real(PREC) :: dist2
      
     call util_dihangle(imp1, imp2, imp3, imp4, dih_angle, co_dih, si_dih, xyz_dih)
 
     dist2 = (xyz_dih(1,imp1) - xyz_dih(1,imp4)) ** 2 &
            +(xyz_dih(2,imp1) - xyz_dih(2,imp4)) ** 2 &
            +(xyz_dih(3,imp1) - xyz_dih(3,imp4)) ** 2
     if (      dih_angle >= inrna%dfhelix_BSSB_lower &
         .and. dih_angle <= inrna%dfhelix_BSSB_upper &
         .and. dist2 <= inrna%dfcontact_st**2        ) then
        irna_st = irna_st + 1
        if (irna_st > MXRNAST) then
           error_message = 'Error: number of base stack is larger than MXRNAST'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        irna_st2mp(1, irna_st) = imp1
        irna_st2mp(2, irna_st) = imp4
        irna_st2unit(1, irna_st) = iunit
        irna_st2unit(2, irna_st) = iunit
        rna_st_nat2(irna_st) = dist2
        rna_st_nat(irna_st)  = sqrt(dist2)
        factor_rna_st(irna_st) = inmisc%factor_local_unit(iunit, iunit)
        coef_rna_st(irna_st)   = factor_rna_st(irna_st) * inrna%cst1210
        coef_rna_st_fD(irna_st)= factor_rna_st(irna_st) * inrna%cstmorse_D
        coef_rna_st_a(irna_st) = factor_rna_st(irna_st) * inrna%cstmorse_a
     endif
  end subroutine nat_stack_rna

  subroutine nat_dih_ligand(idih, imp1, imp2, imp3, imp4, xyz_dih)

    ! -------------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3, imp4
    integer, intent(inout) :: idih
    real(PREC), intent(in) :: xyz_dih(3, MXMP)

    ! -------------------------------------------------------------------
    ! local variables
    real(PREC) :: dih_angle, si_dih, co_dih

    ! -------------------------------------------------------------------
    idih = idih + 1
    idih2mp(1, idih) = imp1
    idih2mp(2, idih) = imp2
    idih2mp(3, idih) = imp3
    idih2mp(4, idih) = imp4
    call util_dihangle(imp1, imp2, imp3, imp4, dih_angle, co_dih, si_dih, xyz_dih)
    dih_nat(idih) = dih_angle
    dih_sin_nat(idih) = si_dih
    dih_cos_nat(idih) = co_dih
    factor_dih(idih) = inmisc%factor_local_unit(iunit, iunit)
    coef_dih(1, idih) = factor_dih(idih) * inligand%cdih
    coef_dih(2, idih) = 0.0e0_PREC
  end subroutine nat_dih_ligand

end subroutine setp_native_dih
