! setp_native_bangle
!> @brief This subroutine is to calculate the native bond angle.

! ************************************************************************
subroutine setp_native_bangle(xyz_mp_init, xyz_dna) 

  use const_maxsize
  use const_physical
  use const_index
  use var_setp, only : inpro, inmisc, indna, indna2, inlip, &
                       inrna, indtrna13, indtrna15, inarna, &
                       inligand, insimu
  use var_struct, only : nunit_all, lunit2mp, iclass_unit, &
                         cmp2seq, nba, iba2mp, ba_nat, factor_ba, coef_ba, &
                         imp2type, iba2type, iunit2ba, ires_mp, cmp2atom, &
                         imp2unit, &
                         aicg13_nat, factor_aicg13   !aicg2
                         !fba_para_x, fba_para_y, fba_para_y2

!#ifdef MPI_PAR
!  use mpiconst
!#endif
  
  implicit none    

  ! -----------------------------------------------------------------
  real(PREC), intent(in) :: xyz_mp_init(SPACE_DIM, MXMP)
  real(PREC), intent(in) :: xyz_dna(SPACE_DIM, MXMP)
  ! intent(out) :: nba, iba2mp, ba_nat, factor_ba

  ! -----------------------------------------------------------------
  ! function
!  integer :: ifunc_seq2id
!  real(PREC) :: rfunc_propensity

  ! -----------------------------------------------------------------
  ! local variables
  integer :: imp, impmod, imp1, imp2, imp3, iba, lmp
  integer :: impmod_P, impmod_S, impmod_B
  integer :: iunit
!  integer :: idel, ini, las
!  integer :: id, ip
! integer :: itype, id_mp, ier
!  integer :: isumba(10)
!  real(PREC) :: sumba(10)
!  real(PREC) :: ave, sum
!  real(PREC) :: pre(3), pre_ba(MXBA)
!  character :: char3*3
  character(CARRAY_MSG_ERROR) :: error_message
  character(4) :: cmp

  ! -----------------------------------------------------------------
  ! calc native bond angle
  iba = 0
!  isumba(1:10) = 0
!  sumba(1:10) = 0.0
  do iunit = 1, nunit_all

!     if (mod(inmisc%itype_nlocal(iunit, iunit), INTERACT%ENM) == 0) cycle
     if (inmisc%flag_nlocal_unit(iunit, iunit, INTERACT%ENM)) cycle

     lmp = lunit2mp(1, iunit)
     iunit2ba(1, iunit) = iba + 1
     
     if(iclass_unit(iunit) == CLASS%PRO) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_GO) .or. &
             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG1) .or. &
             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .or. &
             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS) .or. &
             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_FLP)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit) - 2
              imp1 = imp
              imp2 = imp + 1
              imp3 = imp + 2
              call nat_bangle(iba, imp1, imp2, imp3, xyz_mp_init)
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%LIP) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_LIP_BROWN) .or. inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_LIP_NOGU)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit) - 2
              impmod = mod(imp - lmp, inlip%num_lip_total)
              if(impmod < inlip%num_lip_total - 2) then
                 imp1 = imp
                 imp2 = imp + 1
                 imp3 = imp + 2
                 impmod = impmod + 1
!              call nat_bangle_lipid(iba, imp1, imp2, imp3, xyz_mp_rep(:,:,1))
                 call nat_bangle_lipid(iba, imp1, imp2, imp3)
              end if
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%DNA) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_BDNA)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit) - 2
!              impmod = mod(imp - lmp, 3)
!              if(impmod == 0) then
              if(imp2type(imp) == MPTYPE%DNA_SUGAR) then
                 ! b-S(3')-P A:108.38 T:112.71 G:108.12 C:112.39
                 imp1 = imp + 1 
                 imp2 = imp
                 imp3 = imp + 2
                 call nat_bangle_dna(iba, imp1, imp2, imp3, xyz_dna)
                 if(cmp2seq(imp1) == ' DA') then
!                 sumba(3) = sumba(3) + ba_nat(iba)
!                 isumba(3) = isumba(3) + 1
                 else if(cmp2seq(imp1) == ' DT') then
!                 sumba(4) = sumba(4) + ba_nat(iba)
!                 isumba(4) = isumba(4) + 1
                 else if(cmp2seq(imp1) == ' DG') then
!                 sumba(5) = sumba(5) + ba_nat(iba)
!                 isumba(5) = isumba(5) + 1
                 else if(cmp2seq(imp1) == ' DC') then
!                 sumba(6) = sumba(6) + ba_nat(iba)
!                 isumba(6) = isumba(6) + 1
                 else
                    error_message = 'Error: not DNA sequence in setp_native_bangle'
                    call util_error(ERROR%STOP_ALL, error_message)
                 end if
                 
                 if(imp - 1 < lunit2mp(1, iunit)) cycle
                 ! P-(5')S(3')-P 120.15
                 imp1 = imp - 1 
                 imp2 = imp
                 imp3 = imp + 2
                 call nat_bangle_dna(iba, imp1, imp2, imp3, xyz_dna)
!              sumba(1) = sumba(1) + ba_nat(iba)
!              isumba(1) = isumba(1) + 1
                 
!              else if(impmod == 1) then
              else if(imp2type(imp) == MPTYPE%DNA_BASE) then
                 if(imp - 1 < lunit2mp(1, iunit)) cycle
                 ! S(3')-P-(5')S 94.49
                 imp1 = imp - 1
                 imp2 = imp + 1
                 imp3 = imp + 2
                 call nat_bangle_dna(iba, imp1, imp2, imp3, xyz_dna)
!              sumba(2) = sumba(2) + ba_nat(iba)
!              isumba(2) = isumba(2) + 1

              else
                 ! P-(5')S-b A:113.13 T:102.79 G:113.52 C:103.49
                 imp1 = imp
                 imp2 = imp + 1
                 imp3 = imp + 2
                 call nat_bangle_dna(iba, imp1, imp2, imp3, xyz_dna)
                 if(cmp2seq(imp3) == ' DA') then
!                 sumba(7) = sumba(7) + ba_nat(iba)
!                 isumba(7) = isumba(7) + 1
                 else if(cmp2seq(imp3) == ' DT') then
!                 sumba(8) = sumba(8) + ba_nat(iba)
!                 isumba(8) = isumba(8) + 1
                 else if(cmp2seq(imp3) == ' DG') then
!                 sumba(9) = sumba(9) + ba_nat(iba)
!                 isumba(9) = isumba(9) + 1
                 else if(cmp2seq(imp3) == ' DC') then
!                 sumba(10) = sumba(10) + ba_nat(iba)
!                 isumba(10) = isumba(10) + 1
                 else
                    error_message = 'Error: not DNA sequence in setp_native_bangle'
                    call util_error(ERROR%STOP_ALL, error_message)
                 end if
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
           case (MPTYPE%DNA2_SUGAR) ! chain starts with S
              impmod_S = 0
              impmod_B = 1
              impmod_P = 2
           case (MPTYPE%DNA2_BASE) ! chain must NOT start with B
              error_message = 'Error: DNA sequence is broken in setp_native_bangle'
              call util_error(ERROR%STOP_ALL, error_message)
           case default
              error_message = 'Error: not DNA sequence in setp_native_bangle'
              call util_error(ERROR%STOP_ALL, error_message)
           end select

           do imp = lmp + 2, lunit2mp(2, iunit)
              impmod = mod(imp - lmp, 3)
              
              if(impmod == impmod_S) then
                 ! S-P-S
                 imp1 = imp - 3
                 imp2 = imp - 1
                 imp3 = imp
                 call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_SPS, xyz_mp_init)
                 
              else if(impmod == impmod_B) then

                 ! P-S-B
                 imp1 = imp - 2
                 imp2 = imp - 1
                 imp3 = imp
                 if (cmp2seq(imp3) == 'DA ') then
                    call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_PSA, xyz_mp_init)
                 else if (cmp2seq(imp3) == 'DT ') then
                    call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_PST, xyz_mp_init)
                 else if (cmp2seq(imp3) == 'DC ') then
                    call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_PSC, xyz_mp_init)
                 else if (cmp2seq(imp3) == 'DG ') then
                    call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_PSG, xyz_mp_init)
                 else
                    error_message = 'Error: (1) logical defect in setp_native_bangle, invalid cmp2seq(imp3)'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif
                 
              else if (impmod == impmod_P) then
                 ! B-S-P
                 imp1 = imp - 1 
                 imp2 = imp - 2
                 imp3 = imp 
                 if (cmp2seq(imp1) == 'DA ') then
                    call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_ASP, xyz_mp_init)
                 else if (cmp2seq(imp1) == 'DT ') then
                    call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_TSP, xyz_mp_init)
                 else if (cmp2seq(imp1) == 'DC ') then
                    call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_CSP, xyz_mp_init)
                 else if (cmp2seq(imp1) == 'DG ') then
                    call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_GSP, xyz_mp_init)
                 else
                    error_message = 'Error: (2) logical defect in setp_native_bangle, invalid cmp2seq(imp3)'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif

                 
                 if (imp < lunit2mp(1, iunit) + 3) cycle
                 ! P-S-P
                 imp1 = imp - 3
                 imp2 = imp - 2
                 imp3 = imp 
                 call nat_bangle_dna2(iba, imp1, imp2, imp3, BATYPE%DNA2_PSP ,xyz_mp_init)
                 
              else 
                 error_message = 'Error: (3) logical defect in setp_native_bangle'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end do
        end if


     else if(iclass_unit(iunit) == CLASS%RNA) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_GO) .OR. &
           inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_DTRNA)) then
           select case (imp2type(lmp))
           case (MPTYPE%RNA_PHOS) ! chain starts with P
              impmod_P = 0
              impmod_S = 1
              impmod_B = 2
           case (MPTYPE%RNA_SUGAR) ! chain starts with S
              impmod_S = 0
              impmod_B = 1
              impmod_P = 2
           case (MPTYPE%RNA_BASE) ! chain must NOT start with B
              error_message = 'Error: RNA sequence is broken in setp_native_bangle'
              call util_error(ERROR%STOP_ALL, error_message)
           case default
              error_message = 'Error: not RNA sequence in setp_native_bangle'
              call util_error(ERROR%STOP_ALL, error_message)
           end select

           do imp = lmp + 2, lunit2mp(2, iunit)
              impmod = mod(imp - lmp, 3)
              
              if(impmod == impmod_S) then
                 ! S(3')-P-(5')S
                 imp1 = imp - 3
                 imp2 = imp - 1
                 imp3 = imp
                 call nat_bangle_rna(iba, iunit, imp1, imp2, imp3, cmp, BATYPE%RNA_SPS, xyz_mp_init)
                 
              else if(impmod == impmod_B) then
!              if (imp-3 >= lmp) then
!                 ! B-S-B (stack)
!                 imp1 = imp - 3
!                 imp2 = imp - 1
!                 imp3 = imp
!                 call nat_bangle_rna(iba, imp1, imp2, imp3, BATYPE%RNA_BSB, xyz_mp_init)
!              endif

                 ! P-(5')S-b
                 imp1 = imp - 2
                 imp2 = imp - 1
                 imp3 = imp
                 cmp = cmp2atom(imp3)
                 if (cmp == ' Ab ' .OR. cmp == ' Gb ' .OR. cmp == ' Rb ') then
                    call nat_bangle_rna(iba, iunit, imp1, imp2, imp3, cmp, BATYPE%RNA_PSR, xyz_mp_init)
                 else if (cmp == ' Ub ' .OR. cmp == ' Cb ' .OR. cmp == ' Yb ') then
                    call nat_bangle_rna(iba, iunit, imp1, imp2, imp3, cmp, BATYPE%RNA_PSY, xyz_mp_init)
                 else
                    error_message = 'Error: logical defect in setp_native_bangle, invalid cmp2atom(imp3)'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif
                 
              else if (impmod == impmod_P) then
                 ! b-S(3')-P
                 imp1 = imp - 1 
                 imp2 = imp - 2
                 imp3 = imp 
                 cmp = cmp2atom(imp1)
                 if (cmp == ' Ab ' .OR. cmp == ' Gb ' .OR. cmp == ' Rb ') then
                    call nat_bangle_rna(iba, iunit, imp1, imp2, imp3, cmp, BATYPE%RNA_RSP, xyz_mp_init)
                 else if (cmp == ' Ub ' .OR. cmp == ' Cb ' .OR. cmp == ' Yb ') then
                    call nat_bangle_rna(iba, iunit, imp1, imp2, imp3, cmp, BATYPE%RNA_YSP, xyz_mp_init)
                 else
                    error_message = 'Error: logical defect in setp_native_bangle, invalid cmp2atom(imp1)'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif
                 
                 if (imp < lunit2mp(1, iunit) + 3) cycle
                 ! P-(5')S(3')-P
                 imp1 = imp - 3
                 imp2 = imp - 2
                 imp3 = imp 
                 call nat_bangle_rna(iba, iunit, imp1, imp2, imp3, cmp, BATYPE%RNA_PSP ,xyz_mp_init)
                 
              else 
                 error_message = 'Error: logical defect in setp_native_bangle'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%LIG) then
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_RIGID_LIG)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit) - 2
              imp1 = imp
              imp2 = imp + 1
              imp3 = imp + 2
              if(ires_mp(imp1) == ires_mp(imp3)) then
                 call nat_bangle_ligand(iba, imp1, imp2, imp3, xyz_mp_init)
              end if
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%ION) then

     else
        write(error_message,*) 'Error: unit',iunit,&
             'has undefined class in setp_native_bangle'
        call util_error(ERROR%STOP_ALL, error_message)

     end if

     iunit2ba(2, iunit) = iba

  end do
  nba = iba
!  if(iclass_unit(1) == CLASS%DNA) then
!     sumba(1:10) = 180.0/3.1416*sumba(1:10)/isumba(1:10)
!     write (*, *) "P-(5')S(3')-P(120.15):", sumba(1)
!     write (*, *) "S(3')-P-(5')S(94.49):", sumba(2)
!     write (*, *) "Ab-S(3')-P(108.38):", sumba(3)
!     write (*, *) "Tb-S(3')-P(112.72):", sumba(4)
!     write (*, *) "Gb-S(3')-P(108.12):", sumba(5)
!     write (*, *) "Cb-S(3')-P(112.39):", sumba(6)
!     write (*, *) "P-(5')S-Ab(113.13):", sumba(7)
!     write (*, *) "P-(5')S-Tb(102.79):", sumba(8)
!     write (*, *) "P-(5')S-Gb(113.52):", sumba(9)
!     write (*, *) "P-(5')S-Cb(103.49):", sumba(10)
!  end if

  ! -----------------------------------------------------------------
  ! bond angle coefficient
!  if(imodel(2) == 2) then 
!     sum = 0.0e0_PREC
!     do iba = 1, nba
!        if(iclass_mp(iba2mp(1, iba)) /= CLASS%PRO) cycle
!        do i = 1, 3
!           imp1 = iba2mp(i, iba)
!           char3 = cmp2seq(imp1)
!           id_mp = ifunc_seq2id(char3)        
!           itype = istype_mp(imp1)         
!           pre(i) = rfunc_propensity(id_mp, itype)
!        end do
!            
!        pre_ba(iba) = (pre(1) * 0.5e0_PREC) * pre(2) * (pre(3) * 0.5e0_PREC)
!        pre_ba(iba) = 1.0e0_PREC / &
!             (-(log10(pre(1)) + log10(pre(2)) + log10(pre(3))))              
!            
!        pre_ba(iba) = 1.0e0_PREC / &
!             (-(log10(pre(1)) + log10(pre(2)) + log10(pre(3))))              
!
!        sum = sum + pre_ba(iba)
!     end do
!
!     ave = sum / real(nba, PREC)
!     do iba = 1, nba
!        if(iclass_mp(iba2mp(1, iba)) /= CLASS%PRO) cycle
!        factor_ba(iba) = factor_ba(iba) * pre_ba(iba) * (1 / ave)
!     end do
!  end if

  !------------------------------------------------------------------
  ! aicg2
!  if (inmisc%force_flag(INTERACT%AICG2) .OR. inmisc%force_flag_local(LINTERACT%L_AICG2)) then
  if (inmisc%force_flag_local(LINTERACT%L_AICG2) .or. &
      inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
     do iba = 1, nba
        iunit = imp2unit(iba2mp(1, iba))
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .or. &
           inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then
           imp1 = iba2mp(1, iba)
           imp2 = iba2mp(2, iba)
           imp3 = iba2mp(3, iba)
           aicg13_nat(iba) = sqrt((xyz_mp_init(1, imp1) - xyz_mp_init(1, imp3))**2  &
                                + (xyz_mp_init(2, imp1) - xyz_mp_init(2, imp3))**2  &
                                + (xyz_mp_init(3, imp1) - xyz_mp_init(3, imp3))**2)
           factor_aicg13(iba) = factor_ba(iba)
           factor_ba(iba) = 0.0e0_PREC
           coef_ba(1, iba) = 0.0e0_PREC
           coef_ba(2, iba) = 0.0e0_PREC
        end if
     end do
  end if

  ! -----------------------------------------------------------------
  ! del bond angel interaction
!  if(inmisc%ndel_ba_dih > 0) then
!     do idel = 1, inmisc%ndel_ba_dih
!        ini = inmisc%idel_ba_dih(1, idel)
!        las = inmisc%idel_ba_dih(2, idel)
!
!        do iba = 1, nba
!           imp1 = iba2mp(1, iba)
!           imp2 = iba2mp(2, iba)
!           imp3 = iba2mp(3, iba)
!           if(imp1 <= las .and. imp3 >= ini) then
!              factor_ba(iba) = 0.0e0_PREC
!              coef_ba(1, iba) = 0.0e0_PREC
!              coef_ba(2, iba) = 0.0e0_PREC
              !nfba = nfba + 1
              !ifba2mp(1, nfba) = imp1
              !ifba2mp(2, nfba) = imp2
              !ifba2mp(3, nfba) = imp3
!
!           end if
!        end do
!
!     end do
!  end if

  !------------------------------------------------------------------
  ! AICG
!  if (inmisc%force_flag(INTERACT%AICG2) .OR. inmisc%force_flag_local(LINTERACT%L_AICG2)) then
!     do iba = 1, nba
!        iunit = imp2unit(iba2mp(1, iba))
!        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2))then
!           factor_aicg13(iba) = factor_ba(iba)
!           factor_ba(iba) = 0.0e0_PREC
!           coef_ba(1, iba) = 0.0e0_PREC
!           coef_ba(2, iba) = 0.0e0_PREC
!        end if
!     end do
!  end if

  ! allocate flexible bond angle potential parameter array
  !allocate( fba_para_x(10, nfba), stat=ier)
  !if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  !fba_para_x(:,:) = 0
  !allocate( fba_para_y(10, nfba), stat=ier)
  !if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  !fba_para_y(:,:) = 0
  !allocate( fba_para_y2(10, nfba), stat=ier)
  !if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  !fba_para_y2(:,:) = 0

  !if (inmisc%i_add_int == 1) then
     ! set parameter of flexible local potential
  !   do iba = 1, nfba
  !      id = ifunc_seq2id( cmp2seq( ifba2mp(2, iba) ) )
  !      do ip = 1, 10
  !         fba_para_x(ip, iba)  = inflp%ang_para_x(ip)
  !         fba_para_y(ip, iba)  = inflp%ang_para_y(id, ip)
  !         fba_para_y2(ip, iba) = inflp%ang_para_y2(id, ip)
  !      end do
  !   end do
     
     !inflp%coeff = BOLTZC * insimu%tempk

!#ifdef MPI_PAR
!     call MPI_Bcast(inflp, inflp%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
!     call MPI_Bcast(fba_para_x, 10*nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
!     call MPI_Bcast(fba_para_y, 10*nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
!     call MPI_Bcast(fba_para_y2, 10*nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
!#endif
     

     
!end if
  
contains

  subroutine nat_bangle(iba, imp1, imp2, imp3, xyz_ba)

    ! ---------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3
    integer, intent(inout) :: iba
    real(PREC), intent(in) :: xyz_ba(3, MXMP)

    ! ---------------------------------------------------------------
    ! local variables
    real(PREC) :: co_theta

    ! ---------------------------------------------------------------
    iba = iba + 1
    iba2type(iba) = BATYPE%PRO
    iba2mp(1, iba) = imp1
    iba2mp(2, iba) = imp2
    iba2mp(3, iba) = imp3
    call util_bondangle(imp1, imp2, imp3, co_theta, xyz_ba)
    ba_nat(iba) = acos(co_theta)
    call check_angle(iba,imp1,imp2,imp3,ba_nat(iba))
    factor_ba(iba) = inmisc%factor_local_unit(iunit, iunit)
    coef_ba(1, iba) = factor_ba(iba) * inpro%cba
    coef_ba(2, iba) = 0.0

  end subroutine nat_bangle

  subroutine nat_bangle_dna(iba, imp1, imp2, imp3, xyz_ba)

    ! ---------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3
    integer, intent(inout) :: iba
    real(PREC), intent(in) :: xyz_ba(3, MXMP)

    ! ---------------------------------------------------------------
    ! local variables
    real(PREC) :: co_theta

    ! ---------------------------------------------------------------
    iba = iba + 1
    iba2mp(1, iba) = imp1
    iba2mp(2, iba) = imp2
    iba2mp(3, iba) = imp3
    call util_bondangle(imp1, imp2, imp3, co_theta, xyz_ba)
    ba_nat(iba) = acos(co_theta)
    call check_angle(iba,imp1,imp2,imp3,ba_nat(iba))
    factor_ba(iba) = inmisc%factor_local_unit(iunit, iunit)
    coef_ba(1, iba) = factor_ba(iba) * indna%cba_dna
    coef_ba(2, iba) = 0.0

  end subroutine nat_bangle_dna

  subroutine nat_bangle_dna2(iba, imp1, imp2, imp3, i_type_ba, xyz_ba)

    ! ---------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3
    integer, intent(inout) :: iba
    integer, intent(in) :: i_type_ba
    real(PREC), intent(in) :: xyz_ba(SPACE_DIM, MXMP)

    ! ---------------------------------------------------------------
    ! local variables
    real(PREC) :: co_theta

    ! ---------------------------------------------------------------
    iba = iba + 1
    iba2type(iba) = i_type_ba
    iba2mp(1, iba) = imp1
    iba2mp(2, iba) = imp2
    iba2mp(3, iba) = imp3
    factor_ba(iba) = inmisc%factor_local_unit(iunit, iunit)

    select case (i_type_ba)
    case (BATYPE%DNA2_SPS)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_SPS
       ba_nat(iba) = indna2%nba_SPS * F_PI / 180.0e0_PREC
    case (BATYPE%DNA2_PSP)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_PSP
       ba_nat(iba) = indna2%nba_PSP * F_PI / 180.0e0_PREC
    case (BATYPE%DNA2_PSA)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_PSA
       ba_nat(iba) = indna2%nba_PSA * F_PI / 180.0e0_PREC
    case (BATYPE%DNA2_PST)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_PST
       ba_nat(iba) = indna2%nba_PST * F_PI / 180.0e0_PREC
    case (BATYPE%DNA2_PSC)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_PSC
       ba_nat(iba) = indna2%nba_PSC * F_PI / 180.0e0_PREC
    case (BATYPE%DNA2_PSG)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_PSG
       ba_nat(iba) = indna2%nba_PSG * F_PI / 180.0e0_PREC
    case (BATYPE%DNA2_ASP)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_ASP
       ba_nat(iba) = indna2%nba_ASP * F_PI / 180.0e0_PREC
    case (BATYPE%DNA2_TSP)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_TSP
       ba_nat(iba) = indna2%nba_TSP * F_PI / 180.0e0_PREC
    case (BATYPE%DNA2_CSP)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_CSP
       ba_nat(iba) = indna2%nba_CSP * F_PI / 180.0e0_PREC
    case (BATYPE%DNA2_GSP)
       coef_ba(1, iba) = factor_ba(iba) * indna2%cba_GSP
       ba_nat(iba) = indna2%nba_GSP * F_PI / 180.0e0_PREC
    case default
       error_message = 'Error: (4) logical defect in setp_native_bangle'
       call util_error(ERROR%STOP_ALL, error_message)
    end select

    coef_ba(2, iba) = 0.0e0_PREC

  end subroutine nat_bangle_dna2

  subroutine nat_bangle_rna(iba, iunit, imp1, imp2, imp3, cmp, i_type_ba, xyz_ba)

    ! ---------------------------------------------------------------
    integer, intent(in) :: iunit, imp1, imp2, imp3
    integer, intent(inout) :: iba
    integer, intent(in) :: i_type_ba
    real(PREC), intent(in) :: xyz_ba(SPACE_DIM, MXMP)
    character(4), intent(in) :: cmp

    ! ---------------------------------------------------------------
    ! local variables
    real(PREC) :: ba
    real(PREC) :: co_theta

    ! ---------------------------------------------------------------
    iba = iba + 1
    iba2type(iba) = i_type_ba
    iba2mp(1, iba) = imp1
    iba2mp(2, iba) = imp2
    iba2mp(3, iba) = imp3
    call util_bondangle(imp1, imp2, imp3, co_theta, xyz_ba)
    ba = acos(co_theta)
    call check_angle(iba,imp1,imp2,imp3,ba_nat(iba))
    factor_ba(iba) = inmisc%factor_local_unit(iunit, iunit)
  
    if (inmisc%flag_local_unit(iunit,iunit,LINTERACT%L_DTRNA)) then ! (Denesyuk's model)
       if (inmisc%i_dtrna_model == 2013) then
          select case (i_type_ba)
          case (BATYPE%RNA_RSP)
             coef_ba(1, iba) = factor_ba(iba) * indtrna13%ba_BSP
             if (cmp == ' Ab ') then
                ba_nat(iba) = inarna%angl_ASP
             else if (cmp == ' Gb ') then
                ba_nat(iba) = inarna%angl_GSP
             else
                error_message = 'Error: logical defect in setp_native_bangle'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          case (BATYPE%RNA_YSP)
             coef_ba(1, iba) = factor_ba(iba) * indtrna13%ba_BSP
             if (cmp == ' Ub ') then
                ba_nat(iba) = inarna%angl_USP
             else if (cmp == ' Cb ') then
                ba_nat(iba) = inarna%angl_CSP
             else
                error_message = 'Error: logical defect in setp_native_bangle'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          case (BATYPE%RNA_PSP)
             coef_ba(1, iba) = factor_ba(iba) * indtrna13%ba_PSP
             ba_nat(iba) = inarna%angl_PSP
          case (BATYPE%RNA_SPS)
             coef_ba(1, iba) = factor_ba(iba) * indtrna13%ba_SPS
             ba_nat(iba) = inarna%angl_SPS
          case (BATYPE%RNA_PSR)
             coef_ba(1, iba) = factor_ba(iba) * indtrna13%ba_PSB
             if (cmp == ' Ab ') then
                ba_nat(iba) = inarna%angl_PSA
             else if (cmp == ' Gb ') then
                ba_nat(iba) = inarna%angl_PSG
             else
                error_message = 'Error: logical defect in setp_native_bangle'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          case (BATYPE%RNA_PSY)
             coef_ba(1, iba) = factor_ba(iba) * indtrna13%ba_PSB
             if (cmp == ' Ub ') then
                ba_nat(iba) = inarna%angl_PSU
             else if (cmp == ' Cb ') then
                ba_nat(iba) = inarna%angl_PSC
             else
                error_message = 'Error: logical defect in setp_native_bangle'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          !case (BATYPE%RNA_BSB)
          !    coef_ba(1, iba) = factor_ba(iba) * indtrna13%cba_BSB
          case default
             error_message = 'Error: logical defect in setp_native_bangle'
             call util_error(ERROR%STOP_ALL, error_message)
          end select
   
          coef_ba(2, iba) = 0.0e0_PREC
   
       else if (inmisc%i_dtrna_model == 2015) then
          select case (i_type_ba)
          case (BATYPE%RNA_RSP)
             coef_ba(1, iba) = factor_ba(iba) * indtrna15%ba_BSP
             if (cmp == ' Ab ') then
                ba_nat(iba) = inarna%angl_ASP
             else if (cmp == ' Gb ') then
                ba_nat(iba) = inarna%angl_GSP
             else
                error_message = 'Error: logical defect in setp_native_bangle'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          case (BATYPE%RNA_YSP)
             coef_ba(1, iba) = factor_ba(iba) * indtrna15%ba_BSP
             if (cmp == ' Ub ') then
                ba_nat(iba) = inarna%angl_USP
             else if (cmp == ' Cb ') then
                ba_nat(iba) = inarna%angl_CSP
             else
                error_message = 'Error: logical defect in setp_native_bangle'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          case (BATYPE%RNA_PSP)
             coef_ba(1, iba) = factor_ba(iba) * indtrna15%ba_PSP
             ba_nat(iba) = inarna%angl_PSP
          case (BATYPE%RNA_SPS)
             coef_ba(1, iba) = factor_ba(iba) * indtrna15%ba_SPS
             ba_nat(iba) = inarna%angl_SPS
          case (BATYPE%RNA_PSR)
             coef_ba(1, iba) = factor_ba(iba) * indtrna15%ba_PSB
             if (cmp == ' Ab ') then
                ba_nat(iba) = inarna%angl_PSA
             else if (cmp == ' Gb ') then
                ba_nat(iba) = inarna%angl_PSG
             else
                error_message = 'Error: logical defect in setp_native_bangle'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          case (BATYPE%RNA_PSY)
             coef_ba(1, iba) = factor_ba(iba) * indtrna15%ba_PSB
             if (cmp == ' Ub ') then
                ba_nat(iba) = inarna%angl_PSU
             else if (cmp == ' Cb ') then
                ba_nat(iba) = inarna%angl_PSC
             else
                error_message = 'Error: logical defect in setp_native_bangle'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          case default
             error_message = 'Error: logical defect in setp_native_bangle'
             call util_error(ERROR%STOP_ALL, error_message)
          end select
   
          coef_ba(2, iba) = 0.0e0_PREC
       endif
    
    else  ! RNA model (Hori 2012)
       select case (i_type_ba)
       case (BATYPE%RNA_RSP)
          coef_ba(1, iba) = factor_ba(iba) * inrna%cba_RSP
       case (BATYPE%RNA_YSP)
          coef_ba(1, iba) = factor_ba(iba) * inrna%cba_YSP
       case (BATYPE%RNA_PSP)
          coef_ba(1, iba) = factor_ba(iba) * inrna%cba_PSP
       case (BATYPE%RNA_SPS)
          coef_ba(1, iba) = factor_ba(iba) * inrna%cba_SPS
       case (BATYPE%RNA_PSR)
          coef_ba(1, iba) = factor_ba(iba) * inrna%cba_PSR
       case (BATYPE%RNA_PSY)
          coef_ba(1, iba) = factor_ba(iba) * inrna%cba_PSY
       !case (BATYPE%RNA_BSB)
       !    coef_ba(1, iba) = factor_ba(iba) * inrna%cba_BSB
       case default
          error_message = 'Error: logical defect in setp_native_bangle'
          call util_error(ERROR%STOP_ALL, error_message)
       end select

       ba_nat(iba) = ba
       coef_ba(2, iba) = 0.0e0_PREC
    endif

  end subroutine nat_bangle_rna

!  subroutine nat_bangle_lipid(iba, imp1, imp2, imp3, xyz_ba)
  subroutine nat_bangle_lipid(iba, imp1, imp2, imp3)

    ! ---------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3
    integer, intent(inout) :: iba
!    real(PREC), intent(in) :: xyz_ba(3, MXMP)

    ! ---------------------------------------------------------------
    iba = iba + 1
    iba2mp(1, iba) = imp1
    iba2mp(2, iba) = imp2
    iba2mp(3, iba) = imp3
    ba_nat(iba) = 0.0
    factor_ba(iba) = inmisc%factor_local_unit(iunit, iunit)
    coef_ba(1, iba) = 0.0
    coef_ba(2, iba) = factor_ba(iba) * inlip%cba_lipid

  end subroutine nat_bangle_lipid

  subroutine nat_bangle_ligand(iba, imp1, imp2, imp3, xyz_ba)

    ! ---------------------------------------------------------------
    integer, intent(in) :: imp1, imp2, imp3
    integer, intent(inout) :: iba
    real(PREC), intent(in) :: xyz_ba(3, MXMP)

    ! ---------------------------------------------------------------
    ! local variables
    real(PREC) :: co_theta

    ! ---------------------------------------------------------------
    iba = iba + 1
    iba2mp(1, iba) = imp1
    iba2mp(2, iba) = imp2
    iba2mp(3, iba) = imp3
    call util_bondangle(imp1, imp2, imp3, co_theta, xyz_ba)
    ba_nat(iba) = acos(co_theta)
    call check_angle(iba,imp1,imp2,imp3,ba_nat(iba))
    factor_ba(iba) = inmisc%factor_local_unit(iunit, iunit)
    coef_ba(1, iba) = factor_ba(iba) * inligand%cba
    coef_ba(2, iba) = 0.0

  end subroutine nat_bangle_ligand

  subroutine check_angle(iba,imp1,imp2,imp3,theta)
    use const_index
    use const_maxsize
    use const_physical
    integer, intent(in)    :: iba,imp1, imp2, imp3
    real(PREC), intent(in) :: theta

    if (theta > WARN_ANGLE) then
       write(error_message,'(a,i6,a,i6,a,i6,a,i6,a,f6.2,a)') &
          'Warning: A bond angle is unnaturally large. iba=',&
          iba,', (imp1,imp2,imp3)=(',imp1,',',imp2,',',imp3,'),  theta=',theta/F_PI*180.0,'[deg]'
       call util_error(ERROR%WARN_ALL, error_message)
    endif
 endsubroutine check_angle

end subroutine setp_native_bangle
