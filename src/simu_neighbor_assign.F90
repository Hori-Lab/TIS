! simu_neighbor_assign
!> @brief Assignment and sorting of neighborling list for each energy type

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine simu_neighbor_assign(irep, ineigh2mp, lmp2neigh)
  
  use if_neighbor
  use const_maxsize
  use const_index
  use const_physical
  use var_inp,    only : inperi
  use var_setp,   only : inpro, inlip, inmisc, inrna, indtrna13, indtrna15
  use var_struct, only : nunit_real, nmp_real, lunit2mp, iontype_mp, &
                         pxyz_mp_rep, &
                         cmp2seq, imp2unit, lmp2con, icon2mp, coef_go, &
                         lmp2stack, istack2mp, &
                         ipnl2mp, lpnl, itail2mp, ltail2mp, imp2type, &
                         iclass_unit, ires_mp, nmp_all, &
                         lmp2morse, lmp2rna_bp, lmp2rna_st, &
                         imorse2mp, irna_bp2mp, irna_st2mp, &
                         coef_morse_a, coef_morse_fD, &
                         coef_rna_bp, coef_rna_bp_a, coef_rna_bp_fD, &
                         coef_rna_st, coef_rna_st_a, coef_rna_st_fD, & 
                        cmp2atom   !sasa
  use time
  use mpiconst

  implicit none

  ! -------------------------------------------------------------------
  integer, intent(in) :: irep
  integer, intent(in) :: lmp2neigh((nmp_l+nthreads-1)/nthreads  ,0:nthreads-1)
  integer, intent(in) :: ineigh2mp(MXMPNEIGHBOR*nmp_all/nthreads,0:nthreads-1)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: n
  integer :: klen, ksta, kend
  integer :: inum, imp, jmp, kmp, impmod, jmpmod, jimp
  integer :: imp1, imp2, imirror
  integer :: iunit, junit
  integer :: isep_nlocal
  integer :: isep_nlocal_rna
  integer :: icon, istack
  integer :: ipnl, npnl
  integer :: lcore, lint
  integer :: istart, isearch, isearch2, isearch_morse
  integer :: isearch_rna_bp, isearch_rna_st
  integer :: i_exvol, i_lip_brown, i_lip_nogu, i_dna, i_dna2, i_ion_hyd, i_ion_exv
  integer :: i_sasa, i_exv_wca, i_exv_dt15 !sasa
  integer :: ipnl2mp_l  (3, MXMPNEIGHBOR*nmp_all)
  integer :: ipnl2mp_pre(3, MXMPNEIGHBOR*nmp_all)
  integer :: npnl_lall(0:npar_mpi-1)
  !integer :: ii, iz
  real(PREC) :: vx(3)
  integer :: bp_dna2
  logical :: flg_nlocal_dna2, flg_bp_dna2

  type calc_type
     integer :: GO
     integer :: MORSE
     integer :: EXV
     integer :: DNA
     integer :: DNA2
     integer :: LIP_BROWN
     integer :: LIP_NOGU
     integer :: LIP_SOLV
     integer :: PAIR_RNA
     integer :: STACK_RNA
     integer :: ION_HYD
     integer :: ION_EXV
     integer :: AICG1
     integer :: AICG2
     integer :: SASA   !sasa
     integer :: EXV_WCA
     integer :: EXV_DT15
     integer :: MAX
  endtype calc_type
  type(calc_type), parameter :: CALC = calc_type(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,17)
  integer :: icalc(CALC%MAX, nunit_real, nunit_real)

  character(CARRAY_MSG_ERROR) :: error_message

  integer :: imp_l
  !integer :: lip_sta, lip_end
#ifdef SHARE_NEIGH_PNL
  integer :: ipnl2mp_l_sort(3, MXMPNEIGHBOR*nmp_all)
  integer :: disp(0:npar_mpi-1)
  integer :: count(0:npar_mpi-1)
  integer :: npnl_l
#endif

  ! -------------------------------------------------------------------
  isep_nlocal  = inpro%n_sep_nlocal

  ipnl     = 0
  isearch  = 1
  isearch_morse  = 1
  isearch_rna_bp = 1
  isearch_rna_st = 1
  isearch2 = 1

  lcore = inlip%num_lip_core
  lint = lcore + inlip%num_lip_int

  ! --------------------------------------------------------------------
  ! calc icalc
  icalc(1:CALC%MAX, 1:nunit_real, 1:nunit_real) = 0

  do iunit = 1, nunit_real
     do junit = iunit, nunit_real
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%GO)) then
           icalc(CALC%GO, iunit, junit) = 1
        end if
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%MORSE)) then
           icalc(CALC%MORSE, iunit, junit) = 1
        end if
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV)) then
           icalc(CALC%EXV, iunit, junit) = 1
        end if
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%DNA)) then
           icalc(CALC%DNA, iunit, junit) = 1
        end if
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%DNA2) .OR. &
           inmisc%flag_nlocal_unit(iunit, junit, INTERACT%DNA2C)) then
           icalc(CALC%DNA2, iunit, junit) = 1
        end if
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%LIP_BROWN)) then
           icalc(CALC%LIP_BROWN, iunit, junit) = 1
        end if
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%LIP_NOGU)) then
           icalc(CALC%LIP_NOGU, iunit, junit) = 1
        end if
        if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%LIP_SOLV)) then
           icalc(CALC%LIP_SOLV, iunit, junit) = 1
        endif
        if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%PAIR_RNA)) then
           icalc(CALC%PAIR_RNA, iunit, junit) = 1
        endif
        if ((iunit == junit) .AND. (iclass_unit(iunit) == CLASS%RNA)) then
           icalc(CALC%STACK_RNA, iunit, junit) = 1
        endif
        if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%ION_HYD)) then
           icalc(CALC%ION_HYD, iunit, junit) = 1
        endif
        if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%ION_EXV)) then
           icalc(CALC%ION_EXV, iunit, junit) = 1
        endif
        if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG1)) then  !AICG
           icalc(CALC%AICG1, iunit, junit) = 1
        endif
        if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG2)) then  !AICG
           icalc(CALC%AICG2, iunit, junit) = 1
        endif
        if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%SASA)) then  !sasa
           icalc(CALC%SASA, iunit, junit) = 1
        endif
        if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV_WCA)) then
           icalc(CALC%EXV_WCA, iunit, junit) = 1
        endif
        if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV_DT15)) then
           icalc(CALC%EXV_DT15, iunit, junit) = 1
        endif

     end do
  end do


  ! --------------------------------------------------------------------
  if( imp_l2g(1) /= 1 ) then
    isearch  = lmp2con  (imp_l2g(1)-1) + 1
    isearch2 = lmp2stack(imp_l2g(1)-1) + 1
    isearch_morse  = lmp2morse(imp_l2g(1)-1) + 1
    if (inmisc%class_flag(CLASS%RNA)) then
       isearch_rna_bp = lmp2rna_bp(imp_l2g(1)-1) + 1
       isearch_rna_st = lmp2rna_st(imp_l2g(1)-1) + 1
    endif
  end if

!!$omp parallel
!!$omp do private(klen,ksta,kend,imp,iunit,jmp,junit,kmp,isearch, &
!!$omp&           
  do n = 0, nthreads-1
  klen=(nmp_l-1+nthreads)/nthreads
  ksta=1+klen*n
  kend=min(ksta+klen-1,nmp_l)

  istart = 1
  do imp_l = ksta, kend
     imp = imp_l2g(imp_l)
     iunit = imp2unit(imp)

     !write(*, *) 'simu_neighbor_assign: lmp2neigh', lmp2neigh(imp_l-ksta+1,n)
     loop_lneigh: do inum = istart, lmp2neigh(imp_l-ksta+1,n)
        
        jmp = ineigh2mp(inum,n)
        junit = imp2unit(jmp)
        
        i_exvol = 1
        i_sasa  = 1 !sasa
        ! -----------------------------------------------------------------
        ! go
!        if(icalc(CALC%GO, iunit, junit) == 1) then
        if(icalc(CALC%GO, iunit, junit) == 1 .OR. &
           icalc(CALC%AICG1, iunit, junit) == 1 .OR. &
           icalc(CALC%AICG2, iunit, junit) == 1) then ! AICG

           do icon = isearch, lmp2con(imp)
              kmp = icon2mp(2, icon)

              if(jmp < kmp) exit

              if(jmp == kmp) then
                 isearch = icon + 1
!                 cycle loop_lneigh
                 if(coef_go(icon) > ZERO_JUDGE) then
                    i_exvol = 0
                 end if
              end if
           end do
        end if

        ! -----------------------------------------------------------------
        ! morse
        if(icalc(CALC%MORSE, iunit, junit) == 1) then
           do icon = isearch_morse, lmp2morse(imp)
              kmp = imorse2mp(2, icon)

              if(jmp < kmp) exit

              if(jmp == kmp) then
                 isearch_morse = icon + 1
!                 cycle loop_lneigh
                 if(coef_morse_a(icon) > ZERO_JUDGE .and. coef_morse_fD(icon) > ZERO_JUDGE) then
                    i_exvol = 0
                 end if
              end if
           end do
        end if

        ! -----------------------------------------------------------------
        ! RNA base pairing
        if(icalc(CALC%PAIR_RNA, iunit, junit) == 1) then
           do icon = isearch_rna_bp, lmp2rna_bp(imp)
              kmp = irna_bp2mp(2, icon)

              if(jmp < kmp) exit

              if(jmp == kmp) then
                 isearch_rna_bp = icon + 1
!                 cycle loop_lneigh
                 if(coef_rna_bp(icon) > ZERO_JUDGE .or. (coef_rna_bp_a(icon) > ZERO_JUDGE .and. coef_rna_bp_fD(icon) > ZERO_JUDGE)) then
                    i_exvol = 0
                 end if
                 i_exvol = 0
              end if
           end do
        end if

        ! -----------------------------------------------------------------
        ! RNA base stacking
        if(icalc(CALC%STACK_RNA, iunit, junit) == 1) then
           do icon = isearch_rna_st, lmp2rna_st(imp)
              kmp = irna_st2mp(2, icon)

              if(jmp < kmp) exit

              if(jmp == kmp) then
                 isearch_rna_st = icon + 1
!                 cycle loop_lneigh
                 if(coef_rna_st(icon) > ZERO_JUDGE .or. (coef_rna_st_a(icon) > ZERO_JUDGE .and. coef_rna_st_fD(icon) > ZERO_JUDGE)) then
                    i_exvol = 0
                 end if
              end if
           end do
        end if

        ! -----------------------------------------------------------------
        ! exvol
        if(icalc(CALC%EXV, iunit, junit) == 1) then
           if(iunit == junit) then
              if (iclass_unit(iunit) == CLASS%RNA) then
                 select case (imp2type(imp))
                 case (MPTYPE%RNA_PHOS) !P
                    isep_nlocal_rna = inrna%n_sep_nlocal_P
                 case (MPTYPE%RNA_SUGAR)!S
                    isep_nlocal_rna = inrna%n_sep_nlocal_S
                 case (MPTYPE%RNA_BASE) !B 
                    isep_nlocal_rna = inrna%n_sep_nlocal_B
                 case default 
                    error_message = 'Error: logical defect in simu_neighbor_assign'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endselect
                 
                 if (jmp < imp + isep_nlocal_rna) then
!                    cycle loop_lneigh
                    i_exvol = 0
                 end if

              else if (iclass_unit(iunit) == CLASS%LIG) then !explig
                 if(ires_mp(imp) == ires_mp(jmp)) then
!                    cycle loop_lneigh
                    i_exvol = 0
                 end if

              else
                 if(jmp < imp + isep_nlocal) then
!                    cycle loop_lneigh
                    i_exvol = 0
                 end if

              endif
           end if ! (iunit==junit)

           if(i_exvol == 1) then
              ipnl = ipnl + 1
              ipnl2mp_l(1, ipnl) = imp 
              ipnl2mp_l(2, ipnl) = jmp
              ipnl2mp_l(3, ipnl) = E_TYPE%EXV
!             cycle loop_lneigh
           end if
        end if
        
        ! -----------------------------------------------------------------
        ! DNA2-DNA2 (for 3SPN.2 model)
        i_dna2 = 1
        
        if (icalc(CALC%DNA2, iunit, junit) == 1) then

           flg_nlocal_dna2 = .false.
           flg_bp_dna2     = .false.
           
           if (iunit == junit) then
              ! check if the particle pair is in the non-local location
              if (abs(jmp - imp) > 4) flg_nlocal_dna2 = .true.
              ! check if the particle pair is W.C. base pair
              if (imp2type(imp) == MPTYPE%DNA2_BASE .and. imp2type(jmp) == MPTYPE%DNA2_BASE) then
                 call seq2bp(cmp2seq(imp), cmp2seq(jmp), bp_dna2, flg_bp_dna2)
              end if
           else
              ! check if the particle pair is in the non-local location (always!)
              flg_nlocal_dna2 = .true.
              ! check if the particle pair is W.C. base pair
              if (imp2type(imp) == MPTYPE%DNA2_BASE .and. imp2type(jmp) == MPTYPE%DNA2_BASE) then
                 call seq2bp(cmp2seq(imp), cmp2seq(jmp), bp_dna2, flg_bp_dna2)
              end if
           end if

           if (flg_nlocal_dna2) then

              if (.not. flg_bp_dna2) then
                 ! If the particle pair is in the non-local location
                 ! and does not participate in the base paring (stacking) interaction,
                 ! we impose excluded volume interaction
                 ipnl = ipnl + 1 
                 ipnl2mp_l(1, ipnl) = imp
                 ipnl2mp_l(2, ipnl) = jmp
                 ipnl2mp_l(3, ipnl) = E_TYPE%EXV_DNA2
              end if
              
           end if
           
        end if

        ! -----------------------------------------------------------------
        ! excluded volume (DTRNA)
        if(icalc(CALC%EXV_WCA, iunit, junit) == 1) then
           i_exv_wca = 1 
           if(iunit == junit) then
              if (iclass_unit(iunit) == CLASS%RNA) then
                 select case (imp2type(imp))
                 case (MPTYPE%RNA_PHOS) !P
                    isep_nlocal_rna = indtrna13%n_sep_nlocal_P
                 case (MPTYPE%RNA_SUGAR)!S
                    isep_nlocal_rna = indtrna13%n_sep_nlocal_S
                 case (MPTYPE%RNA_BASE) !B 
                    isep_nlocal_rna = indtrna13%n_sep_nlocal_B
                 case default 
                    error_message = 'Error: logical defect in simu_neighbor_assign'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endselect
                 
                 if (jmp < imp + isep_nlocal_rna) then
                    i_exv_wca = 0
                 end if

              else
                 !if(jmp < imp + isep_nlocal) then
                 !   i_exvol = 0
                 !end if
                 error_message = 'Error: EXV_WCA is only for RNA currently, in simu_neighbor_assign'
                 call util_error(ERROR%STOP_ALL, error_message)
              endif
           end if ! (iunit==junit)

           if(i_exv_wca == 1) then
              ipnl = ipnl + 1
              ipnl2mp_l(1, ipnl) = imp 
              ipnl2mp_l(2, ipnl) = jmp
              ipnl2mp_l(3, ipnl) = E_TYPE%EXV_WCA
           end if
        end if

        ! -----------------------------------------------------------------
        ! excluded volume (DTRNA)
        if(icalc(CALC%EXV_DT15, iunit, junit) == 1) then
           i_exv_dt15 = 0 
           if(iunit == junit) then
              i_exv_dt15 = 1 
              if (iclass_unit(iunit) == CLASS%RNA) then
                 select case (imp2type(imp))
                 case (MPTYPE%RNA_PHOS) !P
                    isep_nlocal_rna = indtrna15%n_sep_nlocal_P
                 case (MPTYPE%RNA_SUGAR)!S
                    isep_nlocal_rna = indtrna15%n_sep_nlocal_S
                 case (MPTYPE%RNA_BASE) !B 
                    isep_nlocal_rna = indtrna15%n_sep_nlocal_B
                 case default 
                    error_message = 'Error: logical defect in simu_neighbor_assign'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endselect
                 
                 if (jmp < imp + isep_nlocal_rna) then
                    i_exv_dt15 = 0
                 end if

              else if (iclass_unit(iunit) == CLASS%ION) then
                 continue

              else
                 !if(jmp < imp + isep_nlocal) then
                 !   i_exvol = 0
                 !end if
                 error_message = 'Error: EXV_DT15 is only for RNA currently, in simu_neighbor_assign'
                 call util_error(ERROR%STOP_ALL, error_message)
              endif

           elseif ((iclass_unit(iunit) == CLASS%RNA .and. iclass_unit(junit) == CLASS%ION) .OR. &
                   (iclass_unit(iunit) == CLASS%ION .and. iclass_unit(junit) == CLASS%RNA)) then

              i_exv_dt15 = 1 

           else
              error_message = 'Error: EXV_DT15 is only for RNA currently, in simu_neighbor_assign'
              call util_error(ERROR%STOP_ALL, error_message)
                
           end if ! (iunit==junit)

           if(i_exv_dt15 == 1) then
              ipnl = ipnl + 1
              ipnl2mp_l(1, ipnl) = imp 
              ipnl2mp_l(2, ipnl) = jmp
              ipnl2mp_l(3, ipnl) = E_TYPE%EXV_DT15
           end if
        end if

        ! -----------------------------------------------------------------
        ! DNA-DNA
        i_dna = 1
        if(icalc(CALC%DNA, iunit, junit) == 1) then

           ! base stacking
           if(iunit == junit) then
              do istack = isearch2, lmp2stack(imp)
                 kmp = istack2mp(2, istack)
                 if(jmp < kmp) exit

                 if(jmp == kmp) then
                    isearch2 = istack + 1
!                    cycle loop_lneigh
                    i_dna = 0
                 end if
              end do
           end if

!           impmod = mod(imp - lunit2mp(1, iunit), 3)
!           jmpmod = mod(jmp - lunit2mp(1, junit), 3)

           if(i_dna == 1) then
!              if(impmod == 1 .and. jmpmod == 1) then
              if(imp2type(imp) == MPTYPE%DNA_BASE .and. imp2type(jmp) == MPTYPE%DNA_BASE) then
                 ! base pair
                 if((cmp2seq(imp) == ' DA' .and. cmp2seq(jmp) == ' DT') .or. &
                      (cmp2seq(imp) == ' DT' .and. cmp2seq(jmp) == ' DA')) then
                    ipnl = ipnl + 1
                    ipnl2mp_l(1, ipnl) = imp 
                    ipnl2mp_l(2, ipnl) = jmp
                    ipnl2mp_l(3, ipnl) = E_TYPE%BP_AT
                    !                 cycle loop_lneigh
                    i_dna = 0
                    
                 else if((cmp2seq(imp) == ' DG' .and. cmp2seq(jmp) == ' DC') .or. &
                      (cmp2seq(imp) == ' DC' .and. cmp2seq(jmp) == ' DG')) then
                    ipnl = ipnl + 1
                    ipnl2mp_l(1, ipnl) = imp
                    ipnl2mp_l(2, ipnl) = jmp
                    ipnl2mp_l(3, ipnl) = E_TYPE%BP_GC
                    !                 cycle loop_lneigh
                    i_dna = 0

                    ! mismatch base pair
                 else
                    ipnl = ipnl + 1
                    ipnl2mp_l(1, ipnl) = imp 
                    ipnl2mp_l(2, ipnl) = jmp
                    ipnl2mp_l(3, ipnl) = E_TYPE%MBP
                    !                 cycle loop_lneigh
                    i_dna = 0
                 end if
              end if

              ! exvol dna; skip local part
              if(iunit == junit) then
                 jimp = jmp - imp
                 if(jimp <= 2) then
!                    if(impmod == 0) then
                    if(imp2type(imp) == MPTYPE%DNA_SUGAR) then
                       !                    cycle loop_lneigh
                       i_dna = 0
!                    else if(impmod == 1) then
                    else if(imp2type(imp) == MPTYPE%DNA_BASE) then
                       if(jimp <= 0) then
                          !                       cycle loop_lneigh
                          i_dna = 0
                       end if
                    else
                       if(jimp <= 1) then
                          !                    cycle loop_lneigh
                          i_dna = 0
                       end if
                    end if
                 end if
              end if
              
              ! skip phosphate-phosphate interaction when using explicit ion model
              if(icalc(CALC%ION_HYD, iunit, junit) == 1 .and. iontype_mp(imp) /= 0 .and. iontype_mp(imp) /= 0) then
                 i_dna = 0
              end if

              if(i_dna == 1) then
                 ipnl = ipnl + 1 
                 ipnl2mp_l(1, ipnl) = imp
                 ipnl2mp_l(2, ipnl) = jmp
                 ipnl2mp_l(3, ipnl) = E_TYPE%EXV_DNA
                 !              cycle loop_lneigh
              end if
           end if
        end if

        ! -----------------------------------------------------------------
        ! lipid (Brown)
        i_lip_brown = 1
        if(icalc(CALC%LIP_BROWN, iunit, junit) == 1) then

           impmod = mod(imp - lunit2mp(1, iunit), inlip%num_lip_total)
           jmpmod = mod(jmp - lunit2mp(1, junit), inlip%num_lip_total)
           if(iunit == junit) then
              if(jmp - imp <= inlip%num_lip_total - 1 .and. impmod < jmpmod) then
                 if(jmpmod - impmod <= 3) then
!                    cycle loop_lneigh
                    i_lip_brown = 0
                 end if
              end if
           end if
           
           if(i_lip_brown == 1) then
              if(impmod < lcore .or. jmpmod < lcore) then
                 ipnl = ipnl + 1 
                 ipnl2mp_l(1, ipnl) = imp 
                 ipnl2mp_l(2, ipnl) = jmp
                 ipnl2mp_l(3, ipnl) = E_TYPE%CORE
                 
              else if(impmod < lint .and. jmpmod < lint) then
                 ipnl = ipnl + 1 
                 ipnl2mp_l(1, ipnl) = imp 
                 ipnl2mp_l(2, ipnl) = jmp
                 ipnl2mp_l(3, ipnl) = E_TYPE%INT
              else
                 ipnl = ipnl + 1 
                 ipnl2mp_l(1, ipnl) = imp 
                 ipnl2mp_l(2, ipnl) = jmp
                 ipnl2mp_l(3, ipnl) = E_TYPE%TAIL
              end if
           end if
        end if

        ! -----------------------------------------------------------------
        ! lipid (Noguchi without solvation)

        i_lip_nogu = 1
        if(icalc(CALC%LIP_NOGU, iunit, junit) == 1) then

           impmod = mod(imp - lunit2mp(1, iunit), inlip%num_lip_total)
           jmpmod = mod(jmp - lunit2mp(1, junit), inlip%num_lip_total)
           if(iunit == junit) then
              if(jmp - imp <= inlip%num_lip_total - 1 .and. impmod < jmpmod) then
                 if(jmpmod - impmod <= 3) then
!                    cycle loop_lneigh
                    i_lip_nogu = 0
                 end if
              end if
           end if

           if(i_lip_nogu == 1) then
              ipnl = ipnl + 1 
              ipnl2mp_l(1, ipnl) = imp 
              ipnl2mp_l(2, ipnl) = jmp
              ipnl2mp_l(3, ipnl) = E_TYPE%CORE_NOGU
           end if
        end if
        
        ! -----------------------------------------------------------------
        ! solvation (lipid-lipid, lipid-protein)
!        if(icalc(CALC%LIP_SOLV, iunit, junit) == 1) then
!
!           impmod = mod(imp - lunit2mp(1, iunit), inlip%num_lip_total)
!           jmpmod = mod(jmp - lunit2mp(1, junit), inlip%num_lip_total)
!           if(iunit == junit) then
!              if(jmp - imp <= inlip%num_lip_total - 1 .and. impmod < jmpmod) then
!                 if(jmpmod - impmod <= 3) then
!                    cycle loop_lneigh
!                 end if
!              end if
!           end if
!
!           if(impmod >= lcore .and. jmpmod >= lcore) then
!              ipnl = ipnl + 1 
!              ipnl2mp_l(1, ipnl) = imp 
!              ipnl2mp_l(2, ipnl) = jmp
!              ipnl2mp_l(3, ipnl) = E_TYPE%TAIL_NOGU
!           end if
!
!        end if
        
        ! -----------------------------------------------------------------
        ! Ion
        if(icalc(CALC%ION_HYD, iunit, junit) == 1 .or. icalc(CALC%ION_EXV, iunit, junit) == 1) then

           i_ion_hyd = 0
           i_ion_exv = 0

           if(iontype_mp(imp) /= 0 .and. iontype_mp(jmp) /= 0 .and. icalc(CALC%ION_HYD, iunit, junit) == 1) then
              i_ion_hyd = 1

           else if(iontype_mp(imp) /= 0 .or. iontype_mp(jmp) /= 0) then

              if((iontype_mp(imp) == 0 .and. iontype_mp(jmp) == IONTYPE%P) .or. (iontype_mp(imp) == IONTYPE%P .and. iontype_mp(jmp) == 0)) then
                 ! skip interaction between phosphate and non-ion molecules
              else
                 i_ion_exv = 1
              end if

           end if

           if(i_ion_hyd == 1) then
              ipnl = ipnl + 1 
              ipnl2mp_l(1, ipnl) = imp 
              ipnl2mp_l(2, ipnl) = jmp
              ipnl2mp_l(3, ipnl) = E_TYPE%HYD_ION
           else if(i_ion_exv == 1) then
              ipnl = ipnl + 1 
              ipnl2mp_l(1, ipnl) = imp 
              ipnl2mp_l(2, ipnl) = jmp
              ipnl2mp_l(3, ipnl) = E_TYPE%EXV_ION
           end if
           
        end if
!sasa
        ! -----------------------------------------------------------------
        ! SASA
        if(icalc(CALC%SASA, iunit, junit) == 1) then
           if (iclass_unit(iunit) /= CLASS%PRO .and. &
               cmp2atom(imp) /= ' P  ' .and. cmp2atom(imp) /= ' O  ') i_sasa = 0
           if (iclass_unit(junit) /= CLASS%PRO .and. &
               cmp2atom(jmp) /= ' P  ' .and. cmp2atom(jmp) /= ' O  ') i_sasa = 0
           if(i_sasa == 1) then
              ipnl = ipnl + 1
              ipnl2mp_l(1, ipnl) = imp
              ipnl2mp_l(2, ipnl) = jmp
              ipnl2mp_l(3, ipnl) = E_TYPE%SASA
           end if
        end if
        ! ------------------------------------------------------------------         
     end do loop_lneigh

     istart   = lmp2neigh(imp_l-ksta+1,n) + 1 
     isearch  = lmp2con  (imp_l2g(min(imp_l+1,nmp_l))-1) + 1
     isearch2 = lmp2stack(imp_l2g(min(imp_l+1,nmp_l))-1) + 1
     isearch_morse = lmp2morse(imp_l2g(min(imp_l+1,nmp_l))-1) + 1
     if (inmisc%class_flag(CLASS%RNA)) then
        isearch_rna_bp= lmp2rna_bp(imp_l2g(min(imp_l+1,nmp_l))-1) + 1
        isearch_rna_st= lmp2rna_st(imp_l2g(min(imp_l+1,nmp_l))-1) + 1
     endif
  end do
  end do


  ! --------------------------------------------------------------------
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_PNL
  npnl_l = ipnl

  call simu_neighbor_sort(irep, npnl_l, ipnl2mp_l, ipnl2mp_l_sort)

  TIME_S( tmc_neighbor )

  call mpi_allgather(npnl_l   ,1,MPI_INTEGER, &
                     npnl_lall,1,MPI_INTEGER, &
                     mpi_comm_local,ierr)

  npnl = sum( npnl_lall(0:npar_mpi-1) )

  disp (0) = 0
  count(0) = 3*npnl_lall(0)
  do n = 1, npar_mpi-1
    disp (n) = disp(n-1) + 3*npnl_lall(n-1)
    count(n) = 3*npnl_lall(n)
  end do

  call mpi_allgatherv( ipnl2mp_l_sort,3*npnl_l  ,MPI_INTEGER, &
                       ipnl2mp_pre  ,count, disp,MPI_INTEGER, &
                       mpi_comm_local,ierr )

  TIME_E( tmc_neighbor )

  call simu_neighbor_sort(irep, npnl, ipnl2mp_pre, npnl_lall=npnl_lall)
#else
  npnl = ipnl
  npnl_lall(:) = 0  ! NIS aza
  npnl_lall(0) = npnl

  call simu_neighbor_sort(irep, npnl, ipnl2mp_l, ipnl2mp_pre)
  call simu_neighbor_sort(irep, npnl, ipnl2mp_pre, npnl_lall=npnl_lall)
#endif

#else
  npnl = ipnl
  npnl_lall(:) = 0  ! NIS aza
  npnl_lall(0) = npnl

  call simu_neighbor_sort(irep, npnl, ipnl2mp_l, ipnl2mp_pre)
  call simu_neighbor_sort(irep, npnl, ipnl2mp_pre, npnl_lall=npnl_lall)

#endif


  if(inperi%i_periodic == 1) then
!     call util_pbindex(ipnl2mp, npnl, irep)
     do ipnl = 1, npnl
        imp1 = ipnl2mp(1, ipnl, irep)
        imp2 = ipnl2mp(2, ipnl, irep)
 
        vx(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        
!        do ix = 1, 3
!           if(vx(ix) > inperi%psizeh(ix)) then
!           vx(ix) = vx(ix) - inperi%psize(ix)
!              imi(ix) = 1
!           else if(vx(ix) < -inperi%psizeh(ix)) then
!              !           vx(ix) = vx(ix) + inperi%psize(ix)
!              imi(ix) = 2
!           else
!              imi(ix) = 0
!           end if
!        end do
!        imirror = 9*imi(1) + 3*imi(2) + imi(3)

        call util_pbneighbor(vx, imirror)

        ipnl2mp(3, ipnl, irep) = imirror
     end do

  end if


  ! -----------------------------------------------------------------
  ! lipid (Noguchi)
!  if (inmisc%class_flag(CLASS%LIP)) then
!     iz = 0
!     do iunit = 1, nunit_real
!        if(icalc(CALC%LIP_SOLV, iunit, iunit) /= 1) cycle
!
!#if defined(MPI_PAR2) && defined(MPI_PAR3)
!        lip_sta = lunit2mp(1, iunit)
!        lip_end = lunit2mp(2, iunit)
!        klen=(lip_end-lip_sta+nprocs)/nprocs
!        ksta=lip_sta+klen*myrank
!        kend=min(ksta+klen-1,lip_end)
!
!        do ii = ksta, kend
!#else
!        do ii = lunit2mp(1, iunit), lunit2mp(2, iunit)
!#endif
!           ltail2mp(1, ii, irep) = iz + 1
!   
!           do ipnl = lpnl(1, E_TYPE%TAIL_NOGU, irep), lpnl(2, E_TYPE%TAIL_NOGU, irep)
!              if(ipnl2mp(1, ipnl, irep) == ii) then
!                 iz = iz + 1
!                 itail2mp(iz, irep) = ipnl2mp(2, ipnl, irep)
!              else if(ipnl2mp(2, ipnl, irep) == ii) then
!                 iz = iz + 1
!                 itail2mp(iz, irep) = ipnl2mp(1, ipnl, irep)
!              end if
!           end do
!           
!           ltail2mp(2, ii, irep) = iz
!        end do
!     end do
!  endif

end subroutine simu_neighbor_assign
