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
  use var_setp,   only : inpro, inmisc, inrna, indtrna13, indtrna15
  use var_struct, only : nunit_real, iontype_mp, pxyz_mp_rep, &
                         imp2unit, lmp2con, icon2mp, coef_go, iexv2mp, imp2type, &
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
  integer :: inum, imp, jmp, kmp
  integer :: imp1, imp2, imirror
  integer :: iunit, junit
  integer :: isep_nlocal
  integer :: isep_nlocal_rna
  integer :: icon, iexv, nexv
  integer :: istart, isearch, isearch_morse
  integer :: isearch_rna_bp, isearch_rna_st
  integer :: i_exvol, i_ion_hyd, i_ion_exv
  integer :: i_sasa, i_exv_wca, i_exv_dt15 !sasa
  integer :: iexv2mp_l  (3, MXMPNEIGHBOR*nmp_all)
  integer :: iexv2mp_pre(3, MXMPNEIGHBOR*nmp_all)
  integer :: nexv_lall(0:npar_mpi-1)
  !integer :: ii, iz
  real(PREC) :: vx(3)

  type calc_type
     integer :: GO
     integer :: MORSE
     integer :: EXV
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
  type(calc_type), parameter :: CALC = calc_type(1,2,3,4,5,6,7,8,9,10,11,12,12)
  integer :: icalc(CALC%MAX, nunit_real, nunit_real)

  character(CARRAY_MSG_ERROR) :: error_message

  integer :: imp_l
#ifdef SHARE_NEIGH_PNL
  integer :: iexv2mp_l_sort(3, MXMPNEIGHBOR*nmp_all)
  integer :: disp(0:npar_mpi-1)
  integer :: count(0:npar_mpi-1)
  integer :: nexv_l
#endif

  ! -------------------------------------------------------------------
  isep_nlocal  = inpro%n_sep_nlocal

  iexv     = 0
  isearch  = 1
  isearch_morse  = 1
  isearch_rna_bp = 1
  isearch_rna_st = 1

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
              iexv = iexv + 1
              iexv2mp_l(1, iexv) = imp 
              iexv2mp_l(2, iexv) = jmp
              iexv2mp_l(3, iexv) = E_TYPE%EXV
!             cycle loop_lneigh
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
              iexv = iexv + 1
              iexv2mp_l(1, iexv) = imp 
              iexv2mp_l(2, iexv) = jmp
              iexv2mp_l(3, iexv) = E_TYPE%EXV_WCA
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
              iexv = iexv + 1
              iexv2mp_l(1, iexv) = imp 
              iexv2mp_l(2, iexv) = jmp
              iexv2mp_l(3, iexv) = E_TYPE%EXV_DT15
           end if
        end if
        
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
              iexv = iexv + 1 
              iexv2mp_l(1, iexv) = imp 
              iexv2mp_l(2, iexv) = jmp
              iexv2mp_l(3, iexv) = E_TYPE%HYD_ION
           else if(i_ion_exv == 1) then
              iexv = iexv + 1 
              iexv2mp_l(1, iexv) = imp 
              iexv2mp_l(2, iexv) = jmp
              iexv2mp_l(3, iexv) = E_TYPE%EXV_ION
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
              iexv = iexv + 1
              iexv2mp_l(1, iexv) = imp
              iexv2mp_l(2, iexv) = jmp
              iexv2mp_l(3, iexv) = E_TYPE%SASA
           end if
        end if
        ! ------------------------------------------------------------------         
     end do loop_lneigh

     istart   = lmp2neigh(imp_l-ksta+1,n) + 1 
     isearch  = lmp2con  (imp_l2g(min(imp_l+1,nmp_l))-1) + 1
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
  nexv_l = iexv

  call simu_neighbor_sort(irep, nexv_l, iexv2mp_l, iexv2mp_l_sort)

  TIME_S( tmc_neighbor )

  call mpi_allgather(nexv_l   ,1,MPI_INTEGER, &
                     nexv_lall,1,MPI_INTEGER, &
                     mpi_comm_local,ierr)

  nexv = sum( nexv_lall(0:npar_mpi-1) )

  disp (0) = 0
  count(0) = 3*nexv_lall(0)
  do n = 1, npar_mpi-1
    disp (n) = disp(n-1) + 3*nexv_lall(n-1)
    count(n) = 3*nexv_lall(n)
  end do

  call mpi_allgatherv( iexv2mp_l_sort,3*nexv_l  ,MPI_INTEGER, &
                       iexv2mp_pre  ,count, disp,MPI_INTEGER, &
                       mpi_comm_local,ierr )

  TIME_E( tmc_neighbor )

  call simu_neighbor_sort(irep, nexv, iexv2mp_pre, nexv_lall=nexv_lall)
#else
  nexv = iexv
  nexv_lall(:) = 0  ! NIS aza
  nexv_lall(0) = nexv

  call simu_neighbor_sort(irep, nexv, iexv2mp_l, iexv2mp_pre)
  call simu_neighbor_sort(irep, nexv, iexv2mp_pre, nexv_lall=nexv_lall)
#endif

#else
  nexv = iexv
  nexv_lall(:) = 0  ! NIS aza
  nexv_lall(0) = nexv

  call simu_neighbor_sort(irep, nexv, iexv2mp_l, iexv2mp_pre)
  call simu_neighbor_sort(irep, nexv, iexv2mp_pre, nexv_lall=nexv_lall)

#endif


  if(inperi%i_periodic == 1) then
!     call util_pbindex(iexv2mp, nexv, irep)
     do iexv = 1, nexv
        imp1 = iexv2mp(1, iexv, irep)
        imp2 = iexv2mp(2, iexv, irep)
 
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

        iexv2mp(3, iexv, irep) = imirror
     end do

  end if

end subroutine simu_neighbor_assign
