! simu_energy
!> @brief The subroutine to calculate all the energy terms

!#define MEM_ALLOC

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! ************************************************************************
! subroutine for the energy
! ************************************************************************
subroutine simu_energy(irep,          &
                       velo_mp,       &
                       pnlet,         &
                       pnle_unit)

  use if_energy
  use const_maxsize
  use const_index
  use var_inp,     only : outfile
  use var_setp,    only : inmisc, inflp, inele
  use var_struct,  only : nunit_all, ncon, nmorse, nrna_bp
  use var_mgo,     only : inmgo

  use time
  use mpiconst

  implicit none

! --------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(in)  :: velo_mp(:,:)      ! (SPACE_DIM, nmp_real)
  real(PREC), intent(out) :: pnlet(:)          ! (E_TYPE%MAX)
  real(PREC), intent(out) :: pnle_unit(:,:,:)  ! (nunit_all, nunit_all, E_TYPE%MAX)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: ier
  integer :: i, iunit, junit
  real(PREC) :: sume
#ifdef MEM_ALLOC
  integer, allocatable :: now_con(:, :)  ! (ncon)
  integer, allocatable :: now_morse(:, :)
  integer, allocatable :: now_rna_bp(:, :)
  integer, allocatable :: now_allcon(:, :)
#else
  integer :: now_con(2, ncon)
  integer :: now_morse(2, nmorse)
  integer :: now_rna_bp(2, nrna_bp)
  integer :: now_allcon(2, ncon+nmorse+nrna_bp)
#endif
  character(CARRAY_MSG_ERROR),parameter :: msg_er_allocate = &
     'failed in memory allocation at mloop_simulator, PROGRAM STOP'
  character(CARRAY_MSG_ERROR),parameter :: msg_er_deallocate = &
     'failed in memory deallocation at mloop_simulator, PROGRAM STOP'

  integer :: n, tn
  real(PREC), allocatable :: pnlet_l(:,:)         !(E_TYPE%MAX, 0:nthreads-1)
  real(PREC), allocatable :: pnle_unit_l(:,:,:,:) !(nunit_all,nunit_all,E_TYPE%MAX,0:nthreads-1)
  integer, allocatable :: now_allcon_l(:, :)      !(2, ncon+nmorse+nrna_bp)


! ------------------------------------------------------------------------
! zero clear
  pnlet(:)         = 0.0e0_PREC
  pnle_unit(:,:,:) = 0.0e0_PREC

#ifdef MEM_ALLOC
  allocate(now_con(2, ncon), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  allocate(now_morse(2, nmorse), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  allocate(now_rna_bp(2, nrna_bp, stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  allocate(now_allcon(2, ncon+nmorse+nrna_bp), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
#endif

  now_con(:, :) = 0
  now_morse(:, :) = 0
  now_rna_bp(:, :) = 0
  now_allcon(:, :) = 0

! --------------------------------------------------------------------
  allocate( pnlet_l(E_TYPE%MAX, 0:nthreads-1), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  allocate( pnle_unit_l(nunit_all, nunit_all, E_TYPE%MAX, 0:nthreads-1), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  allocate( now_allcon_l(2, ncon+nmorse+nrna_bp), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)

!$omp parallel private(tn)
  tn = 0
!$  tn = omp_get_thread_num()

  pnlet_l(:,tn) = 0.0e0_PREC
  pnle_unit_l(:,:,:,tn) = 0.0e0_PREC

  TIME_S( tm_energy_velo) 
  call simu_energy_velo(velo_mp, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  TIME_E( tm_energy_velo) 

  TIME_S( tm_energy_bond) 
  call simu_energy_bond  (irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  TIME_E( tm_energy_bond) 

  TIME_S( tm_energy_bangle)
  call simu_energy_bangle(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  TIME_E( tm_energy_bangle)

  if (inmisc%force_flag_local(LINTERACT%L_AICG2) .or. &
      inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
     call simu_energy_aicg13_gauss(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  end if

  !if (inmisc%i_add_int == 1) then
  if (inflp%i_flp == 1 .or. inmisc%force_flag_local(LINTERACT%L_FLP)) then
     call simu_energy_fbangle(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  endif

  TIME_S( tm_energy_dih)
  call simu_energy_dih   (irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  TIME_E( tm_energy_dih)

  if (inmisc%force_flag_local(LINTERACT%L_AICG2)) then
     call simu_energy_aicg14_gauss(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  else if (inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
     call simu_energy_dih_gauss(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  end if

  !if (inmisc%i_add_int == 1) then
  if (inflp%i_flp == 1 .or. inmisc%force_flag_local(LINTERACT%L_FLP)) then
     call simu_energy_fdih(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  endif

  TIME_S( tm_energy_dih_harmonic) 
  call simu_energy_dih_harmonic   (irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  TIME_E( tm_energy_dih_harmonic) 

  call simu_energy_rna_stack(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))

  if (inmisc%i_dtrna_model == 2013) then
     call simu_energy_dtrna_stack(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     call simu_energy_dtrna_hbond13(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  else if (inmisc%i_dtrna_model == 2015) then
     call simu_energy_dtrna_stack_nlocal(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     call simu_energy_dtrna_stack(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     call simu_energy_dtrna_hbond15(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     call simu_energy_exv_dt15(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  endif

  if(inmgo%i_multi_mgo >= 1) then
     TIME_S( tm_energy_nlocal_mgo) 
     call simu_energy_nlocal_mgo(irep, now_con, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     call simu_energy_nlocal_rna_bp(irep, now_rna_bp, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     TIME_E( tm_energy_nlocal_mgo) 
  else if (inmisc%force_flag(INTERACT%ENM)) then
     TIME_S( tm_energy_enm) 
     call simu_energy_enm(irep, now_con, pnlet_l(:,tn), pnle_unit_l(:,:,:,tn))
     TIME_E( tm_energy_enm) 
  else
     TIME_S( tm_energy_nlocal_go) 
     call simu_energy_nlocal_go(irep, now_con, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     call simu_energy_nlocal_morse(irep, now_morse, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     call simu_energy_nlocal_rna_bp(irep, now_rna_bp, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     TIME_E( tm_energy_nlocal_go) 
  end if

!$omp barrier

!$omp master
  TIME_S( tm_energy_orderpara) 
  ! sum total
  now_allcon_l(1:2, 1:ncon) = now_con(1:2, 1:ncon)
  now_allcon_l(1:2, ncon+1:ncon+nmorse) = now_morse(1:2, 1:nmorse)
  now_allcon_l(1:2, ncon+nmorse+1:ncon+nmorse+nrna_bp) = now_rna_bp(1:2, 1:nrna_bp)
#ifdef MPI_PAR3
  call mpi_allreduce(now_allcon_l, now_allcon, &
       2*(ncon+nmorse+nrna_bp), MPI_INTEGER, &
       MPI_SUM, mpi_comm_local, ierr)
#else
  now_allcon(:,:) = now_allcon_l(:,:)
#endif

  call simu_energy_orderpara(irep, now_allcon)
  TIME_E( tm_energy_orderpara)
!$omp end master

  TIME_S( tm_energy_pnl) 
  if (inmisc%i_residue_exv_radii == 0) then
     call simu_energy_pnl (irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  else if (inmisc%i_residue_exv_radii == 1) then
     call simu_energy_pnl_restype(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  endif
  if (inmisc%force_flag(INTERACT%EXV_WCA)) then 
     call simu_energy_exv_wca (irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
  endif
  TIME_E( tm_energy_pnl) 

  TIME_S( tm_energy_ele)
  if (inmisc%force_flag(INTERACT%ELE)) then

     if (inele%i_function_form == 0) then     ! Debye-Huckel (default)

        if(inele%i_calc_method == 0) then         ! neighboring list (default)
           call simu_energy_ele(irep, pnlet_l(:,tn), pnle_unit_l(:,:,:,tn))
        else if(inele%i_calc_method == 1) then    ! neighboring list for K computer
           call simu_energy_ele2(irep, pnlet_l(:,tn), pnle_unit_l(:,:,:,tn))
        else                                      ! direct calculation for K computer
           call simu_energy_ele3(irep, pnlet_l(:,tn), pnle_unit_l(:,:,:,tn))
        end if

     elseif (inele%i_function_form == 1) then ! Coulomb potential
        call simu_energy_ele_coulomb(irep, pnlet_l(:,tn), pnle_unit_l(:,:,:,tn))
     endif
  endif
  TIME_E( tm_energy_ele)

  if (inmisc%class_flag(CLASS%ION)) then
     TIME_S( tm_energy_pnl) 
     call simu_energy_ion(irep, pnle_unit_l(:,:,:,tn), pnlet_l(:,tn))
     TIME_E( tm_energy_pnl) 
  endif

  if (inmisc%force_flag(INTERACT%HP)) then
     TIME_S( tm_energy_hp) 
     call simu_energy_hp(irep, pnlet_l(:,tn), pnle_unit_l(:,:,:,tn))
     TIME_E( tm_energy_hp) 
  endif

!sasa
  if (inmisc%force_flag(INTERACT%SASA)) then
     TIME_S( tm_energy_sasa)
     call simu_energy_sasa(irep, pnlet_l(:,tn))
     TIME_E( tm_energy_sasa)
  endif

!$omp end parallel

  TIME_S( tmc_energy )
  do n = 1, nthreads-1
    pnle_unit_l(1:nunit_all,1:nunit_all,1:E_TYPE%MAX,0) = &
    pnle_unit_l(1:nunit_all,1:nunit_all,1:E_TYPE%MAX,0) + &
    pnle_unit_l(1:nunit_all,1:nunit_all,1:E_TYPE%MAX,n)
    pnlet_l(1:E_TYPE%MAX,0) = &
    pnlet_l(1:E_TYPE%MAX,0) + &
    pnlet_l(1:E_TYPE%MAX,n)
  end do

  TIME_E( tmc_energy )

  TIME_S( tmc_energy )
#ifdef MPI_PAR3
  call mpi_allreduce(pnle_unit_l, pnle_unit, &
                     nunit_all*nunit_all*E_TYPE%MAX, PREC_MPI, &
                     MPI_SUM, mpi_comm_local, ierr)
  call mpi_allreduce(pnlet_l, pnlet, &
                     E_TYPE%MAX, PREC_MPI, &
                     MPI_SUM, mpi_comm_local, ierr)
#else
  pnle_unit(:,:,:) = pnle_unit_l(:,:,:,0)
  pnlet(:) = pnlet_l(:,0)
#endif
  TIME_E( tmc_energy )
 
  if(inmgo%i_multi_mgo >= 1) then
     TIME_S( tm_energy_mgo)
     call simu_energy_mgo(pnle_unit, pnlet)
     TIME_E( tm_energy_mgo)
  end if

  deallocate( pnlet_l,     stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)
  deallocate( pnle_unit_l, stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)
  deallocate( now_allcon_l,stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)


  if(inmisc%i_in_box == 1) then
     call simu_energy_box(irep, pnle_unit, pnlet)
  end if

  if(inmisc%i_in_cap == 1) then
     call simu_energy_cap(irep, pnle_unit, pnlet)
  end if

  if(inmisc%i_cylinder == 1) then
     call simu_energy_cylinder(irep, pnle_unit, pnlet)
  end if

  if(inmisc%i_bridge == 1) then
     call simu_energy_bridge(irep, pnle_unit, pnlet)
  end if

  if(inmisc%i_pulling == 1) then
     call simu_energy_pulling(irep, pnle_unit, pnlet)
  end if

  if(inmisc%i_anchor == 1) then
     call simu_energy_anchor(irep, pnle_unit, pnlet)
  end if

  if(inmisc%i_rest1d == 1) then
     call simu_energy_rest1d(irep, pnle_unit, pnlet)
  end if

  if(inmisc%i_window == 1) then
     call simu_energy_window(irep, pnle_unit, pnlet)
  end if
  
  if(inmisc%i_winz == 1) then
     call simu_energy_winz(irep, pnle_unit, pnlet)
  end if

  ! implicit ligand model
  if(inmisc%i_implig == 1) then
     call simu_energy_implig(irep, pnle_unit, pnlet, IMPLIGENERGY_TYPE%FOR_NON_MC)
     ! calculate implicit_ligand binding energy based on [state of implicit-ligand &  structure].
  end if
  
!  pnle_unit(1:nunit_all, 1:nunit_all, 1:E_TYPE%MAX) = pnle_unit_l(1:nunit_all, 1:nunit_all, 1:E_TYPE%MAX, 1)
!  pnlet(1:E_TYPE%MAX) = pnlet_l(1:E_TYPE%MAX, 1)

  ! --------------------------------------------------------------------
  ! calc energy of each interaction unit
  ! --------------------------------------------------------------------

  TIME_S( tm_energy_unit)

  if(inmisc%i_output_energy_style == 0) then
     do iunit = 1, nunit_all
        do junit = 1, iunit - 1
           pnle_unit(iunit, iunit, 1:E_TYPE%MAX)                     &
                =  pnle_unit(iunit, iunit, 1:E_TYPE%MAX)             &
                + 0.5e0_PREC * pnle_unit(junit, iunit, 1:E_TYPE%MAX)
        end do
        do junit = iunit + 1, nunit_all
           pnle_unit(iunit, iunit, 1:E_TYPE%MAX)                     &
                =  pnle_unit(iunit, iunit, 1:E_TYPE%MAX)             &
                + 0.5e0_PREC * pnle_unit(iunit, junit, 1:E_TYPE%MAX)
        end do
     end do
  else
     do iunit = 1, nunit_all
        do junit = iunit + 1, nunit_all
           do i = 1, E_TYPE%MAX
              sume = pnle_unit(iunit, junit, i) &
                   + pnle_unit(junit, iunit, i)
              pnle_unit(iunit, junit, i) = sume
              pnle_unit(junit, iunit, i) = sume
           end do
        end do
     end do

  end if
  do i = 1, E_TYPE%MAX
     if(i == E_TYPE%TOTAL .or. i == E_TYPE%VELO) cycle
     pnle_unit(1:nunit_all, 1:nunit_all, E_TYPE%TOTAL)                 &
          =  pnle_unit(1:nunit_all, 1:nunit_all, E_TYPE%TOTAL)  &
          + pnle_unit(1:nunit_all, 1:nunit_all, i)
  end do
  TIME_E( tm_energy_unit)

  ! --------------------------------------------------------------------
  ! calc total energy
  ! --------------------------------------------------------------------
  TIME_S( tm_energy_total)
  if(inmgo%i_multi_mgo >= 1) then
     do i = 1, E_TYPE%MAX
        if(i == E_TYPE%TOTAL  .or. i == E_TYPE%VELO    .or. i == E_TYPE%BOND .or. &
           i == E_TYPE%BANGLE .or. i == E_TYPE%DIHE    .or. i == E_TYPE%GO   .or. &
           i == E_TYPE%MORSE  .or. i == E_TYPE%DIHE_HARMONIC) cycle
        pnlet(E_TYPE%TOTAL) = pnlet(E_TYPE%TOTAL) + pnlet(i)
     end do
  else
     do i = 1, E_TYPE%MAX
        if(i == E_TYPE%TOTAL .or. i == E_TYPE%VELO) cycle
        pnlet(E_TYPE%TOTAL) = pnlet(E_TYPE%TOTAL) + pnlet(i)
     end do
  end if
  TIME_E( tm_energy_total)

#ifdef MEM_ALLOC
  deallocate( now_con )
  deallocate( now_morse )
  deallocate( now_rna_bp )
  deallocate( now_allcon )
#endif

end subroutine simu_energy
