!energy_sumup
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
subroutine energy_sumup(irep,          &
                        velo_mp,       &
                        energy,        &
                        energy_unit )

  use if_energy
  use const_maxsize
  use const_index
  use var_setp,    only : inmisc, inele !,inflp
  use var_struct,  only : nunit_all, ncon, nLJ, ndtrna_st, ndtrna_tst, ndtrna_hb, nLJ1210, nHPS
  use var_simu,    only : st_status, ene_st, ene_tst, ene_hb
  use time
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(in)  :: velo_mp(:,:)      ! (SDIM, nmp_real)
  real(PREC), intent(out) :: energy(:)          ! (E_TYPE%MAX)
  real(PREC), intent(out) :: energy_unit(:,:,:)  ! (nunit_all, nunit_all, E_TYPE%MAX)

  integer :: ier
  integer :: i, iunit, junit
  real(PREC) :: sume
#ifdef MEM_ALLOC
  integer, allocatable :: now_con(:, :)  ! (ncon)
  integer, allocatable :: now_LJ(:, :)
  integer, allocatable :: now_HPS(:, :)
!   integer, allocatable :: now_LJ1210(:, :)
!  integer, allocatable :: now_morse(:, :)
!  integer, allocatable :: now_rna_bp(:, :)
  integer, allocatable :: now_allcon(:, :)
#else
  integer :: now_con(2, ncon)
  integer :: now_LJ(2, nLJ)
  integer :: now_HPS(2, nHPS)
!   integer :: now_LJ1210(2, nLJ1210)
!  integer :: now_morse(2, nmorse)
!  integer :: now_rna_bp(2, nrna_bp)
  !integer :: now_allcon(2, ncon+nmorse+nrna_bp+nLJ)
  integer :: now_allcon(2, ncon+nLJ)

#endif
  character(CARRAY_MSG_ERROR),parameter :: msg_er_allocate = &
     'failed in memory allocation at mloop_simulator, PROGRAM STOP'
  character(CARRAY_MSG_ERROR),parameter :: msg_er_deallocate = &
     'failed in memory deallocation at mloop_simulator, PROGRAM STOP'

  integer :: n, tn
  real(PREC), allocatable :: energy_l(:,:)          !(E_TYPE%MAX, 0:nthreads-1)
  real(PREC), allocatable :: energy_unit_l(:,:,:,:) !(nunit_all,nunit_all,E_TYPE%MAX,0:nthreads-1)
  integer,    allocatable :: now_allcon_l(:, :)     !(2, ncon+nmorse+nrna_bp+nLJ)
  real(PREC), allocatable :: ene_st_l(:,:)          !(1:ndtrna_st, 0:nthreads-1)
  real(PREC), allocatable :: ene_tst_l(:,:)         !(1:ndtrna_st, 0:nthreads-1)
  logical,    allocatable :: st_status_l(:,:)       !(1:ndtrna_st, 0:nthreads-1)
  real(PREC), allocatable :: ene_hb_l(:,:)          !(1:ndtrna_hb, 0:nthreads-1)

#ifdef _DEBUG
  write(6,*) '######## start energy_sumup'
#endif

! ------------------------------------------------------------------------
! zero clear
  energy(:)         = 0.0e0_PREC
  energy_unit(:,:,:) = 0.0e0_PREC

  if (inmisc%i_dtrna_model /= 0) then
     ene_st(:,irep)   = 0.0e0_PREC
     ene_tst(:,irep)  = 0.0e0_PREC

     if (inmisc%i_dtrna_model == 2013 .or. inmisc%i_dtrna_model == 2018) then
        ene_hb(:,irep) = 0.0e0_PREC
     endif
  endif

#ifdef MEM_ALLOC
  allocate(now_con(2, ncon), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  allocate(now_LJ(2, nLJ), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
   allocate(now_HPS(2, nHPS), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
!   allocate(now_LJ1210(2, nLJ1210), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
!  allocate(now_morse(2, nmorse), stat=ier)
!  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
!  allocate(now_rna_bp(2, nrna_bp, stat=ier)
!  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
!  allocate(now_allcon(2, ncon+nmorse+nrna_bp+nLJ), stat=ier)
  allocate(now_allcon(2, ncon+nLJ), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
#endif

  now_con(:, :) = 0
  now_LJ(:, :) = 0
  now_HPS(:, :) = 0
!   now_LJ1210(:, :) = 0
!  now_morse(:, :) = 0
!  now_rna_bp(:, :) = 0
  now_allcon(:, :) = 0

! --------------------------------------------------------------------
  allocate( energy_l(E_TYPE%MAX, 0:nthreads-1), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  allocate( energy_unit_l(nunit_all, nunit_all, E_TYPE%MAX, 0:nthreads-1), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  !allocate( now_allcon_l(2, ncon+nmorse+nrna_bp+nLJ), stat=ier)
  allocate( now_allcon_l(2, ncon+nLJ), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)

  if (inmisc%i_dtrna_model /= 0) then
     allocate( ene_st_l(1:ndtrna_st, 0:nthreads-1), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
     ene_st_l(:,:) = 0.0e0_PREC
     allocate( ene_tst_l(1:ndtrna_tst, 0:nthreads-1), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
     ene_tst_l(:,:) = 0.0e0_PREC

     if (inmisc%i_dtrna_model == 2013 .or. inmisc%i_dtrna_model == 2018) then
        allocate( ene_hb_l(1:ndtrna_hb, 0:nthreads-1), stat=ier)
        if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
        ene_hb_l(:,:) = 0.0e0_PREC
     endif

     if (inmisc%i_dtrna_model == 2015 .or. &
         inmisc%i_dtrna_model == 2019 ) then
        allocate( st_status_l(1:ndtrna_st, 0:nthreads-1), stat=ier)
        if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
        st_status_l(:,:) = .True.
     endif
  endif

!$omp parallel private(tn)
  tn = 0
!$  tn = omp_get_thread_num()

  energy_l(:,tn) = 0.0e0_PREC
  energy_unit_l(:,:,:,tn) = 0.0e0_PREC

!$omp master
  TIME_S( tm_energy_kinetic)
!$omp end master

  call energy_kinetic(velo_mp, energy_unit_l(:,:,:,tn), energy_l(:,tn))

!$omp master
  TIME_E( tm_energy_kinetic)
  TIME_S( tm_energy_local)
!$omp end master

  call energy_bond  (irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  call energy_fene  (irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  call energy_rouse (irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))

  call energy_bangle(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))

!$omp master
  TIME_E( tm_energy_local)
!$omp end master

!  !if (inmisc%i_add_int == 1) then
!  if (inflp%i_flp == 1 .or. inmisc%force_flag_local(LINTERACT%L_FLP)) then
!     call energy_fbangle(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
!  endif

!  TIME_S( tm_energy_dih)
!  call energy_dih   (irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
!  TIME_E( tm_energy_dih)

!  !if (inmisc%i_add_int == 1) then
!  if (inflp%i_flp == 1 .or. inmisc%force_flag_local(LINTERACT%L_FLP)) then
!     call energy_fdih(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
!  endif

!  TIME_S( tm_energy_dih_harmonic) 
!  call energy_dih_harmonic   (irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
!  TIME_E( tm_energy_dih_harmonic) 

!  call energy_rna_stack(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))

  if (inmisc%i_dtrna_model == 2013) then
     call energy_dtrna_stack(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn), ene_st_l(:,tn))
     call energy_dtrna_hbond13(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn), ene_hb_l(:,tn))

  else if (inmisc%i_dtrna_model == 2018) then
     call energy_dtrna_stack(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn), ene_st_l(:,tn))
     call energy_dtrna_hbond13(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn), ene_hb_l(:,tn))
     call energy_exv_dt15(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))

  else if (inmisc%i_dtrna_model == 2015 .or.&
           inmisc%i_dtrna_model == 2019 ) then
     call energy_dtrna_stack_nlocal(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn), ene_tst_l(:,tn), st_status_l(:,tn))

!############### Collect all information into st_status_l(:,0) within a process
!$omp barrier
!$omp critical
        st_status_l(1:ndtrna_st,0) = st_status_l(1:ndtrna_st,0) .and. st_status_l(1:ndtrna_st,tn)
!$omp end critical
!$omp barrier

!############### Collect and distribute all information in ALL processes to st_status
!$omp master
#ifdef MPI_PAR3
     call mpi_allreduce(st_status_l(:,0), st_status(:,irep), ndtrna_st,&
                        MPI_LOGICAL, MPI_LAND, mpi_comm_local, ierr)
#else
     st_status(:,irep) = st_status_l(:,0)
#endif
!$omp end master
!$omp barrier

     call energy_dtrna_stack(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn), ene_st_l(:,tn))
     call energy_dtrna_hbond15(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
        !!!! Do not use "if (local_rank_mpi == 0) for this routine,
        !!!! even though only local_rank_mpi=0 process is relevant to additional force.
        !!!! Other processes also have to execute this routine
        !!!! to keep consistency of random number generator (mts).
     call energy_exv_dt15(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  endif

!  if(inmgo%i_multi_mgo >= 1) then
!     TIME_S( tm_energy_nlocal_mgo) 
!     call energy_nlocal_mgo(irep, now_con, energy_unit_l(:,:,:,tn), energy_l(:,tn))
!     call energy_nlocal_rna_bp(irep, now_rna_bp, energy_unit_l(:,:,:,tn), energy_l(:,tn))
!     TIME_E( tm_energy_nlocal_mgo) 
!  else if (inmisc%force_flag(INTERACT%ENM)) then
!     TIME_S( tm_energy_enm) 
!     call energy_enm(irep, now_con, energy_l(:,tn), energy_unit_l(:,:,:,tn))
!     TIME_E( tm_energy_enm) 
!  else
!$omp master
     TIME_S( tm_energy_nlocal_go) 
!$omp end master
     call energy_LJ(irep, now_LJ, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  if(inmisc%force_flag(INTERACT%LJ1210)) then
     call energy_LJ_1210(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  end if
  if(inmisc%force_flag(INTERACT%HPS)) then
     call energy_HPS(irep, now_HPS, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  end if
     call energy_wca(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
!     call energy_nlocal_go(irep, now_con, energy_unit_l(:,:,:,tn), energy_l(:,tn))
!     call energy_nlocal_morse(irep, now_morse, energy_unit_l(:,:,:,tn), energy_l(:,tn))
!     call energy_nlocal_rna_bp(irep, now_rna_bp, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  if (inmisc%force_flag(INTERACT%CON_GAUSS)) then 
     call energy_con_gauss(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  endif
!$omp master
     TIME_E( tm_energy_nlocal_go) 
!$omp end master
!  end if

!$omp barrier

!$omp master
  TIME_S( tm_energy_orderpara) 
  ! sum total
  now_allcon_l(1:2, 1:ncon) = now_con(1:2, 1:ncon)
!  now_allcon_l(1:2, ncon+1:ncon+nmorse) = now_morse(1:2, 1:nmorse)
!  now_allcon_l(1:2, ncon+nmorse+1:ncon+nmorse+nrna_bp) = now_rna_bp(1:2, 1:nrna_bp)
!  now_allcon_l(1:2, ncon+nmorse+nrna_bp+1:ncon+nmorse+nrna_bp+nLJ) = now_LJ(1:2, 1:nLJ)
  now_allcon_l(1:2, ncon+1:ncon+nLJ) = now_LJ(1:2, 1:nLJ)
#ifdef MPI_PAR3
  call mpi_allreduce(now_allcon_l, now_allcon, &
       !2*(ncon+nmorse+nrna_bp+nLJ), MPI_INTEGER, &
       2*(ncon+nLJ), MPI_INTEGER, &
       MPI_SUM, mpi_comm_local, ierr)
#else
  now_allcon(:,:) = now_allcon_l(:,:)
#endif

  call energy_orderpara(irep, now_allcon)
  TIME_E( tm_energy_orderpara)
!$omp end master

!$omp master
  TIME_S( tm_energy_exv) 
!$omp end master
  if (inmisc%i_residuenergy_radii == 0) then
     call energy_exv_rep12 (irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
     call energy_exv_rep6 (irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  else if (inmisc%i_residuenergy_radii == 1) then
     call energy_exv_restype(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  endif
  if (inmisc%force_flag(INTERACT%EXV_WCA)) then 
     call energy_exv_wca (irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  endif
  if (inmisc%force_flag(INTERACT%EXV_GAUSS)) then 
     call energy_exv_gauss (irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  endif

  if (inmisc%i_BBR == 1) then
     call energy_BBR(irep, energy_unit_l(:,:,:,tn), energy_l(:,tn))
  endif
!$omp master
  TIME_E( tm_energy_exv) 

  TIME_S( tm_energy_ele)
!$omp end master
  if (inmisc%force_flag(INTERACT%ELE)) then

     if (inele%i_function_form == 0) then     ! Debye-Huckel (default)

        call energy_ele_DH(irep, energy_l(:,tn), energy_unit_l(:,:,:,tn))
        !if(inele%i_calc_method == 0) then         ! neighboring list (default)
        !   call energy_ele(irep, energy_l(:,tn), energy_unit_l(:,:,:,tn))
        !else if(inele%i_calc_method == 1) then    ! neighboring list for K computer
        !   call energy_ele2(irep, energy_l(:,tn), energy_unit_l(:,:,:,tn))
        !else                                      ! direct calculation for K computer
        !   call energy_ele3(irep, energy_l(:,tn), energy_unit_l(:,:,:,tn))
        !end if

     elseif (inele%i_function_form == 1) then ! Coulomb potential
        call energy_ele_coulomb(irep, energy_l(:,tn), energy_unit_l(:,:,:,tn))

     elseif (inele%i_function_form == 2) then ! Coulomb potential (Ewald)
        call energy_ele_coulomb_ewld(irep, energy_l(:,tn), energy_unit_l(:,:,:,tn))

     elseif (inele%i_function_form == 3) then ! Coulomb potential (Brute-force to check Ewald results)
        call energy_ele_coulomb_brute(irep, energy_l(:,tn), energy_unit_l(:,:,:,tn))
     endif
  endif
!$omp master
  TIME_E( tm_energy_ele)
!$omp end master

!  if (inmisc%force_flag(INTERACT%HP)) then
!     TIME_S( tm_energy_hp) 
!     call energy_hp(irep, energy_l(:,tn), energy_unit_l(:,:,:,tn))
!     TIME_E( tm_energy_hp) 
!  endif

!sasa
!  if (inmisc%force_flag(INTERACT%SASA)) then
!     TIME_S( tm_energy_sasa)
!     call energy_sasa(irep, energy_l(:,tn))
!     TIME_E( tm_energy_sasa)
!  endif

!$omp end parallel

  TIME_S( tmc_energy )
  do n = 1, nthreads-1
    energy_unit_l(1:nunit_all,1:nunit_all,1:E_TYPE%MAX,0) = &
    energy_unit_l(1:nunit_all,1:nunit_all,1:E_TYPE%MAX,0) + &
    energy_unit_l(1:nunit_all,1:nunit_all,1:E_TYPE%MAX,n)
    energy_l(1:E_TYPE%MAX,0) = &
    energy_l(1:E_TYPE%MAX,0) + &
    energy_l(1:E_TYPE%MAX,n)
  end do

  TIME_E( tmc_energy )

  TIME_S( tmc_energy )
#ifdef MPI_PAR3
  call mpi_allreduce(energy_unit_l, energy_unit, &
                     nunit_all*nunit_all*E_TYPE%MAX, PREC_MPI, &
                     MPI_SUM, mpi_comm_local, ierr)
  call mpi_allreduce(energy_l, energy, &
                     E_TYPE%MAX, PREC_MPI, &
                     MPI_SUM, mpi_comm_local, ierr)
#else
  energy_unit(:,:,:) = energy_unit_l(:,:,:,0)
  energy(:) = energy_l(:,0)
#endif

  deallocate( energy_l,     stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)
  deallocate( energy_unit_l, stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)
  deallocate( now_allcon_l,stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)


  if (inmisc%i_dtrna_model /= 0) then

     do n = 1, nthreads-1
        ene_st_l(:,0) = ene_st_l(:,0) + ene_st_l(:,n)
        ene_tst_l(:,0) = ene_tst_l(:,0) + ene_tst_l(:,n)
     end do
     if (inmisc%i_dtrna_model == 2013 .or. inmisc%i_dtrna_model == 2018) then
        do n = 1, nthreads-1
           ene_hb_l(:,0) = ene_hb_l(:,0) + ene_hb_l(:,n)
        enddo
     endif

#ifdef MPI_PAR3
     if (ndtrna_st > 0 ) then
        call mpi_allreduce(ene_st_l, ene_st(1,irep), &
                           ndtrna_st, PREC_MPI, &
                           MPI_SUM, mpi_comm_local, ierr)
     endif

     if (ndtrna_tst > 0 ) then
        call mpi_allreduce(ene_tst_l, ene_tst(1,irep), &
                           ndtrna_tst, PREC_MPI, &
                           MPI_SUM, mpi_comm_local, ierr)
     endif

     if ((inmisc%i_dtrna_model == 2013 .or. inmisc%i_dtrna_model == 2018) .and. ndtrna_hb > 0) then
        call mpi_allreduce(ene_hb_l, ene_hb(1,irep), &
                           ndtrna_hb, PREC_MPI, &
                           MPI_SUM, mpi_comm_local, ierr)
     endif
#else
     ene_st(:,irep) = ene_st_l(:,0)
     ene_tst(:,irep) = ene_tst_l(:,0)
     if (inmisc%i_dtrna_model == 2013 .or. inmisc%i_dtrna_model == 2018) then
        ene_hb(:,irep) = ene_hb_l(:,0)
     endif
#endif
     deallocate( ene_st_l,stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)
     deallocate( ene_tst_l,stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)
     if (inmisc%i_dtrna_model == 2013 .or. inmisc%i_dtrna_model == 2018) then
        deallocate( ene_hb_l,stat=ier)
        if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)
     endif
  endif

  TIME_E( tmc_energy )
 
!  if(inmgo%i_multi_mgo >= 1) then
!     TIME_S( tm_energy_mgo)
!     call energy_mgo(energy_unit, energy)
!     TIME_E( tm_energy_mgo)
!  end if

!  if(inmisc%i_in_box == 1) then
!     call energy_box(irep, energy_unit, energy)
!  end if

!  if(inmisc%i_in_cap == 1) then
!     call energy_cap(irep, energy_unit, energy)
!  end if

!  if(inmisc%i_cylinder == 1) then
!     call energy_cylinder(irep, energy_unit, energy)
!  end if

  if(inmisc%i_bridge == 1) then
     call energy_bridge(irep, energy_unit, energy)
  end if

  if(inmisc%i_pulling == 1) then
     call energy_pulling(irep, energy_unit, energy)
  end if

  if(inmisc%i_anchor == 1) then
     call energy_anchor(irep, energy_unit, energy)
  end if

  if(inmisc%i_rest1d == 1) then
     call energy_rest1d(irep, energy_unit, energy)
  end if

  if(inmisc%i_window == 1) then
     call energy_window(irep, energy_unit, energy)
  end if
  
!   if(inmisc%force_flag(INTERACT%LJ1210)) then
!      call energy_LJ_1210(irep, energy_unit, energy)
!   end if

!  if(inmisc%i_winz == 1) then
!     call energy_winz(irep, energy_unit, energy)
!  end if

!  ! implicit ligand model
!  if(inmisc%i_implig == 1) then
!     call energy_implig(irep, energy_unit, energy, IMPLIGENERGY_TYPE%FOR_NON_MC)
!     ! calculate implicit_ligand binding energy based on [state of implicit-ligand &  structure].
!  end if
  
!  energy_unit(1:nunit_all, 1:nunit_all, 1:E_TYPE%MAX) = energy_unit_l(1:nunit_all, 1:nunit_all, 1:E_TYPE%MAX, 1)
!  energy(1:E_TYPE%MAX) = energy_l(1:E_TYPE%MAX, 1)

  ! --------------------------------------------------------------------
  ! calc energy of each interaction unit
  ! --------------------------------------------------------------------

  TIME_S( tm_energy_unit)

  if(inmisc%i_output_energy_style == 0) then
     do iunit = 1, nunit_all
        do junit = 1, iunit - 1
           energy_unit(iunit, iunit, 1:E_TYPE%MAX)                     &
                =  energy_unit(iunit, iunit, 1:E_TYPE%MAX)             &
                + 0.5e0_PREC * energy_unit(junit, iunit, 1:E_TYPE%MAX)
        end do
        do junit = iunit + 1, nunit_all
           energy_unit(iunit, iunit, 1:E_TYPE%MAX)                     &
                =  energy_unit(iunit, iunit, 1:E_TYPE%MAX)             &
                + 0.5e0_PREC * energy_unit(iunit, junit, 1:E_TYPE%MAX)
        end do
     end do
  else
     do iunit = 1, nunit_all
        do junit = iunit + 1, nunit_all
           do i = 1, E_TYPE%MAX
              sume = energy_unit(iunit, junit, i) &
                   + energy_unit(junit, iunit, i)
              energy_unit(iunit, junit, i) = sume
              energy_unit(junit, iunit, i) = sume
           end do
        end do
     end do

  end if
  do i = 1, E_TYPE%MAX
     if(i == E_TYPE%TOTAL .or. i == E_TYPE%VELO) cycle
     energy_unit(1:nunit_all, 1:nunit_all, E_TYPE%TOTAL)                 &
          =  energy_unit(1:nunit_all, 1:nunit_all, E_TYPE%TOTAL)  &
          + energy_unit(1:nunit_all, 1:nunit_all, i)
  end do
  TIME_E( tm_energy_unit)

  ! --------------------------------------------------------------------
  ! calc total energy
  ! --------------------------------------------------------------------
  TIME_S( tm_energy_total)
!  if(inmgo%i_multi_mgo >= 1) then
!     do i = 1, E_TYPE%MAX
!        if(i == E_TYPE%TOTAL  .or. i == E_TYPE%VELO    .or. i == E_TYPE%BOND .or. &
!           i == E_TYPE%BANGLE .or. i == E_TYPE%DIHE    .or. i == E_TYPE%GO   .or. &
!           i == E_TYPE%MORSE  .or. i == E_TYPE%DIHE_HARMONIC) cycle
!        energy(E_TYPE%TOTAL) = energy(E_TYPE%TOTAL) + energy(i)
!     end do
!  else
     do i = 1, E_TYPE%MAX
        if(i == E_TYPE%TOTAL .or. i == E_TYPE%VELO) cycle
        energy(E_TYPE%TOTAL) = energy(E_TYPE%TOTAL) + energy(i)
     end do
!  end if
  TIME_E( tm_energy_total)

#ifdef MEM_ALLOC
  deallocate( now_con )
  deallocate( now_LJ )
  deallocate( now_HPS )
!  deallocate( now_LJ1210 )
!  deallocate( now_morse )
!  deallocate( now_rna_bp )
  deallocate( now_allcon )
#endif

#ifdef _DEBUG
  write(6,*) '######## end energy_sumup'
#endif

end subroutine energy_sumup
