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
                       e_exv,         &
                       e_exv_unit)

  use if_energy
  use const_maxsize
  use const_index
  use var_setp,    only : inmisc, inflp, inele
  use var_struct,  only : nunit_all, ncon, nmorse, nrna_bp
  use var_mgo,     only : inmgo

  use time
  use mpiconst

  implicit none

! --------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(in)  :: velo_mp(:,:)      ! (SPACE_DIM, nmp_real)
  real(PREC), intent(out) :: e_exv(:)          ! (E_TYPE%MAX)
  real(PREC), intent(out) :: e_exv_unit(:,:,:)  ! (nunit_all, nunit_all, E_TYPE%MAX)

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
  real(PREC), allocatable :: e_exv_l(:,:)         !(E_TYPE%MAX, 0:nthreads-1)
  real(PREC), allocatable :: e_exv_unit_l(:,:,:,:) !(nunit_all,nunit_all,E_TYPE%MAX,0:nthreads-1)
  integer, allocatable :: now_allcon_l(:, :)      !(2, ncon+nmorse+nrna_bp)


! ------------------------------------------------------------------------
! zero clear
  e_exv(:)         = 0.0e0_PREC
  e_exv_unit(:,:,:) = 0.0e0_PREC

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
  allocate( e_exv_l(E_TYPE%MAX, 0:nthreads-1), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  allocate( e_exv_unit_l(nunit_all, nunit_all, E_TYPE%MAX, 0:nthreads-1), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)
  allocate( now_allcon_l(2, ncon+nmorse+nrna_bp), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)

!$omp parallel private(tn)
  tn = 0
!$  tn = omp_get_thread_num()

  e_exv_l(:,tn) = 0.0e0_PREC
  e_exv_unit_l(:,:,:,tn) = 0.0e0_PREC

  TIME_S( tm_energy_velo) 
  call simu_energy_velo(velo_mp, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  TIME_E( tm_energy_velo) 

  TIME_S( tm_energy_bond) 
  call simu_energy_bond  (irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  TIME_E( tm_energy_bond) 

  TIME_S( tm_energy_bangle)
  call simu_energy_bangle(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  TIME_E( tm_energy_bangle)

  if (inmisc%force_flag_local(LINTERACT%L_AICG2) .or. &
      inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
     call simu_energy_aicg13_gauss(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  end if

  !if (inmisc%i_add_int == 1) then
  if (inflp%i_flp == 1 .or. inmisc%force_flag_local(LINTERACT%L_FLP)) then
     call simu_energy_fbangle(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  endif

  TIME_S( tm_energy_dih)
  call simu_energy_dih   (irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  TIME_E( tm_energy_dih)

  if (inmisc%force_flag_local(LINTERACT%L_AICG2)) then
     call simu_energy_aicg14_gauss(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  else if (inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
     call simu_energy_dih_gauss(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  end if

  !if (inmisc%i_add_int == 1) then
  if (inflp%i_flp == 1 .or. inmisc%force_flag_local(LINTERACT%L_FLP)) then
     call simu_energy_fdih(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  endif

  TIME_S( tm_energy_dih_harmonic) 
  call simu_energy_dih_harmonic   (irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  TIME_E( tm_energy_dih_harmonic) 

  call simu_energy_rna_stack(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))

  if (inmisc%i_dtrna_model == 2013) then
     call simu_energy_dtrna_stack(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
     call simu_energy_dtrna_hbond13(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  else if (inmisc%i_dtrna_model == 2015) then
     call simu_energy_dtrna_stack_nlocal(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
     call simu_energy_dtrna_stack(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
     call simu_energy_dtrna_hbond15(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
     call simu_energy_exv_dt15(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  endif

  if(inmgo%i_multi_mgo >= 1) then
     TIME_S( tm_energy_nlocal_mgo) 
     call simu_energy_nlocal_mgo(irep, now_con, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
     call simu_energy_nlocal_rna_bp(irep, now_rna_bp, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
     TIME_E( tm_energy_nlocal_mgo) 
  else if (inmisc%force_flag(INTERACT%ENM)) then
     TIME_S( tm_energy_enm) 
     call simu_energy_enm(irep, now_con, e_exv_l(:,tn), e_exv_unit_l(:,:,:,tn))
     TIME_E( tm_energy_enm) 
  else
     TIME_S( tm_energy_nlocal_go) 
     call simu_energy_nlocal_go(irep, now_con, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
     call simu_energy_nlocal_morse(irep, now_morse, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
     call simu_energy_nlocal_rna_bp(irep, now_rna_bp, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
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

  TIME_S( tm_energy_exv) 
  if (inmisc%i_residue_exv_radii == 0) then
     call simu_energy_exv_rep12 (irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  else if (inmisc%i_residue_exv_radii == 1) then
     call simu_energy_exv_restype(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  endif
  if (inmisc%force_flag(INTERACT%EXV_WCA)) then 
     call simu_energy_exv_wca (irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
  endif
  TIME_E( tm_energy_exv) 

  TIME_S( tm_energy_ele)
  if (inmisc%force_flag(INTERACT%ELE)) then

     if (inele%i_function_form == 0) then     ! Debye-Huckel (default)

        if(inele%i_calc_method == 0) then         ! neighboring list (default)
           call simu_energy_ele(irep, e_exv_l(:,tn), e_exv_unit_l(:,:,:,tn))
        else if(inele%i_calc_method == 1) then    ! neighboring list for K computer
           call simu_energy_ele2(irep, e_exv_l(:,tn), e_exv_unit_l(:,:,:,tn))
        else                                      ! direct calculation for K computer
           call simu_energy_ele3(irep, e_exv_l(:,tn), e_exv_unit_l(:,:,:,tn))
        end if

     elseif (inele%i_function_form == 1) then ! Coulomb potential
        call simu_energy_ele_coulomb(irep, e_exv_l(:,tn), e_exv_unit_l(:,:,:,tn))
     endif
  endif
  TIME_E( tm_energy_ele)

  if (inmisc%class_flag(CLASS%ION)) then
     TIME_S( tm_energy_exv) 
     call simu_energy_ion(irep, e_exv_unit_l(:,:,:,tn), e_exv_l(:,tn))
     TIME_E( tm_energy_exv) 
  endif

  if (inmisc%force_flag(INTERACT%HP)) then
     TIME_S( tm_energy_hp) 
     call simu_energy_hp(irep, e_exv_l(:,tn), e_exv_unit_l(:,:,:,tn))
     TIME_E( tm_energy_hp) 
  endif

!sasa
  if (inmisc%force_flag(INTERACT%SASA)) then
     TIME_S( tm_energy_sasa)
     call simu_energy_sasa(irep, e_exv_l(:,tn))
     TIME_E( tm_energy_sasa)
  endif

!$omp end parallel

  TIME_S( tmc_energy )
  do n = 1, nthreads-1
    e_exv_unit_l(1:nunit_all,1:nunit_all,1:E_TYPE%MAX,0) = &
    e_exv_unit_l(1:nunit_all,1:nunit_all,1:E_TYPE%MAX,0) + &
    e_exv_unit_l(1:nunit_all,1:nunit_all,1:E_TYPE%MAX,n)
    e_exv_l(1:E_TYPE%MAX,0) = &
    e_exv_l(1:E_TYPE%MAX,0) + &
    e_exv_l(1:E_TYPE%MAX,n)
  end do

  TIME_E( tmc_energy )

  TIME_S( tmc_energy )
#ifdef MPI_PAR3
  call mpi_allreduce(e_exv_unit_l, e_exv_unit, &
                     nunit_all*nunit_all*E_TYPE%MAX, PREC_MPI, &
                     MPI_SUM, mpi_comm_local, ierr)
  call mpi_allreduce(e_exv_l, e_exv, &
                     E_TYPE%MAX, PREC_MPI, &
                     MPI_SUM, mpi_comm_local, ierr)
#else
  e_exv_unit(:,:,:) = e_exv_unit_l(:,:,:,0)
  e_exv(:) = e_exv_l(:,0)
#endif
  TIME_E( tmc_energy )
 
  if(inmgo%i_multi_mgo >= 1) then
     TIME_S( tm_energy_mgo)
     call simu_energy_mgo(e_exv_unit, e_exv)
     TIME_E( tm_energy_mgo)
  end if

  deallocate( e_exv_l,     stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)
  deallocate( e_exv_unit_l, stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)
  deallocate( now_allcon_l,stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, msg_er_deallocate)


  if(inmisc%i_in_box == 1) then
     call simu_energy_box(irep, e_exv_unit, e_exv)
  end if

  if(inmisc%i_in_cap == 1) then
     call simu_energy_cap(irep, e_exv_unit, e_exv)
  end if

  if(inmisc%i_cylinder == 1) then
     call simu_energy_cylinder(irep, e_exv_unit, e_exv)
  end if

  if(inmisc%i_bridge == 1) then
     call simu_energy_bridge(irep, e_exv_unit, e_exv)
  end if

  if(inmisc%i_pulling == 1) then
     call simu_energy_pulling(irep, e_exv_unit, e_exv)
  end if

  if(inmisc%i_anchor == 1) then
     call simu_energy_anchor(irep, e_exv_unit, e_exv)
  end if

  if(inmisc%i_rest1d == 1) then
     call simu_energy_rest1d(irep, e_exv_unit, e_exv)
  end if

  if(inmisc%i_window == 1) then
     call simu_energy_window(irep, e_exv_unit, e_exv)
  end if
  
  if(inmisc%i_winz == 1) then
     call simu_energy_winz(irep, e_exv_unit, e_exv)
  end if

  ! implicit ligand model
  if(inmisc%i_implig == 1) then
     call simu_energy_implig(irep, e_exv_unit, e_exv, IMPLIGENERGY_TYPE%FOR_NON_MC)
     ! calculate implicit_ligand binding energy based on [state of implicit-ligand &  structure].
  end if
  
!  e_exv_unit(1:nunit_all, 1:nunit_all, 1:E_TYPE%MAX) = e_exv_unit_l(1:nunit_all, 1:nunit_all, 1:E_TYPE%MAX, 1)
!  e_exv(1:E_TYPE%MAX) = e_exv_l(1:E_TYPE%MAX, 1)

  ! --------------------------------------------------------------------
  ! calc energy of each interaction unit
  ! --------------------------------------------------------------------

  TIME_S( tm_energy_unit)

  if(inmisc%i_output_energy_style == 0) then
     do iunit = 1, nunit_all
        do junit = 1, iunit - 1
           e_exv_unit(iunit, iunit, 1:E_TYPE%MAX)                     &
                =  e_exv_unit(iunit, iunit, 1:E_TYPE%MAX)             &
                + 0.5e0_PREC * e_exv_unit(junit, iunit, 1:E_TYPE%MAX)
        end do
        do junit = iunit + 1, nunit_all
           e_exv_unit(iunit, iunit, 1:E_TYPE%MAX)                     &
                =  e_exv_unit(iunit, iunit, 1:E_TYPE%MAX)             &
                + 0.5e0_PREC * e_exv_unit(iunit, junit, 1:E_TYPE%MAX)
        end do
     end do
  else
     do iunit = 1, nunit_all
        do junit = iunit + 1, nunit_all
           do i = 1, E_TYPE%MAX
              sume = e_exv_unit(iunit, junit, i) &
                   + e_exv_unit(junit, iunit, i)
              e_exv_unit(iunit, junit, i) = sume
              e_exv_unit(junit, iunit, i) = sume
           end do
        end do
     end do

  end if
  do i = 1, E_TYPE%MAX
     if(i == E_TYPE%TOTAL .or. i == E_TYPE%VELO) cycle
     e_exv_unit(1:nunit_all, 1:nunit_all, E_TYPE%TOTAL)                 &
          =  e_exv_unit(1:nunit_all, 1:nunit_all, E_TYPE%TOTAL)  &
          + e_exv_unit(1:nunit_all, 1:nunit_all, i)
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
        e_exv(E_TYPE%TOTAL) = e_exv(E_TYPE%TOTAL) + e_exv(i)
     end do
  else
     do i = 1, E_TYPE%MAX
        if(i == E_TYPE%TOTAL .or. i == E_TYPE%VELO) cycle
        e_exv(E_TYPE%TOTAL) = e_exv(E_TYPE%TOTAL) + e_exv(i)
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
