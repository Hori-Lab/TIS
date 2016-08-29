! var_simu
!> @brief Module defining the global variables used in mloop_simulator.F90

module var_simu

  use const_maxsize
  use const_physical
  use const_index
  implicit none

  ! time step
  integer, save        :: istep_sim, mstep_sim
  integer, save        :: imstep, mstep
  integer(L_INT), save :: istep, ntstep, nstep_opt_temp , ntstep_max
  integer(L_INT), save :: ibefore_time = 0

  !
  integer       :: n_exchange, max_exchange
  integer       :: iopt_stage

  ! Neighbor list
  real(PREC), allocatable,save  :: dxyz_mp(:,:,:)  ! (SDIM, nmp_real, n_replica_mpi)

  ! physical variables
  real(PREC), save :: tempk
  real(PREC), save :: tstep, tsteph, tstep2
  real(PREC), save :: accelaf(SDIM)
  real(PREC), allocatable,save  :: velo_mp(:,:,:)  ! (SDIM, nmp_real, n_replica_mpi)
  real(PREC), allocatable,save  :: accel_mp(:,:,:) ! (SDIM, nmp_real, n_replica_mpi)
  real(PREC), allocatable, save :: force_mp(:,:)   ! (SDIM, nmp_all)
  real(PREC), allocatable, save :: rcmass_mp(:)    ! (nmp_all)

!  ! mcanonical
!  real(PREC), save :: e_md, fac_mmc, em_mid, em_depth, em_sigma
!  real(PREC), save :: energy_muca(E_TYPE%MAX)
!  real(PREC), allocatable, save :: energy_unit_muca(:,:,:) ! (nunit_all, nunit_all, E_TYPE%MAX)
  
  ! Langevin
  integer, parameter          :: nLAN_CONST = 4
  real(PREC), save              :: r_force(SDIM)
  real(PREC), save              :: tstep_fric_h, ulconst1, ulconst2
  real(PREC), allocatable, save :: rlan_const(:,:,:) ! (nLAN_CONST, mp, replica)

  ! Brownian with hydrodynamic interaction
  real(PREC), allocatable, save :: diffuse_tensor(:,:)
  real(PREC), allocatable, save :: random_tensor(:,:)

  ! Nose-Hoover  
  ! Now, Replica is not available
  integer, save    :: ics, jcs, ncs
  real(PREC), save :: velo_yojou(MXCS), evcs(MXCS)
  real(PREC), save :: xyz_cs(MXCS), velo_cs(MXCS), cmass_cs(MXCS)

  ! energy
  real(PREC), allocatable, save :: energy(:,:)          ! (E_TYPE%MAX, replica)
  real(PREC), allocatable, save :: energy_unit(:,:,:,:)  ! (unit, unit, E_TYPE%MAX, replica)
  real(PREC), allocatable, save :: qscore(:)           ! (replica)
  real(PREC), allocatable, save :: qscore_unit(:,:,:)  ! (unit, unit, replica)
  real(PREC), allocatable, save :: rg(:)               ! (replica)
  real(PREC), allocatable, save :: rg_unit(:,:)        ! (unit, replica)
  real(PREC), allocatable, save :: rmsd(:)             ! (replica)
  real(PREC), allocatable, save :: rmsd_unit(:,:)      ! (unit, replica)
  real(PREC), allocatable, save :: replica_energy(:,:) ! (2, replica)
#ifdef MPI_PAR
  real(PREC), allocatable, save :: replica_energy_l(:,:) ! (2, replica)
#endif

  ! sasa
  real(PREC), allocatable, save :: sasa(:)          ! (E_TYPE%MAX, replica)

  ! DTRNA15
  integer, allocatable, save  :: hbsite_excess(:)     ! (1:nhbsite)  ! used only in force_dtrna_hbond15
  real(PREC), allocatable, save :: hb_energy(:,:)     ! (1:ndtrna_hb, REPLICA)
  logical, allocatable, save :: hb_status(:,:)        ! (1:ndtrna_hb, REPLICA)
  logical, save :: flg_hb_energy
  logical, allocatable, save :: st_status(:,:)        ! (1:ndtrna_st, REPLICA)


  logical, save :: flg_ppr_release(MXBRIDGE) = .false.

  ! Widom
  integer(L_INT), save :: widom_iw
  real(PREC), save :: widom_chp
  logical, save :: widom_flg_exv_inf

  ! Ewald
  integer, save :: ewld_f_n
  real(PREC), allocatable, save :: ewld_h(:,:)      ! Reciplocal lattice vector
  real(PREC), allocatable, save :: ewld_f_rlv(:,:)  ! Reciplocal lattice vector (*coef)
  real(PREC), allocatable, save :: ewld_f_coef(:)
  real(PREC), save :: ewld_s_coef    ! Self interaction
  real(PREC), save :: ewld_s_sum     ! Self interaction

end module var_simu
