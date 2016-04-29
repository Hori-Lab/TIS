! simu_checkforce
!> @brief This subroutine is to calculate the force from the gradient of the energy &
!>        and compare with the value of force which is calculated by ``force_sumup" subroutine.&
!>        This subroutine is especially for debug (checking).

! ***********************************************************************
subroutine simu_checkforce()
      
  use if_mloop
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : outfile
  use var_setp,    only : inele, insimu
  use var_struct,  only : nmp_real, xyz_mp_rep, pxyz_mp_rep, nunit_all
  use var_mgo,     only : inmgo
  use var_replica, only : n_replica_all, n_replica_mpi
  use var_simu,    only : flg_hb_energy
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! -----------------------------------------------------------------
  ! function
  real(PREC) :: rfunc_boxmuller

  ! --------------------------------------------------------------------
  ! local varables
  integer :: irep
  integer :: lunout
  integer :: i, imp, idimn
  integer :: ier = 0
  real(PREC) :: af_ene, xyz_save, tempk, ddrand
  real(PREC) :: force_mp(SPACE_DIM, MXMP), ff(SPACE_DIM), zure(3)
  character(CARRAY_MSG_ERROR) :: error_message
  ! -----------------------------------------------------------------
  real(PREC), allocatable :: velo_mp(:,:,:)     ! (SPACE_DIM, MXMP, REPLICA=1)
  real(PREC), allocatable :: energy(:,:)         ! (E_TYPE%MAX, REPLICA=1)
  real(PREC), allocatable :: energy_unit(:,:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX, REPLICA=1)
  real(PREC), allocatable :: replica_energy(:,:)! (2, REPLICA=1)

  ! -----------------------------------------------------------------
  ! constant parameter
  integer,    parameter :: IDX_REPLICA = 1
  real(PREC), parameter :: small = 0.0001e0_PREC
  ! TO DEVELOPER : you must select appropriate value about small.

  ! --------------------------------------------------------------------
  ! Allocation
  error_message = 'failed in memory allocation at mloop_simulator, PROGRAM STOP'
  ! velo_mp
  allocate( velo_mp(SPACE_DIM, nmp_real, n_replica_all), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL,error_message)
  ! energy
  allocate( energy(E_TYPE%MAX, n_replica_all), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL,error_message)
  ! energy_unit
  allocate( energy_unit(nunit_all, nunit_all, E_TYPE%MAX, n_replica_all), stat=ier) 
  if (ier/=0) call util_error(ERROR%STOP_ALL,error_message)
  ! replica_energy
  allocate( replica_energy(2, n_replica_all), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  

  ! --------------------------------------------------------------------
  lunout = outfile%data
  tempk  = insimu%tempk
  velo_mp(:,:,:) = 0.0e0_PREC

  ! -----------------------------------------------------------------
  ! initial structure
  ! -----------------------------------------------------------------
  call simu_initial()
  
  ! set inele%coef, inele%cdist
  call simu_para2(tempk, inele%ionic_strength)

  do imp = 1, nmp_real
     do idimn = 1, 3
        ddrand = small * rfunc_boxmuller(IDX_REPLICA-1, 0)
        xyz_mp_rep(idimn, imp, IDX_REPLICA) = xyz_mp_rep(idimn, imp, IDX_REPLICA) &
                                            + ddrand
        pxyz_mp_rep(idimn, imp, IDX_REPLICA) = pxyz_mp_rep(idimn, imp, IDX_REPLICA) &
                                            + ddrand
    end do
  end do

  do irep = 1, n_replica_mpi
     call simu_copyxyz(irep)
     call simu_neighbor(irep)
  enddo
  
  call energy_allrep(energy_unit, energy, &
                   velo_mp, replica_energy, .false., tempk)

  call force_sumup(force_mp, IDX_REPLICA)
  flg_hb_energy = .False.

  do imp = 1, nmp_real
     do idimn = 1, 3
        xyz_save = xyz_mp_rep(idimn, imp, IDX_REPLICA)
        xyz_mp_rep(idimn, imp, IDX_REPLICA) = xyz_mp_rep(idimn, imp, IDX_REPLICA) + small
        pxyz_mp_rep(idimn, imp, IDX_REPLICA) = pxyz_mp_rep(idimn, imp, IDX_REPLICA) + small
        call simu_copyxyz(IDX_REPLICA)
   
        call energy_allrep(energy_unit, energy, &
                         velo_mp, replica_energy, .false., tempk)
            
        af_ene = energy(E_TYPE%TOTAL, IDX_REPLICA)
        xyz_mp_rep(idimn, imp, IDX_REPLICA) = xyz_mp_rep(idimn, imp, IDX_REPLICA) - 2.0_PREC * small
        pxyz_mp_rep(idimn, imp, IDX_REPLICA) = pxyz_mp_rep(idimn, imp, IDX_REPLICA) - 2.0_PREC * small
        call simu_copyxyz(IDX_REPLICA)
   
        call energy_allrep(energy_unit, energy, &
                         velo_mp, replica_energy, .false., tempk)
   
        ff(idimn) = -(af_ene - energy(E_TYPE%TOTAL, IDX_REPLICA)) * 0.5_PREC / small
        xyz_mp_rep(idimn, imp, IDX_REPLICA) = xyz_save
        pxyz_mp_rep(idimn, imp, IDX_REPLICA) = xyz_save
     end do
 
#ifdef MPI_PAR
     if (myrank == 0) then
#endif
        write (lunout, "('imp=', i4, 2x, 'f_finite', 3f10.4, 2x, 'f_analysis=', 3f10.4)") &
               imp, (ff(i), i = 1, 3), (force_mp(i, imp), i = 1, 3)

        zure(1:3) = ff(1:3) - force_mp(1:3, imp)
        
        write (lunout, "('f_finite - f_analysis', i4, 3f10.4)") imp, (zure(i), i = 1, 3)

        if(abs(zure(1)) > small .or. abs(zure(2)) > small .or. abs(zure(3)) > small) then
           write (*, "('f_finite - f_analysis', i4, 3f10.4)") &
                  imp, (zure(i), i = 1, 3)
        end if
#ifdef MPI_PAR
     end if
#endif
  end do
  

  ! Deallocation
  if (allocated(velo_mp)    ) deallocate(velo_mp,     stat=ier)
  if (ier/=0                ) call util_error(ERROR%STOP_ALL,error_message)
  if (allocated(energy)      ) deallocate(energy,       stat=ier)
  if (ier/=0                ) call util_error(ERROR%STOP_ALL,error_message)
  if (allocated(energy_unit)  ) deallocate(energy_unit,   stat=ier)
  if (ier/=0                ) call util_error(ERROR%STOP_ALL,error_message)
  if (allocated(replica_energy)) deallocate(replica_energy, stat=ier)
  if (ier/=0                ) call util_error(ERROR%STOP_ALL,error_message)

  error_message = 'PROGRAM STOP'
  call util_error(ERROR%STOP_ALL, error_message)

end subroutine simu_checkforce
