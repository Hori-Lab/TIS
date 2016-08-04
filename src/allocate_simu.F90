subroutine allocate_simu()

  use const_maxsize
  use const_physical
  use const_index
  use var_io,     only : i_simulate_type
  use var_struct,  only : nmp_all, nmp_real, nunit_all
  use var_replica, only : n_replica_all, n_replica_mpi
  use var_simu,    only : dxyz_mp, velo_mp, accel_mp, force_mp, rcmass_mp, &
                          energy_unit_muca, rlan_const, nLAN_CONST, &
                          energy, energy_unit, qscore, qscore_unit, &
                          rg, rg_unit, rmsd, rmsd_unit, &
#ifdef MPI_PAR
                          replica_energy_l, &
#endif
                          replica_energy, &
                          sasa, &
                          diffuse_tensor, random_tensor

  implicit none
  integer :: ier = 0
  character(CARRAY_MSG_ERROR) :: error_message
! **********************************************************************
#ifdef _DEBUG
  write(*,*) '##### start allocate_simu'
#endif
  
  call check()

  error_message = 'failed in memory allocation at mloop_simulator, PROGRAM STOP'
  ! dxyz_mp
  allocate( dxyz_mp(SDIM, nmp_real, n_replica_mpi), stat=ier)
  if(ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

  ! velo_mp
  allocate( velo_mp(SDIM, nmp_real, n_replica_mpi), stat=ier)
  if(ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

  ! accel_mp
  allocate( accel_mp(SDIM, nmp_real, n_replica_mpi), stat=ier)
  if(ier /= 0) call util_error(ERROR%STOP_ALL, error_message)

  ! force_mp
  allocate( force_mp(SDIM, nmp_all), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
  
  ! rcmass_mp
  allocate( rcmass_mp(nmp_all), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
  
  ! energy_unit_muca
  allocate( energy_unit_muca(nunit_all, nunit_all, E_TYPE%MAX), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
  
  ! rlan_const 
  allocate( rlan_const(nLAN_CONST, nmp_real, n_replica_mpi), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
  
  ! energy
  allocate( energy(E_TYPE%MAX, n_replica_mpi), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  
  ! energy_unit
  allocate( energy_unit(nunit_all, nunit_all, E_TYPE%MAX, n_replica_mpi), stat=ier) 
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  
  ! qscore
  allocate( qscore(n_replica_mpi), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  
  ! qscore_unit
  allocate( qscore_unit(nunit_all, nunit_all, n_replica_mpi), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  
  ! rg
  allocate( rg(n_replica_mpi), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  
  ! rg_unit
  allocate( rg_unit(nunit_all, n_replica_mpi), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  
  ! rmsd
  allocate( rmsd(n_replica_mpi), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  
  ! rmsd_unit
  allocate( rmsd_unit(nunit_all, n_replica_mpi), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  
  ! replica_energy
  allocate( replica_energy(2, n_replica_all), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  
#ifdef MPI_PAR
  allocate(replica_energy_l(2, n_replica_all), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL, error_message)
#endif

  !sasa
  allocate( sasa(nmp_all), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)

  ! Tensor used for Brownian dynamics with hydrodynamic interaction
  if(i_simulate_type == SIM%BROWNIAN_HI) then
     allocate(diffuse_tensor(1:3*nmp_real, 1:3*nmp_real), stat=ier)
     if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)

     allocate(random_tensor(1:3*nmp_real, 1:3*nmp_real), stat=ier)
     if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  endif

#ifdef _DEBUG
  write(*,*) '##### end allocate_simu'
#endif

! **********************************************************************
contains


! **********************************************************************
  subroutine check()

    if ( allocated(dxyz_mp) .OR.& 
         allocated(velo_mp) .OR.&
         allocated(accel_mp) .OR.&
         allocated(force_mp) .OR.&
         allocated(rcmass_mp) .OR.&
         allocated(energy_unit_muca) .OR.&
         allocated(rlan_const) .OR.&
         allocated(energy) .OR.&
         allocated(energy_unit) .OR.&
         allocated(qscore) .OR.&
         allocated(qscore_unit) .OR.&
         allocated(rg) .OR.&
         allocated(rg_unit) .OR.&
         allocated(rmsd) .OR.&
         allocated(rmsd_unit) .OR.&
         allocated(replica_energy) .OR.&
#ifdef MPI_PAR
         allocated(replica_energy_l) .OR.&
#endif
         allocated(sasa) .OR.&
         allocated(diffuse_tensor) .OR. &
         allocated(random_tensor) .OR. &
         .FALSE.) then
       error_message = 'defect at allocate_simu::check, PROGRAM STOP'
       call util_error(ERROR%STOP_ALL,error_message)
    endif
  end subroutine check

end subroutine allocate_simu


! **********************************************************************
subroutine deallocate_simu()

  use var_simu
  implicit none

  if (allocated(dxyz_mp)) deallocate(dxyz_mp)
  if (allocated(velo_mp)) deallocate(velo_mp)
  if (allocated(accel_mp)) deallocate(accel_mp)
  if (allocated(force_mp)) deallocate(force_mp)
  if (allocated(rcmass_mp)) deallocate(rcmass_mp)
  if (allocated(energy_unit_muca)) deallocate(energy_unit_muca)
  if (allocated(rlan_const)) deallocate(rlan_const)
  if (allocated(energy)) deallocate(energy)
  if (allocated(energy_unit)) deallocate(energy_unit)
  if (allocated(qscore)) deallocate(qscore)
  if (allocated(qscore_unit)) deallocate(qscore_unit)
  if (allocated(rg)) deallocate(rg)
  if (allocated(rg_unit)) deallocate(rg_unit)
  if (allocated(rmsd)) deallocate(rmsd)
  if (allocated(rmsd_unit)) deallocate(rmsd_unit)
  if (allocated(replica_energy)) deallocate(replica_energy)
#ifdef MPI_PAR
  if (allocated(replica_energy_l)) deallocate(replica_energy_l)
#endif
  if (allocated(sasa)) deallocate(sasa)
  if (allocated(diffuse_tensor)) deallocate(diffuse_tensor)
  if (allocated(random_tensor)) deallocate(random_tensor)

end subroutine deallocate_simu
