subroutine simu_calc_energy_dcd(istep_write)
  use const_maxsize
  use const_physical
  use const_index
  use if_mloop
  use if_write
  use var_inp,     only : ifile_out_opt
  use var_replica, only : n_replica_mpi
  use var_simu,    only : ibefore_time, tempk, velo_mp, e_exv, e_exv_unit, &
                          rg, rg_unit, rmsd, rmsd_unit, replica_energy
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  integer(L_INT), intent(in) :: istep_write

  integer :: irep
  logical, save :: flg_header = .true.


  ! -----------------------------------------------------------------
  ! Read a snapshot from DCD
  ! -----------------------------------------------------------------
  call read_dcd()

  ! -----------------------------------------------------------------
  ! calc neighbour list
  ! -----------------------------------------------------------------
  do irep=1, n_replica_mpi
     call simu_neighbor(irep)
  enddo
  
  ! calc energy
  ! ## Here, replica_energy is dummy.
  call energy_allrep(e_exv_unit, e_exv, &
       velo_mp, replica_energy, .false., tempk)
  call simu_radiusg(rg_unit, rg)
  call simu_rmsd(rmsd_unit, rmsd)
  
#ifdef MPI_PAR
  if (local_rank_mpi == 0) then
#endif

  call write_tseries(ibefore_time, istep_write, &
                     rg_unit, rg, rmsd_unit, rmsd, &
                     e_exv_unit, e_exv, tempk, flg_header)
  if (ifile_out_opt == 1) then
     ! something to write to opt file

     !write(outfile%opt,'(a1,1x,i10,1x,f15.6,1x,f15.6,1x,f15.6)') &
     !                   '#',istep_write, inmisc%rest1d_s(1), &
     !                   e_exv(E_TYPE%REST1D,irep), e_exv(E_TYPE%TOTAL,irep)
  endif
  flg_header = .false.
#ifdef MPI_PAR
  end if
#endif
  
!  flg_step_each_replica(1:n_replica_mpi)  = .false. 

endsubroutine simu_calc_energy_dcd