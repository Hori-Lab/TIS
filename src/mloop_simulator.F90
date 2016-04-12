!mloop_simulator
!> @brief Main time-propagation loop. This subroutine is called   &
!>        by "main_loop".

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! **********************************************************************
! i_run_mode: define bacis run mode 
!           = 1 : Debug Mode, Check the consistence between force and energy
!           = 2 : Constant temperature simulation 
!           = 3 : Simulated annealing (require "<<<< annealing" field)
!           = 4 : Auto-search of T_f (require "<<<< searching_tf" field)
!           = 5 : Energy calculation at single point
!           = 6 : Replica exchange method
!           = 7 : Fluctuation matching method
!
! i_simulate_type: define dynamics 
!           = 0 Newtonian dynamics (velocity Verlet) with the constant energy
!           = 1 Langevin dynamics (recommended)
!           = 2 Newtonian dynamics (velocity Verlet) with Berendsen thermostat
! **********************************************************************
subroutine mloop_simulator()
  
  use const_maxsize
  use const_physical
  use const_index
  use if_mloop
  use if_write
  use if_energy
  use var_inp,     only : i_run_mode, outfile
  use var_setp,    only : inele, insimu
  use var_replica, only : n_replica_mpi
  use var_simu,    only : imstep, mstep, istep, ntstep, ibefore_time, &
                          tempk, qscore, velo_mp
  use time, only : tm_main_loop, time_s, time_e, time_write, time_initialize, &
                   tm_tinte, tm_tinte_post
  use mpiconst

  implicit none

  ! -----------------------------------------------------------------
  ! local variables
  logical :: flg_step_each_replica(n_replica_mpi)
  logical :: flg_exit_loop_mstep


  ! -----------------------------------------------------------------
#ifdef _DEBUG
  write(6,*) 'mloop_simulator: START'
#endif
#ifdef MPI_PAR
  call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif

  call simu_logicalcheck()

  ! -------------------------------------------------------------------
  ! check the force and energy: for debug
  ! -------------------------------------------------------------------
  if(i_run_mode == RUN%CHECK_FORCE) then
     call simu_checkforce()
  end if

  ! -----------------------------------------------------------------
  ! setting parameters and reading parameters from input file
  ! -----------------------------------------------------------------
  call simu_para()

  ! -----------------------------------------------------------------
  ! initial setting for time integral
  ! -----------------------------------------------------------------
#ifdef _DEBUG
  write(6,*) 'mloop_simulator: initial setting for time integral'
#endif

  loop_mstep: do imstep = 1, mstep  ! only for searching Tf

     ! -----------------------------------------------------------------
     ! initialize structure and velocity
     ! -----------------------------------------------------------------
     call simu_initial()

     ! -----------------------------------------------------------------
     ! set parameters depending on temperature or ion strength
     ! -----------------------------------------------------------------
     call simu_para2(tempk, inele%ionic_strength)

     ! -----------------------------------------------------------------
     ! Case of energy calculation for DCD trajectory
     ! -----------------------------------------------------------------
     if (i_run_mode == RUN%ENERGY_DCD) then
        velo_mp(:,:,:) = 0.0e0_PREC
        !call simu_calc_energy_dcd(0)
        do istep = insimu%i_tstep_init-1, ntstep, insimu%n_step_save
           call simu_calc_energy_dcd(istep)
        enddo
        exit loop_mstep
     endif

     ! -----------------------------------------------------------------
     ! preparing time integral
     ! -----------------------------------------------------------------
     call simu_tintegral_pre(flg_step_each_replica)

     ! #################################################################
     !  time integral
     ! #################################################################
#ifdef _DEBUG
     write(6,*) 'mloop_simulator: time integral'
#endif 
#ifdef TIME
     call time_initialize()
#endif
     TIME_S( tm_main_loop )

     do istep = insimu%i_tstep_init, ntstep
        flg_exit_loop_mstep = .FALSE.

        TIME_S( tm_tinte )
        if (i_run_mode == RUN%EMIN) then
           ! -----------------------------------------------------------------
           ! 1-step of energy minimization
           ! -----------------------------------------------------------------
           call simu_minimize(flg_exit_loop_mstep)
        else
           ! -----------------------------------------------------------------
           ! time integral
           ! -----------------------------------------------------------------
           call simu_tintegral(flg_step_each_replica)
        endif
        TIME_E( tm_tinte )

        ! -----------------------------------------------------------------
        ! treatment after time integral
        ! -----------------------------------------------------------------
        TIME_S( tm_tinte_post )
        call simu_tintegral_post(flg_step_each_replica, flg_exit_loop_mstep)
        TIME_E( tm_tinte_post )
        if(flg_exit_loop_mstep) then
           TIME_E( tm_main_loop )
           exit loop_mstep
        endif

     end do ! for ntstep
     TIME_E( tm_main_loop )

     if(i_run_mode == RUN%SEARCH_TF) then
        call simu_searchingtf(istep, ntstep, qscore(1), tempk)
     end if

  end do loop_mstep

  ibefore_time = ibefore_time + ntstep

#ifdef TIME
  call time_write(6)
  if (myrank == 0) then
     call time_write(outfile%data)
  endif
#endif

#ifdef _DEBUG
  write(6,*) 'mloop_simulator: END'
#endif

end subroutine mloop_simulator
