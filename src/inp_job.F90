! inp_job
!> @brief This subroutine is to read and set the parameters for job-control.

! ********************************************************
! subroutine for reading parameters of job-control
! ********************************************************

! NOTICE: i_run_mode has been already read in setp_replica.

subroutine inp_job()

  use const_maxsize
  use const_index
  use var_io,   only : infile, outfile,  i_run_mode, i_simulate_type, &
                       i_initial_state, i_initial_velo, flg_rst
  use var_setp, only : inperi
  use mpiconst

  implicit none

  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  ! --------------------------------------------------------------------
  ! default 
!  i_run_mode = -1   ! i_run_mode is already read in inp_replica
  i_simulate_type = -1
  i_initial_state = INISTAT%VOID
  i_initial_velo  = INIVELO%MAXWELL
  inperi%i_periodic = 0


  ! --------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'job_cntl        ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)   
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "job_cntl" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        cvalue = 'i_simulate_type'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             i_simulate_type, cvalue)
        
        cvalue = 'i_initial_state'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             i_initial_state, cvalue)
           
        cvalue = 'i_initial_velo'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             i_initial_velo, cvalue)
        
        cvalue = 'i_periodic'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inperi%i_periodic, cvalue)
     end do
  end do
  

  ! For the case of restart
  if (flg_rst) then
     i_initial_state = INISTAT%RST
     i_initial_velo = INIVELO%RST
  endif
  
  if(i_simulate_type == -1) then
     error_message = 'Error: cannot find "i_simulate_type" in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(i_initial_state == INISTAT%VOID) then
     error_message = 'Error: cannot find "i_initial_state" in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  
  ! -----------------------------------------------------------------
  ! The explanation of i_run_mode
  write (lunout, *) 'i_run_mode = ', i_run_mode
  if(i_run_mode == RUN%CHECK_FORCE) then
     write (lunout, *) 'Debug mode'
     write (lunout, *) 'Check the consitence between force and energy'
  else if(i_run_mode == RUN%CONST_TEMP) then
     write (lunout, *) 'Constant temperature simulation'
  else if(i_run_mode == RUN%SA) then
     write (lunout, *) 'Simulated annealing (require "<<<< annealing" field)'
  else if(i_run_mode == RUN%ENERGY_CALC) then
     write (lunout, *) 'Energy calculation at single point'
  else if(i_run_mode == RUN%REPLICA) then
     write (lunout, *) 'Replica exchange method (require "<<<< replica" field)'
  else if(i_run_mode == RUN%FMAT) then              !fmat
     write (lunout, *) 'Fluctuation matching'    !fmat
  else if(i_run_mode == RUN%ENERGY_DCD) then
     write (lunout, *) 'Energy calculation for DCD trajectory'
  else if(i_run_mode == RUN%EMIN) then
     write (lunout, *) 'Energy minimization'
  else if(i_run_mode == RUN%WIDOM) then
     write (lunout, *) 'Widom method to calculate chemical potential(s)'
  else
     error_message = 'Error: invalid value about i_run_mode'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  

  ! -----------------------------------------------------------------
  ! The way of i_simulate_type 
  write (lunout, *) 'i_simulate_type = ', i_simulate_type
  if(i_simulate_type == SIM%CONST_ENERGY) then
     write (lunout, *) 'Newtonian dynamics (velocity Verlet) with the constant energy'
  else if(i_simulate_type == SIM%LANGEVIN) then
     write (lunout, *) 'Langevin dynamics (recommended)'
  else if(i_simulate_type == SIM%BERENDSEN) then
     write (lunout, *) 'Newtonian dynamics (velocity Verlet) with Berendsen thermostat'
  else if(i_simulate_type == SIM%NOSEHOOVER) then
     write (lunout, *) 'Newtonian dynamics (volocity Verlet) with Nose-Hoover thermostat'
  else if(i_simulate_type == SIM%MPC) then
     write (lunout, *) 'MPC dynamics'
  else if(i_simulate_type == SIM%BROWNIAN) then
     write (lunout, *) 'Brownian dynamics'
  else if(i_simulate_type == SIM%BROWNIAN_HI) then
     write (lunout, *) 'Brownian dynamics with hydrodynamic interaction'
  else if(i_simulate_type == SIM%PS_BROWNIAN) then
     write (lunout, *) 'Brownian dynamics (time scale: ps)'
  else if(i_simulate_type == SIM%ND_LANGEVIN) then
     write (lunout, *) 'Langevin dynamics implemented by ND'
  else if(i_simulate_type == SIM%LANGEVIN_GJF) then
     write (lunout, *) 'Langevin dynamics by GJF algorithm'
  else
     error_message = 'Error: invalid value about i_simulate_type'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  

  ! -----------------------------------------------------------------
  ! The explanation of i_initial_state
  ! i_initial_state : initial structure 
  ! -----------------------------------------------------------------
  write (lunout, *) 'i_initial_state = ', i_initial_state
  if(i_initial_state == INISTAT%RANDOM) then
     write(lunout, *) 'from random configuration.'
  else if(i_initial_state == INISTAT%NATIVE) then
     write(lunout, *) 'from native configuration.'
  else if(i_initial_state == INISTAT%INPUT) then
     write(lunout, *) 'from configuration given in the input.'
  else if(i_initial_state == INISTAT%CG) then
     write(lunout, *) 'from configuration given in the input with cafemol(CG) style.'
  else if(i_initial_state == INISTAT%RST) then
     write(lunout, *) 'from configuration given in the restart file.'
  else
     error_message = 'Error: invalid value about i_initial_state'
  end if
  
  ! -----------------------------------------------------------------
  ! The explanation of i_initial_velo
  ! i_initial_velo : initial velocities
  ! -----------------------------------------------------------------
  write (lunout, *) 'i_initial_velo = ', i_initial_velo
  if(i_initial_velo == INIVELO%MAXWELL) then
     write(lunout, *) 'Maxwell-Boltzmann distribution using random numbers.'
  else if(i_initial_velo == INIVELO%CARD) then
     write(lunout, *) 'from .velo file (CARD style).'
  else if(i_initial_velo == INIVELO%RST) then
     write(lunout, *) 'from restart file.'
  else
     error_message = 'Error: invalid value about i_initial_state'
  end if
  
  ! -----------------------------------------------------------------
  ! using periodic boundary condition, inperi%i_periodic
  if(inperi%i_periodic == 0) then
  else if(inperi%i_periodic == 1) then
     write(lunout, *) 'using periodic boundary condition: i_periodic = 1'
  else
     error_message = 'Error: invalid value for i_periodic'
     call util_error(ERROR%STOP_ALL, error_message)
  end if


#ifdef MPI_PAR
  end if

  call MPI_Bcast(i_simulate_type, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(i_initial_state, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(i_initial_velo,  1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inperi%i_periodic, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

#endif

!  write (*, *) myrank, inperi%i_periodic
!  stop

end subroutine inp_job
