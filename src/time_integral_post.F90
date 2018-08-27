!time_integral_post
!> @brief

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
subroutine time_integral_post(flg_step_each_replica, flg_exit_loop_mstep)
  
  use const_maxsize
  use const_physical
  use const_index
  use if_mloop
  use if_write
  use if_energy
  use var_io,     only : i_run_mode, i_simulate_type, flg_file_out, outfile
  use var_setp,    only : insimu, inmisc, inele, inwidom !, inann
  use var_struct,  only : nmp_real, cmass_mp, fric_mp, xyz_mp_rep
  use var_replica, only : inrep, rep2val, rep2step, flg_rep, &
                          n_replica_mpi, irep2grep, exchange_step
  use var_simu,    only : istep, ntstep, nstep_opt_temp, ibefore_time, &
                          n_exchange, iopt_stage, &
                          tstep, tempk, velo_mp, rlan_const, &
                          energy, energy_unit, rg, rg_unit, rmsd, rmsd_unit, &
                          ene_st, ene_tst, ene_hb,&
                          replica_energy
  use time, only : tm_energy, tm_radiusg_rmsd, &
                   tm_output, tm_replica, & 
                   time_s, time_e, tm_others !, tm_step_adj, &
!  use var_fmat,    only : infmat
#ifdef MPI_PAR
  use mpiconst
  use var_simu, only : replica_energy_l
  use time , only : tmc_replica
  use var_replica, only : n_replica_all
#endif

  implicit none

  logical, intent(inout) :: flg_step_each_replica(n_replica_mpi)
  logical, intent(inout) :: flg_exit_loop_mstep

  logical :: flg_step_save, flg_step_rep_exc, flg_step_rep_save, flg_step_rep_opt
  logical :: flg_step_rst, flg_step_widom
  logical :: flg_CCX_cross
  integer :: imp1, imp2
  integer :: imp, irep, grep
  real(PREC) :: v21(3), dee
  character(CARRAY_MSG_ERROR) :: error_message
#ifdef _DUMP_COMMON
  integer :: lundump = outfile%dump
#endif

  ! -----------------------------------------------------------------
  TIME_S( tm_others )
  flg_CCX_cross     = .false.
  flg_step_save     = .false. 
  flg_step_rep_exc  = .false.
  flg_step_rep_save = .false.
  flg_step_rep_opt  = .false.
  flg_step_rst      = .false.
  flg_step_widom    = .false.

  if (mod(istep, insimu%n_step_save)   == 0) flg_step_save = .true.

  if (i_run_mode == RUN%REPLICA) then
     if (n_replica_mpi == count(flg_step_each_replica))    flg_step_rep_exc  = .true.
     if (mod(istep, inrep%n_step_save) == 0)               flg_step_rep_save = .true.
     if (inrep%flg_opt_temp .and. istep == nstep_opt_temp) flg_step_rep_opt  = .true.
  endif
  
  if (flg_file_out%rst) then
     if (mod(istep, insimu%n_step_rst) == 0) then
        flg_step_rst = .true.
     endif
  endif

  if (i_run_mode == RUN%WIDOM) then
     if (istep >= inwidom%n_step_skip .and. mod(istep, inwidom%n_step_interval) == 0) then
        flg_step_widom = .true.
     endif
  endif
  
  TIME_E( tm_others )

  ! --------------------------------------------------------------
  ! Check chain crossing (CCX)
  ! --------------------------------------------------------------
  if (inmisc%flg_CCX) then
     call check_chaincrossing(flg_CCX_cross)
      
     if (flg_CCX_cross) then
        flg_step_save     = .true.
        flg_step_rep_save = .true.
        !flg_step_rst      = .true.
     endif
  endif

  ! --------------------------------------------------------------
  ! energy calculation
  ! --------------------------------------------------------------
  if (    flg_step_save      &     ! to save
     .OR. istep == 1       ) then  ! to save 1st step
     
     TIME_S( tm_energy )
     replica_energy(:,:) = 0.0e0_PREC
     call energy_allrep(energy_unit, energy, velo_mp, replica_energy, flg_step_rep_exc, tempk)
     TIME_E( tm_energy )
     
     TIME_S( tm_radiusg_rmsd)
     call simu_radiusg(rg_unit, rg)
     call simu_rmsd(rmsd_unit, rmsd)
     TIME_E( tm_radiusg_rmsd )

  else if (flg_step_rep_exc) then

     replica_energy(:,:) = 0.0e0_PREC

     TIME_S( tm_energy )
     call energy_allrep(energy_unit, energy, velo_mp, replica_energy, flg_step_rep_exc, tempk)
     TIME_E( tm_energy )

  endif

  if (flg_step_widom) then
     call widom()
  endif
  
  ! --------------------------------------------------------------
  ! write data
  ! --------------------------------------------------------------
  TIME_S( tm_output )
  if(flg_step_save .OR. istep == 1 ) then
     
     if(flg_step_save) then 
        if(insimu%i_com_zeroing == 1 .or. insimu%i_com_zeroing == 2) then
           call simu_xyz_adjst()
        end if
        call write_traject_file(ibefore_time, istep, tempk, velo_mp)
     end if
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     call write_tseries(ibefore_time, istep, &
                        rg_unit, rg,                  &
                        rmsd_unit, rmsd,              &
                        energy_unit, energy, tempk)

     if (flg_file_out%st .or. flg_file_out%tst .or. flg_file_out%stall .or. flg_file_out%tstall) then
        call write_stack(ene_st, ene_tst, tempk)
     endif

     if (flg_file_out%hb .or. flg_file_out%hball) then
        call write_hbond(ene_hb, tempk)
     endif

     if (flg_file_out%opt) then
        ! something to write to opt file
     endif
     if (flg_file_out%ee) then
        irep = 1
        imp1 = 1
        imp2 = nmp_real
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        dee = sqrt(dot_product(v21,v21))
        write(outfile%ee(irep), '(i12,1x,f7.2,3(1xf7.2))') istep, dee, v21(1),v21(2),v21(3)
     endif

#ifdef MPI_PAR
     end if
     ! call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
  end if
  TIME_E( tm_output )
  
  TIME_S( tm_others )

  ! --------------------------------------------------------------
  ! write RECORD file
  ! --------------------------------------------------------------
  if(istep == 0 .OR. istep == ntstep) then
     if(insimu%i_com_zeroing == 1 .OR. insimu%i_com_zeroing == 2) then
        call simu_xyz_adjst()
     end if
     call write_record_file(istep)
  end if
  
  ! --------------------------------------------------------------
  ! annealing
  ! --------------------------------------------------------------
!  if(i_run_mode == RUN%SA .AND. mod(istep, (ntstep/(inann%n_time_change+1))) == 0) then
!     call simu_anneal(istep, ntstep, tempk)
!     
!     call simu_para2(tempk, inele%ionic_strength)
!     
!     if(i_simulate_type == SIM%LANGEVIN) then
!        do imp = 1, nmp_real
!           rlan_const(1, imp, 1) = sqrt(2.0e0_PREC * fric_mp(imp) * BOLTZ_KCAL_MOL * tempk &
!                / (tstep * cmass_mp(imp) )        )   
!           !NOTE: SA should be only 1 replica.
!        end do
!     elseif (i_simulate_type == SIM%BROWNIAN .or. i_simulate_type == SIM%PS_BROWNIAN) then
!        do imp = 1, nmp_real
!           rlan_const(2, imp, irep) = sqrt( 2.0e0_PREC * BOLTZ_KCAL_MOL * tempk * tstep / fric_mp(imp))
!           !NOTE: SA should be only 1 replica.
!        end do
!     end if
!  end if
  
  ! ------------------------
  ! Replica Exchange
  ! ------------------------
  TIME_E( tm_others )
  TIME_S( tm_replica )
  if (i_run_mode == RUN%REPLICA) then
     
     ! write trajectory
     if (flg_file_out%rep .AND. (istep == insimu%i_tstep_init .OR. flg_step_rep_save)) then
        call write_rep_traject(istep)
     endif
     
     if (flg_step_rep_exc) then
        
        !              write(6,*) ' << Replica Exchange istep >> ', istep, n_exchange
        
        flg_step_each_replica(1:n_replica_mpi)  = .false. 
#ifdef MPI_PAR
        TIME_S( tmc_replica )
        replica_energy_l(:,:) = replica_energy(:,:)
        call mpi_allreduce(replica_energy_l, replica_energy, 2*n_replica_all, PREC_MPI, &
             MPI_SUM, mpi_comm_rep, ierr)
        
        TIME_E( tmc_replica )
#endif
        
        n_exchange = n_exchange + 1
!              write(6,*) ' << Replica Exchange istep >> ', istep, n_exchange
        call simu_replica_exchange(velo_mp, replica_energy, tempk)
! DBGs
!             do irep=1, n_replica_all
!               write(6,*) irep,replica_energy(1,irep)
!             enddo
! DBGe
        if (flg_rep(REPTYPE%TEMP)) then
           ! update rlan_const (depend on temperature)
           if(i_simulate_type == SIM%LANGEVIN) then
              
              do irep = 1, n_replica_mpi
                 grep  = irep2grep(irep)
                 tempk = rep2val(grep, REPTYPE%TEMP)
                 do imp   = 1, nmp_real
                    rlan_const(1, imp, irep)  &
                         = sqrt( 2.0e0_PREC * fric_mp(imp) * BOLTZ_KCAL_MOL * tempk & 
                                     / (tstep * cmass_mp(imp)) )
                 enddo
              enddo

           else if(i_simulate_type == SIM%ND_LANGEVIN) then
              
              do irep = 1, n_replica_mpi
                 grep  = irep2grep(irep)
                 tempk = rep2val(grep, REPTYPE%TEMP)
                 do imp   = 1, nmp_real
                    rlan_const(3, imp, irep) = sqrt(2.0e0_PREC * fric_mp(imp) * BOLTZ_KCAL_MOL * tempk / tstep)
                 enddo
              enddo

           endif ! LANGEVIN or ND_LANGEVIN
        endif
        
        if (flg_rep(REPTYPE%TEMP) .OR. flg_rep(REPTYPE%ION) .OR. &
            flg_rep(REPTYPE%WIND) .OR. flg_rep(REPTYPE%PULL)) then

           call simu_para2(tempk, inele%ionic_strength)

        end if
        
        !if (inrep%i_loadbalance >= 1) then
        !   if (mod(n_exchange, inrep%n_adjust_interval) == 0) then
        !      
        !      TIME_S(tm_step_adj)
        !      call step_adjustment(istep, n_exchange, inrep%i_loadbalance)
        !      TIME_E(tm_step_adj)
        !
        !   end if
        !endif
        
        do irep = 1, n_replica_mpi
           grep = irep2grep(irep)
           exchange_step(grep) = istep + rep2step(grep)
!#ifdef _DEBUG
!                write(6,'(a,4i5)') ' - irep,grep,exchange_step(grep),rep2step(grep) - ' ,&
!                                       irep,grep,exchange_step(grep),rep2step(grep)
!#endif
        enddo

     endif ! flg_step_rep_exc
     
     if (flg_step_rep_opt) then
        call simu_replica_opt_temp(iopt_stage)
        call write_rep_table()
        nstep_opt_temp = nstep_opt_temp + inrep%n_step_opt_temp
        
        if (iopt_stage >= inrep%n_stage_opt_temp) then
           flg_exit_loop_mstep = .TRUE.
           return
        endif
     endif
     
  endif ! Replica exchange
  TIME_E( tm_replica )
  TIME_S( tm_others )
  
!  ! --------------------------------------------------------------
!  ! fmat
!  ! --------------------------------------------------------------
!  if(i_run_mode == RUN%FMAT .AND. mod(istep, infmat%n_step_interval) == 0) then
!     
!     call simu_fmat()
!     
!     if (mod(istep, infmat%n_step_save) == 0) then
!     endif
!     
!  end if
  
#ifdef _DUMP_COMMON
  if (istep <= 5) then
#endif 
#ifdef _DUMP
  call dump_var_all()
#endif
#ifdef _DUMP_REPLICA
  call dump_var_replica(lundump)
#endif
#ifdef _DUMP
#include   "dump_mloop_simulator.F90"
#endif
#ifdef _DUMP_COMMON
  endif
#endif

  TIME_E( tm_others )
  TIME_S( tm_output )
#ifdef MPI_PAR
  if (local_rank_mpi == 0) then
#endif
  if (flg_step_rst) then
     call write_rst()
  endif
#ifdef MPI_PAR
  endif
#endif

  TIME_E( tm_output )
  TIME_S( tm_others )
  !if (i_run_mode == RUN%REPLICA) then
  !   if (inrep%i_loadbalance >= 1) then
  !      if (n_exchange == max_exchange) then
  !         flg_exit_loop_mstep = .TRUE.
  !         return
  !      endif
  !   endif
  !endif

  if (inmisc%nbrid_ppr > 0) then
     call simu_ppr()
  endif
  TIME_E( tm_others )

  if (flg_CCX_cross) then
     write(error_message,*) 'Stopped because a chain crossing detected at istep=', istep
     call util_error(ERROR%STOP_ALL, error_message)
  endif

end subroutine time_integral_post
