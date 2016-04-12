!simu_tingegral_pre
!> @brief

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

subroutine simu_tintegral_pre(flg_step_each_replica)

  use const_maxsize
  use const_physical
  use const_index
  use if_mloop
  use if_write
  use if_energy
  use var_inp,     only : i_run_mode, i_simulate_type, flg_rst
  use var_setp,    only : inpara, insimu, inmmc, inmisc
  use var_struct,  only : nmp_real, cmass_mp, fric_mp, grp, xyz_mp_rep
  use var_replica, only : rep2val, flg_rep, &
                          n_replica_all, n_replica_mpi, irep2grep
  use var_simu,    only : istep_sim, &
                          ntstep, nstep_opt_temp, ibefore_time, &
                          tstep, tstep2, tsteph, tempk, accelaf, &
                          accel_mp, velo_mp, force_mp, rcmass_mp, &
                          cmass_cs, &
                          e_md, fac_mmc, em_mid, em_depth, em_sigma, &
                          pnlet_muca, pnle_unit_muca, &
                          rlan_const, &
                          tstep_fric_h, ulconst1, ulconst2, &
                          ics, jcs, ncs, velo_yojou, evcs, xyz_cs, velo_cs, &
                          pnlet, pnle_unit, &
                          rg, rg_unit, rmsd, rmsd_unit, &
                          replica_energy
  use time, only : tm_random, tm_muca, time_s, time_e

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! -----------------------------------------------------------------
  logical, intent(inout) :: flg_step_each_replica(n_replica_mpi)

  ! --------------------------------------------------------------------
  ! function
  real(PREC) :: rfunc_boxmuller

  ! -----------------------------------------------------------------
  ! local variables
  integer :: imp, irep, grep, idimn, istream
  integer(L_INT) :: istep_dummy
  real(PREC) :: r_force(1:SPACE_DIM)

  real(PREC) :: r_boxmuller(SPACE_DIM, nmp_real, n_replica_mpi)
  character(CARRAY_MSG_ERROR) error_message

  integer :: ianc, igrp

  ! -----------------------------------------------------------------
  ! calc neighbour list
  ! -----------------------------------------------------------------
  do irep=1, n_replica_mpi
     call simu_neighbor(irep)
  enddo
  


  ! -------------------------------------
  ! prepare random numbers for Langevin
  ! -------------------------------------
  if (i_simulate_type == SIM%LANGEVIN .AND.  &
       (istep_sim == 1 .OR. inmisc%i_reset_struct == 1)) then
     TIME_S( tm_random)
     do irep = 1, n_replica_mpi
        istream = irep
        do imp= 1, nmp_real
           do idimn = 1, SPACE_DIM
              r_boxmuller(idimn, imp, irep) = rfunc_boxmuller(istream, 0)
           enddo
        enddo
     enddo
     TIME_E( tm_random)
  endif
  
  ! -------------------
  !  loop for REPLICAs
  ! -------------------
  do irep = 1, n_replica_mpi
     
     grep = irep2grep(irep)
     
     ! correcting velocity for removing translation and rotation motion
     if(insimu%i_no_trans_rot == 1) then
        call simu_velo_adjst(velo_mp,irep)
     end if
     
     if (flg_rep(REPTYPE%TEMP)) then
        tempk = rep2val(grep, REPTYPE%TEMP)
     endif
#ifdef _DEBUG
     write(6,*) 'mloop_simulator: tempk = ',tempk
#endif

     ! Set parameter for anchor of COM
     if (inmisc%nanc_com_ini > 0) then
        do ianc = 1, inmisc%nanc_com_ini
           igrp = inmisc%ianc_com_ini2grp(ianc)
           
           inmisc%anc_com_ini_xyz(:, ianc, irep) = 0.0
           do imp = 1, grp%nmp(igrp)
              inmisc%anc_com_ini_xyz(:, ianc, irep) = &
                   inmisc%anc_com_ini_xyz(:, ianc, irep) + &
                   xyz_mp_rep(:, grp%implist(imp, igrp), irep) * &
                   grp%mass_fract(imp, igrp)
           end do
        end do
     endif
     
     ! calc force and acceleration
     call simu_force(force_mp, irep)

     !mcanonical
     ! multicanonical algorithm --------------------
     ! based on Gosavi et al. JMB,2006,357,986
     TIME_S( tm_muca )
     if(inmmc%i_modified_muca == 1)then
        call simu_energy(irep, velo_mp(:,:,irep), pnlet_muca, pnle_unit_muca)
        e_md = pnlet_muca(E_TYPE%TOTAL)
        fac_mmc = 1 + em_depth * (e_md - em_mid) / (em_sigma*em_sigma) * &
             exp(-(e_md - em_mid)**2 / (2.0e0_PREC*em_sigma*em_sigma))
        do imp = 1, nmp_real
           force_mp(1:3, imp) = force_mp(1:3, imp) * fac_mmc
        end do
     endif
     TIME_E( tm_muca )
     !----------------------------------------------

#ifdef _DEBUG
     do imp=1, nmp_real
        write(6,'(2i5,1p3d15.7)'),irep,imp,force_mp(1,imp),force_mp(2,imp),force_mp(3,imp)
     enddo
#endif
        
     ! Langevin
     if(i_simulate_type == SIM%LANGEVIN) then
        do imp = 1, nmp_real
           tstep_fric_h = 0.5e0_PREC * tstep * fric_mp(imp)
           ulconst1 = 1.0e0_PREC - tstep_fric_h
           ulconst2 = 1.0e0_PREC - tstep_fric_h + tstep_fric_h**2
           
           rlan_const(1, imp, irep)  = sqrt( 2.0e0_PREC * fric_mp(imp) * BOLTZC * tempk &
                / (tstep * cmass_mp(imp) )    )
           rlan_const(2, imp, irep) = ulconst1 * ulconst2
           rlan_const(3, imp, irep) = tsteph * ulconst1
           rlan_const(4, imp, irep) = tstep * ulconst1
        end do
        
        if (flg_rst) then
           call read_rst(RSTBLK%ACCEL)
        else if (istep_sim == 1 .OR. inmisc%i_reset_struct == 1) then
           do imp = 1, nmp_real
              r_force(1:3) = rlan_const(1, imp, irep) *  r_boxmuller(1:3, imp, irep)
              accel_mp(1:3, imp, irep) = force_mp(1:3, imp) * rcmass_mp(imp) + r_force(1:3)
           end do
        endif
        
#ifdef _DEBUG
        do imp=1, nmp_real
           write(6,'(2i5,1pd15.7)'),irep,imp,accel_mp(1,imp,irep)
        enddo
#endif
     ! Berendsen or Constant Energy
     else if(i_simulate_type == SIM%BERENDSEN .OR. i_simulate_type == SIM%CONST_ENERGY .or. i_simulate_type == SIM%MPC) then
        if (flg_rst) then
           call read_rst(RSTBLK%ACCEL)
        else if (istep_sim == 1 .OR. inmisc%i_reset_struct == 1) then
           do imp = 1, nmp_real
              accel_mp(1:3, imp, irep) = force_mp(1:3, imp) * rcmass_mp(imp)
           end do
           if((i_simulate_type /= SIM%MPC)) then
              call simu_velo_settemp(velo_mp, irep, tempk)
           end if
        endif
        
     ! Nose-Hoover
     else if(i_simulate_type == SIM%NOSEHOOVER) then
        if (flg_rst) then
           call read_rst(RSTBLK%ACCEL)
        else if (istep_sim == 1 .OR. inmisc%i_reset_struct == 1) then
           do imp = 1, nmp_real
              accel_mp(1:3, imp, irep) = force_mp(1:3, imp) * rcmass_mp(imp)
           end do
        end if
           
        ncs = MXCS
        cmass_cs(1:ncs) = inpara%csmass_per*nmp_real
        xyz_cs(1:ncs) = 0.0e0_PREC
        velo_cs(1:ncs) = sqrt(BOLTZC * tempk / cmass_cs(1:ncs))
        call simu_velo_nosehoover(velo_mp, irep, tempk, velo_yojou(1))
        do ics = 2, ncs
           velo_yojou(ics) = cmass_cs(ics-1) * velo_cs(ics-1)**2 - BOLTZC * tempk
        end do
     endif
     
  enddo ! irep -----------------------------------------------
  
  ! calc energy
  ! ## Here, replica_energy is dummy.
  call simu_energy_allrep(pnle_unit, pnlet, &
       velo_mp, replica_energy, .false., tempk)
  call simu_radiusg(rg_unit, rg)
  call simu_rmsd(rmsd_unit, rmsd)
  
#ifdef MPI_PAR
  if (local_rank_mpi == 0) then
#endif

  ! writing initial energy and tag for t-series
  if (flg_rst) then
     istep_dummy = insimu%i_tstep_init - 1
  else
     istep_dummy = 0
  endif
  call write_tseries(ibefore_time, istep_dummy, &
                     rg_unit, rg, rmsd_unit, rmsd, &
                     pnle_unit, pnlet, tempk, .true.)
  
  if(i_run_mode == RUN%ENERGY_CALC) then
     write(error_message, *) 'PROGRAM STOP'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
#ifdef MPI_PAR
  end if
#endif
  
  if(insimu%i_com_zeroing == 1 .or. insimu%i_com_zeroing == 2) then
     call simu_xyz_adjst()
  end if
  
  if(istep_sim == insimu%i_step_sim_init) then
     call write_traject_file(ibefore_time, istep_dummy, tempk, velo_mp)
  endif
  
  call write_record_file(istep_dummy, velo_mp)
  
  flg_step_each_replica(1:n_replica_mpi)  = .false. 

  if (inmisc%nbrid_ppr > 0) then
     call simu_ppr()
  endif

  if (inmisc%nanc_com_ini > 0) then
     do irep = 1, n_replica_mpi
        do ianc = 1, inmisc%nanc_com_ini
           igrp = inmisc%ianc_com_ini2grp(ianc)
           
           inmisc%anc_com_ini_xyz(:, ianc, irep) = 0.0
           do imp = 1, grp%nmp(igrp)
              inmisc%anc_com_ini_xyz(:, ianc, irep) = &
                   inmisc%anc_com_ini_xyz(:, ianc, irep) + &
                   xyz_mp_rep(:, grp%implist(imp, igrp), irep) * &
                   grp%mass_fract(imp, igrp)
           end do
        end do
     end do
  endif

end subroutine simu_tintegral_pre
