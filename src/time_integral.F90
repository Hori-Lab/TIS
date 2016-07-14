!time_integral
!> @brief

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

  
subroutine time_integral(flg_step_each_replica)

  use const_maxsize
  use const_physical
  use const_index
  use if_mloop
  use if_write
  use if_energy
  use var_io,     only : i_run_mode, i_simulate_type
  use var_setp,    only : insimu, ifix_mp, inmmc
  use var_struct,  only : nmp_real, xyz_mp_rep, pxyz_mp_rep
  use var_replica, only : inrep, rep2val, rep2step, flg_rep, n_replica_mpi, exchange_step
  use var_simu,    only : istep, n_exchange, tstep, tstep2, tsteph, tempk, accelaf, &
                          accel_mp, velo_mp, force_mp, rcmass_mp, cmass_cs, &
                          e_md, fac_mmc, em_mid, em_depth, em_sigma, &
                          energy_muca, energy_unit_muca, rlan_const, &
                          ics, jcs, ncs, velo_yojou, evcs, xyz_cs, velo_cs, &
                          diffuse_tensor, random_tensor
  use time, only : tm_lap, tm_random, tmc_random, tm_muca, &
                   tm_neighbor, tm_update, tm_copyxyz, tm_force, &
                   time_s, time_e
  use mpiconst

  implicit none

  ! -----------------------------------------------------------------
  logical, intent(inout) :: flg_step_each_replica(n_replica_mpi)

  ! -----------------------------------------------------------------
  integer    :: i,k,imp, irep, grep
  real(PREC) :: r_force(1:SDIM), dxyz(1:3)
  real(PREC) :: r_boxmuller(SDIM, nmp_real, n_replica_mpi)
  real(PREC) :: random_vector(3*nmp_real) ! BROWNIAN_HI
  real(PREC) :: force_vector(3*nmp_real) ! BROWNIAN_HI

  ! --------------------------------------------------------------
  ! calc neighbour list
  ! --------------------------------------------------------------
  if(mod(istep, insimu%n_step_neighbor) == 1 .OR. istep == insimu%i_tstep_init) then  
     do irep=1, n_replica_mpi
        TIME_S( tm_neighbor )
        call neighbor(irep)
        TIME_E( tm_neighbor )
     enddo
  end if
  
  !if (inrep%i_loadbalance >= 1) then
  !   if(istep == insimu%i_tstep_init) then  
  !      call step_adjustment(istep, n_exchange, inrep%i_loadbalance)
  !      TIME_S( tm_lap )
  !   endif
  !endif

  ! -------------------------------------
  ! prepare random numbers for Langevin
  ! -------------------------------------
  r_boxmuller(:,:,:) = 0.0
  if (i_simulate_type == SIM%LANGEVIN .OR. &
      i_simulate_type == SIM%BROWNIAN .OR. i_simulate_type == SIM%BROWNIAN_HI) then
     !call get_random_number(r_boxmuller)
     call get_random_number()
  end if
  
  ! -------------------
  !  loop for REPLICAs
  ! -------------------
  do irep = 1, n_replica_mpi

     if (.not. flg_step_each_replica(irep)) then
        
!        grep = irep2grep(irep)
        
        if (flg_rep(REPTYPE%TEMP)) then
           tempk = rep2val(grep, REPTYPE%TEMP)
        endif
#ifdef _DEBUG
        write(6,*) 'mloop_simulator: tempk = ',tempk
#endif
        ! --------------------------------------------------------------
        ! move atoms
        ! --------------------------------------------------------------
        
        ! Langevin
        if(i_simulate_type == SIM%LANGEVIN) then
           
           TIME_S( tm_update )
           do imp = 1, nmp_real
              if(ifix_mp(imp) == 1) cycle
              
              ! xyz(t+h) update coordinates
              dxyz(1:3) = rlan_const(4, imp, irep) * velo_mp(1:3, imp, irep) &
                   + tstep2 * accel_mp(1:3, imp, irep)
              xyz_mp_rep(1:3, imp, irep) = xyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
              pxyz_mp_rep(1:3, imp, irep) = pxyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
           enddo
           TIME_E( tm_update )
           
           TIME_S( tm_copyxyz )
           call simu_copyxyz(irep)
           TIME_E( tm_copyxyz )
   
           TIME_S( tm_force )
           call force_sumup(force_mp, irep)
           TIME_E( tm_force )
           
           !mcanonical
           ! multicanonical algorithm --------------------
           ! based on Gosavi et al. JMB,2006,357,986
           TIME_S( tm_muca )
           if(inmmc%i_modified_muca == 1)then
              call energy_sumup(irep, velo_mp(:,:,irep), energy_muca, energy_unit_muca)
              e_md = energy_muca(E_TYPE%TOTAL)
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
           TIME_S( tm_update ) 
           do imp = 1, nmp_real
              if(ifix_mp(imp) == 1) cycle
              
              ! R(t+h)
              r_force(1:3) = rlan_const(1, imp, irep) *  r_boxmuller(1:3, imp, irep)
              ! a(t+h) temporary
              accelaf(1:3) = force_mp(1:3, imp) * rcmass_mp(imp) + r_force(1:3)
              
              ! v(h+h) update velocity
              velo_mp(1:3, imp, irep) = rlan_const(2, imp, irep) * velo_mp(1:3, imp, irep) &
                   + rlan_const(3, imp, irep) * (accel_mp(1:3, imp, irep) + accelaf(1:3))
              
              ! a(t+h) update acceleration
              accel_mp(1:3, imp, irep) = accelaf(1:3)
           end do
           TIME_E( tm_update )
           
           
           ! correcting velocity for removing translation and rotation motion
           if((insimu%i_no_trans_rot == 1) .and. (mod(istep, 200) == 1)) then
              call simu_velo_adjst(velo_mp,irep)
           end if
           
        ! Berendsen
        else if(i_simulate_type == SIM%BERENDSEN .or. i_simulate_type == SIM%CONST_ENERGY) then
           TIME_S( tm_update )
           ! xyz(t+h) update coordinates
           do imp = 1, nmp_real
              if(ifix_mp(imp) == 1) cycle
              dxyz(1:3) = tstep * velo_mp(1:3, imp, irep)    &
                   + tstep2 * accel_mp(1:3, imp, irep)
              xyz_mp_rep(1:3, imp, irep) = xyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
              pxyz_mp_rep(1:3, imp, irep) = pxyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
           end do
           TIME_E( tm_update )
           
           TIME_S( tm_copyxyz )
           call simu_copyxyz(irep)
           TIME_E( tm_copyxyz )
           
           TIME_S( tm_force )
           call force_sumup(force_mp, irep)
           TIME_E( tm_force )
           
           !mcanonical
           ! multicanonical algorithm --------------------
           ! based on Gosavi et al. JMB,2006,357,986
           TIME_S( tm_muca )
           if(inmmc%i_modified_muca == 1)then
              call energy_sumup(irep, velo_mp(:,:,irep), energy_muca, energy_unit_muca)
              e_md = energy_muca(E_TYPE%TOTAL)
              fac_mmc = 1 + em_depth * (e_md - em_mid) / (em_sigma*em_sigma) * &
                   exp(-(e_md - em_mid)**2 / (2.0e0_PREC*em_sigma*em_sigma))
              do imp = 1, nmp_real
                 force_mp(1:3, imp) = force_mp(1:3, imp) * fac_mmc
              end do
           endif
           TIME_E( tm_muca )
           !----------------------------------------------
           
           TIME_S( tm_update )
           do imp = 1, nmp_real
              if(ifix_mp(imp) == 1) cycle
              
              ! a(t+h) temporary
              accelaf(1:3) = force_mp(1:3, imp) * rcmass_mp(imp)
              
              ! v(t+h) update velocity
              velo_mp(1:3, imp, irep) = velo_mp(1:3, imp, irep)             &
                   + tsteph * (accel_mp(1:3, imp, irep) + accelaf(1:3))
              
              ! a(t+h) update acceleration
              accel_mp(1:3, imp, irep) = accelaf(1:3)
           end do
           
           ! correcting velocity for removing translation and rotation motion
           if((insimu%i_no_trans_rot == 1) .and. (mod(istep, 200) == 1)) then
              call simu_velo_adjst(velo_mp, irep)
           end if
           
           ! set temperature
           if(i_simulate_type == SIM%BERENDSEN) then
              call simu_velo_settemp(velo_mp, irep, tempk)
           end if
           TIME_E( tm_update )
           
#ifdef _DEBUG
           do imp=1, nmp_real
              write(6,'(2i5,1p3d15.7)'),irep,imp,force_mp(1,imp),force_mp(2,imp),force_mp(3,imp)
           enddo
#endif
   
        ! Nose-Hoover
        else if(i_simulate_type == SIM%NOSEHOOVER) then
           TIME_S( tm_update )
           velo_cs(ncs) = velo_cs(ncs) + tsteph * velo_yojou(ncs) / cmass_cs(ncs)
           xyz_cs(ncs) = xyz_cs(ncs) + tstep * velo_cs(ncs)
           evcs(ncs) = exp(-tsteph * velo_cs(ncs))
           do ics = 1, ncs-1
              jcs = ncs - ics
              velo_cs(jcs) = evcs(jcs+1) * velo_cs(jcs) + tsteph * velo_yojou(jcs) / cmass_cs(jcs)
              xyz_cs(jcs) = xyz_cs(jcs) + tstep * velo_cs(jcs)
              evcs(jcs) = exp(-tsteph * velo_cs(jcs))
           end do
           
           do imp = 1, nmp_real
              if(ifix_mp(imp) == 1) cycle
              velo_mp(1:3, imp, irep) = evcs(1) * velo_mp(1:3, imp, irep) &
                   + tsteph * accel_mp(1:3, imp, irep)
              dxyz(1:3) = tstep * velo_mp(1:3, imp, irep)
              xyz_mp_rep(1:3, imp, irep) = xyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
              pxyz_mp_rep(1:3, imp, irep) = pxyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
           end do
           TIME_E( tm_update )
           
           TIME_S( tm_copyxyz )
           call simu_copyxyz(irep)
           TIME_E( tm_copyxyz )
           
           TIME_S( tm_force )
           call force_sumup(force_mp, irep)
           TIME_E( tm_force )
           
           !mcanonical
           ! multicanonical algorithm --------------------
           ! based on Gosavi et al. JMB,2006,357,986
           TIME_S( tm_muca )
           if(inmmc%i_modified_muca == 1)then
              call energy_sumup(irep, velo_mp(:,:,irep), energy_muca, energy_unit_muca)
              e_md = energy_muca(E_TYPE%TOTAL)
              fac_mmc = 1 + em_depth * (e_md - em_mid) / (em_sigma*em_sigma) * &
                   exp(-(e_md - em_mid)**2 / (2.0e0_PREC*em_sigma*em_sigma))
              do imp = 1, nmp_real
                 force_mp(1:3, imp) = force_mp(1:3, imp) * fac_mmc
              end do
           endif
           TIME_E( tm_muca )
           !----------------------------------------------
           
           TIME_S( tm_update )
           do imp = 1, nmp_real
              if(ifix_mp(imp) == 1) cycle
              accel_mp(1:3, imp, irep) = force_mp(1:3, imp) * rcmass_mp(imp)
              velo_mp(1:3, imp, irep) = evcs(1) &
                   * (  velo_mp(1:3, imp, irep) &
                   + tsteph * accel_mp(1:3, imp, irep) )
           end do
           
           call simu_velo_nosehoover(velo_mp, irep, tempk, velo_yojou(1))
           
           do ics = 1, ncs - 1
              velo_cs(ics) = evcs(ics + 1) * (velo_cs(ics) + tsteph * velo_yojou(ics) / cmass_cs(ics))
              velo_yojou(ics + 1) = cmass_cs(ics) * velo_cs(ics)**2 - BOLTZ_KCAL_MOL * tempk
           end do
           velo_cs(ncs) = velo_cs(ncs) + tsteph * velo_yojou(ncs) / cmass_cs(ncs)
           
           ! correcting velocity for removing translation and rotation motion
           if ((insimu%i_no_trans_rot == 1) .and. (mod(istep, 200) == 1)) then
              call simu_velo_adjst(velo_mp, irep)
           end if
           
           TIME_E( tm_update )
           
        else if(i_simulate_type == SIM%BROWNIAN) then
           
           TIME_S( tm_update )
           do imp = 1, nmp_real
              if(ifix_mp(imp) == 1) cycle
              
              dxyz(1:3) =  rlan_const(1, imp, irep) * force_mp(1:3, imp)  &
                         + rlan_const(2, imp, irep) * r_boxmuller(1:3, imp, irep)

              xyz_mp_rep(1:3, imp, irep) = xyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
              pxyz_mp_rep(1:3, imp, irep) = pxyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
           enddo
           TIME_E( tm_update )
           
           TIME_S( tm_copyxyz )
           call simu_copyxyz(irep)
           TIME_E( tm_copyxyz )
   
           TIME_S( tm_force )
           call force_sumup(force_mp, irep)
           TIME_E( tm_force )
           
        else if(i_simulate_type == SIM%BROWNIAN_HI) then
           
           call simu_hydro_tensors(irep,tempk)

           random_vector = reshape(r_boxmuller(:,:,irep), (/ 3*nmp_real /) )
           force_vector = reshape(force_mp(:, 1:nmp_real), (/ 3*nmp_real /) )

           TIME_S( tm_update )
           do imp = 1, nmp_real
              if(ifix_mp(imp) == 1) cycle

              do i=1,3
                 k = 3*(imp-1) + i
                 dxyz(i) = rlan_const(1,imp,irep)  &
                              * dot_product( diffuse_tensor(:,k), force_vector(:) ) &
                           + rlan_const(2,imp,irep)  &
                              * dot_product( random_tensor(k, 1:k), random_vector(1:k) )
              enddo

              xyz_mp_rep(1:3, imp, irep) = xyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
              pxyz_mp_rep(1:3, imp, irep) = pxyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
           enddo
           TIME_E( tm_update )
           
           TIME_S( tm_copyxyz )
           call simu_copyxyz(irep)
           TIME_E( tm_copyxyz )
   
           TIME_S( tm_force )
           call force_sumup(force_mp, irep)
           TIME_E( tm_force )
           !----------------------------------------------
#ifdef _DEBUG
           do imp=1, nmp_real
              write(6,'(2i5,1p3d15.7)'),irep,imp,force_mp(1,imp),force_mp(2,imp),force_mp(3,imp)
           enddo
#endif
        end if
        
        if (i_run_mode == RUN%REPLICA .and. inrep%flg_exchange) then
           ! write(6,*) ' -- ', istep, exchange_step(irep)
           if (istep == exchange_step(grep)) then
              flg_step_each_replica(irep) = .true.
           endif
        end if

     end if   !  if (.not. flg_step_each_replica(irep)) then
  enddo
  ! irep --------------------------------------------------------- 


contains

  subroutine get_random_number()
  !subroutine get_random_number(r_boxmuller)

    use mt_stream
    use mt_kind_defs
#ifdef MPI_PAR
    use var_setp, only : mts
#endif
    implicit none

    ! --------------------------------------------------------------------
    !real(PREC), intent(out) :: r_boxmuller(SDIM, nmp_real, n_replica_mpi)
    
    ! --------------------------------------------------------------------
    ! function
    real(PREC) :: rfunc_boxmuller
    
    ! --------------------------------------------------------------------
    integer :: irep, idimn, istream
#ifdef MPI_PAR
    real(PREC) :: vx, vy, r2, rf
    integer :: klen, ksta, kend, tn
    real(PREC) :: r_boxmuller_l(SDIM, nmp_real, n_replica_mpi)

    integer(INT32) :: umask,lmask,n,m,is
    integer(INT32) :: k,nm,n1
    real(REAL64)   :: a,b,rx,ry
    integer(INT32) :: ia,ib
#endif
    
    TIME_S( tm_random)
    
    if(insimu%i_rand_type == 0) then
       do irep = 1, n_replica_mpi
          istream = irep
          do imp= 1, nmp_real
             do idimn = 1, SDIM
                r_boxmuller(idimn, imp, irep) = rfunc_boxmuller(istream, 0)
             end do
          end do
       end do
       
    else
       
#ifdef MPI_PAR

       r_boxmuller_l(:, :, :) = 0.0

       do irep = 1, n_replica_mpi
          if(insimu%i_rand_type == 1) then
             klen=(nmp_real-1+npar_mpi)/npar_mpi
             ksta=1+klen*local_rank_mpi
             kend=min(ksta+klen-1,nmp_real)
          else
             ksta=1
             kend=nmp_real
          end if
!$omp parallel private(tn,istream)
          tn = 0
!$ tn = omp_get_thread_num()
          if(insimu%i_rand_type == 1) then
             istream = irep + local_rank_mpi*n_replica_mpi
          else
             istream = irep
          end if
!$omp do private(idimn,vx,vy,r2,rf,&
!$omp&           n,m,lmask,umask,nm,n1,k,is,&
!$omp&           ia,ib,a,b,rx,ry)
          do imp = ksta, kend, 2
!          do imp = ksta, kend
             do idimn = 1, 3
                do
#ifndef SUB_COPY
                   vx = 2.0e0_PREC * genrand_double1(mts(istream, tn)) - 1.0e0_PREC
                   vy = 2.0e0_PREC * genrand_double1(mts(istream, tn)) - 1.0e0_PREC

#else               
                   if (mts(istream, tn)%i >= mts(istream, tn)%nn-3) then
                      n = mts(istream, tn)%nn
                      m = mts(istream, tn)%mm
                      lmask = mts(istream, tn)%lmask
                      umask = mts(istream, tn)%umask
                      nm = n - m
                      n1 = n - 1
                      do k=0,nm-1
                         is = IOR(IAND(mts(istream, tn)%state(k),umask),IAND(mts(istream, tn)%state(k+1),lmask))
                         mts(istream, tn)%state(k) = IEOR(IEOR(mts(istream, tn)%state(k+m),ISHFT(is,-1)),mts(istream, tn)%mag(IAND(is,1)))
                      enddo
                      do k=nm,n1-1
                         is = IOR(IAND(mts(istream, tn)%state(k),umask),IAND(mts(istream, tn)%state(k+1),lmask))
                         mts(istream, tn)%state(k) = IEOR(IEOR(mts(istream, tn)%state(k+m-n),ISHFT(is,-1)),mts(istream, tn)%mag(IAND(is,1)))
                      enddo
                      is = IOR(IAND(mts(istream, tn)%state(n-1),umask),IAND(mts(istream, tn)%state(0),lmask))
                      mts(istream, tn)%state(n-1) = IEOR(IEOR(mts(istream, tn)%state(m-1),ISHFT(is,-1)),mts(istream, tn)%mag(IAND(is,1)))
                      mts(istream, tn)%i = 0
                   endif
                   
                   is = mts(istream, tn)%state(mts(istream, tn)%i)
                   mts(istream, tn)%i = mts(istream, tn)%i + 1
                   is = IEOR(is,ISHFT(is,-mts(istream, tn)%shift0))
                   is = IEOR(is,IAND(ISHFT(is, mts(istream, tn)%shiftB),mts(istream, tn)%maskB))
                   is = IEOR(is,IAND(ISHFT(is, mts(istream, tn)%shiftC),mts(istream, tn)%maskC))
                   is = IEOR(is,ISHFT(is,-mts(istream, tn)%shift1))
                   ia = is
                   
                   is = mts(istream, tn)%state(mts(istream, tn)%i)
                   mts(istream, tn)%i = mts(istream, tn)%i + 1
                   is = IEOR(is,ISHFT(is,-mts(istream, tn)%shift0))
                   is = IEOR(is,IAND(ISHFT(is, mts(istream, tn)%shiftB),mts(istream, tn)%maskB))
                   is = IEOR(is,IAND(ISHFT(is, mts(istream, tn)%shiftC),mts(istream, tn)%maskC))
                   is = IEOR(is,ISHFT(is,-mts(istream, tn)%shift1))
                   ib = is
                   
                   ia = ISHFT(ia,-5)             ! ia in [0,2^27-1]
                   ib = ISHFT(ib,-6)             ! ib in [0,2^26-1]
                   a = REAL(ia,kind=KIND(rx))
                   b = REAL(ib,kind=KIND(rx))
                   !===============================
                   ! ( a*2^26 + b ) in [0,2^53-1]
                   ! r = ( a*2^26 + b )/(2^53-1)
                   !===============================
                   rx = (a*67108864.0_REAL64 + b)*(1.0_REAL64/9007199254740991.0_REAL64)
                   is = mts(istream, tn)%state(mts(istream, tn)%i)
                   mts(istream, tn)%i = mts(istream, tn)%i + 1
                   is = IEOR(is,ISHFT(is,-mts(istream, tn)%shift0))
                   is = IEOR(is,IAND(ISHFT(is, mts(istream, tn)%shiftB),mts(istream, tn)%maskB))
                   is = IEOR(is,IAND(ISHFT(is, mts(istream, tn)%shiftC),mts(istream, tn)%maskC))
                   is = IEOR(is,ISHFT(is,-mts(istream, tn)%shift1))
                   ia = is
                   
                   is = mts(istream, tn)%state(mts(istream, tn)%i)
                   mts(istream, tn)%i = mts(istream, tn)%i + 1
                   is = IEOR(is,ISHFT(is,-mts(istream, tn)%shift0))
                   is = IEOR(is,IAND(ISHFT(is, mts(istream, tn)%shiftB),mts(istream, tn)%maskB))
                   is = IEOR(is,IAND(ISHFT(is, mts(istream, tn)%shiftC),mts(istream, tn)%maskC))
                   is = IEOR(is,ISHFT(is,-mts(istream, tn)%shift1))
                   ib = is
                   
                   ia = ISHFT(ia,-5)             ! ia in [0,2^27-1]
                   ib = ISHFT(ib,-6)             ! ib in [0,2^26-1]
                   a = REAL(ia,kind=KIND(ry))
                   b = REAL(ib,kind=KIND(ry))
                   !===============================
                   ! ( a*2^26 + b ) in [0,2^53-1]
                   ! r = ( a*2^26 + b )/(2^53-1)
                   !===============================
                   ry = (a*67108864.0_REAL64 + b)*(1.0_REAL64/9007199254740991.0_REAL64)
            
                   vx = 2.0e0_PREC * rx - 1.0e0_PREC
                   vy = 2.0e0_PREC * ry - 1.0e0_PREC             
                   
#endif                   
                   r2 = vx*vx + vy*vy
                   if(r2 < 1.0e0_PREC .and. r2 > 0.0e0_PREC) exit
                end do
                
                rf = sqrt(-2.0e0_PREC * log(r2) / r2)
                r_boxmuller_l(idimn, imp, irep) = vy * rf
                
                if(imp < kend) then
                   r_boxmuller_l(idimn, imp+1, irep) = vx * rf
                end if
             end do
          end do

!$omp end do
!$omp end parallel
       end do

       TIME_S( tmc_random)
       if(insimu%i_rand_type == 1) then
          call mpi_allreduce(r_boxmuller_l, r_boxmuller, &
               SDIM*nmp_real*n_replica_mpi, PREC_MPI, &
               MPI_SUM, mpi_comm_local, ierr)
       else
          r_boxmuller(:,:,:) = r_boxmuller_l(:,:,:)
       end if
       TIME_E( tmc_random)

#else
       do irep = 1, n_replica_mpi
          istream = irep
          do imp = 1, nmp_real
             do idimn = 1, SDIM
                r_boxmuller(idimn, imp, irep) = rfunc_boxmuller(istream, 0)
             end do
          end do
       end do
#endif
    end if
    TIME_E( tm_random)

  end subroutine get_random_number

end subroutine time_integral
