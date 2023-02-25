! simu_xyz_adjst
!> @brief Correct coordinates to move the system to the origin

subroutine simu_xyz_adjst()

  use const_maxsize
  use const_index
  use const_physical
  use var_io,      only : i_simulate_type, flg_file_out, outfile
  use var_setp,    only : insimu, inperi
  use var_struct,  only : nmp_real, xyz_mp_rep, pxyz_mp_rep, iclass_mp
  use var_replica, only : n_replica_mpi, irep2grep
  use var_simu,    only : accel_mp, dxyz_mp, istep, dxyz_mp
  use mpiconst

  implicit none

  integer :: i, imp, irep, n_mp, grep
  real(PREC) :: xyzg(SDIM), sumxyz(SDIM)
  real(PREC) :: pxyzg(SDIM), psumxyz(SDIM)
  real(PREC) :: theta
  real(PREC) :: s_cos(SDIM), s_sin(SDIM)


  if (inperi%i_periodic == 1) then

     do irep = 1, n_replica_mpi
   
        s_cos(1:3) = 0.0e0_PREC
        s_sin(1:3) = 0.0e0_PREC
        n_mp = 0
        do imp = 1, nmp_real
           if (iclass_mp(imp) == CLASS%ION) then
               cycle
           endif
           do i = 1, SDIM
              theta = (pxyz_mp_rep(i,imp,irep) + inperi%psizeh(i)) / inperi%psize(i) * 2 * F_PI
              s_cos(i) = s_cos(i) + cos(theta)
              s_sin(i) = s_sin(i) + sin(theta)
           enddo
           n_mp = n_mp + 1
        enddo
   
        do i = 1, SDIM
           s_cos(i) = s_cos(i) / real(n_mp)
           s_sin(i) = s_sin(i) / real(n_mp)
           theta = atan2(-s_sin(i),-s_cos(i)) + F_PI
           pxyzg(i) = 0.5 * inperi%psize(i) * theta / F_PI - inperi%psizeh(i)
        enddo
        
        do imp = 1, nmp_real
           pxyz_mp_rep(:, imp, irep) = pxyz_mp_rep(:, imp, irep) - pxyzg(:)
        end do

        !### Wrap
        if (flg_file_out%neigh) then
#ifdef MPI_PAR
        if (local_rank_mpi == 0)then
#endif
           grep = irep2grep(irep)
           write(outfile%neigh(grep), '(i10,1x,i5,1x,f4.1,1x,f4.1,1x,f4.1)',advance='no') istep, grep, 0.0, 0.0, 0.0
#ifdef MPI_PAR
        endif
#endif
        endif
        call neighbor(irep)   ! util_periodic will be called in neighbor for PBC wrapping.
        dxyz_mp(:,:,irep) = 0.0e0_PREC
     enddo

  else

     do irep = 1, n_replica_mpi
   
        sumxyz(1:3) = 0.0e0_PREC
        psumxyz(1:3) = 0.0e0_PREC
        do imp = 1, nmp_real
           sumxyz(1:3) = sumxyz(1:3) + xyz_mp_rep(1:3, imp, irep)
           psumxyz(1:3) = psumxyz(1:3) + pxyz_mp_rep(1:3, imp, irep)
        end do
      
        xyzg(1:3) = sumxyz(1:3) / real(nmp_real, PREC)
        pxyzg(1:3) = psumxyz(1:3) / real(nmp_real, PREC)
        
        do imp = 1, nmp_real
           xyz_mp_rep(1:3, imp, irep) = xyz_mp_rep(1:3, imp, irep) - xyzg(1:3)
           pxyz_mp_rep(1:3, imp, irep) = pxyz_mp_rep(1:3, imp, irep) - pxyzg(1:3)

           ! In ND_LANGEVIN, accel_mp stores the previous coordinates (xyz_mp_rep).
           if(i_simulate_type == SIM%ND_LANGEVIN) then
              accel_mp(1:3, imp, irep) = accel_mp(1:3, imp, irep) - xyzg(1:3)
           endif
        end do
        
        if(insimu%i_com_zeroing_ini == 2 .or. insimu%i_com_zeroing == 2) then
           xyz_mp_rep(1:3, 1:nmp_real, irep) = xyz_mp_rep(1:3, 1:nmp_real, irep) + 4500.0
           pxyz_mp_rep(1:3, 1:nmp_real, irep) = pxyz_mp_rep(1:3, 1:nmp_real, irep) + 4500.0

           ! In ND_LANGEVIN, accel_mp stores the previous coordinates (xyz_mp_rep).
           if(i_simulate_type == SIM%ND_LANGEVIN) then
              accel_mp(1:3, 1:nmp_real, irep) = accel_mp(1:3, 1:nmp_real, irep) + 4500.0
           endif
        end if
      
     enddo
  endif

end subroutine simu_xyz_adjst
