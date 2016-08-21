! simu_xyz_adjst
!> @brief Correct coordinates to move the system to the origin

subroutine simu_xyz_adjst()

  use const_maxsize
  use const_index
  use const_physical
  use var_setp,    only : insimu, inperi
  use var_struct,  only : nmp_real, xyz_mp_rep, pxyz_mp_rep, imp2type
  use var_replica, only : n_replica_mpi
  use mpiconst

  implicit none

  integer :: i, imp, irep, n_mp
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
           if (imp2type(imp) == MPTYPE%ION_MG .or. imp2type(imp) == MPTYPE%ION_NA .or.&
               imp2type(imp) == MPTYPE%ION_K  .or. imp2type(imp) == MPTYPE%ION_CL) then
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
        call util_periodic(irep)
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
        end do
      
        if(insimu%i_com_zeroing_ini == 2 .or. insimu%i_com_zeroing == 2) then
           xyz_mp_rep(1:3, 1:nmp_real, irep) = xyz_mp_rep(1:3, 1:nmp_real, irep) + 4500.0
           pxyz_mp_rep(1:3, 1:nmp_real, irep) = pxyz_mp_rep(1:3, 1:nmp_real, irep) + 4500.0
        end if
      
     enddo
  endif

end subroutine simu_xyz_adjst
