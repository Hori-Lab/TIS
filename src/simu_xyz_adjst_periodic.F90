! simu_xyz_adjst
!> @brief Correct coordinates to move the system to the origin

subroutine simu_xyz_adjst_periodic()

  use const_maxsize
  use const_physical
  use var_setp,    only : inperi
  use var_struct,  only : nmp_real, pxyz_mp_rep
  use var_replica, only : n_replica_mpi
  use mpiconst

  implicit none

  integer :: imp, irep, n_mp
  real(PREC) :: theta
  real(PREC) :: pxyzg(SDIM)
  real(PREC) :: s_cos(SDIM), s_sin(SDIM)

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

end subroutine simu_xyz_adjst_periodic
