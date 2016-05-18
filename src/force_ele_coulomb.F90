! force_ele_coulomb
!> @brief This subroutine calculates the force of electrostatic interaction.

subroutine force_ele_coulomb(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inele, inperi
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, lele, iele2mp, coef_ele, nmp_all
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SDIM, nmp_all)

  integer :: imp1, imp2
  integer :: ksta, kend
  integer :: iele
  integer :: imirror
  real(PREC) :: dist1, dist2
  real(PREC) :: dv_dr, cutoff2
  real(PREC) :: v21(SDIM), for(SDIM)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! --------------------------------------------------------------------
#ifdef _DEBUG
  write(*,*) '#### start force_ele'
#endif

  ! for speed up
  cutoff2 = inele%cutoff_ele ** 2

#ifdef _DEBUG
  write(*,*) 'lele(irep), ',lele(irep)
#endif
#ifdef MPI_PAR
#ifdef SHARE_NEIGH
  klen=(lele(irep)-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,lele(irep))
#else
  ksta = 1
  kend = lele(irep)
#endif
#ifdef MPI_DEBUG
  print *,"exv2_6       = ", kend-ksta+1
#endif

#else
  ksta = 1
  kend = lele(irep)
#endif

!$omp do private(imp1,imp2,v21,dist2,dist1,dv_dr,for,imirror)
  do iele=ksta, kend
     imp1 = iele2mp(1, iele, irep)
     imp2 = iele2mp(2, iele, irep)

     if(inperi%i_periodic == 0) then
        v21(1:SDIM) = xyz_mp_rep(1:SDIM, imp2, irep) - xyz_mp_rep(1:SDIM, imp1, irep)
     else
        imirror = iele2mp(3, iele, irep)
        v21(1:SDIM) = pxyz_mp_rep(1:SDIM, imp2, irep) - pxyz_mp_rep(1:SDIM, imp1, irep) + inperi%d_mirror(1:SDIM, imirror)
     end if

     dist2 = dot_product(v21,v21)
     if(dist2 > cutoff2) cycle
   
     ! -----------------------------------------------------------------
     dist1 = sqrt(dist2)
     dv_dr = coef_ele(iele, irep) / (dist1 * dist2)
     
     for(1:SDIM) = dv_dr * v21(1:SDIM)
     force_mp(1:SDIM, imp1) = force_mp(1:SDIM, imp1) - for(1:SDIM)
     force_mp(1:SDIM, imp2) = force_mp(1:SDIM, imp2) + for(1:SDIM)
  end do
!$omp end do nowait

end subroutine force_ele_coulomb
