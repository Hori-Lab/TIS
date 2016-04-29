! energy_ele_coulomb
!> @brief Calculate the energy of electrostatic interaction 

subroutine energy_ele_coulomb(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inmisc, inele
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, lele, iele2mp, coef_ele
  use mpiconst

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit, iele, imirror
  real(PREC) :: dist1, dist2
  real(PREC) :: ene, cutoff2
  real(PREC) :: v21(SPACE_DIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  cutoff2 = inele%cutoff_ele ** 2

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH
  klen=(lele(irep)-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,lele(irep))
#else
  ksta = 1
  kend = lele(irep)
#endif
#else
  ksta = 1
  kend = lele(irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,dist1,&
!$omp&           ene,iunit,junit,imirror)
  do iele=ksta, kend

     imp1 = iele2mp(1, iele, irep)
     imp2 = iele2mp(2, iele, irep)
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = iele2mp(3, iele, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
     dist2 = dot_product(v21,v21)
     if(dist2 > cutoff2) cycle
        
     ! --------------------------------------------------------------------
     dist1 = sqrt(dist2)
     
     ene = coef_ele(iele,irep)/dist1

     energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene
     
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + ene
  end do
!$omp end do nowait

end subroutine energy_ele_coulomb
