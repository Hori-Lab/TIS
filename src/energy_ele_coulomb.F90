! energy_ele_coulomb
!> @brief Calculate the energy of electrostatic interaction 

subroutine energy_ele_coulomb(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inele, inperi
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
  real(PREC) :: v21(SDIM)
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
!$omp do private(imp1,imp2,v21,dist2,dist1,ene,iunit,junit,imirror)
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


subroutine energy_ele_coulomb_tp(irep, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inele, inperi
  use var_struct,  only : pxyz_mp_rep, nmp_real, lmp2charge, coef_charge, &
                          ntp, xyz_tp, charge_tp
  use var_replica, only : irep2grep
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(E_TYPE%MAX)

  integer :: itp1, itp2, imp2, imirror, grep
  real(PREC) :: q1, dist2, ene
  real(PREC) :: cutoff2
  real(PREC) :: v21(SDIM), vx(SDIM)

  cutoff2 = inele%cutoff_ele ** 2
  grep = irep2grep(irep)

  ene = 0.0

!$omp do private(itp1,imp2,itp2,v21,vx,dist2,imirror)
  do itp1 = 1, ntp

     q1 = charge_tp(itp1)

     do imp2 = 1, nmp_real
   
        !if(inperi%i_periodic == 0) then
        !   v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_tp(1:3, itp1)
        !else
           vx(1:3) = pxyz_mp_rep(1:3, imp2, irep) - xyz_tp(1:3, itp1)
           call util_pbneighbor(vx, imirror)
           v21(1:3) = vx(1:3) + inperi%d_mirror(1:3, imirror)
        !end if
        
        dist2 = dot_product(v21,v21)
        if(dist2 > cutoff2) cycle
           
        ene = ene + q1 * coef_charge( lmp2charge(imp2),grep ) / sqrt(dist2)
     end do

     do itp2 = itp1+1, ntp
        !if(inperi%i_periodic == 0) then
        !   v21(1:3) = xyz_tp(1:3, itp2) - xyz_tp(1:3, itp1)
        !else
           vx(1:3) = xyz_tp(1:3, itp2) - xyz_tp(1:3, itp1)
           call util_pbneighbor(vx, imirror)
           v21(1:3) = vx(1:3) + inperi%d_mirror(1:3, imirror)
        !end if
        
        dist2 = dot_product(v21,v21)
        if(dist2 > cutoff2) cycle
           
        ene = ene + q1 * charge_tp(itp2) / sqrt(dist2)
     enddo
  enddo
!$omp end do nowait

  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + inele%coef(grep) * ene

end subroutine energy_ele_coulomb_tp
