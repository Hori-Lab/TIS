subroutine  simu_energy_ele(irep, pnlet, pnle_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inmisc, inele
  use var_struct,  only : nunit_real, lunit2mp, imp2unit, xyz_mp_rep, &
                          lele, iele2mp, coef_ele
  use var_replica, only : inrep, irep2grep
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: pnlet(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: pnle_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  ! local variables
  integer :: imp1, imp2, iunit, junit, grep, iele1
  integer :: imp3, imp4
  integer :: ielem
  real(PREC) :: dist1, dist2
  real(PREC) :: pnl, rcdist, cutoff2
  real(PREC) :: dist11, dist12, dist21, dist22
  real(PREC) :: v21(SPACE_DIM)
  real(PREC) :: v22(SPACE_DIM)
  real(PREC) :: pnlet_l
#ifdef MPI_PAR3
  integer :: klen, ksta, kend
#endif

  if (inmisc%force_flag(INTERACT%ELE)) then

     grep = irep2grep(irep)
     cutoff2 = (inele%cutoff_ele * inele%cdist(grep))**2
     rcdist = 1.0 / inele%cdist(grep)
#ifdef MPI_PAR3

#ifdef SHARE_NEIGH
     klen=(lele(irep)-1+npar_mpi)/npar_mpi
     ksta=1+klen*local_rank_mpi
     kend=min(ksta+klen-1,lele(irep))

     ielem = mod(kend-ksta+1, 2)
     do iele1=ksta, ksta+ielem-1
#else
     do iele1 = 1, lele(irep)
#endif

#else
     ielem = mod(lele(irep), 2)
     do iele1 = 1, ielem
#endif
        imp1 = iele2mp(1, iele1, irep)
        imp2 = iele2mp(2, iele1, irep)
        
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
        if(dist2 > cutoff2) cycle
           
        ! --------------------------------------------------------------------
        dist1 = sqrt(dist2)
        pnl = coef_ele(iele1,irep) / dist1 * exp(-dist1 * rcdist)
   
        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%ELE) = pnlet(E_TYPE%ELE) + pnl
   
        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%ELE) = pnle_unit(iunit, junit, E_TYPE%ELE) + pnl
     end do

     pnlet_l = pnlet(E_TYPE%ELE)
#ifdef MPI_PAR3

#ifdef SHARE_NEIGH
     do iele1=ksta+ielem, kend, 2
#else
     do iele1=ielem+1, lele(irep), 2
#endif

#else
     do iele1 =ielem+1, lele(irep), 2
#endif
        imp1 = iele2mp(1, iele1, irep)
        imp2 = iele2mp(2, iele1, irep)

        imp3 = iele2mp(1, iele1+1, irep)
        imp4 = iele2mp(2, iele1+1, irep)
        
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        v22(1:3) = xyz_mp_rep(1:3, imp4, irep) - xyz_mp_rep(1:3, imp3, irep)

        dist21 = v21(1)**2 + v21(2)**2 + v21(3)**2
        dist22 = v22(1)**2 + v22(2)**2 + v22(3)**2

        if(dist21 > cutoff2) goto 200
           
        ! --------------------------------------------------------------------
        dist11 = sqrt(dist21)
        pnl = coef_ele(iele1,irep) / dist11 * exp(-dist11 * rcdist)
   
        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet_l = pnlet_l + pnl
   
        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%ELE) = pnle_unit(iunit, junit, E_TYPE%ELE) + pnl

 200    continue

        if(dist22 > cutoff2) goto 210
           
        ! --------------------------------------------------------------------
        dist12 = sqrt(dist22)
        pnl = coef_ele(iele1+1,irep) / dist12 * exp(-dist12 * rcdist)
   
        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet_l = pnlet_l + pnl
   
        iunit = imp2unit(imp3)
        junit = imp2unit(imp4)
        pnle_unit(iunit, junit, E_TYPE%ELE) = pnle_unit(iunit, junit, E_TYPE%ELE) + pnl

 210    continue

     end do
     pnlet(E_TYPE%ELE) = pnlet_l

  endif  ! (inmisc%force_flag(INTERACT%ELE))

end subroutine simu_energy_ele
