! simu_energy_ele
!> @brief Calculate the energy of electrostatic interaction 

subroutine simu_energy_ele_coulomb(irep, pnlet, pnle_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inmisc, inele, inion, insimu
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, &
                          lele, iele2mp, coef_ele, iontype_mp, imp2type
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
  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit, iele1, imirror
  real(PREC) :: dist1, dist2
  real(PREC) :: pnl, cutoff2
  real(PREC) :: rk
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
!$omp&           pnl,iunit,junit,imirror)
  do iele1=ksta, kend

     imp1 = iele2mp(1, iele1, irep)
     imp2 = iele2mp(2, iele1, irep)
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = iele2mp(3, iele1, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
     dist2 = dot_product(v21,v21)
     if(dist2 > cutoff2) cycle
        
     ! --------------------------------------------------------------------
     dist1 = sqrt(dist2)
     
     if (inmisc%i_temp_independent == 0) then
        pnl = coef_ele(iele1,irep)/dist1

     else if (inmisc%i_temp_independent == 1) then
     !   pnl = coef_ele(iele1,irep) * (1.0e0_PREC / dist1 - 0.5e0_PREC * rcdist) &
     !        * exp(-dist1*rcdist)
        write(*,*) 'error in simu_energy_ele_coulomb'
        stop

     else if (inmisc%i_temp_independent == 2) then
     !   rk = dist1 * rcdist
     !   pnl = coef_ele(iele1,irep) / dist1 * exp(-rk)    &
     !        * ( - (1.0e0_PREC + 0.5e0_PREC*rk) * inele%diele_dTcoef  )
        write(*,*) 'error in simu_energy_ele_coulomb'
        stop
     endif

     pnlet(E_TYPE%ELE) = pnlet(E_TYPE%ELE) + pnl
     
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     pnle_unit(iunit, junit, E_TYPE%ELE) = pnle_unit(iunit, junit, E_TYPE%ELE) + pnl
  end do
!$omp end do nowait

end subroutine simu_energy_ele_coulomb
