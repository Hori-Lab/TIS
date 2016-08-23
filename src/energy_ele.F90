! energy_ele
!> @brief Calculate the energy of electrostatic interaction 

subroutine energy_ele(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inmisc, inele, inperi
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, lele, iele2mp, coef_ele
  use var_replica, only : irep2grep
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit, grep, iele1, imirror
  real(PREC) :: dist1, dist2, ew, ek, rk
  real(PREC) :: ene, rcdist, cutoff2
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  grep = irep2grep(irep)
  cutoff2 = (inele%cutoff_ele * inele%cdist(grep))**2
  rcdist = 1.0 / inele%cdist(grep)

  ew = inele%diele_water
  ek = inele%diele

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
!$omp do private(imp1,imp2,v21,dist2,dist1, &
!$omp&           ene,iunit,junit,imirror)
  do iele1=ksta, kend

     imp1 = iele2mp(1, iele1, irep)
     imp2 = iele2mp(2, iele1, irep)
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = iele2mp(3, iele1, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
     ! v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
     if(dist2 > cutoff2) cycle
        
     ! --------------------------------------------------------------------
     dist1 = sqrt(dist2)
     
     if (inmisc%i_temp_independent == 0) then
        ene = coef_ele(iele1,irep)/dist1*exp(-dist1*rcdist)

     else if (inmisc%i_temp_independent == 1) then
        ene = coef_ele(iele1,irep) * (1.0e0_PREC / dist1 - 0.5e0_PREC * rcdist) &
             * exp(-dist1*rcdist)

     else if (inmisc%i_temp_independent == 2) then
        rk = dist1 * rcdist
        ene = coef_ele(iele1,irep) / dist1 * exp(-rk)    &
             * ( - (1.0e0_PREC + 0.5e0_PREC*rk) * inele%diele_dTcoef  )

     else
        ene = 0.0  ! to suppress compiler message
     endif
     
     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene
     
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + ene
  end do
!$omp end do nowait

end subroutine energy_ele
