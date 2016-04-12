!simu_energy_fbangle
!> @brief Calculates the energy related to flexible bond-angle term. &
!>        Values are added into "pnlet(E_TYPE%BANGLE)" and           &
!>        "pnle_unit(,,E_TYPE%BANGLE)"

subroutine simu_energy_fbangle(irep, pnle_unit, pnlet)

  use const_maxsize
  use const_index
  use const_physical
  use var_struct, only : nfba, ifba2mp, xyz_mp_rep, imp2unit, nunit_all, fba_ener_corr, &
                         fba_max_th, fba_min_th, fba_max_th_ener, fba_min_th_ener
  use var_setp,   only : inflp

#ifdef MPI_PAR3
  use mpiconst
#endif
  
  implicit none

  ! ---------------------------------------------------------------------------
  integer, intent(in) :: irep
  real(PREC), intent(inout) :: pnlet(E_TYPE%MAX) 
  real(PREC), intent(inout) :: pnle_unit(nunit_all, nunit_all, E_TYPE%MAX) 
  
  ! ---------------------------------------------------------------------------
  ! local variables
  integer    :: ksta, kend
  integer    :: imp(3)
  integer    :: ifba, i
  integer    :: iunit, junit
  real(PREC) :: co_theta, th
  real(PREC) :: efull
#ifdef MPI_PAR3
  integer    :: klen
#endif

#ifdef MPI_PAR3
  klen=(nfba-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nfba)
#else
  ksta = 1
  kend = nfba
#endif
!$omp do private(i, imp, co_theta, th, efull, iunit,junit)
  do ifba = ksta, kend

     do i = 1, 3
        imp(i) = ifba2mp(i, ifba)
     end do

     call util_bondangle( imp(1), imp(2), imp(3), co_theta, xyz_mp_rep(:, :, irep) )
     
     if (co_theta > 1.0e0_PREC) then
        co_theta = 1.0e0_PREC
     else if (co_theta < -1.0) then
        co_theta = -1.0e0_PREC
     end if
     th = acos(co_theta)

     ! Calculate energy
     if ( th < fba_min_th(ifba) ) then
        efull = FBA_MIN_ANG_FORCE * th + fba_min_th_ener(ifba) - FBA_MIN_ANG_FORCE * fba_min_th(ifba)
     else if ( th > fba_max_th(ifba) ) then
        efull = FBA_MAX_ANG_FORCE * th + fba_max_th_ener(ifba) - FBA_MAX_ANG_FORCE * fba_max_th(ifba)
     else
        call splint(th, ifba, efull)
     end if
     
     efull = inflp%k_ang * ( efull - fba_ener_corr(ifba) )

     ! -------------------------------------------------------------------
     ! sum of the energy
     pnlet(E_TYPE%BANGLE) = pnlet(E_TYPE%BANGLE) + efull

     iunit = imp2unit( imp(1) )
     junit = imp2unit( imp(2) )
     
     pnle_unit(iunit, junit, E_TYPE%BANGLE) = pnle_unit(iunit, junit, E_TYPE%BANGLE) + efull
     
  end do
!$omp end do nowait
  
end subroutine simu_energy_fbangle
  
  
subroutine splint(theta, ifba, energy)
  
  use const_maxsize
  use const_physical
  use var_struct, only : fba_para_x, fba_para_y, fba_para_y2 
    
  implicit none
  
  ! ----------------------------------------------------------------------
  real(PREC), intent(in) :: theta
  integer, intent(in) :: ifba
  real(PREC), intent(inout) :: energy
  
  ! ----------------------------------------------------------------------
  ! local variables
  integer :: n, klo, khi, k
  real(PREC) :: h, a, b
  ! ----------------------------------------------------------------------

  n = size(fba_para_x(:,1)) ! Get the number of bin
  
  klo = 1
  khi = n 

  do while (khi - klo > 1) 
     k = (khi+klo) / 2
     if ( fba_para_x(k, ifba) > theta) then
        khi = k;
     else
        klo = k;
     end if
  end do

  h = fba_para_x(khi, ifba) - fba_para_x(klo, ifba)

  a = ( fba_para_x(khi, ifba) - theta ) / h
  b = ( theta - fba_para_x(klo, ifba) ) / h

  energy = a * fba_para_y(klo, ifba) + b * fba_para_y(khi, ifba) &
       + ( (a**3 - a) * fba_para_y2(klo, ifba) &
       +   (b**3 - b) * fba_para_y2(khi, ifba) ) * h * h / 6.0e0_PREC

end subroutine splint

