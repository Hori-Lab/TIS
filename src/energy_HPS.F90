! energy_HPS

! ************************************************************************
! formula of LJ
! eLJ = coef * {(x0/x)**12 -2*(x0/x)**6}
! ************************************************************************
subroutine energy_HPS(irep, now_HPS, energy_unit, energy)

  use if_util
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inperi, inprotrna
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, &
                          nHPS, iHPS2mp, coef_HPS, HPS_nat2, lambda
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  integer,    intent(out)   :: now_HPS(:,:)
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  integer :: imp1, imp2, iunit, junit
  integer :: ksta, kend
  integer :: iHPS
  real(PREC) :: rcut_off2 !, rcut_off2_pro, rcut_off2_rna
  real(PREC) :: rjudge_contact, rjudge
  real(PREC) :: roverdist2, roverdist6, roverdist12
  real(PREC) :: inv_dist2, inv_dist6
  real(PREC) :: dist2, efull
  real(PREC) :: v21(SDIM)
  real(PREC) :: sigma_6th, sig_over_dist_6th
  real(PREC) :: ehps_tmp
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
!!$omp master

  ! --------------------------------------------------------------------
  ! set parameter 
  rcut_off2 = 1.0e0_PREC / inprotrna%cutoff_HPS**2
  rjudge_contact = 1.2e0_PREC**2

  ! --------------------------------------------------------------------
  ! zero clear
!  now_con(:,:) = 0
      
  ! --------------------------------------------------------------------
#ifdef MPI_PAR3
   klen=(nHPS-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,nHPS)
#else
   ksta = 1
   kend = nHPS
#endif
!$omp do private(imp1, imp2, iunit, junit, v21, dist2, roverdist2, rjudge, efull)
   do iHPS=ksta,kend
   
     imp1 = iHPS2mp(1, iHPS)
     imp2 = iHPS2mp(2, iHPS)

#ifdef _DEBUG_NLOCAL
     write(*,*) 'HPS ', imp1, imp2
#endif
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21)
     end if
     
     sigma_6th = HPS_nat2(iHPS) * HPS_nat2(iHPS) * HPS_nat2(iHPS) 

     dist2 = dot_product(v21,v21)

     roverdist2 = HPS_nat2(iHPS) / dist2
     roverdist6 = roverdist2 * roverdist2 * roverdist2
     roverdist12 = roverdist6 * roverdist6

     if(coef_HPS(iHPS) >= ZERO_JUDGE) then
        now_HPS(2, iHPS) = 1
     else
        now_HPS(2, iHPS) = 0
     end if

     if(roverdist2 < rcut_off2) cycle

     rjudge = HPS_nat2(iHPS) * rjudge_contact
     !  judging contact 
     if(dist2 < rjudge) then
        now_HPS(1, iHPS) = 1
     else
        now_HPS(1, iHPS) = 0
     end if

     inv_dist2       = 1.0_PREC / dist2
     inv_dist6       = inv_dist2 * inv_dist2 * inv_dist2

     sig_over_dist_6th = sigma_6th * inv_dist6

     ehps_tmp = 4.0_PREC * coef_HPS(iHPS) * (roverdist12 - roverdist6)

   !   lambda = HPS_lambda_half(iHPS) + HPS_lambda_half(iHPS)

     ! this judgement criteria is separating the repulsive and attractive part, r < 2^(1/6) * sigma
     if (sig_over_dist_6th >= 0.5_PREC) then
        efull = ehps_tmp + (1.0_PREC - lambda(iHPS)) * coef_HPS(iHPS)
     else
        efull = ehps_tmp + lambda(iHPS) * coef_HPS(iHPS)
     end if

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%HPS) = energy(E_TYPE%HPS) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%HPS) = energy_unit(iunit, junit, E_TYPE%HPS) + efull
  end do
!$omp end do nowait
!!$omp end master

end subroutine energy_HPS
