! energy_ele_DH
!> @brief Calculate the energy of electrostatic interaction 

subroutine energy_ele_DH(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inmisc, inele, inperi, inpmf
  use var_struct,  only : imp2unit, pxyz_mp_rep, lele, iele2mp, coef_ele
  use var_replica, only : irep2grep
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit, grep, iele1, imir, ipmf
  real(PREC) :: dist1, dist2, ew, ek, rk
  real(PREC) :: ene, rcdist, cutoff2
  real(PREC) :: v21(SDIM)
  real(PREC) :: Rmin, Rbininv
#ifdef MPI_PAR3
  integer :: klen
#endif

  grep = irep2grep(irep)
  if (inele%i_DH_cutoff_type == 0) then
     cutoff2 = (inele%cutoff_ele * inele%cdist(grep))**2
  else
     cutoff2 = inele%cutoff_ele**2
  endif
  rcdist = 1.0 / inele%cdist(grep)

  ew = inele%diele_water
  ek = inele%diele

  if (inmisc%i_dtrna_model == 2019) then
     Rmin = inpmf%Rmin(PMFTYPE%MG_P)
     Rbininv = 1.0 / inpmf%Rbin(PMFTYPE%MG_P)
  endif

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
!$omp&           ene,iunit,junit,imir)
  do iele1=ksta, kend

     imp1 = iele2mp(1, iele1, irep)
     imp2 = iele2mp(2, iele1, irep)
     imir = iele2mp(3, iele1, irep)
     ipmf = iele2mp(4, iele1, irep)

     v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imir)

     dist2 = dot_product(v21,v21)
     if(dist2 > cutoff2) cycle
        
     ! --------------------------------------------------------------------
     dist1 = sqrt(dist2)
     
     ! PMF + DH (semiexplicit model)
     if (ipmf > 0 .and. dist1 <= inpmf%Rmax(ipmf)) then 

        ene = energy_pmfdh(dist1, grep)

     ! DH
     else
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

     endif
     
     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene
     
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + ene
  end do
!$omp end do nowait

contains
   
   real(PREC) function energy_pmfdh(r, grep)
      use var_simu, only : pmfdh_energy, pmfdh_force
      implicit none
      real(PREC),intent(in) :: r
      integer, intent(in) :: grep
      integer :: ibin

      ! Rmin = 2.6
      ! Rbin = 0.025
      !
      ! r = 2.624 ==> bin 1
      ! r = 2.625 ==> bin 2
      ! r = 2.626 ==> bin 2
      ! r = 2.649 ==> bin 2
      ! r = 2.650 ==> bin 3

      ! If r = 2.624
      ! (r - Rmin) = 0.024
      ! (r - Rmin) * Rbininv = 0.96  ==floor==>  0
      ! floor((r - Rmin) * Rbininv) + 1 = 1

      ! If r = 2.625
      ! (r - Rmin) = 0.025
      ! (r - Rmin) * Rbininv = 1.0  ==floor==>  1
      ! floor((r - Rmin) * Rbininv) + 1 = 2

      ! If r = 2.626
      ! (r - Rmin) = 0.026
      ! (r - Rmin) * Rbininv = 1.04  ==floor==>  1
      ! floor((r - Rmin) * Rbininv) + 1 = 2

      ! If r = 2.649
      ! (r - Rmin) = 0.049
      ! (r - Rmin) * Rbininv = 1.96 ==floor==>  1
      ! floor((r - Rmin) * Rbininv) + 1 = 2

      ! If r = 2.650
      ! (r - Rmin) = 0.050
      ! (r - Rmin) * Rbininv = 2.0 ==floor==>  2
      ! floor((r - Rmin) * Rbininv) + 1 = 3

      ibin = floor( (r - Rmin) * Rbininv) + 1

      energy_pmfdh = pmfdh_energy(ibin, grep, PMFTYPE%MG_P) &
                   - pmfdh_force(ibin,grep,PMFTYPE%MG_P) * (r - (Rmin + real(ibin-1,kind=PREC)/Rbininv))

   endfunction energy_pmfdh

end subroutine energy_ele_DH
