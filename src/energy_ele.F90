! energy_ele
!> @brief Calculate the energy of electrostatic interaction 

subroutine energy_ele(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inmisc, inele, inion
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, lele, iele2mp, coef_ele, iontype_mp, imp2type
  use var_replica, only : irep2grep
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit, grep, iele1, imirror
  integer :: itype1, itype2
  real(PREC) :: dist1, dist2, ew, ek, rk
  real(PREC) :: exv, rcdist, cutoff2, xtanh, rsig, rek_corr
  real(PREC) :: v21(SPACE_DIM)
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
!$omp do private(imp1,imp2,v21,dist2,dist1,itype1,itype2, &
!$omp&           rsig,xtanh,rek_corr,exv,iunit,junit,imirror)
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
     
     itype1 = iontype_mp(imp1)
     itype2 = iontype_mp(imp2)
     
     if((.not. inmisc%class_flag(CLASS%ION)) .or. itype1 <= 0 .or. itype1 > IONTYPE%MAX_ALL .or. itype2 <= 0 .or. itype2 > IONTYPE%MAX_ALL) then

        if (inmisc%i_temp_independent == 0) then
           exv = coef_ele(iele1,irep)/dist1*exp(-dist1*rcdist)

        else if (inmisc%i_temp_independent == 1) then
           exv = coef_ele(iele1,irep) * (1.0e0_PREC / dist1 - 0.5e0_PREC * rcdist) &
                * exp(-dist1*rcdist)

        else if (inmisc%i_temp_independent == 2) then
           rk = dist1 * rcdist
           exv = coef_ele(iele1,irep) / dist1 * exp(-rk)    &
                * ( - (1.0e0_PREC + 0.5e0_PREC*rk) * inele%diele_dTcoef  )

        endif

     else
        rsig = 1.0/inion%csigmame(itype1, itype2)
        xtanh = tanh((dist1-inion%cdistme(itype1, itype2))*rsig)
        rek_corr = 1.0/(0.5*(ew + 5.2) + 0.5*(ew - 5.2)*xtanh)
        
        exv = ek*rek_corr*coef_ele(iele1, irep)/dist1*exp(-dist1*rcdist) 
     end if
     
     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + exv
     
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + exv
  end do
!$omp end do nowait

end subroutine energy_ele
