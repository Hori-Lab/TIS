! simu_force_ele
!> @brief This subroutine calculates the force of electrostatic interaction.

subroutine simu_force_ele(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : outfile, inperi
  use var_setp,   only : inmisc, inele, inion
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, &
                         lele, iele2mp, coef_ele, nmp_all, iontype_mp, imp2type
  use var_replica,only : irep2grep
  use var_simu, only   : istep
  use mpiconst

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

  ! --------------------------------------------------------------------
  ! local variables                                                             
  integer :: imp1, imp2
  integer :: ksta, kend
  integer :: grep, iele1, itype1, itype2
  integer :: imptype1, imptype2
  integer :: imirror
  real(PREC) :: dist1, dist2, rdist1, xtanh, rsig, rek_corr
  real(PREC) :: dvdw_dr, rcdist, cutoff2
  real(PREC) :: ek, ek_simu
  real(PREC) :: v21(3), for(3)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! --------------------------------------------------------------------
#ifdef _DEBUG
#ifdef MPI_PAR
  call MPI_barrier(MPI_COMM_local,ierr)
#endif
  write(*,*) '#### start simu_force_ele'
  call flush(6)
#endif

  ! for speed up
  grep = irep2grep(irep)
  cutoff2 = (inele%cutoff_ele * inele%cdist(grep))**2
  rcdist = 1.0e0_PREC / inele%cdist(grep)

  ek_simu = inele%diele
  ek = inele%diele_water

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
  print *,"pnl2_6       = ", kend-ksta+1
#endif

#else
  ksta = 1
  kend = lele(irep)
#endif
  write(*,*) 'ksta,kend=',ksta,kend

!$omp do private(imp1,imp2,v21,dist2,dist1,rdist1,itype1,itype2, imptype1, imptype2,&
!$omp&           rsig,xtanh,rek_corr,dvdw_dr,for,imirror)
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
   
     ! -----------------------------------------------------------------
     dist1 = sqrt(dist2)
     rdist1 = 1.0 / dist1
     
     itype1 = iontype_mp(imp1)
     itype2 = iontype_mp(imp2)
     
     if((.not. inmisc%class_flag(CLASS%ION)) .or. itype1 <= 0 .or. itype1 > IONTYPE%MAX_ALL .or. itype2 <= 0 .or. itype2 > IONTYPE%MAX_ALL) then
        dvdw_dr = coef_ele(iele1, irep) * rdist1 * rdist1 * &
             (rdist1 + rcdist) * &
             exp(-dist1 * rcdist)
     else
        rsig = 1.0/inion%csigmame(itype1, itype2)
        xtanh = tanh((dist1-inion%cdistme(itype1, itype2))*rsig)
        rek_corr = 1.0/(0.5*(ek + 5.2) + 0.5*(ek - 5.2)*xtanh)
        
        dvdw_dr = ek_simu*rek_corr*coef_ele(iele1, irep)*rdist1*rdist1 &
             *(rdist1 + rcdist &
             + rek_corr*rsig*0.5*(ek - 5.2)*(1.0 - xtanh**2)) * &
             exp(-dist1 * rcdist)
     end if
     
     if(dvdw_dr > DE_MAX) then
        ! write (*, *) "electrostatic interaction", istep, imp1, imp2, dist1, dvdw_dr, DE_MAX
        dvdw_dr = DE_MAX
     end if
     ! if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
        
     imptype1 = imp2type(imp1)
     imptype2 = imp2type(imp2)

     for(1:3) = dvdw_dr * v21(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
  end do
!$omp end do nowait
#ifdef _DEBUG
  write(*,*) '#### end simu_force_ele'
  call flush(6)
#ifdef MPI_PAR
  call MPI_barrier(MPI_COMM_local,ierr)
#endif
#endif

end subroutine simu_force_ele
