! simu_energy_ion
!> @brief Calculate all the energy terms related to ion, &
!>        including LJ, hydration term, excluded volume term

subroutine  simu_energy_ion(irep, e_exv_unit, e_exv)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inmisc, inion
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, &
                          lexv, iexv2mp, iontype_mp
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout)   :: e_exv_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)
  real(PREC), intent(inout)   :: e_exv(:)         ! (E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit, itype1, itype2
  integer :: iexv, imirror
  real(PREC) :: dist1, dist2, cdist2, ddist1, ddist2
  real(PREC) :: roverdist2, roverdist6, roverdist12
  real(PREC) :: exv, coef, cutoff2
  real(PREC) :: v21(SPACE_DIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ------------------------------------------------------------------------

  if (inmisc%class_flag(CLASS%ION)) then

     ! ------------------------------------------------------------------------
     ! LJ and hydration interaction of ion
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
     klen=(lexv(2,E_TYPE%HYD_ION,irep)-lexv(1,E_TYPE%HYD_ION,irep)+npar_mpi)/npar_mpi
     ksta=lexv(1,E_TYPE%HYD_ION,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lexv(2,E_TYPE%HYD_ION,irep))
#else
     ksta = lexv(1, E_TYPE%HYD_ION, irep)
     kend = lexv(2, E_TYPE%HYD_ION, irep)
#endif
#else
     ksta = lexv(1, E_TYPE%HYD_ION, irep)
     kend = lexv(2, E_TYPE%HYD_ION, irep)
#endif
!$omp do private(imp1,imp2,itype1,itype2,v21,dist2,dist1,&
!$omp&           ddist1,ddist2,roverdist2, &
!$omp&           roverdist6,roverdist12,exv,&
!$omp&           iunit,junit,imirror)
     do iexv=ksta, kend

        imp1 = iexv2mp(1, iexv, irep)
        imp2 = iexv2mp(2, iexv, irep)
        itype1 = iontype_mp(imp1)
        itype2 = iontype_mp(imp2)
     
        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           imirror = iexv2mp(3, iexv, irep)
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

        if(dist2 > inion%cutofflj2(itype1, itype2)) cycle
                                                                
        ! --------------------------------------------------------------------
        roverdist2 = inion%cdistlj2(itype1, itype2) / dist2
        roverdist6 = roverdist2 * roverdist2 * roverdist2
        roverdist12 = roverdist6 * roverdist6

        dist1 = sqrt(dist2)
        ddist1 = dist1 - inion%cdistmh1(itype1, itype2)
        ddist2 = dist1 - inion%cdistmh2(itype1, itype2)

        exv = inion%clj_energy(itype1, itype2) * (roverdist12 - roverdist6) &
             + inion%cmh1_energy(itype1, itype2) * exp(-inion%rsigmamh1(itype1, itype2) * ddist1**2) &
             + inion%cmh2_energy(itype1, itype2) * exp(-inion%rsigmamh2(itype1, itype2) * ddist2**2)

        ! --------------------------------------------------------------------
        ! sum of the energy
        e_exv(E_TYPE%HYD_ION) = e_exv(E_TYPE%HYD_ION) + exv

        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        e_exv_unit(iunit, junit, E_TYPE%HYD_ION) = e_exv_unit(iunit, junit, E_TYPE%HYD_ION) + exv

     end do
!$omp end do nowait

     ! ------------------------------------------------------------------------
     ! exvol ion
     ! for speed up
     cdist2 = inion%cdist_exv_ion**2
     cutoff2 = cdist2*inion%cutoff_exv_ion**2
     coef = 4.0e0_PREC * inion%cexv_ion

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
     klen=(lexv(2,E_TYPE%EXV_ION,irep)-lexv(1,E_TYPE%EXV_ION,irep)+npar_mpi)/npar_mpi
     ksta=lexv(1,E_TYPE%EXV_ION,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lexv(2,E_TYPE%EXV_ION,irep))
#else
     ksta = lexv(1, E_TYPE%EXV_ION, irep)
     kend = lexv(2, E_TYPE%EXV_ION, irep)
#endif
#else
     ksta = lexv(1, E_TYPE%EXV_ION, irep)
     kend = lexv(2, E_TYPE%EXV_ION, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2, &
!$omp&           roverdist6,roverdist12,exv,iunit,junit,imirror)
     do iexv=ksta, kend

        imp1 = iexv2mp(1, iexv, irep)
        imp2 = iexv2mp(2, iexv, irep)
     
        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           imirror = iexv2mp(3, iexv, irep)
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if
        
!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

        if(dist2 > cutoff2) cycle
                                                                
        ! --------------------------------------------------------------------
        roverdist2 = cdist2 / dist2
        roverdist6 = roverdist2 * roverdist2 * roverdist2
        roverdist12 = roverdist6 * roverdist6
        exv = coef * (roverdist12 - roverdist6 + 0.25e0_PREC)

        ! --------------------------------------------------------------------
        ! sum of the energy
        e_exv(E_TYPE%EXV_ION) = e_exv(E_TYPE%EXV_ION) + exv

        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        e_exv_unit(iunit, junit, E_TYPE%EXV_ION) = e_exv_unit(iunit, junit, E_TYPE%EXV_ION) + exv
     end do
!$omp end do nowait

  endif  ! inmisc%class_flag(CLASS%ION)

end subroutine simu_energy_ion
