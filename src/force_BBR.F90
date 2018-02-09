subroutine force_BBR(irep, force_mp)

  use if_util
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : indtrna15, inperi
  use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep, nbbr, ibbr2mp
  use mpiconst
  use var_simu, only : istep

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SDIM, nmp_all)

  integer :: ksta, kend
  integer :: imp1, imp2, jmp1, jmp2
  integer :: ibbr
  real(PREC) :: dist2
  real(PREC) :: coef, cdist2, cutoff2
  real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist14
  real(PREC) :: dvdw_dr
  real(PREC) :: vi(SDIM), vj(SDIM), vij(SDIM), for(SDIM)
#ifdef MPI_PAR
  integer :: klen
#endif
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  cutoff2 = indtrna15%bbr_cutoff ** 2
  cdist2 = indtrna15%bbr_dist ** 2
  coef = 6.0e0_PREC * indtrna15%bbr_eps / cdist2

#ifdef MPI_PAR
  klen=(nbbr(irep)-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nbbr(irep))
#else
  ksta = 1
  kend = nbbr(irep)
#endif

!$omp do private(imp1,imp2,jmp1,jmp2,vij,vi,vj,dist2,roverdist2,roverdist4,roverdist8,roverdist14,dvdw_dr,for)
  do ibbr=ksta, kend

     imp1 = ibbr2mp(1, ibbr, irep)
     imp2 = ibbr2mp(2, ibbr, irep)
     jmp1 = ibbr2mp(3, ibbr, irep)
     jmp2 = ibbr2mp(4, ibbr, irep)

     if(inperi%i_periodic == 0) then
        !   0.5 * (imp1 + imp2) - 0.5 * (jmp1 + jmp2)
        ! = 0.5 * {(imp1 + imp2) - (jmp1 + jmp2)}
        vi(1:3) = xyz_mp_rep(1:3, imp1, irep) + xyz_mp_rep(1:3,imp2,irep)
        vj(1:3) = xyz_mp_rep(1:3, jmp1, irep) + xyz_mp_rep(1:3,jmp2,irep)
        vij(1:3) = 0.5 * (vi(1:3) - vj(1:3))
     else
        !   0.5 * (imp1 + imp2) - 0.5 * (jmp1 + jmp2)
        !  (Because we need to use util_pbneighbor)
        vi(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(vi)
        vi(1:3) = xyz_mp_rep(1:3, imp1, irep) + 0.5 * vi(1:3)
        call util_pbneighbor(vi)

        vj(1:3) = xyz_mp_rep(1:3, jmp2, irep) - xyz_mp_rep(1:3, jmp1, irep)
        call util_pbneighbor(vj)
        vj(1:3) = xyz_mp_rep(1:3, jmp1, irep) + 0.5 * vj(1:3)
        call util_pbneighbor(vj)

        vij(1:3) = vi(1:3) - vj(1:3)
        call util_pbneighbor(vij)
     end if
     
     dist2 = dot_product(vij,vij)

     if(dist2 > cutoff2) cycle
     
     ! -----------------------------------------------------------------
     roverdist2 = cdist2 / dist2
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist14 = roverdist2 * roverdist4 * roverdist8
     dvdw_dr = coef * roverdist14
     !write(*,*) 'ForceBBR', imp1, imp2, jmp1, jmp2, sqrt(dist2), dvdw_dr
     if(dvdw_dr > 500.0) then
        write(error_message,*) 'BBR force > 500.0', imp1, imp2, jmp1, jmp2, sqrt(dist2), dvdw_dr
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     for(1:3) = dvdw_dr * vij(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, jmp1) = force_mp(1:3, jmp1) - for(1:3)
     force_mp(1:3, jmp2) = force_mp(1:3, jmp2) - for(1:3)
  end do
!$omp end do nowait

end subroutine force_BBR
