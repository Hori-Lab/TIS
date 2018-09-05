!force_exv_wca
!> @brief Calculates the force related to excluded volume

subroutine force_exv_dt15(irep, force_mp)

  use const_maxsize, only : PREC
  use const_physical,only : DE_MAX, SDIM
  use const_index,   only : CLASS, E_TYPE, ERROR, MPTYPE
  use var_setp,   only : indtrna15, inperi
  use var_struct, only : nmp_all, pxyz_mp_rep, imp2type, &
                         lexv, iexv2mp, iclass_mp, exv_radius_mp, exv_epsilon_mp
  use var_simu, only : istep
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SDIM, nmp_all)

  integer :: ksta, kend
  integer :: imp1, imp2, iexv, imirror
  real(PREC) :: dist, dr, dr2, dij, a, a2, coef, dv_dr
  real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist14
  real(PREC) :: v21(SDIM), for(SDIM)
#ifdef SHARE_NEIGH_PNL
  integer :: klen
#endif
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------

  a = indtrna15%exv_adjust
  a2 = a*a

#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lexv(2,E_TYPE%EXV_DT15,irep)-lexv(1,E_TYPE%EXV_DT15,irep)+npar_mpi)/npar_mpi
  ksta=lexv(1,E_TYPE%EXV_DT15,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lexv(2,E_TYPE%EXV_DT15,irep))
#else
  ksta = lexv(1, E_TYPE%EXV_DT15, irep)
  kend = lexv(2, E_TYPE%EXV_DT15, irep)
#endif
#ifdef MPI_DEBUG
  print *,"exv_dt15      = ", kend-ksta+1
#endif
#else
  ksta = lexv(1, E_TYPE%EXV_DT15, irep)
  kend = lexv(2, E_TYPE%EXV_DT15, irep)
#endif
!$omp do private(imp1,imp2,v21,dij,dist,coef,dr,dr2,&
!$omp&           roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist14,dv_dr,for,imirror)
  do iexv=ksta, kend

     imp1 = iexv2mp(1, iexv, irep)
     imp2 = iexv2mp(2, iexv, irep)
     imirror = iexv2mp(3, iexv, irep)

     if (imp2type(imp1) == MPTYPE%RNA_PHOS .AND. imp2type(imp2) == MPTYPE%RNA_SUGAR) then
        dij  = indtrna15%exv_dist_PS
     else if (imp2type(imp1) == MPTYPE%RNA_SUGAR .AND. imp2type(imp2) == MPTYPE%RNA_PHOS) then
        dij  = indtrna15%exv_dist_PS
     else if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
        dij  = indtrna15%exv_dist
     else
        dij  = exv_radius_mp(imp1)  + exv_radius_mp(imp2)
     endif

     !if(inperi%i_periodic == 0) then
     !   v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     !else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     !end if

     dist = sqrt(dot_product(v21,v21))

     if (dist > dij) cycle

     coef = exv_epsilon_mp(imp1) * exv_epsilon_mp(imp2)
     dr = dist + a - dij
     dr2 = dr * dr

     roverdist2 = a2 / dr2
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist14 = roverdist2 * roverdist4 * roverdist8

     dv_dr = abs(12.0e0_PREC * coef * (roverdist14 - roverdist8) * dr / a2 / dist)

     !if (dv_dr < 0) then
     !   dv_dr = DE_MAX
     !else if (dv_dr > DE_MAX) then
     if (dv_dr > DE_MAX) then
        write(error_message,*) 'force_exv_dt15 > DE_MAX', istep, imp1, imp2, dist, dv_dr
        call util_error(ERROR%WARN_ALL, error_message)
        dv_dr = DE_MAX
     end if

     for(1:3) = dv_dr * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
  end do
!$omp end do nowait

end subroutine force_exv_dt15
