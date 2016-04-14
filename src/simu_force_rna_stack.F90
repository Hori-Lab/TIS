!simu_energy_rna_stack
!> @brief Calculates the force related to stacking interaction of  &
!>        RNA particles.

! ************************************************************************
! formula of go1210
! ego = coef_go * {5*(go_nat/x)**12 -6*(go_nat/x)**10}
! coef_go = icon_dummy_mgo * factor_go * cgo1210
!
! parameter list
! cgo1210: constant of go energy
! factor_go: value of amino acid specifity (ex. 1.0, MJ)
! go_nat: distance of native contact 
! ***********************************************************************
subroutine simu_force_rna_stack(irep, force_mp)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inrna, inmisc
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, nrna_st, irna_st2mp, &
                         coef_rna_st, coef_rna_st_a, coef_rna_st_fD,   &
                         rna_st_nat, rna_st_nat2, nmp_all
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: force_mp(SPACE_DIM, nmp_all)

  integer :: imp1, imp2, ist, imirror
  real(PREC) :: dist, dist2, ex, dgo_dr
  real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist12, roverdist14
  real(PREC) :: v21(SPACE_DIM), for(SPACE_DIm)
  character(CARRAY_MSG_ERROR) :: error_message
  integer :: ksta, kend
#ifdef MPI_PAR
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  if (.not. inmisc%class_flag(CLASS%RNA)) then
     return
  endif

#ifdef MPI_PAR
  klen=(nrna_st-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nrna_st)
#ifdef MPI_DEBUG
  print *,"rna stack    = ", kend-ksta+1
#endif
#else
  ksta = 1
  kend = nrna_st
#endif

  if (inrna%i_potential_st == POTTYPE%LJ1210) then

!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,roverdist14,dgo_dr,for,imirror)
     do ist = ksta, kend

        imp1 = irna_st2mp(1, ist)
        imp2 = irna_st2mp(2, ist)

        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
           call util_pbneighbor(v21, imirror)
        end if
     
!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

        roverdist2 = rna_st_nat2(ist) / dist2
     
        roverdist4 = roverdist2 * roverdist2
        roverdist8 = roverdist4 * roverdist4
        roverdist12 = roverdist4 * roverdist8
        roverdist14 = roverdist12 * roverdist2
            
        dgo_dr = 60.0e0_PREC * coef_rna_st(ist) / rna_st_nat2(ist) * (roverdist14 - roverdist12)
     
        if(dgo_dr > DE_MAX) then
           dgo_dr = DE_MAX
        end if

        for(1:3) = dgo_dr * v21(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
        force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     end do
!$omp end do nowait

  else if (inrna%i_potential_st == POTTYPE%MORSE) then

!$omp do private(imp1,imp2,v21,dist,ex,dgo_dr,for,imirror)
     do ist = ksta, kend

        imp1 = irna_st2mp(1, ist)
        imp2 = irna_st2mp(2, ist)
   
        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
           call util_pbneighbor(v21, imirror)
        end if
     
!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist = sqrt(v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3))
     
        !ex = exp(-a * (r - r0))
        ex = exp(coef_rna_st_a(ist) * (rna_st_nat(ist) - dist))
        !force(1:3) = v21(1:3) / r * 2.0e0_PREC * De * ex * (1.0e0_PREC - ex)
        dgo_dr = coef_rna_st_fD(ist) * coef_rna_st_a(ist) * ex * (1.0e0_PREC - ex) / dist
     
        if(dgo_dr > DE_MAX) then
           dgo_dr = DE_MAX
        end if

        for(1:3) = dgo_dr * v21(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
        force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     end do
!$omp end do nowait

  else
     error_message = 'Error: invalid value for i_potential_st in simu_force_rna_stack'
     call util_error(ERROR%STOP_ALL, error_message)

  endif

end subroutine simu_force_rna_stack
