! simu_force_nlocal_rna_bp
!> @brief Calculate force of RNA base pair

! ***********************************************************************
subroutine simu_force_nlocal_rna_bp(irep, force_mp)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inpara, inrna, inpro, inmisc
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, nrna_bp, irna_bp2mp, &
                         coef_rna_bp, coef_rna_bp_a, coef_rna_bp_fD,   &
                         rna_bp_nat, rna_bp_nat2,  &
                         iclass_mp, nunit_all, nmp_all

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: force_mp(3, nmp_all)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: imp1, imp2
  integer :: ksta, kend
  integer :: ibp, imirror
  real(PREC) :: rcut_off, rcut_off2
  real(PREC) :: dist, dist2
  real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist12, roverdist14
  real(PREC) :: ex, dgo_dr
  real(PREC) :: v21(3), for(3)
  character(CARRAY_MSG_ERROR) :: error_message
#ifdef MPI_PAR
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  if (.not. inmisc%force_flag(INTERACT%PAIR_RNA)) then
     return
  endif

#ifdef MPI_PAR
  klen=(nrna_bp-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nrna_bp)
#ifdef MPI_DEBUG
  print *,"nlocal_go    = ", kend-ksta+1
#endif
#else
  ksta = 1
  kend = nrna_bp
#endif

  if (inrna%i_potential_bp == POTTYPE%LJ1210) then
     rcut_off2 = 1.0e0_PREC / inrna%cutoff_bp**2
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,roverdist14,dgo_dr,for,imirror)
     do ibp=ksta,kend
        imp1 = irna_bp2mp(1, ibp)
        imp2 = irna_bp2mp(2, ibp)

        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
           call util_pbneighbor(v21, imirror)
        end if
     
!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

        roverdist2 = rna_bp_nat2(ibp) / dist2
        if(roverdist2 < rcut_off2) cycle
        
        roverdist4 = roverdist2 * roverdist2
        roverdist8 = roverdist4 * roverdist4
        roverdist12 = roverdist4 * roverdist8
        roverdist14 = roverdist12 * roverdist2
            
        dgo_dr = 60.0e0_PREC * coef_rna_bp(ibp) / rna_bp_nat2(ibp) * (roverdist14 - roverdist12)
     
        if(dgo_dr > DE_MAX) then
           dgo_dr = DE_MAX
        end if

        for(1:3) = dgo_dr * v21(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
        force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     end do
!$omp end do nowait

  else if (inrna%i_potential_bp == POTTYPE%MORSE) then
     rcut_off = inrna%cutoff_bp
!$omp do private(imp1,imp2,v21,dist,ex,dgo_dr,for)
     do ibp=ksta,kend
        imp1 = irna_bp2mp(1, ibp)
        imp2 = irna_bp2mp(2, ibp)

        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        dist = sqrt(v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3))

        if(rcut_off*rna_bp_nat(ibp) < dist) cycle
     
        !ex = exp(-a * (r - r0))
        !force(1:3) = v21(1:3) / r * 2.0e0_PREC * De * ex * (1.0e0_PREC - ex)
        ex = exp(coef_rna_bp_a(ibp) * (rna_bp_nat(ibp) - dist))
        dgo_dr = coef_rna_bp_fD(ibp) * coef_rna_bp_a(ibp) * ex * (1.0e0_PREC - ex) / dist
     
        if(dgo_dr > DE_MAX) then
           dgo_dr = DE_MAX
        end if

        for(1:3) = dgo_dr * v21(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
        force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     end do
!$omp end do nowait
  
  else
     error_message = 'Error: invalid value for i_potential_bp in simu_nlocal_rna_bp'
     call util_error(ERROR%STOP_ALL, error_message)

  endif

end subroutine simu_force_nlocal_rna_bp
