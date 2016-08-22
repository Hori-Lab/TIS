!energy_rna_stack
!> @brief Calculates the energy related to stacking interaction of  &
!>        RNA particles.
!>        The values are added into "energy(E_TYPE%STACK_RNA)" and   &
!>        "energy_unit(E_TYPE%STACK_RNA)".

! ************************************************************************
! formula of go1210
! ego = coef_go * {5*(go_nat/x)**12 -6*(go_nat/x)**10}
! coef_go = icon_dummy_mgo * factor_go * cgo1210
!
! parameter list
! cgo1210: constant of go energy
! factor_go: value of amino acid specifity (ex. 1.0, MJ)
! go_nat: distance of native contact 
! ************************************************************************
subroutine energy_rna_stack(irep, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inrna, inmisc, inperi
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, &
                          nrna_st, irna_st2mp, rna_st_nat, rna_st_nat2, &
                          coef_rna_st, coef_rna_st_fD, coef_rna_st_a
  use mpiconst

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp1, imp2, iunit, junit
  integer :: ist, imirror
  real(PREC) :: dist, dist2, efull, ex
  real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist10, roverdist12
  real(PREC) :: v21(SDIM)
  character(CARRAY_MSG_ERROR) :: error_message
  integer :: ksta, kend
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  if (.not. inmisc%class_flag(CLASS%RNA)) then
     return
  endif

#ifdef MPI_PAR3
   klen=(nrna_st-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,nrna_st)
#else
   ksta = 1
   kend = nrna_st
#endif

   if (inrna%i_potential_st == POTTYPE%LJ1210) then
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist10,roverdist12,efull,iunit,junit,imirror)
      do ist=ksta,kend

         imp1 = irna_st2mp(1, ist)
         imp2 = irna_st2mp(2, ist)
        
         if(inperi%i_periodic == 0) then
            v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
         else
            v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
            call util_pbneighbor(v21, imirror)
         end if
     
!         v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

         dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

         roverdist2 = rna_st_nat2(ist) / dist2

         ! calc energy
         roverdist4 = roverdist2 * roverdist2
         roverdist8 = roverdist4 * roverdist4
         roverdist10 = roverdist2 * roverdist8
         roverdist12 = roverdist4 * roverdist8

         efull = coef_rna_st(ist) * (5.0e0_PREC * roverdist12 - 6.0e0_PREC * roverdist10)

         ! --------------------------------------------------------------------
         ! sum of the energy
         energy(E_TYPE%STACK_RNA) = energy(E_TYPE%STACK_RNA) + efull

         iunit = imp2unit(imp1)
         junit = imp2unit(imp2)
         energy_unit(iunit, junit, E_TYPE%STACK_RNA) = energy_unit(iunit, junit, E_TYPE%STACK_RNA) + efull
      end do
!$omp end do nowait


   else if (inrna%i_potential_st == POTTYPE%MORSE) then

!$omp do private(imp1,imp2,v21,dist,ex,efull,iunit,junit,imirror)
      do ist=ksta,kend

         imp1 = irna_st2mp(1, ist)
         imp2 = irna_st2mp(2, ist)
        
         if(inperi%i_periodic == 0) then
            v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
         else
            v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
            call util_pbneighbor(v21, imirror)
         end if
     
!         v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

         dist = sqrt(v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3))

         ! calc energy
         ex = exp(coef_rna_st_a(ist) * (rna_st_nat(ist) - dist))
         efull = coef_rna_st_fD(ist) * (1.0_PREC - ex) ** 2

         ! --------------------------------------------------------------------
         ! sum of the energy
         energy(E_TYPE%STACK_RNA) = energy(E_TYPE%STACK_RNA) + efull

         iunit = imp2unit(imp1)
         junit = imp2unit(imp2)
         energy_unit(iunit, junit, E_TYPE%STACK_RNA) = energy_unit(iunit, junit, E_TYPE%STACK_RNA) + efull
      end do
!$omp end do nowait

   else
      error_message = 'Error: invalid value for i_potential_st in energy_rna_stack'
      call util_error(ERROR%STOP_ALL, error_message)

   endif

end subroutine energy_rna_stack
