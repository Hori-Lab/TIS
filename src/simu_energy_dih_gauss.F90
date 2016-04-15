!simu_energy_dih_gauss
!> @brief Calculates the energy related to 1-4 residues.
!>        Values are added into "e_exv(E_TYPE%DIHE)" and &
!>        "e_exv_unit(,,E_TYPE%DIHE)".

subroutine simu_energy_dih_gauss(irep, e_exv_unit, e_exv)
  
  use const_maxsize
  use const_index
  use const_physical
  use var_struct, only : ndih, idih2mp, imp2unit, nunit_all, xyz_mp_rep, &
                         coef_dih_gauss, wid_dih_gauss, dih_nat, iclass_mp

#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------------
  integer, intent(in) :: irep
  real(PREC), intent(inout) :: e_exv(E_TYPE%MAX) 
  real(PREC), intent(inout) :: e_exv_unit(nunit_all, nunit_all, E_TYPE%MAX)
  
  ! ---------------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: imp1, imp2, imp3, imp4
  integer :: idih, iunit, junit
  real(PREC) :: c11, c12, c13, c22, c23, c33
  real(PREC) :: t1, t3, t4, t3t4
  real(PREC) :: zahyokei 
  real(PREC) :: co_dih
  real(PREC) :: pre1
  real(PREC) :: v21(3), v32(3), v43(3)
  real(PREC) :: efull, dih
#ifdef MPI_PAR3
  integer :: klen
#endif
  ! ----------------------------------------------------------------------

#ifdef MPI_PAR3
  klen=(ndih-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,ndih)
#else
  ksta = 1
  kend = ndih
#endif
!$omp do private(imp1,imp2,imp3,imp4,v21,v32,v43,c11,c22,c33,c12,c13,c23, &
!$omp&           t1,t3,t4,t3t4,co_dih, zahyokei, pre1, dih, &
!$omp&           efull, iunit, junit)
  do idih = ksta, kend

     if(coef_dih_gauss(idih) < ZERO_JUDGE) cycle

     imp1 = idih2mp(1, idih)
     imp2 = idih2mp(2, idih)
     imp3 = idih2mp(3, idih)
     imp4 = idih2mp(4, idih)

     if (iclass_mp( imp1 ) == CLASS%LIG .AND. &
         iclass_mp( imp4 ) == CLASS%LIG) cycle

     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     v32(1:3) = xyz_mp_rep(1:3, imp3, irep) - xyz_mp_rep(1:3, imp2, irep)
     v43(1:3) = xyz_mp_rep(1:3, imp4, irep) - xyz_mp_rep(1:3, imp3, irep)

     c11 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
     c22 = v32(1)*v32(1) + v32(2)*v32(2) + v32(3)*v32(3)
     c33 = v43(1)*v43(1) + v43(2)*v43(2) + v43(3)*v43(3)

     c12 = v21(1)*v32(1) + v21(2)*v32(2) + v21(3)*v32(3)
     c13 = v21(1)*v43(1) + v21(2)*v43(2) + v21(3)*v43(3)
     c23 = v32(1)*v43(1) + v32(2)*v43(2) + v32(3)*v43(3)

     t1 = c12*c23 - c13*c22
     t3 = c11*c22 - c12*c12
     t4 = c22*c33 - c23*c23

     if (t3 < ZERO_JUDGE) t3 = ZERO_JUDGE
     if (t4 < ZERO_JUDGE) t4 = ZERO_JUDGE

     t3t4 = sqrt(t3*t4)
     co_dih      = t1  / t3t4

     ! when co_dih > 1 ,or < -1 
     if(co_dih > 1.0e0_PREC) then
        co_dih = 1.0e0_PREC
     else if(co_dih < -1.0e0_PREC) then
        co_dih = -1.0e0_PREC
     end if

     zahyokei = v21(1) * v32(2) * v43(3) + &
          v32(1) * v43(2) * v21(3) + &
          v43(1) * v21(2) * v32(3) - &
          v43(1) * v32(2) * v21(3) - &
          v21(1) * v43(2) * v32(3) - &
          v32(1) * v21(2) * v43(3)

     if(zahyokei > 0.0e0_PREC) then
        dih = acos(co_dih)
     else
        dih = -acos(co_dih)
     end if

     if((dih - dih_nat(idih)) > F_PI) dih = dih - 2.0e0_PREC * F_PI
     if((dih - dih_nat(idih)) < -F_PI) dih = dih + 2.0e0_PREC * F_PI

     pre1 = exp(-(dih-dih_nat(idih)) * (dih-dih_nat(idih)) /  &
                 (2.0e0_PREC * wid_dih_gauss(idih) * wid_dih_gauss(idih)))
     efull = -1.0e0_PREC * coef_dih_gauss(idih) * pre1
     
     !------------------------------------------------------------------------
     ! sum of the energy
     iunit = imp2unit( imp1 )
     junit = imp2unit( imp4 )

     e_exv(E_TYPE%DIHE) = e_exv(E_TYPE%DIHE) + efull
     e_exv_unit(iunit, junit, E_TYPE%DIHE) = e_exv_unit(iunit, junit, E_TYPE%DIHE) + efull
     
     !end if
  end do
!$omp end do nowait
  
end subroutine simu_energy_dih_gauss
