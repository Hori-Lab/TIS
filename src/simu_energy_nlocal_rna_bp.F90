! simu_energy_nlocal_rna_bp
!> @brief This subroutine is to calculate the nonlocal interaction energy espcesially for RNA base-pairs.

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
subroutine simu_energy_nlocal_rna_bp(irep, now_rna_bp, pnle_unit, pnlet)

   use const_maxsize
   use const_physical
   use const_index
   use var_inp,     only : inperi
   use var_setp,    only : inpro, inrna, inmisc
   use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, iclass_mp, &
                           nrna_bp, irna_bp2mp, rna_bp_nat, rna_bp_nat2, &
                           coef_rna_bp, coef_rna_bp_fD, coef_rna_bp_a
   use var_replica, only : n_replica_mpi
#ifdef MPI_PAR3
   use mpiconst
#endif

   implicit none

   ! --------------------------------------------------------------------
   integer,    intent(in)    :: irep
   integer,    intent(out)   :: now_rna_bp(:,:)
   real(PREC), intent(inout) :: pnlet(:)
   real(PREC), intent(inout) :: pnle_unit(:,:,:)

   ! --------------------------------------------------------------------
   ! local variables
   integer :: imp1, imp2, iunit, junit
   integer :: ibp, imirror
   real(PREC) :: rcut_off, rcut_off2
   real(PREC) :: rjudge
   real(PREC) :: dist, dist2, efull, ex
   real(PREC) :: v21(SPACE_DIM)
   real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist10, roverdist12
   real(PREC), parameter :: rjudge_contact = 1.2e0_PREC
   real(PREC), parameter :: rjudge_contact2= rjudge_contact ** 2
   character(CARRAY_MSG_ERROR) :: error_message
   integer :: ksta, kend
#ifdef MPI_PAR3
   integer :: klen
#endif

   ! --------------------------------------------------------------------
   if (.not. inmisc%force_flag(INTERACT%PAIR_RNA)) then
      return
   endif
      
  ! --------------------------------------------------------------------
!!$omp master

   ! --------------------------------------------------------------------
   ! zero clear
!   now_rna_bp(:,:) = 0
      
   ! --------------------------------------------------------------------
#ifdef MPI_PAR3
   klen=(nrna_bp-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,nrna_bp)
#else
   ksta = 1
   kend = nrna_bp
#endif

   if (inrna%i_potential_bp == POTTYPE%LJ1210) then

      rcut_off2 = 1.0_PREC / inrna%cutoff_bp ** 2

!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist10,roverdist12, &
!$omp&           efull,rjudge,iunit,junit,imirror)
      do ibp=ksta,kend

         imp1 = irna_bp2mp(1, ibp)
         imp2 = irna_bp2mp(2, ibp)
           
         if(inperi%i_periodic == 0) then
            v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
         else
            v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
            call util_pbneighbor(v21, imirror)
         end if
     
!         v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

         dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

         roverdist2 = rna_bp_nat2(ibp) / dist2
     
         ! 1.44 = 1.2 *1.2
         rjudge = rna_bp_nat2(ibp) * rjudge_contact2
         !  judging contact 
         if(dist2 < rjudge) then
            now_rna_bp(1, ibp) = 1
         else
            now_rna_bp(1, ibp) = 0
         end if
         if(coef_rna_bp(ibp) > ZERO_JUDGE) then
            now_rna_bp(2, ibp) = 1
         else
            now_rna_bp(2, ibp) = 0
         end if

         if(roverdist2 < rcut_off2) cycle

         ! calc energy
         roverdist4 = roverdist2 * roverdist2
         roverdist8 = roverdist4 * roverdist4
         roverdist10 = roverdist2 * roverdist8
         roverdist12 = roverdist4 * roverdist8

         efull = coef_rna_bp(ibp) * (5.0e0_PREC * roverdist12 - 6.0e0_PREC * roverdist10)
 
         ! --------------------------------------------------------------------
         ! sum of the energy
         pnlet(E_TYPE%PAIR_RNA) = pnlet(E_TYPE%PAIR_RNA) + efull
    
         iunit = imp2unit(imp1)
         junit = imp2unit(imp2)
         pnle_unit(iunit, junit, E_TYPE%PAIR_RNA) = pnle_unit(iunit, junit, E_TYPE%PAIR_RNA) + efull
      end do
!$omp end do nowait

   else if (inrna%i_potential_bp == POTTYPE%MORSE) then

      rcut_off  =  inrna%cutoff_bp

!$omp do private(imp1,imp2,v21,dist,rjudge,ex,efull,iunit,junit)
      do ibp=ksta,kend

         imp1 = irna_bp2mp(1, ibp)
         imp2 = irna_bp2mp(2, ibp)
        
         v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
         dist = sqrt(v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3))
   
         if(rcut_off*rna_bp_nat(ibp) < dist) cycle
     
         ! rjudge = 1.2 * r0
         rjudge = rna_bp_nat(ibp) * rjudge_contact
         !  judging contact 
         if(dist < rjudge) then
            now_rna_bp(1, ibp) = 1
         else
            now_rna_bp(1, ibp) = 0
         end if
         if(coef_rna_bp_a(ibp) > ZERO_JUDGE .and. coef_rna_bp_fD(ibp) > ZERO_JUDGE) then
            now_rna_bp(2, ibp) = 1
         else
            now_rna_bp(2, ibp) = 0
         end if

         ! calc energy
         ex = exp(coef_rna_bp_a(ibp) * (rna_bp_nat(ibp) - dist))
         efull = coef_rna_bp_fD(ibp) * (1.0_PREC - ex) ** 2

         ! --------------------------------------------------------------------
         ! sum of the energy
         pnlet(E_TYPE%PAIR_RNA) = pnlet(E_TYPE%PAIR_RNA) + efull
 
         iunit = imp2unit(imp1)
         junit = imp2unit(imp2)
         pnle_unit(iunit, junit, E_TYPE%PAIR_RNA) = pnle_unit(iunit, junit, E_TYPE%PAIR_RNA) + efull
      end do
!$omp end do nowait

   else
      error_message = 'Error: invalid value for i_potential_bp in simu_energy_nlocal_rna_bp'
      call util_error(ERROR%STOP_ALL, error_message)

   endif

!!$omp end master

end subroutine simu_energy_nlocal_rna_bp
