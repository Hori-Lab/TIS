! energy_nlocal_morse
!> @brief Calculate non-local contact energy by Morse potential

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
subroutine energy_nlocal_morse(irep, now_morse, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inpro, inrna, inmisc
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, iclass_mp, &
                          nmorse, imorse2mp, morse_nat, &
                          coef_morse_fD, coef_morse_a
  use mpiconst

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  integer,    intent(out)   :: now_morse(:,:)
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  integer :: imp1, imp2, iunit, junit
  integer :: ksta, kend
  integer :: imorse, imirror
  real(PREC) :: rcut_off, rcut_off_pro, rcut_off_rna
  real(PREC) :: rjudge
  real(PREC) :: dist, efull, ex
  real(PREC) :: v21(SDIM)
  real(PREC), parameter :: rjudge_contact = 1.2e0_PREC
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  ! zero clear
  if (.not. inmisc%force_flag(INTERACT%MORSE)) then
     return
  endif
      
  ! --------------------------------------------------------------------
!!$omp master

  ! --------------------------------------------------------------------
  ! zero clear
!  now_morse(:,:) = 0
      
  ! --------------------------------------------------------------------
  ! set parameter 
  rcut_off_pro = inpro%cutoff_go
  if (inmisc%class_flag(CLASS%RNA)) then
     rcut_off_rna =  inrna%cutoff_go
  endif

  ! --------------------------------------------------------------------
#ifdef MPI_PAR3
   klen=(nmorse-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,nmorse)
#else
   ksta = 1
   kend = nmorse
#endif
!$omp do private(imp1,imp2,rcut_off,v21,dist,ex,efull, &
!$omp&           rjudge,iunit,junit,imirror)
   do imorse=ksta,kend
   
     imp1 = imorse2mp(1, imorse)
     imp2 = imorse2mp(2, imorse)
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist = sqrt(v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3))

     ! rjudge = 1.2 * r0
     rjudge = morse_nat(imorse) * rjudge_contact
     !  judging contact 
     if(dist < rjudge) then
        now_morse(1, imorse) = 1
     else
        now_morse(1, imorse) = 0
     end if
     if(coef_morse_a(imorse) > ZERO_JUDGE .and. coef_morse_fD(imorse) > ZERO_JUDGE) then
        now_morse(2, imorse) = 1
     else
        now_morse(2, imorse) = 0
     end if

     if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
        rcut_off = rcut_off_rna
     else 
        rcut_off = rcut_off_pro
     endif
     if(rcut_off*morse_nat(imorse) < dist) cycle
     
     ! calc energy
     ex = exp(coef_morse_a(imorse) * (morse_nat(imorse) - dist))
     efull = coef_morse_fD(imorse) * (1.0_PREC - ex) ** 2

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%MORSE) = energy(E_TYPE%MORSE) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%MORSE) = energy_unit(iunit, junit, E_TYPE%MORSE) + efull
  end do
!$omp end do nowait
!!$omp end master

end subroutine energy_nlocal_morse
