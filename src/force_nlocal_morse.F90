! force_nlocal_morse
!> @brief This subroubine calculates the non-local interaction force by Morse potential.

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
subroutine force_nlocal_morse(irep, force_mp)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inrna, inpro, inmisc, inperi
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, nmorse, imorse2mp, &
                         coef_morse_a, coef_morse_fD, morse_nat, iclass_mp, nmp_all
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: force_mp(3, nmp_all)

  integer :: imp1, imp2
  integer :: ksta, kend
  integer :: imorse, imirror
  real(PREC) :: rcut_off
  real(PREC) :: rcut_off_pro, rcut_off_rna
  real(PREC) :: dist
  real(PREC) :: ex, dgo_dr
  real(PREC) :: v21(3), for(3)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  if (.not. inmisc%force_flag(INTERACT%MORSE)) then
     return
  endif

  rcut_off_pro = inpro%cutoff_go
  if (inmisc%class_flag(CLASS%RNA)) then
     rcut_off_rna = inrna%cutoff_go
  endif

#ifdef MPI_PAR
  klen=(nmorse-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nmorse)
#ifdef MPI_DEBUG
  print *,"nlocal_go    = ", kend-ksta+1
#endif
#else
  ksta = 1
  kend = nmorse
#endif
!$omp do private(imp1,imp2,rcut_off,v21,dist,ex,dgo_dr,for,imirror)
  do imorse=ksta,kend

     imp1 = imorse2mp(1, imorse)
     imp2 = imorse2mp(2, imorse)
     if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
        rcut_off = rcut_off_rna
     else
        rcut_off = rcut_off_pro
     endif

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist = sqrt(v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3))

     if(rcut_off < dist) cycle
     
     !ex = exp(-a * (r - r0))
     ex = exp(coef_morse_a(imorse) * (morse_nat(imorse) - dist))
     !force(1:3) = v21(1:3) / r * 2.0e0_PREC * De * ex * (1.0e0_PREC - ex)
     dgo_dr = coef_morse_fD(imorse) * coef_morse_a(imorse) * ex * (1.0e0_PREC - ex) / dist
     
     if(dgo_dr > DE_MAX) then
!        write (*, *) "go", imp1, imp2, dgo_dr
        dgo_dr = DE_MAX
     end if
!     if(dgo_dr > 5.0e0_PREC) dgo_dr = 5.0e0_PREC

     for(1:3) = dgo_dr * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
  end do
!$omp end do nowait

end subroutine force_nlocal_morse
