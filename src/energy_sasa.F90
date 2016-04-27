!energy_sasa
!> @brief Calculates the energy related to solvent accessible surface area (sasa)
!>        The values are added into "e_exv(ENERGY%SASA)" and  &
!>        "e_exv_unit(ENERGY%SASA)".

subroutine  energy_sasa(irep, e_exv)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi, outfile
  use var_setp,    only : insasa 
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, lexv, iexv2mp, nmp_real, &
                          para_sasa, rad_sasa, surf, connect,cmp2atom
  use var_simu,    only : istep, sasa
  use mpiconst
  
  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: e_exv(:)         ! (E_TYPE%MAX)

  integer :: lunout
  integer :: ksta, kend
  integer :: imp1, imp2, iunit1, iunit2
  integer :: imirror, imp, isasa
  real(PREC) :: dist2(lexv(1, E_TYPE%SASA, irep):lexv(2, E_TYPE%SASA, irep))
  real(PREC) :: dist(lexv(1, E_TYPE%SASA, irep):lexv(2, E_TYPE%SASA, irep))
  real(PREC) :: radsum
  real(PREC) :: v21(SPACE_DIM)
  real(PREC) :: sasa_t, c2, c3ij, c3ji, bij, bji, fij, fji

!$omp master
  lunout = outfile%data

  ksta = lexv(1, E_TYPE%SASA, irep)
  kend = lexv(2, E_TYPE%SASA, irep)

!initialization

!!$omp do
  do imp=1,nmp_real
     if(cmp2atom(imp) == ' CA ' .or. cmp2atom(imp) == ' P  ' .or. cmp2atom(imp) == ' O  ')then
        sasa(imp) = surf(imp)
     end if
  end do
!!$omp end do

!----------------------------------------------------------------------------
!!$omp parallel do shared(sasa,dist2,dist) private(imp1,imp2,iunit1,iunit2,v21,&
!!$omp&                                    radsum,c2,c3ij,c3ji,bij,bji,fij,fji)
  do isasa=ksta, kend

     imp1 = iexv2mp(1, isasa, irep)
     imp2 = iexv2mp(2, isasa, irep)

     if(cmp2atom(imp1) /= ' CA ' .and. cmp2atom(imp1) /= ' P  ' .and. cmp2atom(imp1) /= ' O  ')then
       write(*,*)'*********************************'
       write(*,*)'CG atoms other than "CA, O, and P" not supported in sasa'
       stop
     end if

     if(cmp2atom(imp2) /= ' CA ' .and. cmp2atom(imp2) /= ' P  ' .and. cmp2atom(imp2) /= ' O  ')then
       write(*,*)'*********************************'
       write(*,*)'CG atoms other than "CA, O, and P" not supported in sasa'
       stop
     end if

     iunit1 = imp2unit(imp1)
     iunit2 = imp2unit(imp2)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = iexv2mp(3, isasa, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if

     dist2(isasa) = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
     dist(isasa) = sqrt(dist2(isasa))
     radsum =  rad_sasa(imp1) + rad_sasa(imp2)

     if(dist(isasa) .gt. radsum)cycle

     c2 = radsum - dist(isasa)
     c3ij = (1.0e0_PREC + (rad_sasa(imp2) - rad_sasa(imp1)) / dist(isasa))
     c3ji = (1.0e0_PREC + (rad_sasa(imp1) - rad_sasa(imp2)) / dist(isasa))

     bij = para_sasa(imp1) * c2 * c3ij
     bji = para_sasa(imp2) * c2 * c3ji

     fij = 1.0e0_PREC - connect(imp1-imp2) * bij
     fji = 1.0e0_PREC - connect(imp2-imp1) * bji

     sasa(imp1) = sasa(imp1) * fij
     sasa(imp2) = sasa(imp2) * fji

  end do
!!$omp end parallel do

!!$omp master
  if(istep == 0)then
     if (myrank == 0) then
        write(lunout,'(A45)')'**********************************************'
        write(lunout,'(A35)')'residue No.,  SASA of ith residue'
     endif
  end if
!!$omp end master

!--------------------------------------------------------------------------
  sasa_t=0.
!!$omp do
  do imp = 1, nmp_real
     if(cmp2atom(imp) == ' CA ' .or. cmp2atom(imp) == ' P  ' .or. cmp2atom(imp) == ' O  ')then
        sasa_t = sasa_t + sasa(imp)
        if(istep == 0)then
           if (myrank == 0) then
              write(lunout,'(I10,f12.2)')imp, sasa(imp)
           endif
        end if
        e_exv(E_TYPE%SASA) = e_exv(E_TYPE%SASA) + insasa%coef_surf * sasa(imp)
     end if
  end do
!!$omp end do

!$omp end master

! The sasa energy is not decomposed into subunit at this moment!!!!!
end subroutine energy_sasa
