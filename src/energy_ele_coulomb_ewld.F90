! energy_ele_coulomb
!> @brief Calculate the energy of electrostatic interaction 

subroutine energy_ele_coulomb_ewld(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inele
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, lele, iele2mp, coef_ele, &
                          ncharge, coef_charge, icharge2mp
  use var_simu,    only : ewld_f_n, ewld_f_coef, ewld_f_rlv, ewld_s_coef
  use mpiconst

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit, iele, imirror, ich1, ich2, ig
  real(PREC) :: dist1, dist2, ssin, scos, ene, cutoff2, q1, q2, dp
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  !================================================
  !================== Real space ==================
  !================================================
  cutoff2 = inele%cutoff_ele ** 2
  !cutoff2 = inperi%psizeh(1) ** 2   ! Assuming a cubic box

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH
  klen=(lele(irep)-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,lele(irep))
#else
  ksta = 1
  kend = lele(irep)
#endif
#else
  ksta = 1
  kend = lele(irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,dist1,ene,iunit,junit,imirror)
  do iele=ksta, kend

     imp1 = iele2mp(1, iele, irep)
     imp2 = iele2mp(2, iele, irep)

     imirror = iele2mp(3, iele, irep)
     v21(1:SDIM) = pxyz_mp_rep(1:SDIM, imp2, irep) - pxyz_mp_rep(1:SDIM, imp1, irep) + inperi%d_mirror(1:SDIM, imirror)
     
     dist2 = dot_product(v21,v21)
     if(dist2 > cutoff2) cycle
        
     ! --------------------------------------------------------------------
     dist1 = sqrt(dist2)
     
     ene = coef_ele(iele,irep ) * erfc(inele%ewld_alpha*dist1) / dist1

     energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene
     
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + ene
  end do
!$omp end do nowait

!! In case no neighbor list used
!!$omp do private(ich1,imp1,q1,ich2,imp2,v21,dist2,dist1,ene,iunit,junit,imirror)
!  do ich1 = 1, ncharge
!
!     imp1 = icharge2mp(ich1)
!     q1 = coef_charge(ich1,irep)
!
!     do ich2 = 1, ncharge
!
!        imp2 = icharge2mp(ich2)
!              
!        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
!        call util_pbneighbor(v21, imirror)
!
!        dist2 = dot_product(v21,v21)
!        if(dist2 > cutoff2) cycle
!        dist1 = sqrt(dist2)
!        
!        q2 = coef_charge(ich2,irep)
!        ene = inele%coef(irep) * q1 * q2 * erfc(inele%ewld_alpha*dist1) / dist1
!   
!        energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene
!        
!        iunit = imp2unit(imp1)
!        junit = imp2unit(imp2)
!        energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + ene
!     end do
!  end do
!!$omp end do nowait

  !================================================
  !================= Fourier space ================
  !================================================

  ene = 0.0
  do ig = 1, ewld_f_n
    
     scos = 0.0
     ssin = 0.0

     do ich1 = 1, ncharge
        imp1 = icharge2mp(ich1)
        q1 = coef_charge(ich1, irep)

        dp = dot_product(ewld_f_rlv(:,ig), pxyz_mp_rep(:,imp1, irep))
        scos = scos + q1 * cos(dp)
        ssin = ssin + q1 * sin(dp)
     end do
      
     ene = ene + ewld_f_coef(ig) * (scos ** 2 + ssin ** 2)
  end do

  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + inele%coef(irep) * ene


  !================================================
  !======= Correction for self interaction ========
  !================================================

  ene = 0.0
  do ich1 = 1, ncharge
     q1 = coef_charge(ich1, irep)
     ene = ene + q1**2
  end do

  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + inele%coef(irep) * ewld_s_coef * ene

end subroutine energy_ele_coulomb_ewld


subroutine energy_ele_coulomb_ewld_tp(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inele
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, lele, iele2mp, coef_ele, &
                          ncharge, coef_charge, icharge2mp, nmp_real, lmp2charge, &
                          ntp, xyz_tp, charge_tp
  use var_replica, only : irep2grep
  use var_simu,    only : ewld_f_n, ewld_f_coef, ewld_f_rlv, ewld_s_coef
  use mpiconst

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  integer :: ksta, kend
  integer :: itp1, itp2, imp2, imirror, ig, grep
  real(PREC) :: dist1, dist2, ssin, scos, ene, cutoff2, q1, q2, dp
  real(PREC) :: v21(SDIM), vx(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  grep = irep2grep(irep)

  !================================================
  !================== Real space ==================
  !================================================
  cutoff2 = inele%cutoff_ele ** 2
  !cutoff2 = inperi%psizeh(1) ** 2   ! Assuming a cubic box

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH
  klen=(lele(irep)-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,lele(irep))
#else
  ksta = 1
  kend = lele(irep)
#endif
#else
  ksta = 1
  kend = lele(irep)
#endif
!!$omp do private(imp1,imp2,v21,dist2,dist1,ene,iunit,junit,imirror)
  do itp1 = 1, ntp

     q1 = charge_tp(itp1)

     do imp2 = 1, nmp_real

        vx(1:3) = pxyz_mp_rep(1:3, imp2, irep) - xyz_tp(1:3, itp1)
        call util_pbneighbor(vx, imirror)
        v21(1:3) = vx(1:3) + inperi%d_mirror(1:3, imirror)
        
        dist2 = dot_product(v21,v21)
        if(dist2 > cutoff2) cycle

        dist1 = sqrt(dist2)
        ene = ene + q1 * coef_charge( lmp2charge(imp2), grep) * erfc(inele%ewld_alpha*dist1) / dist1
     end do

     do itp2 = itp1+1, ntp
        vx(1:3) = xyz_tp(1:3, itp2) - xyz_tp(1:3, itp1)
        call util_pbneighbor(vx, imirror)
        v21(1:3) = vx(1:3) + inperi%d_mirror(1:3, imirror)
        
        dist2 = dot_product(v21,v21)
        if(dist2 > cutoff2) cycle
           
        dist1 = sqrt(dist2)
        ene = ene + q1 * charge_tp(itp2) * erfc(inele%ewld_alpha*dist1) / dist1
     end do
     
  end do


  !================================================
  !================= Fourier space ================
  !================================================
  do ig = 1, ewld_f_n
    
     scos = 0.0
     ssin = 0.0

     do itp1 = 1, ntp

        dp = dot_product(ewld_f_rlv(:,ig), xyz_tp(:,itp1))

        q1 = charge_tp(itp1)

        scos = scos + q1 * cos(dp)
        ssin = ssin + q1 * sin(dp)
     end do
      
     ene = ene + ewld_f_coef(ig) * (scos ** 2 + ssin ** 2)
  end do


  !================================================
  !======= Correction for self interaction ========
  !================================================
  ene = 0.0
  do itp1 = 1, ntp
     q1 = charge_tp(itp1)
     ene = ene + ewld_s_coef * q1**2
  end do

  ene = ene * inele%coef(irep)

end subroutine energy_ele_coulomb_ewld_tp

!subroutine energy_ele_coulomb_tp(irep, energy)
!
!  use const_maxsize
!  use const_physical
!  use const_index
!  use var_inp,     only : inperi
!  use var_setp,    only : inele
!  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, nmp_real, lmp2charge, coef_charge, &
!                          ntp, xyz_tp, charge_tp
!  use var_replica, only : irep2grep
!  use mpiconst
!
!  implicit none
!
!  integer,    intent(in)    :: irep
!  real(PREC), intent(out)   :: energy(E_TYPE%MAX)
!
!  integer :: itp1, itp2, imp2, imirror, grep
!  real(PREC) :: dist2, coef1, ene
!  real(PREC) :: cutoff2
!  real(PREC) :: v21(SDIM), vx(SDIM)
!
!  cutoff2 = inele%cutoff_ele ** 2
!  grep = irep2grep(irep)
!
!  ene = 0.0
!
!!$omp do private(itp1,imp2,itp2,v21,vx,dist2,imirror)
!  do itp1 = 1, ntp
!
!     coef1 = charge_tp(itp1) * inele%coef(grep)
!
!     do imp2 = 1, nmp_real
!   
!        !if(inperi%i_periodic == 0) then
!        !   v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_tp(1:3, itp1)
!        !else
!           vx(1:3) = pxyz_mp_rep(1:3, imp2, irep) - xyz_tp(1:3, itp1)
!           call util_pbneighbor(vx, imirror)
!           v21(1:3) = vx(1:3) + inperi%d_mirror(1:3, imirror)
!        !end if
!        
!        dist2 = dot_product(v21,v21)
!        if(dist2 > cutoff2) cycle
!           
!        ene = ene + coef1 * coef_charge( lmp2charge(imp2),grep ) / sqrt(dist2)
!     end do
!
!     do itp2 = itp1+1, ntp
!        !if(inperi%i_periodic == 0) then
!        !   v21(1:3) = xyz_tp(1:3, itp2) - xyz_tp(1:3, itp1)
!        !else
!           vx(1:3) = xyz_tp(1:3, itp2) - xyz_tp(1:3, itp1)
!           call util_pbneighbor(vx, imirror)
!           v21(1:3) = vx(1:3) + inperi%d_mirror(1:3, imirror)
!        !end if
!        
!        dist2 = dot_product(v21,v21)
!        if(dist2 > cutoff2) cycle
!           
!        ene = ene + coef1 * charge_tp(itp2) / sqrt(dist2)
!     end do
!  end do
!!$omp end do nowait
!
!  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene
!
!end subroutine energy_ele_coulomb_tp
