! energy_ele_coulomb
!> @brief Calculate the energy of electrostatic interaction 

subroutine energy_ele_coulomb_ewld(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inele, inperi
  use var_struct,  only : imp2unit, pxyz_mp_rep, lele, iele2mp, coef_ele, &
                          ncharge, coef_charge, icharge2mp
  use var_simu,    only : ewld_f_n, ewld_f_coef, ewld_f_rlv, ewld_s_sum
  use var_replica, only : irep2grep
  use mpiconst

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit, iele, imirror, ich1, ig, grep
  real(PREC) :: dist1, dist2, ssin, scos, ene, cutoff2, q1, dp
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  grep = irep2grep(irep)

  !================================================
  !================== Real space ==================
  !================================================
  cutoff2 = inele%cutoff_ele ** 2
  !cutoff2 = inperi%psizeh(1) ** 2   ! Assuming a cubic box
  !write(*,*) 'coef_ele(1)=',coef_ele(1,irep)
  !write(*,*) 'coef_ele(196)=',coef_ele(196,irep)
  !write(*,*) 'inele%coef=',inele%coef(grep)
  !write(*,*) 'energy(ELE)=',energy(E_TYPE%ELE)

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
  ene = 0.0
!$omp do private(imp1,imp2,v21,dist2,dist1,ene,iunit,junit,imirror)
  do iele=ksta, kend

     imp1 = iele2mp(1, iele, irep)
     imp2 = iele2mp(2, iele, irep)

     imirror = iele2mp(3, iele, irep)
     v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     
     dist2 = dot_product(v21,v21)
     if(dist2 > cutoff2) cycle
        
     ! --------------------------------------------------------------------
     dist1 = sqrt(dist2)
     
     ene = ene +  coef_ele(iele,irep ) * erfc(inele%ewld_alpha*dist1) / dist1
     if (iele < 20) then
         write(*,*) 'real',iele,dist1,coef_ele(iele,irep ) * erfc(inele%ewld_alpha*dist1) / dist1
     endif

     !iunit = imp2unit(imp1)
     !junit = imp2unit(imp2)
     !energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + ene
  end do
  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene
  !write(*,*) 'Real:',energy(E_TYPE%ELE)
  !write(*,*) 'energy(ELE)=',energy(E_TYPE%ELE)
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
!        ene = inele%coef(grep) * q1 * q2 * erfc(inele%ewld_alpha*dist1) / dist1
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
      
     ene = ene + ewld_f_coef(ig) * (scos * scos + ssin * ssin)
     if (ig <= 10 .or. ig > ewld_f_n - 10) then
        write(*,*) 'ewld_f_coef:', ig, ewld_f_coef(ig), ewld_f_coef(ig) * (scos * scos + ssin * ssin)
     endif
  end do
  !write(*,*) 'Fourier:',ene

  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + inele%coef(grep) * ene
  !write(*,*) 'energy(ELE)=',energy(E_TYPE%ELE)


  !================================================
  !======= Correction for self interaction ========
  !================================================

  !ene = 0.0
  !do ich1 = 1, ncharge
  !   q1 = coef_charge(ich1, irep)
  !   ene = ene + q1**2
  !end do
  !write(*,*) 'Self:',inele%coef(grep), ewld_s_coef,ene, - inele%coef(grep) * ewld_s_coef * ene
  !ene = - inele%coef(grep) * ewld_s_coef * ene
  !ene = - inele%coef(grep) * ewld_s_sum

  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + inele%coef(grep) * ewld_s_sum
  !write(*,*) 'energy(ELE)=',energy(E_TYPE%ELE)

end subroutine energy_ele_coulomb_ewld


subroutine energy_ele_coulomb_ewld_tp(irep, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inele, inperi
  use var_struct,  only : pxyz_mp_rep, coef_charge, nmp_real, lmp2charge, &
                          ntp, xyz_tp, charge_tp
  use var_replica, only : irep2grep
  use var_simu,    only : ewld_f_n, ewld_f_coef, ewld_f_rlv, ewld_s_coef
  use mpiconst

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  integer :: itp1, itp2, imp2, imirror, ig, grep
  real(PREC) :: dist1, dist2, ssin, scos, ene, cutoff2, q1, dp
  real(PREC) :: v21(SDIM), vx(SDIM)

  grep = irep2grep(irep)

  !================================================
  !================== Real space ==================
  !================================================
  cutoff2 = inele%cutoff_ele ** 2
  !cutoff2 = inperi%psizeh(1) ** 2   ! Assuming a cubic box

  ene = 0.0
!!$omp do private(q1,imp2,v21,vx,ene,dist2,dist1,ene,imirror)
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
  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene


  !================================================
  !================= Fourier space ================
  !================================================
  ene = 0.0
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
  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene


  !================================================
  !======= Correction for self interaction ========
  !================================================
  ene = 0.0
  do itp1 = 1, ntp
     q1 = charge_tp(itp1)
     ene = ene + q1**2
  end do
  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ewld_s_coef * ene

  ! Multiply coef to the total
  energy(E_TYPE%ELE) = inele%coef(grep) * energy(E_TYPE%ELE)

end subroutine energy_ele_coulomb_ewld_tp