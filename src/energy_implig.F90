! energy_implig
!> @brief Calculate implicit ligand binding energy.
!>     If (iflag_for_mc==IMPLIGENERGY_TYPE%FOR_MC) then &
!>     calculate implicit-ligand energy for MC_procedure (simu_mc_implig) &
!>     based on [ligand-binding site's coordinate].
!>
!>     If (iflag_for_mc==IMPLIGENERGY_TYPE%FOR_NON_MC) then &
!>     calculate implicit-ligand energy for output (NON_MC_procedure) &
!>     based on [istate_implig and ligand-binding site's coordinate].

subroutine energy_implig(irep, energy_unit, energy, iflag_for_mc)

  use const_maxsize 
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inpro
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, imp2unit, nunit_all
  use var_implig, only : inimplig, Etbind_implig, Ebind_implig, &
                         ncon_implig, icon2mp_implig,           &
                         vdwrad_implig, istate_implig

  implicit none

  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep             ! replica number 
  real(PREC), intent(inout) :: energy_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)
  real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX)
  integer,    intent(in)    :: iflag_for_mc     ! flag for Monte Carlo ((simu_mc_implig).
  ! if (iflag_for_mc==IMPLIGENERGY_TYPE%FOR_MC) then 
  !   calculate implicit-ligand energy for MC_procedure (simu_mc_implig).
  ! if (iflag_for_mc==IMPLIGENERGY_TYPE%FOR_NON_MC) then 
  !   calculate implicit-ligand energy for output (NON_MC_procedure).
  ! --------------------------------------------------------------------------

  integer :: ist, iadd, icon, icon2, imp1, imp2
  integer :: iunit, junit, imirror
  real(PREC) :: pre_go, dist2, vdwrad2, vdwrad4, vdwrad8
  real(PREC) :: pre_ist
  real(PREC) :: v21(3)
  real(PREC) :: vdwrad12, vdwrad10, dist4, dist8, dist12, dist10
  real(PREC) :: roverdist10, roverdist12, distoverr

  ! --------------------------------------------------------------------------
! set parameter
  pre_go = inpro%cgo1210
! write(*,*) 'pre_go=',pre_go
! pre = pre_go*inimplig%pre_implig
! write(*,*) 'pre=',pre


! initialization  (especially for_MC)
!       energy(E_TYPE%IMPLIG) = 0.0e0_PREC
!       energy_unit(iunit, junit, E_TYPE%IMPLIG) = 0.0e0_PREC
  if(iflag_for_mc == IMPLIGENERGY_TYPE%FOR_MC) then 
     energy(E_TYPE%IMPLIG) = 0.0e0_PREC
     energy_unit(1:nunit_all, 1:nunit_all, E_TYPE%IMPLIG) = 0.0e0_PREC
  endif
! FOR_NON_MC case, energy and energy_unit are already initialized at simu_energy.

! initialization
! Etbind_implig(MXSITE_IMPLIG, MXREPLICA)
  Etbind_implig(1:inimplig%nsite_implig, irep) = 0.0e0_PREC

  iadd = 0
  do ist = 1,inimplig%nsite_implig
     pre_ist = pre_go*inimplig%pre_implig(ist)

     if (iflag_for_mc == IMPLIGENERGY_TYPE%FOR_MC .OR.        &
         istate_implig(ist, irep) == IMPLIGBOUND_STATE%BOUND) then
        do icon=1, ncon_implig(ist)
           icon2 = icon + iadd
           imp1 = icon2mp_implig(1,icon2)
           imp2 = icon2mp_implig(2,icon2)
           iunit = imp2unit(imp1)
           junit = imp2unit(imp2)
           !! distance calculation for ligand mediated contact
           if(inperi%i_periodic == 0) then
              v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
           else
              v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
              call util_pbneighbor(v21, imirror)
           end if
           dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
           
!           dx = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
!           dy = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
!           dz = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)
!           dist2 = dx**2 + dy**2 + dz**2        
           
           !! LJ12-10 type implicit ligand binding energy
           if (inimplig%itype_ene_implig == IMPLIGENE_FUNCTYPE%FUNC_LJ12_10) then
              vdwrad2 = vdwrad_implig(icon2)**2
              vdwrad4 = vdwrad2**2
              vdwrad8 = vdwrad4**2
              vdwrad12 = vdwrad8 * vdwrad4
              vdwrad10 = vdwrad8 * vdwrad2
              
              dist4 = dist2 * dist2
              dist8 = dist4 * dist4
              dist12 = dist4 * dist8
              dist10 = dist2 * dist8
              
              roverdist10 = vdwrad10/dist10
              roverdist12 = vdwrad12/dist12
              Ebind_implig(icon2) = pre_ist*(5.0e0_PREC*roverdist12 - 6.0e0_PREC*roverdist10)
              
           !! Gaussian type implicit ligand binding energy
           elseif (inimplig%itype_ene_implig == IMPLIGENE_FUNCTYPE%FUNC_GAUSSIAN) then
              distoverr = sqrt(dist2) / vdwrad_implig(icon2)
              Ebind_implig(icon2) = -pre_ist * exp(-(distoverr-1.0e0_PREC)**2 &
                                                   /(2.0e0_PREC*inimplig%gauss_d_implig(ist)**2))
           endif
           Etbind_implig(ist, irep) = Etbind_implig(ist, irep) + Ebind_implig(icon2)
           
           if (iflag_for_mc == IMPLIGENERGY_TYPE%FOR_NON_MC) then
              energy(E_TYPE%IMPLIG) = energy(E_TYPE%IMPLIG) + Ebind_implig(icon2)
              energy_unit(iunit, junit, E_TYPE%IMPLIG) =  &
                   energy_unit(iunit, junit, E_TYPE%IMPLIG) + Ebind_implig(icon2)
           endif
       enddo
     end if

     iadd = iadd + ncon_implig(ist)     
  enddo

  ! for debug
  !  do ist=1,inimplig%nsite_implig
  !     write(*,*)'Binding site:',ist,'Ebind=',Etbind_implig(ist)
  ! enddo
  ! return
end subroutine energy_implig
