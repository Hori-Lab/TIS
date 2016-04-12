! simu_force_implig
!> @brief This subroutine calculates the interaction force between the implicit ligand and the protein.

! calculate implicit ligand binding force
subroutine simu_force_implig(irep, force_mp)

  use const_maxsize
  use const_index   
  use const_physical
  use var_inp,    only : inperi
  use var_setp,   only : inpro
  use var_struct, only : nmp_real, xyz_mp_rep, pxyz_mp_rep, nmp_all
  use var_mgo,    only : inmgo, ishadow2real_mp_mgo
  use var_implig, only : inimplig, ncon_implig, icon2mp_implig, &
                         vdwrad_implig, istate_implig

  implicit none
  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)
  ! ------------------------------------------------------------

  integer :: iadd, ist, icon, icon2, imp1, imp2, jmp1, jmp2
  integer :: imirror
  real(PREC) :: pre_go, dist2, vdwrad2, vdwrad4
  real(PREC) :: pre_ist, pre_LJ_ist
  real(PREC) :: v21(3), for(3)
  real(PREC) :: vdwrad8, vdwrad12, vdwrad14, dist4, dist8, dist12, dist14
  real(PREC) :: roverdist12, roverdist14, dvdw_dr, dist1, distoverr_m1, dvdr
  real(PREC) :: vdwrad1_inv, gauss_2_inv
  ! ------------------------------------------------------------------------

  ! set parameter
  pre_go = inpro%cgo1210
  !  pre = pre_go*inimplig%pre_implig
  !  pre_LJ = 60.0e0_PREC*pre_go*inimplig%pre_implig 

  iadd = 0
  do ist=1,inimplig%nsite_implig
       pre_ist = pre_go*inimplig%pre_implig(ist)
       pre_LJ_ist = 60.0e0_PREC*pre_go*inimplig%pre_implig(ist) 

     !! if ligand is bound
     if(istate_implig(ist, irep) == IMPLIGBOUND_STATE%BOUND)then 
        do icon=1,ncon_implig(ist)
           icon2 = icon + iadd
           imp1 = icon2mp_implig(1,icon2)
           imp2 = icon2mp_implig(2,icon2)
           
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

           !! LJ12-10 type implicit ligand force
           if(inimplig%itype_ene_implig == IMPLIGENE_FUNCTYPE%FUNC_LJ12_10)then 
              vdwrad2=vdwrad_implig(icon2)**2
              vdwrad4 = vdwrad2 * vdwrad2
              vdwrad8 = vdwrad4 * vdwrad4
              vdwrad12 = vdwrad8 * vdwrad4
              vdwrad14 = vdwrad12 * vdwrad2
              
              dist4 = dist2 * dist2
              dist8 = dist4 * dist4
              dist12 = dist4 * dist8
              dist14 = dist12 * dist2
              
              roverdist12 = vdwrad12 /dist12
              roverdist14 = vdwrad14 / dist14               
              dvdw_dr=(pre_LJ_ist/vdwrad2)*(roverdist14 - roverdist12)
           
           !! Gaussian type implicit ligand force
           else if(inimplig%itype_ene_implig == IMPLIGENE_FUNCTYPE%FUNC_GAUSSIAN)then 
              dist1 = sqrt(dist2)

              ! vdwrad1 = vdwrad_implig(icon2)
              ! distoverr = dist1/vdwrad1
              ! dvdr =  pre_ist                                                            &
              !       * exp( -(distoverr-1.0)**2 / (2.0*inimplig%gauss_d_implig(ist)**2) ) &
              !       * (distoverr-1.0)                                                    &
              !       / (vdwrad1*inimplig%gauss_d_implig(ist)**2)              

              vdwrad1_inv = 1.0e0_PREC / vdwrad_implig(icon2)
              distoverr_m1 = dist1 * vdwrad1_inv - 1.0e0_PREC
              gauss_2_inv = 1.0e0_PREC / ( inimplig%gauss_d_implig(ist) ** 2 )
              dvdr =  pre_ist                                                &
                    * exp( -distoverr_m1**2 * 0.5e0_PREC * gauss_2_inv)      &
                    * distoverr_m1                                           &
                    * vdwrad1_inv * gauss_2_inv             
              dvdw_dr = -dvdr/dist1              
           end if
           for(1:3) = dvdw_dr*v21(1:3)
           
           if(inmgo%i_multi_mgo>=1)then
              jmp1 = ishadow2real_mp_mgo(imp1)
              jmp2 = ishadow2real_mp_mgo(imp2)
           else
              jmp1 = imp1
              jmp2 = imp2
           end if
           force_mp(1:3, jmp2) = force_mp(1:3, jmp2) + for(1:3)
           force_mp(1:3, jmp1) = force_mp(1:3, jmp1) - for(1:3)
        end do
     end if

     iadd = iadd + ncon_implig(ist)
  end do
  
end subroutine simu_force_implig
