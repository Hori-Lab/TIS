! ***********************************************************************
subroutine simu_stream_mpc2(velo_mp, irep)

  use const_maxsize
  use const_physical
  use var_setp, only : inpara, insimu, ifix_mp
  use var_struct, only : nmp_real, cmass_mp, xyz_mp_rep
  use var_simu, only : istep, tstep, tempk
  use var_mpc, only : inmpc, cmass_solv_mpc, velo_solv_mpc, xyz_solv_mpc
  implicit none
      
  ! --------------------------------------------------------------------
  real(PREC), intent(inout) :: velo_mp(:,:,:)
  integer,    intent(in) :: irep

  ! --------------------------------------------------------------------
  !local variables
  integer, parameter :: ncs = MXCS
  integer :: idimn, is, imp, n_solv
  integer :: i_count_mp_All
  integer :: i, njiyu
  integer :: ics, jcs
  integer, save :: init_nh_mpc = 0
  real(PREC) :: tstep_colli
  real(PREC) :: P_total_mpc(3), velo_cm_mpc(3), total_mass_mpc
  real(PREC) :: rescale_fact_All, total_kine_E2_All
  real(PREC) :: tke, tsteph
  real(PREC) :: xtrans, vtrans, accum, ev_nh2
  real(PREC) :: velo_yojou(ncs), ev_nh(ncs)
  real(PREC), save :: cmass_nh(ncs), velo_nh(ncs)


  ! --------------------------------------------------------------------
  !!initialization
  if(init_nh_mpc == 0) then
     init_nh_mpc = 1
!     cmass_nh(1:ncs) = 100000.0
     cmass_nh(1:ncs) = inpara%csmass_mpc_per*inmpc%n_all_solv
     velo_nh(1:ncs) = sqrt(BOLTZC * tempk / cmass_nh(1:ncs))
  end if


  ! --------------------------------------------------------------------
  ! total-linear momentum should be zero.
  n_solv = inmpc%n_all_solv
  P_total_mpc(1:3) = 0.0e0_PREC
  total_mass_mpc = 0.0e0_PREC
  
  do is = 1, n_solv
     P_total_mpc(1:3) = P_total_mpc(1:3) &
          + cmass_solv_mpc(is) * velo_solv_mpc(1:3, is)
     total_mass_mpc = total_mass_mpc + cmass_solv_mpc(is)
  end do
  
  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     P_total_mpc(1:3) = P_total_mpc(1:3) &
          + cmass_mp(imp) * velo_mp(1:3, imp, irep)
     total_mass_mpc = total_mass_mpc + cmass_mp(imp)
  end do

  
  velo_cm_mpc(1:3) = P_total_mpc(1:3) / total_mass_mpc
  
  !!velocity correction
  do is = 1, n_solv
     velo_solv_mpc(1:3, is) = velo_solv_mpc(1:3, is) - velo_cm_mpc(1:3)
  end do

  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     velo_mp(1:3, imp, irep) = velo_mp(1:3, imp, irep) - velo_cm_mpc(1:3)
  end do
  
  
  ! --------------------------------------------------------------------
  ! streaming phase
  i_count_mp_All = 0
  total_kine_E2_All = 0.0e0_PREC 
  
  !!<<for mpc solvent particle>>   
  n_solv = inmpc%n_all_solv
  do is = 1, n_solv
     total_kine_E2_All = total_kine_E2_All + cmass_solv_mpc(is) &
          * (velo_solv_mpc(1, is)**2 &
          + velo_solv_mpc(2, is)**2 + velo_solv_mpc(3, is)**2)
  end do
  i_count_mp_All = i_count_mp_All + 3*n_solv
  
  if(insimu%i_no_trans_rot == 1) then
     njiyu = i_count_mp_All - 9
  else
     njiyu = i_count_mp_All - 3
  end if
  

  if(inmpc%i_thermal_mpc == 1) then

     !!------------------------------------------
     !!calculate rescaling factor for solvent particle velocity to reset T
     !!-----------------------------------------
     if(mod(istep, inmpc%nratio_vcorrect_step) == 0) then
        !!-----------------------------------------
        !! velocity correct
        !!-----------------------------------------
        rescale_fact_All = sqrt(BOLTZC*tempk*njiyu/total_kine_E2_All)
        do is = 1, n_solv
           velo_solv_mpc(1:3, is) = velo_solv_mpc(1:3, is) * rescale_fact_All
        end do
     end if
     
     ! -----------------------------------------------------------------
     !! av_xyz_coord_solv(1:3)=0.0
     !! av2_xyz_coord_solv(1:3)=0.0
     !! free-streaming process for mpc solvent particle
     tstep_colli = inmpc%nratio_colli_step*tstep
     do is = 1, n_solv
        xyz_solv_mpc(1:3, is) = xyz_solv_mpc(1:3, is) &
             + velo_solv_mpc(1:3, is)*tstep_colli
     end do

  else
     tsteph = tstep
     tke = total_kine_E2_All
     velo_yojou(1) = tke - njiyu*BOLTZC*tempk
     do ics = 2, ncs
        velo_yojou(ics) = cmass_nh(ics-1) * velo_nh(ics-1)**2 - BOLTZC * tempk
     end do

     vtrans = 1
     xtrans = 0
     accum = 1
     do i = 1, inmpc%nratio_colli_step

        velo_nh(ncs) = velo_nh(ncs) + tsteph * velo_yojou(ncs) / cmass_nh(ncs)
        ev_nh(ncs) = exp(-tsteph * velo_nh(ncs))
        do ics = 1, ncs-1
           jcs = ncs - ics
           velo_nh(jcs) = ev_nh(jcs+1) * velo_nh(jcs) + tsteph * velo_yojou(jcs) / cmass_nh(jcs)
           ev_nh(jcs) = exp(-tsteph * velo_nh(jcs))
        end do

        ev_nh2 = ev_nh(1)*ev_nh(1)
        tke = ev_nh2*ev_nh2*tke
        velo_yojou(1) = tke - njiyu*BOLTZC*tempk

        vtrans = ev_nh2*vtrans
        if(i == 1) then
           xtrans = ev_nh(1)
           accum = ev_nh2
        else
           xtrans = xtrans + ev_nh(1)*accum
           accum = ev_nh2*accum
        end if

        do ics = 1, ncs - 1
           velo_nh(ics) = ev_nh(ics + 1) * (velo_nh(ics) + tsteph * velo_yojou(ics) / cmass_nh(ics))
           velo_yojou(ics + 1) = cmass_nh(ics) * velo_nh(ics)**2 - BOLTZC * tempk
        end do
        velo_nh(ncs) = velo_nh(ncs) + tsteph * velo_yojou(ncs) / cmass_nh(ncs)

     end do

     xtrans = xtrans * tstep
     
     do is = 1, n_solv
        xyz_solv_mpc(1:3, is) = xyz_solv_mpc(1:3, is) + xtrans * velo_solv_mpc(1:3, is)
     end do

     do is = 1, n_solv
        velo_solv_mpc(1:3, is) = velo_solv_mpc(1:3, is) * vtrans
     end do

  end if


end subroutine simu_stream_mpc2
