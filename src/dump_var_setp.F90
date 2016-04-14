subroutine dump_var_setp(lunout)

  use const_maxsize
  use var_setp

  integer, intent(in) :: lunout
  integer :: i, j

  write(lunout,*) '### var_setp'

!  type input_parameter
!     ! general.para
!     real(PREC) :: velo_adjst
!     real(PREC) :: rneighbor_dist
!     real(PREC) :: cmass
!     real(PREC) :: fric_const
!     ! protein.para
!     real(PREC) :: energy_unit_protein
!     real(PREC) :: cbd
!     real(PREC) :: cba
!     real(PREC) :: cdih_1
!     real(PREC) :: cdih_3
!     integer :: n_sep_nlocal
!     integer :: n_sep_contact
!     real(PREC) :: cutoff_go
!     real(PREC) :: cutoff_exvol
!     real(PREC) :: dfcontact
!     real(PREC) :: cgo1210
!     real(PREC) :: cdist_rep12
!     real(PREC) :: crep12
!  end type input_parameter
!  type(input_parameter), save :: inpara

  write(lunout,*) 'inpara % velo_adjst,',inpara%velo_adjst
  write(lunout,*) 'inpara % rneighbor_dist,',inpara%rneighbor_dist
  write(lunout,*) 'inpara % cmass,',inpara%cmass
  write(lunout,*) 'inpara % fric_const,',inpara%fric_const
  write(lunout,*) 'inpro % energy_unit_protein,',inpro%energy_unit_protein
  write(lunout,*) 'inpro % cbd,',inpro%cbd
  write(lunout,*) 'inpro % cba,',inpro%cba
  write(lunout,*) 'inpro % cdih_1,',inpro%cdih_1
  write(lunout,*) 'inpro % cdih_3,',inpro%cdih_3
  write(lunout,*) 'inpro % n_sep_nlocal,',inpro%n_sep_nlocal
  write(lunout,*) 'inpro % n_sep_contact,',inpro%n_sep_contact
  write(lunout,*) 'inpro % cutoff_go,',inpro%cutoff_go
  write(lunout,*) 'inpro % cutoff_exvol,',inpro%cutoff_exvol
  write(lunout,*) 'inpro % dfcontact,',inpro%dfcontact
  write(lunout,*) 'inpro % cgo1210,',inpro%cgo1210
  write(lunout,*) 'inpro % cdist_rep12,',inpro%cdist_rep12
  write(lunout,*) 'inpro % crep12,',inpro%crep12

  write(lunout,*) 'inele % i_diele,', inele%i_diele
  write(lunout,*) 'inele % cutoff_ele,', inele%cutoff_ele
  write(lunout,*) 'inele % ionic_strength,', inele%ionic_strength
  do i = 1, MXREPLICA
     write(lunout,*) 'inele % cdist(',i,'),', inele%cdist(i)
  enddo
  do i = 1, MXREPLICA
     write(lunout,*) 'inele % coef(',i,'),', inele%coef(i)
  enddo

!  type input_ligandparameter
!     ! ligand.para
!     real(PREC) :: energy_unit
!     real(PREC) :: cbd
!     real(PREC) :: cba
!     real(PREC) :: cdih
!     real(PREC) :: cutoff_exvol
!     real(PREC) :: cdist_rep12_lpro
!     real(PREC) :: cdist_rep12_llig
!     real(PREC) :: crep12
!     
!     integer :: sz
!
!  end type input_ligandparameter
!
!  type(input_ligandparameter), save :: inligand
  write(lunout, *) 'inligand % energy_unit', inligand%energy_unit
  write(lunout, *) 'inligand % cbd', inligand%cbd
  write(lunout, *) 'inligand % cba', inligand%cba
  write(lunout, *) 'inligand % cdih', inligand%cdih
  write(lunout, *) 'inligand % cutoff_exvol', inligand%cutoff_exvol
  write(lunout, *) 'inligand % cdist_rep12_lpro', inligand%cdist_rep12_lpro
  write(lunout, *) 'inligand % cdist_rep12_llig', inligand%cdist_rep12_llig
  write(lunout, *) 'inligand % crep12', inligand%crep12
  write(lunout, *) 'inligand % sz', inligand%sz

!  type input_hpparameter
!     ! hydrophobic.para
!     real(PREC) :: rho_min_hp
!     real(PREC) :: coef_rho_hp
!     real(PREC) :: coef_hp
!     real(PREC) :: coefaa_para_hp(21)
!     real(PREC) :: ncoor_para_hp(21)
!     real(PREC) :: ncoormax_para_hp(21)
!     real(PREC) :: cutoffdmin_para_hp(21,21)
!     real(PREC) :: cutoffdmax_para_hp(21,21)
!
!     integer :: sz
!
!  end type input_hpparameter
!
!  type(input_hpparameter), save :: inhp
  write(lunout, *) 'inhp % rho_min_hp', inhp%rho_min_hp
  write(lunout, *) 'inhp % coef_rho_hp', inhp%coef_rho_hp
  write(lunout, *) 'inhp % coef_hp', inhp%coef_hp
  write(lunout,*) 'inhp % coefaa_para_hp(1:21)', inhp%coefaa_para_hp(1:21)
  write(lunout,*) 'inhp % ncoor_para_hp(1:21)', inhp%ncoor_para_hp(1:21)
  write(lunout,*) 'inhp % ncoormax_para_hp(1:21)', inhp%ncoormax_para_hp(1:21)
  do i = 1, 21 
     write(lunout,*) 'inhp % cutoffdmin_para_hp(',i,', 1:21)',&
     inhp%cutoffdmin_para_hp(i, 1:21)
  enddo
  do i = 1, 21 
     write(lunout,*) 'inhp % cutoffdmax_para_hp(',i,', 1:21)',&
     inhp%cutoffdmax_para_hp(i, 1:21)
  enddo
  write(lunout, *) 'inhp% sz', inhp%sz

!  type input_simuparameter
!     integer :: n_step_sim
!     integer :: n_tstep(MXSIM)
!     integer :: n_step_save
!     integer :: n_step_neighbor
!     real(PREC) :: tstep_size
!
!     integer :: i_com_zeroing
!     integer :: i_no_trans_rot
!
!     real(PREC) :: tempk
!     integer :: n_seed
!  end type input_simuparameter
!  type(input_simuparameter), save :: insimu

  write(lunout,*) 'insimu % n_step_sim,', insimu%n_step_sim
  do i = 1, MXSIM
     write(lunout,*) 'insimu % n_tstep(',i,'),', insimu%n_tstep(i)
  enddo
  write(lunout,*) 'insimu % n_step_save,', insimu%n_step_save
  write(lunout,*) 'insimu % n_step_neighbor,', insimu%n_step_neighbor
  write(lunout,*) 'insimu % tstep_size,', insimu%tstep_size
  write(lunout,*) 'insimu % i_com_zeroing,', insimu%i_com_zeroing
  write(lunout,*) 'insimu % i_no_trans_rot,', insimu%i_no_trans_rot
  write(lunout,*) 'insimu % tempk,', insimu%tempk
  write(lunout,*) 'insimu % n_seed,', insimu%n_seed

!  integer, save :: irand
  write(lunout,*) 'irand,',irand

!  type input_annealing
!     integer :: n_time_change
!     real(PREC) :: tempk_init
!     real(PREC) :: tempk_last
!  end type input_annealing
!  type(input_annealing), save :: inann

  write(lunout,*) 'inann % n_time_change,', inann%n_time_change
  write(lunout,*) 'inann % tempk_init,', inann%tempk_init
  write(lunout,*) 'inann % tempk_last,', inann%tempk_last

!  type input_searchingtf
!     real(PREC) :: tempk_upper
!     real(PREC) :: tempk_lower
!  end type input_searchingtf
!  type(input_searchingtf), save :: insear

  write(lunout,*) 'insear % tempk_upper,', insear%tempk_upper
  write(lunout,*) 'insear % tempk_lower,', insear%tempk_lower

!  type input_miscellaneous
!     integer :: i_use_atom_protein
!     integer :: itype_nlocal(MXUNIT, MXUNIT)
!     logical :: force_flag(INTERACT%MAX)
!
!     ! redefine_parameter
!     integer :: i_redef_para
!
!     ! box_interaction
!     integer :: i_in_box
!     real(PREC) :: xbox
!     real(PREC) :: ybox
!     real(PREC) :: zbox
!     real(PREC) :: boxsigma
!
!     ! delete interaction
!     integer :: i_del_int
!     integer :: ndel_lgo
!     integer :: ndel_go
!     integer :: idel_lgo(2, MXDEL_LGO)
!     integer :: idel_go(4, MXDEL_GO)
!
!     ! energy coefficient
!     integer :: i_energy_para
!     real(PREC) :: factor_go_unit(MXUNIT, MXUNIT)
!     real(PREC) :: factor_local_unit(MXUNIT, MXUNIT)
!
!     ! mass and friction
!     integer :: i_redef_mass_fric
!
!     ! neighbordist
!     integer :: i_neigh_dist
!     real(PREC) :: rneighbordist2_unit(MXUNIT, MXUNIT)
!
!     ! bridge
!     integer :: i_bridge
!     integer :: nbrid
!     integer :: ibrid2mp(2, MXBRIDGE)
!     real(PREC) :: coef_brid(MXBRIDGE)
!     real(PREC) :: brid_dist(MXBRIDGE)
!
!     ! pulling
!     integer :: i_pulling
!     integer :: npull
!     integer :: ipull2mp(MXPULLING)
!     real(PREC) :: coef_pull(MXPULLING)
!     real(PREC) :: pull_xyz(3, MXPULLING)
!     real(PREC) :: pu_xyz(3, MXPULLING)
!
!     ! anchor
!     integer :: i_anchor
!     integer :: nanc
!     integer :: ianc2mp(MXANCHOR)
!     real(PREC) :: coef_anc(MXANCHOR)
!     real(PREC) :: anc_dist(MXANCHOR)
!     real(PREC) :: anc_xyz(3, MXANCHOR)
!
!     ! fix
!     integer :: i_fix
!
!  end type input_miscellaneous
!  type(input_miscellaneous), save :: inmisc

  write(lunout,*) 'inmisc % i_use_atom_protein,', inmisc%i_use_atom_protein
!  do j = 1, MXUNIT
!     do i = 1, MXUNIT
!        write(lunout,*) 'inmisc % itype_nlocal(',i,',',j,'),', inmisc%itype_nlocal(i,j)
!     enddo
!  enddo
  do i = 1, INTERACT%MAX
     write(lunout,*) 'inmisc % force_flag(',i,'),', inmisc%force_flag(i)
  enddo
  write(lunout,*) 'inmisc % i_redef_para,', inmisc%i_redef_para
  write(lunout,*) 'inmisc % i_in_box,', inmisc%i_in_box
  write(lunout,*) 'inmisc % xbox,', inmisc%xbox
  write(lunout,*) 'inmisc % ybox,', inmisc%ybox
  write(lunout,*) 'inmisc % zbox,', inmisc%zbox
  write(lunout,*) 'inmisc % boxsigma,', inmisc%boxsigma
  write(lunout,*) 'inmisc % i_in_cap,', inmisc%i_in_cap
  write(lunout,*) 'inmisc % rcap,', inmisc%rcap
  write(lunout,*) 'inmisc % kcap,', inmisc%kcap
  write(lunout,*) 'inmisc % center_cap,', inmisc%center_cap(1), inmisc%center_cap(2), inmisc%center_cap(3)
  write(lunout,*) 'inmisc % boxsigma,', inmisc%boxsigma
  write(lunout,*) 'inmisc % i_del_int,', inmisc%i_del_int
  write(lunout,*) 'inmisc % ndel_lgo,', inmisc%ndel_lgo
  write(lunout,*) 'inmisc % ndel_go,', inmisc%ndel_go
  do j = 1, MXDEL_LGO
     write(lunout,*) 'inmisc % idel_lgo(:,',j,'),', inmisc%idel_lgo(:,j)
  enddo
  do j = 1, MXDEL_GO
     write(lunout,*) 'inmisc % idel_go(:,',j,'),', inmisc%idel_go(:,j)
  enddo
  write(lunout,*) 'inmisc % i_energy_para,', inmisc%i_energy_para
  do j = 1, MXUNIT
     do i = 1, MXUNIT
        write(lunout,*) 'inmisc % factor_go_unit(',i,',',j,'),', inmisc%factor_go_unit(i,j)
     enddo
  enddo
  do j = 1, MXUNIT
     do i = 1, MXUNIT
        write(lunout,*) 'inmisc % factor_local_unit(',i,',',j,'),', inmisc%factor_local_unit(i,j)
     enddo
  enddo
  write(lunout,*) 'inmisc % i_redef_mass_fric,', inmisc%i_redef_mass_fric
  write(lunout,*) 'inmisc % i_neigh_dist,', inmisc%i_neigh_dist
  do j = 1, MXUNIT
     do i = 1, MXUNIT
        write(lunout,*) 'inmisc % rneighbordist2_unit(',i,',',j,'),', inmisc%rneighbordist2_unit(i,j)
     enddo
  enddo
  write(lunout,*) 'inmisc % i_bridge,', inmisc%i_bridge
  write(lunout,*) 'inmisc % nbrid,', inmisc%nbrid
  do j = 1, MXBRIDGE
     write(lunout,*) 'inmisc % ibrid2mp(:,',j,'),', inmisc%ibrid2mp(:,j)
  enddo
  do j = 1, MXBRIDGE
     write(lunout,*) 'inmisc % coef_brid(',j,'),', inmisc%coef_brid(j)
  enddo
  do j = 1, MXBRIDGE
     write(lunout,*) 'inmisc % brid_dist(',j,'),', inmisc%brid_dist(j)
  enddo
  write(lunout,*) 'inmisc % i_pulling,', inmisc%i_pulling
  write(lunout,*) 'inmisc % npull,', inmisc%npull
  write(lunout,*) 'inmisc % ipull2mp(MXPULLING),', inmisc%ipull2mp(MXPULLING)
  write(lunout,*) 'inmisc % coef_pull(MXPULLING),', inmisc%coef_pull(MXPULLING)
  write(lunout,*) 'inmisc % pull_xyz(3,MXPULLING),', inmisc%pull_xyz(3,MXPULLING)
  write(lunout,*) 'inmisc % pu_xyz(3,MXPULLING),', inmisc%pu_xyz(3,MXPULLING)
  write(lunout,*) 'inmisc % i_anchor,', inmisc%i_anchor
  write(lunout,*) 'inmisc % nanc,', inmisc%nanc
  write(lunout,*) 'inmisc % ianc2mp(MXANCHOR),', inmisc%ianc2mp(MXANCHOR)
  write(lunout,*) 'inmisc % coef_anc(MXANCHOR),', inmisc%coef_anc(MXANCHOR)
  write(lunout,*) 'inmisc % anc_dist(MXANCHOR),', inmisc%anc_dist(MXANCHOR)
  write(lunout,*) 'inmisc % anc_xyz(3,MXANCHOR),', inmisc%anc_xyz(3,MXANCHOR)
  write(lunout,*) 'inmisc % i_fix,', inmisc%i_fix

!  integer :: ifix_mp(MXMP)
  do i = 1, MXMP
     write(lunout,*) 'ifix_mp(',i,'),',ifix_mp(i)
  enddo

endsubroutine dump_var_setp

