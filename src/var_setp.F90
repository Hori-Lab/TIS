!var_setp
!> @brief Contains coefficients of energy function and some other  &
!>        parameters almost from parameter files.   

module var_setp

  use const_maxsize
  use const_index
  use mt_stream

  !==========================================
  !> structure for parameters reading from "para_cafemol_gen" field in "general.para" file
  type input_parameter 
     real(PREC) :: velo_adjst !< the coupling parameter of Berendsen thermostat
     real(PREC) :: csmass_per !< mass of thermal particle for Nose-Hoover thermostat
     real(PREC) :: rneighbor_dist !< the cutoff distance to define a neighbor
     real(PREC) :: cmass(0:CHEMICALTYPE%MAX) !< particle mass (cmass(0) is the default value.)
     real(PREC) :: fric_const !< friction constant for Langevin dynamics simulation
     real(PREC) :: radius(0:CHEMICALTYPE%MAX)
     real(PREC) :: viscosity  !< viscosity for Langevin dynamics when i_fric = 1 
     integer    :: sz !< size of the structure
  end type input_parameter
  type(input_parameter), save :: inpara

  !==========================================
  !> structure for parameters reading from "para_cafemol_pro" field in "protein.para" file
  type input_protein_parameter
     real(PREC) :: energy_unit_protein !< the energy unit for protein modeling (kcal/mol)
     real(PREC) :: cbd !< constant coefficient of the energy function for bond length
     real(PREC) :: cba !< constant coefficient of the energy function for bond angle
     real(PREC) :: cdih_1 !< constant coefficient of the 2pi-period term of the energy function for dihedral angle
     real(PREC) :: cdih_3 !< constant coefficient of the 2pi/3-period term of the energy function for dihedral angle
     integer :: n_sep_nlocal !< the minimum number of amino acids that separate a non-local pair
     integer :: n_sep_contact !< the minimum number of amino acids that separate a contact pair
     real(PREC) :: cutoff_go !< truncation distance for computing non-local Go interaction
     real(PREC) :: cutoff_LJ !< truncation distance for computing non-local Go interaction
     real(PREC) :: cutoff_exvol !< truncation distance for computing nonlocal non-native repulsion
     real(PREC) :: dfcontact !< the cutoff distance to define the native contact
     real(PREC) :: cgo1210 !< constant coefficient "go of the energy function for non-local Go interaction
     real(PREC) :: cdist_rep6  !< reference distance d in the non-native repulsive interaction
     real(PREC) :: cdist_rep12 !< reference distance d in the non-native repulsive interaction
     real(PREC) :: crep12 !< constant coefficient "ev in the non-native repulsive interaction
     real(PREC) :: crep6  !< constant coefficient "ev in the non-native repulsive interaction
     integer    :: sz !< size of the structure
  endtype input_protein_parameter
  type(input_protein_parameter), save :: inpro


  !==========================================
  !> structure for parameters reading from "para_cafemol_ion" field in "ion.para" file
  type input_ion_parameter

     integer :: num_na_ion
     integer :: num_k_ion
     integer :: num_cl_ion
     integer :: num_mg_ion
     integer :: num_ion(IONTYPE%MAX_ION)
     character(3) :: char_ion(IONTYPE%MAX_ION)

     ! parameters
     real(PREC) :: energy_unit_ion

     ! exv
     real(PREC) :: cexv_ion
     real(PREC) :: cdist_exv_ion
     real(PREC) :: cutoff_exv_ion

     ! LJ+Hydration
     real(PREC) :: cutoff_lj_ion
     real(PREC) :: cutoff_hyd_ion

     real(PREC) :: clj(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cdistlj(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cdistme(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: csigmame(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cdistmh1(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: csigmamh1(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cmh1(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cdistmh2(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: csigmamh2(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cmh2(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)

     ! LJ force, energy
     real(PREC) :: cdistlj2(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cutofflj2(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: clj_force(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: clj_energy(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)

     ! hydration force, energy
     real(PREC) :: rsigmamh1(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cmh1_force(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cmh1_energy(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: rsigmamh2(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cmh2_force(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)
     real(PREC) :: cmh2_energy(IONTYPE%MAX_ALL, IONTYPE%MAX_ALL)

     integer    :: sz
  end type input_ion_parameter

  type(input_ion_parameter), save :: inion

  !==========================================
  !> structure for parameters reading from "para_cafemol_rna" field in "rna.para" file
  type input_rna_parameter
     real(PREC) :: energy_unit
     integer :: i_use_atom_base
     integer :: i_use_atom_sugar

     real(PREC) :: cbd_PS
     real(PREC) :: cbd_SR
     real(PREC) :: cbd_SY
     real(PREC) :: cbd_SP

     real(PREC) :: cba_PSR
     real(PREC) :: cba_PSY
     real(PREC) :: cba_PSP
     real(PREC) :: cba_RSP
     real(PREC) :: cba_YSP
     real(PREC) :: cba_SPS
!     real(PREC) :: cba_BSB

     real(PREC) :: cdih_1_PSPS
     real(PREC) :: cdih_1_SPSR
     real(PREC) :: cdih_1_SPSY
     real(PREC) :: cdih_1_SPSP
     real(PREC) :: cdih_1_RSPS
     real(PREC) :: cdih_1_YSPS
     real(PREC) :: cdih_3_PSPS
     real(PREC) :: cdih_3_SPSR
     real(PREC) :: cdih_3_SPSY
     real(PREC) :: cdih_3_SPSP
     real(PREC) :: cdih_3_RSPS
     real(PREC) :: cdih_3_YSPS

     real(PREC) :: dfhelix_BSSB_lower
     real(PREC) :: dfhelix_BSSB_upper

     integer :: n_sep_nlocal_P
     integer :: n_sep_nlocal_S
     integer :: n_sep_nlocal_B
     integer :: n_sep_contact_P
     integer :: n_sep_contact_S
     integer :: n_sep_contact_B

     real(PREC) :: cutoff_go
     real(PREC) :: cutoff_bp
     real(PREC) :: cutoff_exvol

     real(PREC) :: dfcontact
     real(PREC) :: dfcontact_pro
     real(PREC) :: cgo1210_pro_P
     real(PREC) :: cgo1210_pro_S
     real(PREC) :: cgo1210_pro_B
     real(PREC) :: cgo1210_P_P
     real(PREC) :: cgo1210_P_S
     real(PREC) :: cgo1210_P_B
     real(PREC) :: cgo1210_S_S
     real(PREC) :: cgo1210_S_B
     real(PREC) :: cgo1210_B_B
     real(PREC) :: cgomorse_D_pro_P
     real(PREC) :: cgomorse_D_pro_S
     real(PREC) :: cgomorse_D_pro_B
     real(PREC) :: cgomorse_D_P_P
     real(PREC) :: cgomorse_D_P_S
     real(PREC) :: cgomorse_D_P_B
     real(PREC) :: cgomorse_D_S_S
     real(PREC) :: cgomorse_D_S_B
     real(PREC) :: cgomorse_D_B_B
     real(PREC) :: cgomorse_a_pro_P
     real(PREC) :: cgomorse_a_pro_S
     real(PREC) :: cgomorse_a_pro_B
     real(PREC) :: cgomorse_a_P_P
     real(PREC) :: cgomorse_a_P_S
     real(PREC) :: cgomorse_a_P_B
     real(PREC) :: cgomorse_a_S_S
     real(PREC) :: cgomorse_a_S_B
     real(PREC) :: cgomorse_a_B_B
     integer    :: i_potential_go

     real(PREC) :: dfcontact_bp
     real(PREC) :: cbp1210_HB2
     real(PREC) :: cbp1210_HB3
     real(PREC) :: cbpmorse_D
     real(PREC) :: cbpmorse_a
     integer    :: i_potential_bp

     integer    :: n_sep_base_stack
     real(PREC) :: dfcontact_st
     real(PREC) :: cst1210
     real(PREC) :: cstmorse_D
     real(PREC) :: cstmorse_a
     integer    :: i_potential_st

     real(PREC) :: cdist_rep12
     real(PREC) :: crep12
     
     integer    :: sz
  endtype input_rna_parameter
  type(input_rna_parameter), save :: inrna

  !==========================================
  !> structure for parameters reading from "para_cafemol_DT_rna" field in "rna.para" file
  type input_dtrna_parameter
     real(PREC) :: energy_unit
     integer :: i_use_atom_base
     integer :: i_use_atom_sugar

     real(PREC) :: bd_PS
     real(PREC) :: bd_SB
     real(PREC) :: bd_SP

     real(PREC) :: ba_PSB
     real(PREC) :: ba_PSP
     real(PREC) :: ba_BSP
     real(PREC) :: ba_SPS

     real(PREC) :: exv_dist
     real(PREC) :: exv_coef   ! 2013
     real(PREC) :: exv_rad(DT15EXV%MAX)  ! 2015
     real(PREC) :: exv_eps(DT15EXV%MAX)  ! 2015
     real(PREC) :: exv_adjust ! 2015
     real(PREC) :: exv_inf    ! 2015 for Widom insertion
     integer :: n_sep_nlocal_P
     integer :: n_sep_nlocal_S
     integer :: n_sep_nlocal_B

     real(PREC) :: st_dist
     real(PREC) :: st_dih
     real(PREC) :: st_h(16)    ! 1-16 = nucleotide types of stack pair 
     real(PREC) :: st_s(16)    !      (NN in const_index.F90)
     real(PREC) :: st_Tm(16)   ! U0 = - h + kB (T - Tm) * s
     real(PREC) :: st_r0(16)   ! 2015

     real(PREC) :: st_nlocal_dist
     real(PREC) :: st_nlocal_angl
     real(PREC) :: st_nlocal_dih
     real(PREC) :: st_nlocal_u0

     real(PREC) :: hb_dist
     real(PREC) :: hb_angl
     real(PREC) :: hb_dih_hbond
     real(PREC) :: hb_dih_chain
     real(PREC) :: hb_u0

     integer    :: sz
  endtype input_dtrna_parameter
  type(input_dtrna_parameter), save :: indtrna15
  type(input_dtrna_parameter), save :: indtrna13


  !==========================================
  !> structure for parameters reading from "A-form_RNA" in "rna.para" file
  type input_aform_rna
     real(PREC) :: bond_SP
     real(PREC) :: bond_PS
     real(PREC) :: bond_SA
     real(PREC) :: bond_SU
     real(PREC) :: bond_SG
     real(PREC) :: bond_SC
     real(PREC) :: angl_PSP
     real(PREC) :: angl_SPS
     real(PREC) :: angl_PSA
     real(PREC) :: angl_PSU
     real(PREC) :: angl_PSG
     real(PREC) :: angl_PSC
     real(PREC) :: angl_ASP
     real(PREC) :: angl_USP
     real(PREC) :: angl_GSP
     real(PREC) :: angl_CSP
     real(PREC) :: dihd_PSPS
     real(PREC) :: dihd_SPSP
     real(PREC) :: stack_dist(16)
     real(PREC) :: hbond_dist_AU
     real(PREC) :: hbond_dist_GC
     real(PREC) :: hbond_angl_SAU
     real(PREC) :: hbond_angl_SUA
     real(PREC) :: hbond_angl_SGC
     real(PREC) :: hbond_angl_SCG
     real(PREC) :: hbond_dihd_SAUS
     real(PREC) :: hbond_dihd_SGCS
     real(PREC) :: hbond_dihd_PSAU
     real(PREC) :: hbond_dihd_PSUA
     real(PREC) :: hbond_dihd_PSGC
     real(PREC) :: hbond_dihd_PSCG
     integer :: sz
  endtype input_aform_rna
  type(input_aform_RNA), save :: inarna

  !==========================================
  !> structure for parameters reading from "ligand.para" file
  type input_ligandparameter
     ! ligand.para
     real(PREC) :: energy_unit
     real(PREC) :: cbd
     real(PREC) :: cba
     real(PREC) :: cdih
     real(PREC) :: cutoff_exvol
     real(PREC) :: cdist_rep12_lpro
     real(PREC) :: cdist_rep12_llig
     real(PREC) :: crep12

     integer :: sz
  end type input_ligandparameter

  type(input_ligandparameter), save :: inligand

  !==========================================
  !> structure for parameters reading from "hydrophobic.para" file
  type input_hpparameter
     ! hydrophobic.para
     real(PREC) :: rho_min_hp
     real(PREC) :: coef_rho_hp
     real(PREC) :: coef_hp
     real(PREC) :: coefaa_para_hp(21)
     integer    :: ncoor_para_hp(21)
     real(PREC) :: ncoormax_para_hp(21)
     real(PREC) :: cutoffdmin_para_hp(21,21)
     real(PREC) :: cutoffdmax_para_hp(21,21)
     logical    :: flag_hp(MXMP)    ! filter

     integer :: sz
  end type input_hpparameter

  type(input_hpparameter), save :: inhp

  !==========================================
  !> structure for parameters reading from "flexible_local.para" file
  type input_flexible_local_parameter
     
     integer :: i_flp ! flag
!     integer :: nflp  ! # of the region
     ! start and end id of the region calculated by flexible local potential
!     integer :: iflp_lgo(2, MXFLEXIBLE)
          
     real(PREC) :: ang_para_x(10)
     real(PREC) :: ang_para_y(20,10)
     real(PREC) :: ang_para_y2(20,10)
     real(PREC) :: dih_para(400,400,7)
     !real(PREC) :: coeff ! boltzmann constant * temperature
     real(PREC) :: k_ang
     real(PREC) :: k_dih

     ! variable related to boundary
     integer :: boundary
     real(PREC) :: boundary_lower_limit
     real(PREC) :: boundary_upper_limit
     real(PREC) :: boundary_cutoff
     real(PREC) :: boundary_max_force

     integer :: sz
     
  end type input_flexible_local_parameter

  type(input_flexible_local_parameter), save :: inflp

  !==========================================
  !> structure for parameters reading from "md_information" field in input file
  type input_simuparameter
     ! step
     integer    :: n_step_sim
     integer    :: i_step_sim_init
     integer(L_INT) :: n_tstep(MXSIM)
     integer(L_INT) :: n_tstep_all
     integer(L_INT) :: i_tstep_init
     integer    :: n_step_save
     integer    :: n_step_rst
     integer    :: n_step_neighbor
     real(PREC) :: tstep_size

     ! offset
     integer    :: i_com_zeroing_ini
     integer    :: i_com_zeroing
     integer    :: i_no_trans_rot

     ! other
     real(PREC) :: tempk
     real(PREC) :: tempk_ref
     integer    :: i_rand_type !< index of type of random number generator
     integer    :: n_seed

     integer    :: sz
  end type input_simuparameter

  type(input_simuparameter), save :: insimu

  !==========================================
  !> structure for parameters reading from "annealing" field in input file
  type input_annealing
     real(PREC) :: tempk_init
     real(PREC) :: tempk_last
     integer    :: n_time_change
     integer    :: sz
  end type input_annealing

  type(input_annealing), save :: inann

  !==========================================
  !> structure for parameters reading from "searching_tf" field in input file
  type input_searchingtf
     real(PREC) :: tempk_upper
     real(PREC) :: tempk_lower
     integer    :: sz
  end type input_searchingtf

  type(input_searchingtf), save :: insear

  !==========================================
  !> structure for parameters reading from various fields in input file
  type input_miscellaneous
     ! flag
     integer    :: i_use_atom_protein
     integer    :: i_residuenergy_radii
     integer    :: i_output_energy_style
     integer    :: i_dtrna_model
     logical    :: flag_local_unit(MXUNIT, MXUNIT, LINTERACT%MAX)
     logical    :: flag_nlocal_unit(MXUNIT, MXUNIT, INTERACT%MAX)
     logical    :: flag_prime(INTERACT%MAX)
     logical    :: force_flag(INTERACT%MAX)
     logical    :: force_flag_local(LINTERACT%MAX) !AICG
     logical    :: class_flag(CLASS%MAX)
     integer    :: i_triple_angle_term
     integer    :: i_reset_struct
     integer    :: i_temp_independent
     logical    :: flg_coef_from_ninfo

     ! redefine_parameter
     integer    :: i_redef_para

     ! box interaction
     integer    :: i_in_box
     real(PREC) :: xbox
     real(PREC) :: ybox
     real(PREC) :: zbox
     real(PREC) :: boxsigma

     ! cap interaction
     integer    :: i_in_cap
     real(PREC) :: rcap
     real(PREC) :: kcap
     real(PREC) :: center_cap(3)

     ! delete interaction
     integer    :: i_del_int
     !integer    :: i_add_int
     integer    :: ndel_lgo
     integer    :: ndel_go
     integer    :: idel_lgo(2, MXDEL_LGO)
     integer    :: idel_go(4, MXDEL_GO)

     ! energy coefficient
     integer    :: i_energy_para
     real(PREC) :: factor_go_unit(MXUNIT, MXUNIT)
     real(PREC) :: factor_local_unit(MXUNIT, MXUNIT)

     ! mass and friction
     integer    :: i_fric
     integer    :: i_mass
     integer    :: i_redef_mass_fric

     ! neighbordist
     integer    :: i_neigh_dist
     integer    :: i_neigh_dynamic
     real(PREC) :: rneighbordist2_unit(MXUNIT, MXUNIT)

     ! bridge
     integer    :: i_bridge
     integer    :: nbrid
     integer    :: nbrid_com
     integer    :: nbrid_ppr
     integer    :: i_lower_bound
     integer    :: ibrid2mp(2, MXBRIDGE)
     integer    :: ibrid_com2grp(2, MXBRIDGE)
     real(PREC) :: coef_brid(MXBRIDGE)
     real(PREC) :: coef_brid_com(MXBRIDGE)
     real(PREC) :: brid_dist(MXBRIDGE)
     real(PREC) :: brid_com_dist(MXBRIDGE)
     integer    :: ibrid_ppr_com(MXBRIDGE)   !< PPR
     integer    :: ibrid_ppr_gid_r(2,MXBRIDGE) !< PPR
     real(PREC) :: brid_ppr_rmin(MXBRIDGE)   !< PPR
     real(PREC) :: brid_ppr_rcut(MXBRIDGE)   !< PPR
     real(PREC) :: brid_ppr_rmax(MXBRIDGE)   !< PPR
     real(PREC) :: brid_ppr_rzero2(MXBRIDGE) !< PPR
     integer    :: ibrid_ppr_cycl(MXBRIDGE)  !< PPR
     integer    :: ibrid_ppr_opt(MXBRIDGE)   !< PPR

     ! pulling
     logical    :: flg_pull_energy
     integer    :: i_pulling
     integer    :: npull
     integer    :: ipull2mp(MXPULLING)
     real(PREC) :: coef_pull(MXPULLING)
     real(PREC) :: pull_xyz(3, MXPULLING)
     real(PREC) :: pu_xyz(3, MXPULLING)
     integer    :: ipull_force_unit
     integer    :: npull_unravel
     integer    :: ipull_unravel2mp(2, MXPULLING)
     real(PREC) :: pull_unravel_xyz(3, MXPULLING, MXREPLICA)

     ! anchor
     integer    :: i_anchor
     integer    :: nanc
     integer    :: nanc_com_ini
     integer    :: ianc2mp(MXANCHOR)
     integer    :: ianc_com_ini2grp(MXANCHOR)
     real(PREC) :: coef_anc(MXANCHOR)
     real(PREC) :: coef_anc_com_ini(MXANCHOR)
     real(PREC) :: anc_dist(MXANCHOR)
     real(PREC) :: anc_com_ini_dist(MXANCHOR)
     real(PREC) :: anc_xyz(3, MXANCHOR)
     real(PREC) :: anc_com_ini_xyz(3, MXANCHOR, MXREPLICA)

     ! rest1d
     integer    :: i_rest1d
     integer    :: nrest1d
     integer    :: irest1d2mp(MXREST1D)
     real(PREC) :: coef_rest1d(MXREST1D)
     integer    :: irest1d_mp_sa(MXREST1D)
     integer    :: irest1d_mp_sb(MXREST1D)
     real(PREC) :: rest1d_s0(MXREST1D)
     real(PREC) :: rest1d_s(MXREST1D)
     
     integer    :: i_rest1d_center
     integer    :: nrest1d_center
     integer    :: nrest1d_center_mp(MXREST1D)
     integer    :: irest1d_center2mp(MXREST1D, MXMP)
     real(PREC) :: coef_rest1d_center(MXREST1D)
     integer    :: irest1d_center_mp_sa(MXREST1D)
     integer    :: irest1d_center_mp_sb(MXREST1D)
     real(PREC) :: rest1d_center_s0(MXREST1D)
     real(PREC) :: rest1d_center_s(MXREST1D)
     integer    :: rest1d_center_init_flag(MXREST1D)
     real(PREC) :: rest1d_center_v(MXREST1D, 3)
     real(PREC) :: rest1d_center_origin(MXREST1D, 3)
     
     ! fix
     integer    :: i_fix

     ! implicit_ligand
     integer    :: i_implig

     ! Window exchange simulation
     integer    :: i_window
     integer    :: i_winz

     ! Tensor type of hydrodynamic interaction in Brownian dynamics
     integer    :: i_hydro_tensor

     ! cylinder
     integer    :: i_cylinder
     real(PREC) :: cylinder_bgn
     real(PREC) :: cylinder_end
     real(PREC) :: cylinder_radi
     real(PREC) :: cylinder_coef
     
     integer    :: sz
  end type input_miscellaneous

  type(input_miscellaneous), save :: inmisc

  !==========================================
  !mcanonical
  type input_modified_muca
     integer    :: i_modified_muca
     real(PREC) :: em_depth
     real(PREC) :: em_sigma
     real(PREC) :: em_mid
     integer    :: sz
  end type input_modified_muca

  type(input_modified_muca), save :: inmmc

  !==========================================
  ! aicg
  type input_aicg_parameter
     real(PREC) :: cbd_aicg
     real(PREC) :: cba_aicg_G
     real(PREC) :: cba_aicg_H
     real(PREC) :: cba_aicg_E
     real(PREC) :: cba_aicg_T
     real(PREC) :: cba_aicg_C
     real(PREC) :: cdih_aicg_G
     real(PREC) :: cdih_aicg_H
     real(PREC) :: cdih_aicg_E
     real(PREC) :: cdih_aicg_T
     real(PREC) :: cdih_aicg_C
     real(PREC) :: ave_caicg
     real(PREC) :: gen_caicg
     real(PREC) :: ecut_low
     real(PREC) :: ecut_up
     integer    :: iflag_scale
     integer    :: sz
  end type input_aicg_parameter

  type(input_aicg_parameter), save :: inaicg

  !==========================================
  ! aicg2
  type input_aicg2_parameter
     real(PREC) :: cbd_aicg2
     real(PREC) :: wid_aicg13
     real(PREC) :: wid_aicg14
     real(PREC) :: wid_dih
     real(PREC) :: ave_caicg2_loc
     real(PREC) :: gen_caicg2_loc
     real(PREC) :: ave_caicg2plus_13
     real(PREC) :: gen_caicg2plus_13
     real(PREC) :: ave_caicg2plus_14
     real(PREC) :: gen_caicg2plus_14
     real(PREC) :: ave_caicg2_nloc
     real(PREC) :: gen_caicg2_nloc
     real(PREC) :: ave_caicg2plus_nloc
     real(PREC) :: gen_caicg2plus_nloc
     real(PREC) :: ecut_low_aicg2
     real(PREC) :: ecut_up_aicg2
     integer    :: iflag_scale_aicg2
     integer    :: sz
  end type input_aicg2_parameter

  type(input_aicg2_parameter), save :: inaicg2

  !==========================================
  ! electrostatic
  type input_electrostatic
     integer    :: i_diele
     integer    :: i_charge
     integer    :: i_function_form
     integer    :: i_calc_method
     real(PREC) :: diele
     real(PREC) :: diele_water
     real(PREC) :: diele_dTcoef
     real(PREC) :: cutoff_ele
     real(PREC) :: ionic_strength
     real(PREC) :: cdist(MXREPLICA)
     real(PREC) :: coef(MXREPLICA)
     real(PREC) :: coef_charge_type(SEQID%MAX)
     real(PREC) :: length_per_unit(SEQID%MAX)
     real(PREC) :: ewld_alpha
     real(PREC) :: ewld_hmax
     logical :: flag_ele(MXMP)
     logical :: flag_charge_change
     integer :: n_charge_change
     integer :: charge_change_imp(MXCHARGECHANGE)
     real(PREC) :: charge_change_value(MXCHARGECHANGE)
     integer    :: sz
  endtype input_electrostatic

  type(input_electrostatic), save :: inele
  
  !==========================================
  ! sasa
  type input_sasa
     real(PREC) :: p_sasa(SEQID%MAX+1)
     real(PREC) :: r_sasa(SEQID%MAX+1)
     real(PREC) :: r_sol
     real(PREC) :: connectivity(4)
     real(PREC) :: coef_surf
     integer    :: sz
  endtype input_sasa

  type(input_sasa), save :: insasa

  !==========================================
  ! excluded volume
  type input_exv
     real(PREC) :: exv_sigma(0:CHEMICALTYPE%MAX)
     real(PREC) :: exv_cutoff !< truncation distance for exv = exv_cutoff * exv_sigma
     real(PREC) :: exv_coef   !< coefficient epsilon for exv
     integer    :: sz
  endtype input_exv

  type(input_exv), save :: inexv

  !==========================================
  ! Window exchange
  type input_window
     integer :: iwind(MXREPLICA)
     integer :: iwinz(MXREPLICA)
  endtype input_window

  type(input_window), save :: inwind

  !==========================================
  ! Widom method (calculate chemical potential)
  type input_widom
     integer :: n_step_skip
     integer :: n_step_interval
     integer :: n_trial
     integer :: n_Mg_add
     integer :: n_Na_add
     integer :: n_K_add
     integer :: n_Cl_add
     integer :: n_max_mp_add
     integer :: sz
  endtype input_widom

  type(input_widom), save :: inwidom

  !==========================================
  
  integer, save :: irand !< variable for random number
!  integer, parameter :: NSTREAM=128
!  type(mt_state), save :: mts(0:MXREPLICA)
  type(mt_state), allocatable, save :: mts(:,:)

  integer, save :: ifix_mp(MXMP) !< 

  ! periodic boundary  
  type periodic_boundary
     integer :: i_periodic
     integer :: n_mirror_index
     real(PREC) :: psize(3)
     real(PREC) :: psizeh(3)
     real(PREC) :: d_mirror(3, 27)
     integer :: sz
  endtype periodic_boundary
  type(periodic_boundary), save :: inperi

  ! elastic network
  type input_enmparameter
     real(PREC) :: cenm !< constant coefficient of the energy function for elastic network model
     real(PREC) :: dfcontact_enm !< cutoff distance to define the native contact for elastic network model
     integer    :: i_enm !< flag variable for elastic network model
     integer    :: sz !< size of the structure
  end type input_enmparameter

  type(input_enmparameter), save :: inenm
end module var_setp
