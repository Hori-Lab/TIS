! var_struct
!> @brief Module defining the global variables that are
!>        mainly related to structure information
!> "mp" is mass point(MP).
!> "unit" is basic unit of beads.
!> For protein and nucleotide, "unit" corresponds to a single chain.
!> "class" is the type of the unit; protein and RNA

module var_struct

  use const_maxsize
  use const_physical
  implicit none

  ! ----------------------------------------------------------------
  integer, save :: nunit_real !< # of units for real chain (not includ units for shadow units)
  integer, save :: nunit_all !< # of all units
  integer, save :: nmp_real !< # of MP for real chain
  integer, save :: nmp_all !< # of all MP
  integer, save :: nres !< # of residues
  integer, save :: lunit2mp(2, MXUNIT) !< initial mp and final mp of unit
  integer, save :: ires_mp(MXMP) !< residue number of mp
  integer, save :: iontype_mp(MXMP) !< type of ion
  integer, save :: iclass_unit(MXUNIT) !< class of unit
  integer, save :: iclass_mp(MXMP) !< class of mp
  real(PREC), allocatable, save :: xyz_mp_rep(:,:,:) !< coordinate of mp
  real(PREC), allocatable, save :: pxyz_mp_rep(:,:,:) !< coordinate of mp
  real(PREC), save :: xyz_ref_mp(SPACE_DIM, MXMP) !< coordinate of reference structure
  real(PREC), save :: cmass_mp(MXMP) !< mass of mp
  real(PREC), save :: exv_radius_mp(MXMP) !< excluded-volume radius (half of sigma) of mp
  real(PREC), save :: exv_epsilon_mp(MXMP)
  real(PREC), save :: fric_mp(MXMP) !< friction constant of mp
  real(PREC), save :: radius(MXMP)  !< radius of mp (to be used in viscosity calculation)
  character(3), save :: cmp2seq(MXMP) !< residue name of mp
  character(4), save :: cmp2atom(MXMP) !< atom name of mp

  integer, save :: imp2unit(MXMP)  !< unit of mp
  integer, save :: istype_mp(MXMP) !< secondary structure of mp
  integer, save :: imp2type(MXMP)

  type group
     integer    :: ngrp   !< number of group
     integer    :: nmp(MXGRP)  !< number of mp belonging to each group
     integer    :: implist(MXMPGRP, MXGRP)  !< list of mp
     real(PREC) :: mass_fract(MXMPGRP,MXGRP) !< Mass fraction of mp in each group
     integer    :: sz
  endtype group
  type(group), save :: grp
     

  ! ----------------------------------------------------------------
  ! protein structure
  !> parameters for bond length potential
  integer,    save :: nbd = 0
  integer,    allocatable, save :: ibd2mp(:,:)       !(2, MXBD)
  integer,    allocatable, save :: ibd2type(:)       !(MXBD)
  real(PREC), allocatable, save :: bd_nat(:)         !(MXBD)
  real(PREC), allocatable, save :: factor_bd(:)      !(MXBD)
  real(PREC), allocatable, save :: coef_bd(:,:)      !(2, MXBD)
  real(PREC), allocatable, save :: correct_bd_mgo(:) !(MXBD)

  !> parameters for bond angle potential
  integer,    save :: nba = 0
  integer,    save :: nfba = 0
  integer,    allocatable, save :: iba2mp(:,:)       !(3, MXBA)
  integer,    allocatable, save :: ifba2mp(:,:)      !(3, MXBA)
  integer,    allocatable, save :: iba2type(:)       !(MXBA)
  integer,    allocatable, save :: iunit2ba(:,:)     !(2, MXUNIT)
  real(PREC), allocatable, save :: ba_nat(:)         !(MXBA)
  real(PREC), allocatable, save :: factor_ba(:)      !(MXBA)
  real(PREC), allocatable, save :: coef_ba(:,:)      !(2, MXBA)
  real(PREC), allocatable, save :: correct_ba_mgo(:) !(MXBA)
  real(PREC), allocatable, save :: fba_para_x(:,:)   !(10, nfba)
  real(PREC), allocatable, save :: fba_para_y(:,:)   !(10, nfba)
  real(PREC), allocatable, save :: fba_para_y2(:,:)  !(10, nfba)
  real(PREC), allocatable, save :: fba_ener_corr(:)   !(nfba)
  real(PREC), allocatable, save :: fba_max_th(:)     !(nfba)
  real(PREC), allocatable, save :: fba_max_th_ener(:)!(nfba)
  real(PREC), allocatable, save :: fba_min_th(:)     !(nfba)
  real(PREC), allocatable, save :: fba_min_th_ener(:)!(nfba)
  integer,    allocatable, save :: ifba2ba(:)      !(MXBA) !flp_mgo

  !> parameters for dihedral angle potential
  integer,    save :: ndih = 0
  integer,    save :: nfdih = 0
  integer,    allocatable, save :: idih2mp(:,:)        !(4, MXDIH)
  integer,    allocatable, save :: ifdih2mp(:,:)       !(4, MXDIH)
  integer,    allocatable, save :: idih2type(:)        !(MXDIH)
  integer,    allocatable, save :: iunit2dih(:,:)      !(2, MXUNIT)
  real(PREC), allocatable, save :: dih_nat(:)          !(MXDIH)
  real(PREC), allocatable, save :: factor_dih(:)       !(MXDIH)
  real(PREC), allocatable, save :: coef_dih(:,:)       !(2, MXDIH)
  real(PREC), allocatable, save :: dih_sin_nat(:)      !(MXDIH)
  real(PREC), allocatable, save :: dih_cos_nat(:)      !(MXDIH)
  real(PREC), allocatable, save :: correct_dih_mgo(:)  !(MXDIH)
  real(PREC), allocatable, save :: fdih_para(:,:)      !(7, nfdih)
  real(PREC), allocatable, save :: fdih_ener_corr(:)   !(nfdih)
  real(PREC), allocatable, save :: fdih_mult_para(:,:) !(7, nres)
  integer,    allocatable, save :: fdih_mult_resID(:)  !(nres)
  integer,    allocatable, save :: ifdih2dih(:)       !(MXDIH) !flp_mgo      

  !> parameters for go (1210) potential
  integer,    save :: ncon = 0
  integer,    allocatable, save :: icon2mp(:,:)      !(2, MXCON)
  integer,    allocatable, save :: icon2type(:)      !(MXCON)
  integer,    allocatable, save :: lmp2con(:)        !(MXMP)
  integer,    allocatable, save :: icon2unit(:,:)    !(2, MXCON)
  integer,    allocatable, save :: icon_dummy_mgo(:) !(MXCON)
  real(PREC), allocatable, save :: go_nat(:)         !(MXCON)
  real(PREC), allocatable, save :: go_nat2(:)        !(MXCON)
  real(PREC), allocatable, save :: factor_go(:)      !(MXCON)
  real(PREC), allocatable, save :: coef_go(:)        !(MXCON)

  !> parameter for go (morse) potential
  integer,    save :: nmorse = 0
  integer,    allocatable, save :: imorse2mp(:,:)      !(2, MXMORSE)
  integer,    allocatable, save :: imorse2type(:)      !(MXMORSE)
  integer,    allocatable, save :: lmp2morse(:)        !(MXMP)
  integer,    allocatable, save :: imorse2unit(:,:)    !(2, MXMORSE)
  integer,    allocatable, save :: imorse_dummy_mgo(:) !(MXMORSE)
  real(PREC), allocatable, save :: morse_nat(:)        !(MXMORSE)
  real(PREC), allocatable, save :: morse_nat2(:)       !(MXMORSE)
  real(PREC), allocatable, save :: factor_morse(:)     !(MXMORSE)
  real(PREC), allocatable, save :: coef_morse_fD(:)    !(MXMORSE)
  real(PREC), allocatable, save :: coef_morse_a(:)     !(MXMORSE)

  !> parameters for aicg13 potential  !aicg2
  real(PREC), allocatable, save :: aicg13_nat(:)        !(MXBA)
  real(PREC), allocatable, save :: coef_aicg13_gauss(:) !(MXBA)
  real(PREC), allocatable, save :: wid_aicg13_gauss(:)  !(MXBA)
  real(PREC), allocatable, save :: factor_aicg13(:)     !(MXBA)

  !> parameters for aicg14 potential  !aicg2
  real(PREC), allocatable, save :: aicg14_nat(:)          !(MXDIH)
  real(PREC), allocatable, save :: coef_aicg14_gauss(:)   !(MXDIH)
  real(PREC), allocatable, save :: wid_aicg14_gauss(:)    !(MXDIH)
  real(PREC), allocatable, save :: coef_dih_gauss(:)   !(MXDIH)
  real(PREC), allocatable, save :: wid_dih_gauss(:)    !(MXDIH)
  real(PREC), allocatable, save :: factor_aicg14(:)       !(MXDIH)

  !> parameters sasa
  !sasa
  real(PREC), allocatable, save :: para_sasa(:)          !(MXMP)
  real(PREC), allocatable, save :: rad_sasa(:)   !(MXMP)
  real(PREC), allocatable, save :: surf(:)   !(MXMP)
  real(PREC), allocatable, save :: connect(:)   !(-MXMP:MXMP) 

  !> parameters for rna_bp potential
  integer,    save :: nrna_bp = 0
  integer,    allocatable, save :: irna_bp2mp(:,:)      !(2, MXRNABP)
  integer,    allocatable, save :: lmp2rna_bp(:)        !(MXMP)
  integer,    allocatable, save :: irna_bp2unit(:,:)    !(2, MXRNABP)
  integer,    allocatable, save :: nhb_bp(:)            !(MXRNABP)
  integer,    allocatable, save :: irna_bp_dummy_mgo(:) !(MXRNABP)
  real(PREC), allocatable, save :: rna_bp_nat(:)        !(MXRNABP)
  real(PREC), allocatable, save :: rna_bp_nat2(:)       !(MXRNABP)
  real(PREC), allocatable, save :: coef_rna_bp(:)       !(MXRNABP) for LJ1210
  real(PREC), allocatable, save :: coef_rna_bp_a(:)     !(MXRNABP) for Morse
  real(PREC), allocatable, save :: coef_rna_bp_fD(:)    !(MXRNABP) for Morse
  real(PREC), allocatable, save :: factor_rna_bp(:)     !(MXRNABP)

  !> parameters for rna_st potential
  integer,    save :: nrna_st = 0
  integer,    allocatable, save :: irna_st2mp(:,:)       !(2, MXRNAST)
  integer,    allocatable, save :: lmp2rna_st(:)         !(MXMP)
  integer,    allocatable, save :: irna_st2unit(:,:)     !(2,MXRNAST)
  integer,    allocatable, save :: irna_st_dummy_mgo(:)  !(MXRNAST)
  real(PREC), allocatable, save :: rna_st_nat(:)         !(MXRNAST)
  real(PREC), allocatable, save :: rna_st_nat2(:)        !(MXRNAST)
  real(PREC), allocatable, save :: coef_rna_st(:)        !(MXRNAST) for LJ1210
  real(PREC), allocatable, save :: coef_rna_st_a(:)      !(MXRNAST) for Morse
  real(PREC), allocatable, save :: coef_rna_st_fD(:)     !(MXRNAST) for Morse
  real(PREC), allocatable, save :: factor_rna_st(:)      !(MXRNAST)

  !> parameters for rna stack angle potential
  integer,    save :: nstangle = 0
  integer,    allocatable, save :: istangle2mp(:,:)  !(3,MXMP)
  real(PREC), allocatable, save :: stangle_nat(:)    !(MXRNASTANGLE)
  real(PREC), allocatable, save :: factor_stangle(:) !(MXRNASTANGLE)
  real(PREC), allocatable, save :: coef_stangle(:,:) !(2, MXRNASTANGLE)

  !> parameters for base-stack potential of DT-RNA model
  integer,    save :: ndtrna_st = 0
  integer,    allocatable, save :: idtrna_st2mp(:,:)     !(7,MXDTRNAST)
  integer,    allocatable, save :: idtrna_st2nn(:)       !(MXDTRNAST)
  real(PREC), allocatable, save :: dtrna_st_nat(:,:)     !(3, MXDTRNAST)
  real(PREC), allocatable, save :: coef_dtrna_st(:,:,:)  !(0:3, MXDTRNAST,n_replica_all)

  !> parameters for hydrogen-bond potential of DT-RNA model
  integer,    save :: ndtrna_hb = 0
  integer,    allocatable, save :: idtrna_hb2mp(:,:)     !(6,MXDTRNAHB)
  real(PREC), allocatable, save :: dtrna_hb_nat(:,:)     !(6, MXDTRNAHB)
  real(PREC), allocatable, save :: coef_dtrna_hb(:,:)    !(0:6, MXDTRNAHB)

  !> parameters for hydrogen-bond potential of DT-RNA model (2015)
  integer,    save :: nhbsite
  integer, allocatable, save :: imp2hbsite(:,:)     ! (1:2, 1:nmp)
  integer, allocatable, save :: nvalence_hbsite(:)  ! (1:nhbsite)
  integer, allocatable, save :: idtrna_hb2hbsite(:,:,:)   ! (1:3, 1:2, 1:ndtrna_hb)
               ! There are at most three hydrogen bonds for each HB interaction (0 if empty)
               ! (1:3,1,ihb) = IDs of "ihbsite" in one side of the interaction "ihb"
               ! (1:3,2,ihb) = IDs of "ihbsite" in the other side of the interaction "ihb"
  integer, allocatable, save :: list_hb_at_hbsite(:,:,:) ! (MXMPHBNEIGHBOR, 1:nhbsite, REPLICA)
  integer, allocatable, save :: num_hb_at_hbsite(:,:)    ! (1:nhbsite, REPLICA)
  integer, allocatable, save :: nhbneigh(:)         ! (REPLICA)
  integer, allocatable, save :: ineigh2hb(:,:)      ! (1:MXMPHBNEIGH*nmp_all/2, REPLICA)
  logical, allocatable, save :: flg_hb_tertiary(:)  ! (1:MXDTRNHB)

  !> parameters for tertiary-base-stack potential of DT-RNA model (2015)
  integer,    save :: ndtrna_tst = 0
  integer,    allocatable, save :: idtrna_tst2mp(:,:)     !(6,MXDTRNATST)
  integer,    allocatable, save :: idtrna_tst2side(:,:)   !(2,MXDTRNATST)
  real(PREC), allocatable, save :: dtrna_tst_nat(:,:)     !(6, MXDTRNATST)
  real(PREC), allocatable, save :: coef_dtrna_tst(:,:)    !(0:6, MXDTRNATST)

  logical,    allocatable, save :: flg_tst_exclusive(:,:) !(2,MXDTRNATST)  
  integer,    allocatable, save :: idtrna_tst2st(:,:)     !(2,MXDTRNATST)

  !> for simu_energy_orderpara
  integer, allocatable, save :: iallcon2unit(:,:)

  !> parameter for qscore
  integer,    save :: ncon_unit(MXUNIT, MXUNIT)
  integer,    save :: nrna_bp_unit(MXUNIT, MXUNIT)
  integer,    save :: nrna_st_unit(MXUNIT, MXUNIT)

  !> parameters for electrostatic
  integer,    save :: ncharge
  integer, allocatable, save :: icharge2mp(:)         ! (MXCHARGE = nmp_all)
  integer, allocatable, save :: lmp2charge(:)         ! (nmp_all) 
  real(PREC), allocatable, save :: coef_charge(:,:)   ! charge values (MXCHARGE=nmp_all, n_replica_all)
  integer, allocatable, save :: lele(:)          ! (REPLICA)
  integer, allocatable, save :: iele2mp(:,:,:)   ! (2(+1), MXMPELE*nmp_all, REPLICA)
  real(PREC), allocatable, save :: coef_ele(:,:) ! (MXMPELE*nmp_all, REPLICA)
  real(PREC), allocatable, save :: xyz_ele_rep(:,:,:) ! (SPACE_DIM, ncharge, REPLICA)
  real(PREC), allocatable, save :: pxyz_ele_rep(:,:,:) ! (SPACE_DIM, ncharge, REPLICA)

  !> parameters for elctrostatic(K computer)
  integer, allocatable, save :: lele_k(:,:)        ! (ncharge, REPLICA)
  integer, allocatable, save :: iele2charge_k(:,:,:) ! (ncharge, ncharge, REPLICA)

  ! ----------------------------------------------------------------
  !> parameters for neighboring (general) list
  integer, allocatable, save :: lexv(:,:,:)   ! (2, E_TYPE%MAX, REPLICA)  !replica
  integer, allocatable, save :: iexv2mp(:,:,:)! (2(+1), MXMPNEIGHBOR*nmp_all, REPLICA)

  ! ----------------------------------------------------------------
  !> parameters for hydrophobic interaction
  integer, save :: nhp                     ! number of HP mass_points
  integer, save :: ihp2mp(MXHP)
  integer, save :: lunit2hp(2, MXUNIT)
  real(PREC), save :: ncoor_hp(MXHP)     ! number of united atoms
  real(PREC), save :: ncoor_max_hp(MXHP) ! the max coordination number 
  real(PREC), save :: coef_aa_hp(MXHP)   ! coefficients
  integer, allocatable, save :: nhpneigh(:) ! neighbor list for hp (REPLICA)
  integer, allocatable, save :: lhp2neigh(:, :, :)   !(2, MXHP, REPLICA)
  integer, allocatable, save :: ineigh2hp(:, :)      !(MXMPHP*nhp, REPLICA)
  real(PREC), allocatable, save :: cutoff_dmin_hp(:, :) !(MXMPHP*nhp, REPLICA) 
  real(PREC), allocatable, save :: cutoff_dmax_hp(:, :) !(MXMPHP*nhp, REPLICA)

contains

  integer function get_icon_type (imp, jmp)
     use const_index
     use const_maxsize
     implicit none
     integer, intent(in) :: imp, jmp
     character(CARRAY_MSG_ERROR) :: error_message

     get_icon_type = -1

     ! protein : protein
     if (imp2type(imp) == MPTYPE%PRO .AND. imp2type(jmp) == MPTYPE%PRO) then
        get_icon_type = CONTYPE%PRO_PRO

     ! PROTEIN : LIGAND
     else if ((iclass_mp(imp) == CLASS%LIG .AND. imp2type(jmp) == MPTYPE%PRO) .OR. &
              (iclass_mp(jmp) == CLASS%LIG .AND. imp2type(imp) == MPTYPE%PRO)) then
        get_icon_type = CONTYPE%PRO_LIG

     ! protein : RNA
     else if (imp2type(imp) == MPTYPE%PRO .OR. imp2type(jmp) == MPTYPE%PRO) then
        if (imp2type(imp) == MPTYPE%RNA_PHOS .OR. imp2type(jmp) == MPTYPE%RNA_PHOS) then
           get_icon_type = CONTYPE%PRO_RP
        else if (imp2type(imp) == MPTYPE%RNA_BASE .OR. imp2type(jmp) == MPTYPE%RNA_BASE) then
           get_icon_type = CONTYPE%PRO_RB
        else if (imp2type(imp) == MPTYPE%RNA_SUGAR .OR. imp2type(jmp) == MPTYPE%RNA_SUGAR) then
           get_icon_type = CONTYPE%PRO_RS
        else
           write(error_message,*) 'Error: logical defect in get_icon_type'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

     ! RNA(P) : RNA(P)
     else if (imp2type(imp) == MPTYPE%RNA_PHOS .AND. imp2type(jmp) == MPTYPE%RNA_PHOS) then
        get_icon_type = CONTYPE%RP_RP

     ! RNA(P) : RNA(S,B)
     else if (imp2type(imp) == MPTYPE%RNA_PHOS .OR. imp2type(jmp) == MPTYPE%RNA_PHOS) then
        if (imp2type(imp) == MPTYPE%RNA_BASE .OR. imp2type(jmp) == MPTYPE%RNA_BASE) then
           get_icon_type = CONTYPE%RP_RB
        else if (imp2type(imp) == MPTYPE%RNA_SUGAR .OR. imp2type(jmp) == MPTYPE%RNA_SUGAR) then
           get_icon_type = CONTYPE%RP_RS
        else
           write(error_message,*) 'Error: logical defect in get_icon_type'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

     ! RNA(S) : RNA(S)
     else if (imp2type(imp) == MPTYPE%RNA_SUGAR .AND. imp2type(jmp) == MPTYPE%RNA_SUGAR) then
        get_icon_type = CONTYPE%RS_RS

     ! RNA(S) : RNA(B)
     else if (imp2type(imp) == MPTYPE%RNA_SUGAR .OR. imp2type(jmp) == MPTYPE%RNA_SUGAR) then
        if (imp2type(imp) == MPTYPE%RNA_BASE .OR. imp2type(jmp) == MPTYPE%RNA_BASE) then
           get_icon_type = CONTYPE%RS_RB
        else
           write(error_message,*) 'Error: logical defect in get_icon_type'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

     ! RNA(B) : RNA(B)
     else if (imp2type(imp) == MPTYPE%RNA_BASE .AND. imp2type(jmp) == MPTYPE%RNA_BASE) then
        get_icon_type = CONTYPE%RB_RB

     else 
        write(error_message,*) 'Error: logical defect in get_icon_type'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  endfunction get_icon_type

  character(1) function rna_base_type(imp)
     use const_index
     implicit none
     integer, intent(in) :: imp
     character(CARRAY_MSG_ERROR) :: error_message
     if (cmp2atom(imp) == ' Ab ' .OR. cmp2atom(imp) == ' Gb ' .OR. &
         cmp2atom(imp) == ' Rb ') then
        rna_base_type = 'R'
     else if (cmp2atom(imp) == ' Ub ' .OR. cmp2atom(imp) == ' Cb ' .OR. &
              cmp2atom(imp) == ' Yb ') then
        rna_base_type = 'Y'
     else
       error_message = 'Error: invalid atom name, in rna_base_type in setp_native_dih'
       call util_error(ERROR%STOP_ALL, error_message)
     endif
  endfunction rna_base_type

end module var_struct
