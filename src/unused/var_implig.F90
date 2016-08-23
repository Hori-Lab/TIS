!var_implig
!> @brief Contains parameters and variables which are related &
!>        to implicit ligand.

module var_implig

  use const_maxsize
!! MXSITE_IMPLIG: max number of (ligand) binding site
!! MXAA_IMPLIG:   max number of residue for each (ligand) binding-site
!! MXUNIT_IMPLIG: max number of interactive unit for each implicit-ligand
!! MXCON_IMPLIG:  max number of ligand-mediated contact for (ligand) binding site

  implicit none

  !> inimplig parameters are given by input file (<<<< implicit_ligand part).
  type input_impligparameter
     integer :: iexe_implig !< flag for implicit ligand procedure (1==execute, 0==No-execute)
     integer :: nsite_implig !< the number of implicit ligand 
     integer :: initial_state_implig(MXSITE_IMPLIG) !< initial state of implicit ligand 
     integer :: istep_implig !< time steps that are related to the rate constants for (implicit) ligand binding and unbinding process.
     integer :: istep_un_implig !< time steps that are related to the rate constants for (implicit) ligand binding and unbinding process.
     integer :: itype_ene_implig !< interaction type for implicit ligand (1==Gaussian type, 0==LJ12-10 type)
     integer :: isep_contact_implig !< (default: 3)
     real(PREC) :: diff_const_implig(MXSITE_IMPLIG) !< the diffusion constant of implicit ligand
     real(PREC) :: react_rad_implig(MXSITE_IMPLIG) !< reaction radius for implicit-ligand (A)
     real(PREC) :: c_implig(MXSITE_IMPLIG) !< concentration of implicit-ligand    
     real(PREC) :: bind_rate_implig(MXSITE_IMPLIG) !< the ligand-binding rate constant k_b 
     real(PREC) :: pre_implig(MXSITE_IMPLIG) !< decide the strength of ligand mediated interaction
     real(PREC) :: gauss_d_implig(MXSITE_IMPLIG) !< decide the typical length (scale) for interaction
                                                 !< between prorein-residues and impicit-ligand (ligand-binding residues). 
                                                 !< [especially for gauss type interaction]   
     real(PREC) :: dfcontact_implig !< the length of difinition for ligand-mediated contact (default: 10A)
     
     integer    :: sz
  end type input_impligparameter
  type(input_impligparameter), save :: inimplig


  !>  <<meaning>>
  !>  isite: implicit-ligand site number (isite=1..inimplig%nsite_implig)
  !>  integer :: naa_site_implig(isite) !! max number of binding site (residues) for each implicit ligand {isite} 
  integer, save :: naa_site_implig(MXSITE_IMPLIG)


  !> these infromation (bindsiteparameters) are given by imput file (<<<< binding_site part).
  type input_implig_bindsitepara
     integer :: nunit_implig(MXSITE_IMPLIG) !< max number of unit for each implicit ligand {isite} 
     integer :: naa_implig(MXUNIT_IMPLIG, MXSITE_IMPLIG) !< max number of binding sites for each {iu, isite} [iu=1..nunit_implig(isite)]
     integer :: idunit_implig(MXUNIT_IMPLIG, MXSITE_IMPLIG) !< unit&state-ID (number) for each {iu, isite}    [iu=1..nunit_implig(isite)] 
     integer :: list_aa_implig(MXAA_IMPLIG, MXSITE_IMPLIG) !< residue-ID (number) for each {iaa, isite}    [iaa=1..naa_site_implig(isite)]
     integer :: sz
  end type input_implig_bindsitepara
  type(input_implig_bindsitepara), save :: inimplig_bindsite

           
  !> <<meaning>>
  !> binding site(residue) pairs for each ligand-mediated contact {icon_implig}
  !> icon2mp_implig(1, icon_implig) = resid_i,  
  !> icon2mp_implig(2, icon_implig) = resid_j,
  !> (resid_i < resid_j), 
  !> ligand-mediated contact {icon_implig} is composed of resi_i and resi_j.
  integer, save :: icon2mp_implig(2,MXCON_IMPLIG)


  integer, save :: ncon_implig(MXSITE_IMPLIG)                 !< max number of ligand-mediated contact for each implicit-ligand {isite}
  real(PREC), save :: vdwrad_implig(MXCON_IMPLIG)             !< native-length for each ligand-mediated contact {icon_implig}
  real(PREC), save :: Ebind_implig(MXCON_IMPLIG)              !< implicit-ligand binding energy for each ligand-mediated contact {icon_implig}
  real(PREC), save :: Etbind_implig(MXSITE_IMPLIG, MXREPLICA) !< implicit-ligand binding energy for each implicit-liagnd site {isite} 
  integer, save :: istate_implig(MXSITE_IMPLIG, MXREPLICA)    !< state (BOUND or UNBOUND) for each implicit-ligand {isite} at time step-t.
  real(PREC), save :: react_rate_implig(MXSITE_IMPLIG)        !< reaction rate for each implicit-ligand site {isite}
  real(PREC), save :: p_implig(MXSITE_IMPLIG)
end module var_implig
