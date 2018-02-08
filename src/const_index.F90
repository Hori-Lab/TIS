! const_index
!> @brief Module for defining index number

module const_index

  implicit none

  ! i_run_mode: define basic simulation mode 
  type run_mode
     integer :: CHECK_FORCE !< 1: Debug Mode, Check the consistence between force and energy
     integer :: CONST_TEMP  !< 2: Constant temperature simulation 
     integer :: SA          !< 3: Simulated annealing (require "<<<< annealing" field)
     integer :: ENERGY_CALC !< 5: Energy calculation at single point
     integer :: REPLICA     !< 6: Replica exchange method
     integer :: FMAT        !< 7: Fluctuation matching method
     integer :: ENERGY_DCD  !< 8: Energy calculation for DCD trajectory
     integer :: EMIN        !< 9: Energy minimization
     integer :: WIDOM       !<10: Widom method to calculate chemical potential(s)
  endtype run_mode
  type(run_mode), parameter :: RUN = run_mode(1,2,3,5,6,7,8,9,10)

  ! i_simulate_type: define dynamics 
  type simu_type
     integer :: CONST_ENERGY !< 0: Newtonian dynamics (velocity Verlet) with the constant energy
     integer :: LANGEVIN     !< 1: Langevin dynamics (recommended)
     integer :: BERENDSEN    !< 2: Newtonian dynamics (velocity Verlet) with Berendsen thermostat
     integer :: NOSEHOOVER   !< 3: Newtonian dynamics (velocity Verlet) with Nose-Hoover thermostat
     integer :: MPC          !< 4: MPC dynamics
     integer :: BROWNIAN     !< 5: Brownian dynamics without hydrodynamic interaction
     integer :: BROWNIAN_HI  !< 6: Brownian dynamics with hydrodynamic interaction
     integer :: PS_BROWNIAN  !< 7: Brownian dynamics without hydrodynamic interaction
  endtype simu_type
  type(simu_type), parameter :: SIM = simu_type(0,1,2,3,4,5,6,7)

  ! i_hydro_tensor: define tensor type of hydrodynamic interaction (only when SIM%BROWNIAN_HI)
  type tensor_type
     integer :: RPY        !< 1: Rotne-Prager-Yamakawa (RPY)
     integer :: RPY_OVER   !< 2: RPY (overlap allowed)
     integer :: ERMAK_OVER !< 3: Ermak-McCammon's modified RPY (overlap allowed)
     integer :: ZUK_RPY    !< 4: Zuk et al's modified RPY (overlap and different radii allowed)
  endtype tensor_type
  type(tensor_type), parameter :: HI_TENSOR = tensor_type(1,2,3,4)

  type initial_state
     integer :: VOID    !< 0: Invalid
     integer :: RANDOM  !< 1: Random chain polymer
     integer :: NATIVE  !< 2: Copy from native-state structure
     integer :: INPUT   !< 3: Read from input file
     integer :: BDNA    !< 4: B-type DNA is constructed automatically
     integer :: LIPID   !< 5: Lipid structure is constructed
     integer :: CG      !< 6: Read from input file with cafemol(CG) style
     integer :: CARD    !< 7: Read from CARD-style file (.crd)
     integer :: RST     !< 8: Read from restart file (.rst)
  endtype initial_state
  type(initial_state), parameter :: INISTAT = initial_state(0,1,2,3,4,5,6,7,8)

  type initial_velo
     integer :: MAXWELL  !< 0: Boltzmann distribution (default)
     integer :: CARD     !< 1: Read from CARD-style file (.velo)
     integer :: RST      !< 2: Read from restart file (.rst)
  endtype initial_velo
  type(initial_velo), parameter :: INIVELO = initial_velo(0,1,2)

  type seq_read_style
     integer :: PDB         !< 1: Taken from PDB structure
     integer :: INPUT_SEQ   !< 2:
     integer :: INPUT_LIPID !< 3:
     integer :: CG          !< 4: Taken from CafeMol(CG) style
  endtype seq_read_style
  type(seq_read_style), parameter :: SEQREAD = seq_read_style(1,2,3,4)

  type native_read_style
     integer :: PDB         !< 1:
     integer :: INFO        !< 2:
     integer :: NO          !< 3:
  endtype native_read_style
  type(native_read_style), parameter :: NATIVEREAD = native_read_style(1,2,3)

  type native_info_style
     integer :: ALL_IN_ONE
     integer :: ONE_BY_ONE
  endtype native_info_style
  type(native_info_style), parameter :: NATIVEINFO = native_info_style(1,2)

  type class_type
     integer :: VOID
     integer :: LIP    !<  1: Lipid
     integer :: PRO    !<  2: Protein
     integer :: RNA    !<  3: RNA
     integer :: LIG    !<  4: Ligand (not implicit)
     integer :: ION    !<  5: Ion
     integer :: MAX    !<  6: Maximum value
  endtype class_type
  type(class_type), parameter :: CLASS = class_type(0,1,2,3,4,5,6)

  type energy_type
     integer :: TOTAL         !<  1: Total energy
     integer :: VELO          !<  2: Kinetic energy
     integer :: BOND          !<  3: Bond length
     integer :: BANGLE        !<  4: Bond angle
     integer :: DIHE          !<  5: Dihedral
     integer :: GO            !<  6: Go
     integer :: EXV12         !<  7: Excluded volume
     integer :: ELE           !<  8: Electrostatic (Debye Huckel)
     !integer :: BOX           !<  9: Box baundary
     integer :: BRIDGE        !< 10: Bridge constraining
     integer :: PULLING       !< 11: Pulling
     integer :: ANCHOR        !< 12: Anchor constraining
     integer :: DIHE_HARMONIC !< 13: Harmonic dihedral angle
     !integer :: HPENE         !< 14: Hydrophobic interaction
     !integer :: IMPLIG        !< 15: Implicit ligand
     !integer :: MORSE         !< 16: Morse-type Go
     !integer :: STACK_RNA     !< 17: RNA base stacking (distance)
     !integer :: PAIR_RNA      !< 18: RNA base pairing
     !integer :: LJ_ION        !< 19: Ion LJ
     integer :: REST1D        !< 20: 1D-restraint
     !integer :: CAP           !< 21: Cap boundary
     integer :: WINDOW        !< 22: Window potential
     !integer :: SASA          !< 23: SASA potential !sasa
     integer :: EXV_WCA       !< 24: Excluded volume with Weeks-Chandler-Andersen potential
     integer :: STACK_DTRNA   !< 25: Denesyuk-Thirumalai RNA base stacking
     integer :: HBOND_DTRNA   !< 26: Denesyuk-Thirumalai RNA hydrogen bond
     !integer :: CYLINDER      !< 27: Cylindrical boundary
     integer :: EXV_DT15      !< 28: excluded volume of DTRNA2015 model
     integer :: TSTACK_DTRNA  !< 29: Denesyuk-Thirumalai RNA base stacking
     integer :: THBOND_DTRNA  !< 30: Denesyuk-Thirumalai RNA hydrogen bond
     integer :: EXV6          !< 31: Excluded volume
     integer :: EXV_GAUSS     !< 32: Excluded volume
     integer :: WCA_REP
     integer :: WCA_ATT
     integer :: BBR
     integer :: MAX           !< Max value
  endtype energy_type
  type(energy_type), parameter :: E_TYPE  &
     != energy_type(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20, &
     !              21,22,23,24,25,26,27,28,29,30,31,31)
     = energy_type(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,&
                   21,22,23,24,25,25)
     
  type local_interaction_type
     integer :: NOTHING     !<  1: no interaction
     integer :: L_GO        !<  2: local Go interaction
     integer :: L_BOND      !<  3: bond potential only
     integer :: L_ANGL      !<  4: angle potential
     integer :: L_RIGID_LIG !<  5: ligand rigid interaction (not yet released)
     integer :: L_DTRNA     !<  6: Denesyuk-Thirumalai RNA model
     integer :: L_FENE      !<  7: FENE
     integer :: L_ROUSE     !<  8: Rouse
     integer :: MAX         !< Maximum value
     !integer :: L_FLP       !<  5: flexible local potential
     !integer :: L_ENM       !<  7: local elastic network model(=NOTHING)
  endtype local_interaction_type
  type(local_interaction_type), parameter :: LINTERACT  & 
     != local_interaction_type(1,2,3,4,5,6,7,8,9,10,11,11)
     = local_interaction_type(1,2,3,4,5,6,7,8,8)


  type interaction_type
     integer :: NOTHING   !<  1: no interaction (default)
     integer :: GO        !<  2: 12-10 Go potential for native contact pairs
     integer :: EXV12     !<  3: (c/r)**12 repulsive interaction
     integer :: EXV6      !<  4: (c/r)**6 repulsive interaction
     integer :: ELE       !<  5: electrostatic interaction (Debye-Huckel form)
     !integer :: ENM       !<  6: elastic network model (protein)
     !integer :: HP        !<  7: hydrophobic interaction
     !integer :: MORSE     !<  8: Morse Go potential (not yet released)
     !integer :: PAIR_RNA  !<  9: RNA-RNA base pair (not yet released)
     !integer :: AICG1     !< 10: AICG1(protein)
     !integer :: AICG2     !< 11: AICG2(protein)
     !integer :: SASA      !< 12: SASA(protein) !sasa
     integer :: DTRNA     !< 13: Denesyuk-Thirumalai RNA model
     integer :: EXV_WCA   !< 14: Excluded volume with Weeks-Chandler-Andersen potential
     integer :: EXV_DT15  !< 15: Excluded volume with Weeks-Chandler-Andersen potential
     integer :: LJ        !< 16: Lenard-Jones
     integer :: EXV_GAUSS !< 17: Excluded volume with Gaussian function
     integer :: CON_GAUSS !< 
     integer :: WCA
     integer :: MAX       !< Maximum value
  endtype interaction_type
  type(interaction_type), parameter :: INTERACT  & 
     != interaction_type(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16)
     = interaction_type(1,2,3,4,5,6,7,8,9,10,11,12,12)


  type error_handling
     ! Treatment
     !     0* STOP: program stops
     !     1* WARN: program continues
     ! Output:
     !     *0 STD:  standard output (console)
     !     *1 FILE: file
     !     *2 ALL:  file & standard output
     integer :: STOP_STD   !<  0: Outputs error message on STD, then stops.
     integer :: STOP_FILE  !<  1: Outputs error message on .data file, then stops.
     integer :: STOP_ALL   !<  2: Outputs error message on both STD and .data file, then stops.
     integer :: WARN_STD   !< 10: Outputs error message on STD, then continues.
     integer :: WARN_FILE  !< 11: Outputs error message on .data file, then continues.
     integer :: WARN_ALL   !< 12: Outputs error message on both STD and .data file, then continues.
  endtype error_handling
  type(error_handling), parameter :: ERROR = error_handling(0,1,2,10,11,12)

  type use_atom_protein
     integer :: CA        !< 0: Alpha carbon
     integer :: CB        !< 1: Beta carbon
     integer :: COM_SIDE  !< 2: Center of mass of side chain
     integer :: MAX       !< 2: Maximum value
  endtype use_atom_protein
  type(use_atom_protein), parameter :: USE_PRO = use_atom_protein(0,1,2,2)

  type use_atom_rna_base
     integer :: COM        !< 0
     integer :: PuN1_PyN3  !< 1
  endtype use_atom_rna_base
  type(use_atom_rna_base), parameter :: USE_RNA_BASE = use_atom_rna_base(0,1)

  type use_atom_rna_sugar
     integer :: COM      !< 0
     integer :: COM_RING !< 1
     integer :: C4       !< 2
  endtype use_atom_rna_sugar
  type(use_atom_rna_sugar), parameter :: USE_RNA_SUGAR = use_atom_rna_sugar(0,1,2)

  type mp_type
     integer :: VOID       !<  0
     integer :: PRO        !<  1: Protein
     integer :: RNA_PHOS   !<  2: Phosphate of RNA
     integer :: RNA_BASE   !<  3: Base of RNA
     integer :: RNA_SUGAR  !<  4: Sugar of RNA
     integer :: ION_MG     !<  5
     integer :: ION_CA2    !<  6
     integer :: ION_NA     !<  7
     integer :: ION_K      !<  8
     integer :: ION_CL     !<  9
     integer :: LIG_X1     !<  10
  endtype mp_type
  type(mp_type), parameter :: MPTYPE = mp_type(0,1,2,3,4,5,6,7,8,9,10)

  type chemical_type
     integer :: UNKNOWN !<  0: Default value is stored here
     integer :: ALA     !<  1
     integer :: ARG     !<  2
     integer :: ASN     !<  3
     integer :: ASP     !<  4
     integer :: CYS     !<  5
     integer :: GLN     !<  6
     integer :: GLU     !<  7
     integer :: GLY     !<  8
     integer :: HIS     !<  9
     integer :: ILE     !< 10
     integer :: LEU     !< 11
     integer :: LYS     !< 12
     integer :: MET     !< 13
     integer :: PHE     !< 14
     integer :: PRO     !< 15
     integer :: SER     !< 16
     integer :: THR     !< 17
     integer :: TRP     !< 18
     integer :: TYR     !< 19
     integer :: VAL     !< 20
     integer :: P       !< 21: Phosphate
     integer :: S       !< 22: Ribose sugar
     integer :: A       !< 23: Adenine
     integer :: G       !< 24: Guanine
     integer :: U       !< 25: Uracyl
     integer :: C       !< 26: Cytosine
     integer :: MG      !< 27: ion
     integer :: CA2     !< 28: ion
     integer :: K       !< 29: ion
     integer :: NA      !< 30: ion
     integer :: CL      !< 31: ion
     integer :: X1      !< 32: ligand (for rod)
     integer :: MAX     !< 32
  endtype chemical_type
  type(chemical_type), parameter :: CHEMICALTYPE = &
       chemical_type( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,&
                     20,21,22,23,24,25,26,27,28,29,30,31,32,32)

  type ion_type
     integer :: VOID      !< 0
     integer :: NA        !< 1: Na+
     integer :: K         !< 2: K+
     integer :: CL        !< 3: Cl- of DNA
     integer :: MG        !< 4: Mg2+ of DNA
     integer :: CA2       !< 5: Ca2+ of DNA
     integer :: MAX_ION   !< 5: Max number of ion
     integer :: P         !< 6: Phosphate of DNA
     integer :: MAX_ALL   !< 6: Max number of ion and phosphate
  endtype ion_type
  type(ion_type), parameter :: IONTYPE = ion_type(0,1,2,3,4,5,5,6,6)

  ! categories for RNA's bond
  type bd_type
     integer :: VOID     !< 0
     integer :: PRO      !< 1
     integer :: RNA_PS   !< 2
     integer :: RNA_SB   !< 3
     integer :: RNA_SR   !< 4
     integer :: RNA_SY   !< 5
     integer :: RNA_SP   !< 6
  endtype bd_type
  type(bd_type), parameter :: BDTYPE = bd_type(0,1,2,3,4,5,6)

  type ba_type
     integer :: VOID     !< 0
     integer :: PRO      !< 1
     integer :: RNA_BSP  !< 2
     integer :: RNA_RSP  !< 3
     integer :: RNA_YSP  !< 4
     integer :: RNA_PSP  !< 5
     integer :: RNA_SPS  !< 6
     integer :: RNA_PSB  !< 7
     integer :: RNA_PSR  !< 8
     integer :: RNA_PSY  !< 9
  endtype ba_type
  type(ba_type), parameter :: BATYPE = ba_type(0,1,2,3,4,5,6,7,8,9)

  type dih_type
     integer :: VOID      !< 0
     integer :: PRO       !< 1
     integer :: RNA_RSPS  !< 2
     integer :: RNA_YSPS  !< 3
     integer :: RNA_PSPS  !< 4
     integer :: RNA_SPSR  !< 5
     integer :: RNA_SPSY  !< 6
     integer :: RNA_SPSP  !< 7
  endtype dih_type
  type(dih_type), parameter :: DIHTYPE = dih_type(0,1,2,3,4,5,6,7)

  type con_type
     integer :: VOID
     integer :: PRO_PRO      !< 1 protein : protein
     integer :: PRO_RP       !< 2 protein :  RNA_P
     integer :: PRO_RS       !< 3 protein :  RNA_S
     integer :: PRO_RB       !< 4 protein :  RNA_B
     integer :: RP_RP        !< 5  RNA_P  :  RNA_P
     integer :: RP_RS        !< 6  RNA_P  :  RNA_S
     integer :: RP_RB        !< 7  RNA_P  :  RNA_B
     integer :: RS_RS        !< 8  RNA_S  :  RNA_S
     integer :: RS_RB        !< 9  RNA_S  :  RNA_B
     integer :: RB_RB        !<10  RNA_B  :  RNA_B
     integer :: RNA_BP       !<11  RNA_B  :  RNA_B pairing (dist<dfcontact2_rna_bp)
     integer :: PRO_LIG      !<12 protein :  LIGAND
  endtype con_type
  type(con_type), parameter :: CONTYPE = con_type(0,1,2,3,4,5,6,7,8,9,10,11,12)

  
!  !#######################################
!  ! implicit ligand
!  !#######################################
!  type implig_ene_function_type
!     integer :: FUNC_LJ12_10     !< 0: LJ12-10 type 
!     integer :: FUNC_GAUSSIAN    !< 1: Gaussian type
!  endtype implig_ene_function_type
!  type(implig_ene_function_type), parameter :: IMPLIGENE_FUNCTYPE = implig_ene_function_type(0,1)
!  ! this constant is used for energy_implig, force_implig.
!
!
!  type implig_energy_type
!     integer :: FOR_MC       !< 0:   
!     integer :: FOR_NON_MC   !< 1:
!  endtype implig_energy_type
!  type(implig_energy_type), parameter :: IMPLIGENERGY_TYPE = implig_energy_type(0,1)
!  ! this constant is used especially for energy_implig.
!  ! energy FOR_NON_MC is calculated based on not only structure but also state of ligand.
!  ! energy FOR_NON_MC is for output (ts_file).
!  ! energy FOR_MC is calculated based on structre.
!  ! energy FOR_MC is for Monte Calro (simu_mc_implig). 
!  
!
!  type implig_bound_state
!     integer :: UN_BOUND  !< 0:
!     integer :: BOUND     !< 1:  
!  endtype implig_bound_state
!  type(implig_bound_state), parameter :: IMPLIGBOUND_STATE = implig_bound_state(0,1)
!  ! this constant is used for energy_implig, force_implig.
  

  !#######################################
  ! replica exchange
  !#######################################

  type replica_type
     integer :: TEMP   !< 1: Temperature
     integer :: ION    !< 2: Ionic strength
     integer :: PULL   !< 3: Pulling force
     integer :: WIND   !< 4: Window
     !integer :: WINZ   !< 5: Window2
     integer :: MAX    !< 5: Maximum value
  endtype replica_type
  !type(replica_type), parameter :: REPTYPE = replica_type(1,2,3,4,5,5)
  type(replica_type), parameter :: REPTYPE = replica_type(1,2,3,4,4)

  character(4), parameter :: CHAR_REPTYPE(REPTYPE%MAX) &
     != (/'temp','ion ','pull','wind', 'winz'/)
     = (/'temp','ion ','pull','wind'/)

  type wind_type
     integer :: IMP    !< 1: ID of mass point I
     integer :: JMP    !< 2: ID of mass point J
     integer :: COEF   !< 1: Coeffecient of string of window
     integer :: LENGTH !< 2: Natural length of string of window
     integer :: MAX    !< 2: Maximum value
  endtype wind_type
  type(wind_type), parameter :: WINDTYPE = wind_type(1,2,1,2,2)
       
  type variable_style
     integer :: VOID        !<   0: Invalid (not used in REM)
     integer :: LINEAR      !<   1: Linear interpolation
     integer :: EXPONENTIAL !<  10: Exponential interpolation
     integer :: EXPLICIT    !< 100: Explicit defined
     integer :: MAX         !< 100: Maximum value
  endtype variable_style
  type(variable_style), parameter :: REPVARSTYLE = variable_style(0,1,10,100,100)

  type potential_type
     integer :: LJ1210    !< 1: 12-10 type Lenard-Jones potential
     integer :: MORSE     !< 2: Morse potential
     integer :: MAX
  endtype potential_type
  type(potential_type), parameter :: POTTYPE = potential_type(1,2,2)

!  type multiscale
!     integer :: HOMO  !< 0: conventional homogeneous Go; (Default)
!     integer :: AUTO  !< 1: multiscale Go; parameters generated by CafeMol
!     integer :: USER  !< 2: multiscale Go; parameters given by users
!  endtype multiscale
!  type(multiscale), parameter :: AICG = multiscale(0,1,2)

!  type fluctuation_matching_type
!     integer :: VOID    !< 0:
!     integer :: HOMO    !< 1:
!     integer :: HETERO  !< 2:
!     integer :: MAX
!  endtype fluctuation_matching_type
!  type(fluctuation_matching_type), parameter :: FMATTYPE &
!                                                = fluctuation_matching_type(0,1,2,2)

  type rst_block
     integer :: STEP
     integer :: XYZ
     integer :: VELO
     integer :: ACCEL
     integer :: REPLICA
     integer :: DTRNA15
     integer :: RANDOM
  endtype rst_block
  type(rst_block), parameter :: RSTBLK = rst_block(1,2,3,4,5,6,7)

  type emin_type
     integer :: VOID
     integer :: SD
     integer :: CG
     integer :: MAX
  endtype emin_type
  type(emin_type), parameter :: EMIN_METHOD = emin_type(0,1,2,2)

  type nearest_neighbor
     integer :: AA  !<  1
     integer :: AU  !<  2
     integer :: AG  !<  3
     integer :: AC  !<  4
     integer :: UA  !<  5
     integer :: UU  !<  6
     integer :: UG  !<  7
     integer :: UC  !<  8
     integer :: GA  !<  9
     integer :: GU  !< 10
     integer :: GG  !< 11
     integer :: GC  !< 12
     integer :: CA  !< 13
     integer :: CU  !< 14
     integer :: CG  !< 15
     integer :: CC  !< 16
  endtype nearest_neighbor
  type(nearest_neighbor), parameter :: NN = nearest_neighbor(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)

  type dtrna15_exv
     integer :: P    !<  1
     integer :: S    !<  2
     integer :: A    !<  3
     integer :: G    !<  4
     integer :: C    !<  5
     integer :: U    !<  6
     integer :: MG2  !<  7
     integer :: CA2  !<  8
     integer :: CL   !<  9
     integer :: K    !< 10 
     integer :: NA   !< 11
     integer :: X1   !< 12
     integer :: MAX
  endtype dtrna15_exv
  type(dtrna15_exv), parameter :: DT15EXV = dtrna15_exv(1,2,3,4,5,6,7,8,9,10,11,12,12)

  type seq2id
     integer :: ALA  !< 1
     integer :: ARG  !< 2
     integer :: ASN  !< 3
     integer :: ASP  !< 4
     integer :: CYS  !< 5
     integer :: GLN  !< 6
     integer :: GLU  !< 7
     integer :: GLY  !< 8
     integer :: HIS  !< 9
     integer :: ILE  !< 10
     integer :: LEU  !< 11
     integer :: LYS  !< 12
     integer :: MET  !< 13
     integer :: PHE  !< 14
     integer :: PRO  !< 15
     integer :: SER  !< 16
     integer :: THR  !< 17
     integer :: TRP  !< 18
     integer :: TYR  !< 19
     integer :: VAL  !< 20
     integer :: OTH  !< 21
     integer :: P    !< 22
     integer :: P2   !< 23
     integer :: MG   !< 24
     integer :: CA2  !< 25
     integer :: K    !< 26
     integer :: NA   !< 27
     integer :: CL   !< 28
     integer :: RA   !< 29
     integer :: RC   !< 20
     integer :: RG   !< 31
     integer :: RT   !< 32
     integer :: RU   !< 33
     integer :: RI   !< 34
     integer :: DA   !< 35
     integer :: DC   !< 36
     integer :: DG   !< 37
     integer :: DT   !< 38
     integer :: DU   !< 39
     integer :: DI   !< 40
     integer :: MAX  !< 40
  endtype seq2id
  type(seq2id), parameter :: SEQID = seq2id(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,&
                             20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,40)
  
endmodule const_index
