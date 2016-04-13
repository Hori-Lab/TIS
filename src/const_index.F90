! const_index
!> @brief Module for defining index number

module const_index

  implicit none

  ! i_run_mode: define basic simulation mode 
  type run_mode
     integer :: CHECK_FORCE !< 1: Debug Mode, Check the consistence between force and energy
     integer :: CONST_TEMP  !< 2: Constant temperature simulation 
     integer :: SA          !< 3: Simulated annealing (require "<<<< annealing" field)
     integer :: SEARCH_TF   !< 4: Auto-search of T_f (require "<<<< searching_tf" field)
     integer :: ENERGY_CALC !< 5: Energy calculation at single point
     integer :: REPLICA     !< 6: Replica exchange method
     integer :: FMAT        !< 7: Fluctuation matching method
     integer :: ENERGY_DCD  !< 8: Energy calculation for DCD trajectory
     integer :: EMIN        !< 9: Energy minimization
  endtype run_mode
  type(run_mode), parameter :: RUN = run_mode(1,2,3,4,5,6,7,8,9)

  ! i_simulate_type: define dynamics 
  type simu_type
     integer :: CONST_ENERGY !< 0: Newtonian dynamics (velocity Verlet) with the constant energy
     integer :: LANGEVIN     !< 1: Langevin dynamics (recommended)
     integer :: BERENDSEN    !< 2: Newtonian dynamics (velocity Verlet) with Berendsen thermostat
     integer :: NOSEHOOVER   !< 3: Newtonian dynamics (velocity Verlet) with Nose-Hoover thermostat
     integer :: MPC          !< 4: MPC dynamics
     integer :: BROWNIAN     !< 5: Brownian dynamics without hydrodynamic interaction
     integer :: BROWNIAN_HI  !< 6: Brownian dynamics with hydrodynamic interaction
  endtype simu_type
  type(simu_type), parameter :: SIM = simu_type(0,1,2,3,4,5,6)

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
     integer :: DNA    !<  2: DNA  (SPN.1)
     integer :: DNA2   !<  3: DNA2 (SPN.2)
     integer :: PRO    !<  4: Protein
     integer :: RNA    !<  5: RNA
     integer :: LIG    !<  6: Ligand (not implicit)
     integer :: ION    !<  7: Ion
     integer :: MAX    !<  8: Maximum value
  endtype class_type
  type(class_type), parameter :: CLASS = class_type(0,1,2,3,4,5,6,7,8)

  type energy_type
     integer :: TOTAL         !<  1: Total energy
     integer :: VELO          !<  2: Kinetic energy
     integer :: BOND          !<  3: Bond length
     integer :: BANGLE        !<  4: Bond angle
     integer :: DIHE          !<  5: Dihedral
     integer :: GO            !<  6: Go
     integer :: EXV           !<  7: Excluded volume
     integer :: STACK_DNA     !<  8: Base stacking for DNA
     integer :: BP_DNA        !<  9: Base paring for DNA
     integer :: BP_AT         !< 10: DNA AT base pair
     integer :: BP_GC         !< 11: DNA GC base pair
     integer :: MBP           !< 12: DNA mismatch base pair
     integer :: EXV_DNA       !< 13: DNA excluded volume
     integer :: EXV_DNA2      !< 14: DNA2 excluded volume 
     integer :: ELE           !< 15: Electrostatic (Debye Huckel)
     integer :: SOLV_DNA      !< 16: DNA solvation
     integer :: CORE          !< 17: Lipid core (Brown)
     integer :: INT           !< 18: Lipid int (Brown)
     integer :: TAIL          !< 19: Lipid tail (Brown)
     integer :: CORE_NOGU     !< 20: Lipid core (NOGUCHI)
     integer :: TAIL_NOGU     !< 21: Lipid tail (NOGUCHI)
     integer :: BOX           !< 22: Box baundary
     integer :: BRIDGE        !< 23: Bridge constraining
     integer :: PULLING       !< 24: Pulling
     integer :: ANCHOR        !< 25: Anchor constraining
     integer :: DIHE_HARMONIC !< 26: Harmonic dihedral angle
     integer :: HPENE         !< 27: Hydrophobic interaction
     integer :: IMPLIG        !< 28: Implicit ligand
     integer :: MORSE         !< 29: Morse-type Go
     integer :: STACK_RNA     !< 30: RNA base stacking (distance)
     integer :: PAIR_RNA      !< 31: RNA base pairing
     integer :: LJ_ION        !< 32: Ion LJ
     integer :: HYD_ION       !< 33: Ion hydration
     integer :: EXV_ION       !< 34: Ion excluded volume
     integer :: REST1D        !< 35: 1D-restraint
     integer :: CAP           !< 36: Cap boundary
     integer :: WINDOW        !< 37: Window potential
     integer :: SASA          !< 38: SASA potential !sasa
     integer :: EXV_WCA       !< 39: Excluded volume with Weeks-Chandler-Andersen potential
     integer :: STACK_DTRNA   !< 40: Denesyuk-Thirumalai RNA base stacking
     integer :: HBOND_DTRNA   !< 41: Denesyuk-Thirumalai RNA hydrogen bond
     integer :: CYLINDER      !< 42: Cylindrical boundary
     integer :: EXV_DT15      !< 43: excluded volume of DTRNA2015 model
     integer :: TSTACK_DTRNA  !< 44: Denesyuk-Thirumalai RNA base stacking
     integer :: THBOND_DTRNA  !< 45: Denesyuk-Thirumalai RNA hydrogen bond
     integer :: MAX           !< 45: Max value
  endtype energy_type
  type(energy_type), parameter :: E_TYPE  &
     = energy_type(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20, &
                   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40, &
                   41,42,43,44,45,45)
     
  type local_interaction_type
     integer :: NOTHING     !<  1: no interaction
     integer :: L_GO        !<  2: local Go interaction
     integer :: L_AICG1     !<  3: local AICG1
     integer :: L_AICG2     !<  4: local AICG2
     integer :: L_FLP       !<  5: flexible local potential
     integer :: L_BOND      !<  6: bond potential only
     integer :: L_ENM       !<  7: local elastic network model(=NOTHING)
     integer :: L_BDNA      !<  8: local DNA interaction (not yet released)
     integer :: L_DNA2      !<  9: local DNA interaction (not yet released)
     integer :: L_LIP_BROWN !< 10: local lipid-lipid(Brown, not yet released)
     integer :: L_LIP_NOGU  !< 11: local lipid-lipid(Noguchi, not yet released)
     integer :: L_RIGID_LIG !< 12: ligand rigid interaction (not yet released)
     integer :: L_AICG2_PLUS!< 13: local AICG2_PLUS
     integer :: L_DTRNA     !< 14: Denesyuk-Thirumalai RNA model
     integer :: L_DNA2C     !< 15: local DNA interaction (3SPN.2C model)
     integer :: MAX         !< 15: Maximum value
  endtype local_interaction_type
  type(local_interaction_type), parameter :: LINTERACT  & 
     = local_interaction_type(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15)


  type prime_interaction_type
     integer :: NOTHING   !<  1: no interaction (default)
     integer :: GO        !<  2: 12-10 Go potential for native contact pairs
     integer :: EXV       !<  3: (c/r)**12 repulsive interaction
     integer :: ELE       !<  7: electrostatic interaction (Debye-Huckel form)
     integer :: DNA       !< 11: DNA-DNA interaction (not yet released)
     integer :: LIP_BROWN !< 13: lipid-lipid (Brown, not yet relesed)
     integer :: LIP_NOGU  !< 17: lipid-lipid (Noguchi, without solvation, not yet released)
     integer :: LIP_SOLV  !< 19: solvation (lipid-lipid, lipid-protein, not yet released)
     integer :: ENM       !< 23: elastic network model (protein)
     integer :: HP        !< 29: hydrophobic interaction
     integer :: MORSE     !< 31: Morse Go potential (not yet released)
     integer :: PAIR_RNA  !< 37: RNA-RNA base pair (not yet released)
     integer :: ION_HYD   !< 41: hydration interaction (ion-ion, ion-phosphate) (not yet released)
     integer :: ION_EXV   !< 43: repulsive interaction (ion-ion) (not yet released)
     integer :: AICG1     !< 47: AICG1(protein)
     integer :: AICG2     !< 53: AICG2(protein)
     integer :: SASA      !< 59: SASA(protein) !sasa
     integer :: DNA2      !< 61: DNA-DNA interaction (not yet released)
     integer :: DNA2C     !< 62: DNA-DNA interaction (3SPN.2C model)
     integer :: DTRNA     !< 67: Denesyuk-Thirumalai RNA model
     integer :: EXV_WCA   !< 71: Excluded volume with Weeks-Chandler-Andersen potential
     integer :: MAX       !< 71: Maximum value
  endtype prime_interaction_type
  type(prime_interaction_type), parameter :: PRIME_INTERACT  & 
       = prime_interaction_type(1,2,3,7,11,13,17,19,23,29,31,37,41,43,47,53,59,&
                                61,62,67,71,71) 


  type interaction_type
     integer :: NOTHING   !<  1: no interaction (default)
     integer :: GO        !<  2: 12-10 Go potential for native contact pairs
     integer :: EXV       !<  3: (c/r)**12 repulsive interaction
     integer :: ELE       !<  4: electrostatic interaction (Debye-Huckel form)
     integer :: DNA       !<  5: DNA-DNA interaction (not yet released)
     integer :: DNA2      !<  6: DNA-DNA interaction for 3SPN.2 model (not yet released)
     integer :: LIP_BROWN !<  7: lipid-lipid (Brown, not yet relesed)
     integer :: LIP_NOGU  !<  8: lipid-lipid (Noguchi, without solvation, not yet released)
     integer :: LIP_SOLV  !<  9: solvation (lipid-lipid, lipid-protein, not yet released)
     integer :: ENM       !< 10: elastic network model (protein)
     integer :: HP        !< 11: hydrophobic interaction
     integer :: MORSE     !< 12: Morse Go potential (not yet released)
     integer :: PAIR_RNA  !< 13: RNA-RNA base pair (not yet released)
     integer :: ION_HYD   !< 14: hydration interaction (ion-ion, ion-phosphate) (not yet released)
     integer :: ION_EXV   !< 15: repulsive interaction (ion-ion) (not yet released)
     integer :: AICG1     !< 16: AICG1(protein)
     integer :: AICG2     !< 17: AICG2(protein)
     integer :: SASA      !< 18: SASA(protein) !sasa
     integer :: DTRNA     !< 19: Denesyuk-Thirumalai RNA model
     integer :: EXV_WCA   !< 20: Excluded volume with Weeks-Chandler-Andersen potential
     integer :: EXV_DT15  !< 21: Excluded volume with Weeks-Chandler-Andersen potential
     integer :: DNA2C     !< 22: DNA-DNA interaction for 3SPN.2C model
     integer :: MAX       !< 22: Maximum value
  endtype interaction_type
  type(interaction_type), parameter :: INTERACT  & 
     = interaction_type(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,22)


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

  type use_atom_dna
     integer :: COM_PS  !< 0
  endtype use_atom_dna
  type(use_atom_dna), parameter :: USE_DNA = use_atom_dna(0)

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
     integer :: DNA_PHOS   !<  2: Phosphate of DNA
     integer :: DNA_BASE   !<  3: Base of DNA
     integer :: DNA_SUGAR  !<  4: Sugar of DNA
     integer :: DNA2_PHOS  !<  5: Phosphate of DNA
     integer :: DNA2_BASE  !<  6: Base of DNA
     integer :: DNA2_SUGAR !<  7: Sugar of DNA
     integer :: RNA_PHOS   !<  8: Phosphate of RNA
     integer :: RNA_BASE   !<  9: Base of RNA
     integer :: RNA_SUGAR  !< 10: Sugar of RNA
     integer :: LIP_CORE   !< 11 
     integer :: LIP_TAIL   !< 12
     integer :: ION_MG     !< 13
     integer :: ION_NA     !< 14
     integer :: ION_K      !< 15
     integer :: ION_CL     !< 16
  endtype mp_type
  type(mp_type), parameter :: MPTYPE = mp_type(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)

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
     integer :: DP      !< 27: Phosphate(DNA)
     integer :: DS      !< 28: Deoxyribose sugar(DNA)
     integer :: DA      !< 29: Adenine(DNA)
     integer :: DG      !< 30: Guanine(DNA)
     integer :: DT      !< 31: Thymine(DNA)
     integer :: DC      !< 32: Cytosine(DNA)
     integer :: MG      !< 33: ion
     integer :: K       !< 34: ion
     integer :: CL      !< 35: ion
     integer :: MAX     !< 35
  endtype chemical_type
  type(chemical_type), parameter :: CHEMICALTYPE = &
       chemical_type( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,&
                     20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,35)

  type ion_type
     integer :: VOID      !< 0
     integer :: NA        !< 1: Na+
     integer :: K         !< 2: K+
     integer :: CL        !< 3: Cl- of DNA
     integer :: MG        !< 4: Mg2+ of DNA
     integer :: MAX_ION   !< 4: Max number of ion
     integer :: P         !< 5: Phosphate of DNA
     integer :: MAX_ALL   !< 5: Max number of ion and phosphate
  endtype ion_type
  type(ion_type), parameter :: IONTYPE = ion_type(0,1,2,3,4,4,5,5)

  ! categories for RNA's bond
  type bd_type
     integer :: VOID     !< 0
     integer :: PRO      !< 1
     integer :: RNA_PS   !< 2
     integer :: RNA_SB   !< 3
     integer :: RNA_SR   !< 4
     integer :: RNA_SY   !< 5
     integer :: RNA_SP   !< 6
     integer :: DNA2_PS  !< 7
     integer :: DNA2_SP  !< 8
     integer :: DNA2_SA  !< 9
     integer :: DNA2_ST  !< 10
     integer :: DNA2_SG  !< 11
     integer :: DNA2_SC  !< 12
  endtype bd_type
  type(bd_type), parameter :: BDTYPE = bd_type(0,1,2,3,4,5,6,7,8,9,10,11,12)

  type ba_type
     integer :: VOID     !< 0
     integer :: PRO      !< 1
     integer :: DNA2_SPS !< 2
     integer :: DNA2_PSP !< 3
     integer :: DNA2_PSA !< 4
     integer :: DNA2_PST !< 5
     integer :: DNA2_PSC !< 6
     integer :: DNA2_PSG !< 7
     integer :: DNA2_ASP !< 8
     integer :: DNA2_TSP !< 9
     integer :: DNA2_CSP !< 10
     integer :: DNA2_GSP !< 11
     integer :: RNA_BSP  !< 12
     integer :: RNA_RSP  !< 13
     integer :: RNA_YSP  !< 14
     integer :: RNA_PSP  !< 15
     integer :: RNA_SPS  !< 16
     integer :: RNA_PSB  !< 17
     integer :: RNA_PSR  !< 18
     integer :: RNA_PSY  !< 19
  endtype ba_type
  type(ba_type), parameter :: BATYPE = ba_type(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)

  type dih_type
     integer :: VOID      !< 0
     integer :: PRO       !< 1
     integer :: DNA2_PSPS !< 2
     integer :: DNA2_SPSP !< 3
     integer :: RNA_RSPS  !< 4
     integer :: RNA_YSPS  !< 5
     integer :: RNA_PSPS  !< 6
     integer :: RNA_SPSR  !< 7
     integer :: RNA_SPSY  !< 8
     integer :: RNA_SPSP  !< 9
     integer :: DNA_PER1  !< 10
     integer :: DNA_PER2  !< 11
     integer :: DNA_PER3  !< 12
  endtype dih_type
  type(dih_type), parameter :: DIHTYPE = dih_type(0,1,2,3,4,5,6,7,8,9,10,11,12)

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
     integer :: DNA_DNA      !<12  DNA    :  DNA
     integer :: PRO_DNA      !<13 protein :  DNA
     integer :: PRO_DNA2     !<14 protein :  DNA
     integer :: PRO_LIG      !<15 protein :  LIGAND
  endtype con_type
  type(con_type), parameter :: CONTYPE = con_type(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)

  !#######################################
  ! DNA
  !#######################################
  type bp_type
     integer :: VOID !< 0
     integer :: AA   !< 1
     integer :: AT   !< 2
     integer :: AC   !< 3
     integer :: AG   !< 4
     integer :: TA   !< 5
     integer :: TT   !< 6
     integer :: TC   !< 7
     integer :: TG   !< 8
     integer :: CA   !< 9
     integer :: CT   !< 10
     integer :: CC   !< 11
     integer :: CG   !< 12
     integer :: GA   !< 13
     integer :: GT   !< 14
     integer :: GC   !< 15
     integer :: GG   !< 16
     integer :: MAX  !< 16
  endtype bp_type
  type(bp_type), parameter :: BPTYPE = bp_type(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16)

    type base_type
     integer :: VOID !< 0
     integer :: A    !< 1
     integer :: T    !< 2
     integer :: C    !< 3
     integer :: G    !< 4
     integer :: S    !< 5
     integer :: P    !< 6
     integer :: MAX  !< 6
  endtype base_type
  type(base_type), parameter :: BASETYPE = base_type(0,1,2,3,4,5,6,6)

  
  !#######################################
  ! implicit ligand
  !#######################################
  type implig_ene_function_type
     integer :: FUNC_LJ12_10     !< 0: LJ12-10 type 
     integer :: FUNC_GAUSSIAN    !< 1: Gaussian type
  endtype implig_ene_function_type
  type(implig_ene_function_type), parameter :: IMPLIGENE_FUNCTYPE = implig_ene_function_type(0,1)
  ! this constant is used for simu_energy_implig, simu_force_implig.


  type implig_energy_type
     integer :: FOR_MC       !< 0:   
     integer :: FOR_NON_MC   !< 1:
  endtype implig_energy_type
  type(implig_energy_type), parameter :: IMPLIGENERGY_TYPE = implig_energy_type(0,1)
  ! this constant is used especially for simu_energy_implig.
  ! energy FOR_NON_MC is calculated based on not only structure but also state of ligand.
  ! energy FOR_NON_MC is for output (ts_file).
  ! energy FOR_MC is calculated based on structre.
  ! energy FOR_MC is for Monte Calro (simu_mc_implig). 
  

  type implig_bound_state
     integer :: UN_BOUND  !< 0:
     integer :: BOUND     !< 1:  
  endtype implig_bound_state
  type(implig_bound_state), parameter :: IMPLIGBOUND_STATE = implig_bound_state(0,1)
  ! this constant is used for simu_energy_implig, simu_force_implig.
  

  !#######################################
  ! replica exchange
  !#######################################

  type replica_type
     integer :: TEMP   !< 1: Temperature
     integer :: ION    !< 2: Ionic strength
     integer :: PULL   !< 3: Pulling force
     integer :: WIND   !< 4: Window
     integer :: WINZ   !< 5: Window2
     integer :: MAX    !< 5: Maximum value
  endtype replica_type
  type(replica_type), parameter :: REPTYPE = replica_type(1,2,3,4,5,5)

  character(4), parameter :: CHAR_REPTYPE(REPTYPE%MAX) &
     = (/'temp','ion ','pull','wind', 'winz'/)

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

  type multiscale
     integer :: HOMO  !< 0: conventional homogeneous Go; (Default)
     integer :: AUTO  !< 1: multiscale Go; parameters generated by CafeMol
     integer :: USER  !< 2: multiscale Go; parameters given by users
  endtype multiscale
  type(multiscale), parameter :: AICG = multiscale(0,1,2)

  type fluctuation_matching_type
     integer :: VOID    !< 0:
     integer :: HOMO    !< 1:
     integer :: HETERO  !< 2:
     integer :: MAX
  endtype fluctuation_matching_type
  type(fluctuation_matching_type), parameter :: FMATTYPE &
                                                = fluctuation_matching_type(0,1,2,2)

  type rst_block
     integer :: STEP
     integer :: XYZ
     integer :: VELO
     integer :: ACCEL
     integer :: REPLICA
  endtype rst_block
  type(rst_block), parameter :: RSTBLK = rst_block(1,2,3,4,5)

  type emin_type
     integer :: VOID
     integer :: SD
     integer :: CG
     integer :: MAX
  endtype emin_type
  type(emin_type), parameter :: EMIN_METHOD = emin_type(0,1,2,2)
  
  type record_file_type
     integer :: CRD
     integer :: VELO
  endtype record_file_type
  type(record_file_type), parameter :: RECORD_FILE = record_file_type(1,2)

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
     integer :: P
     integer :: S
     integer :: A
     integer :: G
     integer :: C
     integer :: U
     integer :: MG2
     integer :: CA2
     integer :: CL
     integer :: K
     integer :: MAX
  endtype dtrna15_exv
  type(dtrna15_exv), parameter :: DT15EXV = dtrna15_exv(1,2,3,4,5,6,7,8,9,10,10)

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
     integer :: K    !< 25
     integer :: CL   !< 26
     integer :: RA   !< 27
     integer :: RC   !< 28
     integer :: RG   !< 29
     integer :: RT   !< 30
     integer :: RU   !< 31
     integer :: RI   !< 32
     integer :: DA   !< 33
     integer :: DC   !< 34
     integer :: DG   !< 35
     integer :: DT   !< 36
     integer :: DU   !< 37
     integer :: DI   !< 38
     integer :: MAX  !< 38
  endtype seq2id
  type(seq2id), parameter :: SEQID = seq2id(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,&
                             20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38, 38)
  
endmodule const_index
