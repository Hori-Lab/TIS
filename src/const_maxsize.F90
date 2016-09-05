!const_maxsize
!> @brief Most array sizes and limit numbers are defined in this module.

!*********************************************************************
!
! maxsize of array
!
!*********************************************************************
! Abbreviations
!   MX  : maximum
!   MP  : mass point
!   BD  : bond
!   BA  : bond angle
!   DIH : dihedral angle
!   ELE : electrostatic
!   SOLV: solvation
!   BP  : base pair
!   ST  : stack
!   CON : contact
!   MGO : multiple-Go model
!   DEL : delete
!*********************************************************************

module const_maxsize

  implicit none

#ifdef PREC_INTEL16
  integer, parameter :: PREC = 16 !< precision of floating-point
#else
  integer, parameter :: PREC = 8  !< precision of floating-point (default, 4 or 8)
#endif
  integer, save      :: S_REAL    !< byte size of default real   (single precision)
  integer, save      :: M_INT     !< byte size of default integer   (medium-size integer)
  integer, save      :: LOGIC     !< byte size of default logical
  integer, parameter :: L_INT = 8 !< byte size of larger integer ("long long int" in C)

  ! Files
  integer, parameter :: MXPDB = 10    !< maximum # of input PDB files
  integer, parameter :: MXINI = MXPDB  !< maximum # of input structure files
  integer, parameter :: MXPARA = 15     !< maximum # of parameter files !changed for exv
  integer, parameter :: FILENAME_DIGIT_REPLICA = 4  !< # of digits added to .ts filename in REM

  ! MD-control
  integer, parameter :: MXSIM = 10     !< maximum # of simulation stages
!  integer, parameter :: MXSEARCHINGTF = 1000  !< maximum # of cycles for simulated searching TF
  integer, parameter :: MX_REP_OPT_STAGE = 10000 !< maximum # of cycles for Feedback-optimized REM
  integer(L_INT), parameter :: MX_NTSTEP = 50000000000

  ! Molecules and interactions
  integer, parameter :: MXATOM_MP = 1     !< maximum # of atoms for each mass point
  integer, parameter :: MXUNIT = 3        !< maximum # of interaction units
  integer, parameter :: MXMP =  1397      !< maximum # of mass points
  integer, parameter :: MXPDBATOM = 10*MXMP  !< maximum # of atom in pdb file

!  integer, parameter :: MXHBOND = 2 * MXMP !< maximum # of hydrogen bonds
  integer, parameter :: MXMPBD = 1
  integer, parameter :: MXBD = MXMPBD * MXMP        !< maximum # of bonds
  integer, parameter :: MXMPFENE = 1
  integer, parameter :: MXFENE = MXMPFENE * MXMP        !< maximum # of bonds
  integer, parameter :: MXMPBA = 2
  integer, parameter :: MXBA = MXMPBA * MXMP    !< maximum # of bond angles
  integer, parameter :: MXMPDIH = 2
  integer, parameter :: MXDIH = MXMPDIH * MXMP   !< maximum # of dihedral angles
  integer, parameter :: MXMPCON = 20       !< maximum # of Go-type contacts per residue 
  integer, parameter :: MXCON = MXMPCON * MXMP  !< maximum # of total contacts
  integer, parameter :: MXMPLJ = 20       !< maximum # of Go-type contacts per residue 
  integer, parameter :: MXLJ = MXMPLJ * MXMP  !< maximum # of total contacts
!  integer, parameter :: MXMPRNABP = 2
!  integer, parameter :: MXRNABP = MXMPRNABP * MXMP  !< maximum # of contacts
!  integer, parameter :: MXMPRNAST = 2
!  integer, parameter :: MXRNAST = MXMPRNAST * MXMP  !< maximum # of contacts
  integer, parameter :: MXDTRNAST = int(real(MXMP) / 3.0)
  integer, parameter :: MXMPDTRNAHB = 2
  !integer, parameter :: MXDTRNAHB = MXMPDTRNAHB * MXMP / 2
  integer, parameter :: MXDTRNAHB = 5986
  integer, parameter :: MXDTRNATST= 27
  integer, parameter :: MXMPMORSE = 1    !< maximum # of contacts
  integer, parameter :: MXMORSE = MXMPMORSE * MXMP    !< maximum # of contacts
!  integer, parameter :: MXRNASTANGLE = (MXMP/3) !< maximum # of stack-angles in RNA
  integer, parameter :: MXCHARGE = MXMP    !< maximum # of chage particle
  integer, parameter :: MXCHARGECHANGE = 100 !< maximum # number of CHARGE_CHANGE lines in input file
  integer, parameter :: MXMPELE = 1397    !< maximum # of electrostatic interactions per residue
!  integer, parameter :: MXMPELE = 1000     !< maximum # of electrostatic interactions per residue
  integer, parameter :: MXMPNEIGHBOR = 500 !< maximum # of neighboring residues per residue
  !integer, parameter :: MXELE = MXMPELE * MXMP
  !integer, parameter :: MXNEIGHBOR = MXMPNEIGHBOR * MXMP

  integer, parameter :: MXMPHBNEIGHBOR = 50

!  ! MPC
!  integer, parameter :: MXGRID_NX_MPC = 50 ! the max grid number x-axis
!  integer, parameter :: MXGRID_NY_MPC = MXGRID_NX_MPC ! the max grid number y-axis
!  integer, parameter :: MXGRID_NZ_MPC = MXGRID_NX_MPC ! the max grid number z-axis
!
!  integer, parameter :: MXGRID_N_MPC = MXGRID_NX_MPC * MXGRID_NY_MPC * MXGRID_NZ_MPC   !the max grid number z-axis
!  integer, parameter :: MX_AV_SOLV_GRID_MPC = 20 ! the max average solvent particle in one grid
!  integer, parameter :: MXSOLV_MPC =MX_AV_SOLV_GRID_MPC * MXGRID_N_MPC ! the max number of solovent particle

  ! Replica exchange method
  integer, parameter :: MXREPLICA = 1     !< maximum # of replicas
  integer, parameter :: MXREPDIM = 1       !< maximum # of dimensions in REM

  ! Optional interaction
  integer, parameter :: MXDEL_LGO = 20 !< maximum # of del-interaction residue groups
  integer, parameter :: MXDEL_GO = 100     !< maximum # of del-interaction residue groups
  integer, parameter :: MXBRIDGE = 20     !< maximum # of bridge residue pairs
  integer, parameter :: MXPULLING = 5000    !< maximum # of pulling residues
  integer, parameter :: MXANCHOR = 600    !< maximum # of anchored residues
  integer, parameter :: MXREST1D = 20     !< maximum # of 1D-restrained residues
  integer, parameter :: MXFIX = MXMP      !< maximum # of fixed residues
!  integer, parameter :: MXFLEXIBLE = 20   !< maximum # of the region calculated by flexible local potential
  integer, parameter :: MXGRP = 10
  integer, parameter :: MXMPGRP = 1000

!  ! Multiple Go model
!  integer, parameter :: MXSYSTEM_MGO = 20 !< maximum # of systems in multiple-Go simulation
  integer, parameter :: MXSTATE_MGO = 3   !< maximum # of states in multiple-Go simulation
!  integer, parameter :: MXACT_MGO = 6000  !< maximum # of interaction types
  !integer, parameter :: MXLCON_MGO = 2 * MXCON / MXACT_MGO !< maximum # of contacts belonging to the interaction type

!  ! Implicit ligand 
!  integer, parameter :: MXSITE_IMPLIG = 5  !< maximum # of ligand-binding sites
!  integer, parameter :: MXAA_IMPLIG = 30   !< maximum # of residues per binding site 
!  integer, parameter :: MXUNIT_IMPLIG = 5  !< maximum # of interactive units per binding site
!  integer, parameter :: MXCON_SITE_IMPLIG = 100 !< maximum # of ligand-mediated contacts per binding site
!  integer, parameter :: MXCON_IMPLIG = MXSITE_IMPLIG * MXCON_SITE_IMPLIG !< maximum # of the total of ligand-mediated contacts

  ! Hydrophobic interactions
!  integer, parameter :: MXHP = MXMP    !< maximum # of hydrophobic residues
!  integer, parameter :: MXMPHP = 500   !< maximum # of hydrophobic neighbors per residue

  ! Test particles used in Widom method
  integer, parameter :: MXTP = 3

  ! character length
  integer, parameter :: CARRAY_MSG_ERROR = 512  !< maximum # of characters in error message
  integer, parameter :: CARRAY_MXFILE = 256  !< maximum # of characters in input filename
  integer, parameter :: CARRAY_MXLINE = 1000 !< maximum # of lines in one input-block
  integer, parameter :: CARRAY_MXCOLM = 256 !< maximum # of characters in one line
  integer, parameter :: CARRAY_MXEQUA = 15 !< maximum # of equations in one line

end module const_maxsize
