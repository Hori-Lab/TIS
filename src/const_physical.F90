!const_physical
!> @brief Physical values and numerical limits are defined here.

module const_physical
  use const_maxsize

  integer,    parameter :: SDIM = 3 !< # of space dimension

  real(PREC), parameter :: F_PI    = 3.14159265358979323846264338e0_PREC !< Circular constant (Pi)
  real(PREC), parameter :: F_2PI   = 2.0 * F_PI

  real(PREC), parameter :: EPSI_0  = 8.854187817e-12_PREC  !< Vacuum permittivity [F/m]
  real(PREC), parameter :: ELE     = 1.6021766208e-19_PREC !< Elementary charge [C]
  real(PREC), parameter :: BOLTZ_J = 1.38064852e-23_PREC   !< Boltzmann constant [J/K]
  real(PREC), parameter :: N_AVO   = 6.022140857e23_PREC   !< Avogadro constant [/mol]
  real(PREC), parameter :: KCAL2JOUL = 4184.0              !< (kcal -> J)  [J/kcal]

  real(PREC), parameter :: JOUL2KCAL = 1.0/KCAL2JOUL   !< (J -> kcal)  [kcal/J]
  real(PREC), parameter :: JOUL2KCAL_MOL  = JOUL2KCAL * N_AVO  !< (J -> kcal/mol)
  real(PREC), parameter :: BOLTZ_KCAL_MOL = BOLTZ_J * JOUL2KCAL_MOL   !< Boltzmann constant [kcal/mol/K]
                                        ! = 0.00198720359
!  real(PREC), parameter :: BOLTZ_KCAL_MOL_ND = 0.0019858775  !< BOLTZ used in Denesyuk

  real(PREC), parameter :: DE_MAX  = 300.0e0_PREC !< limit value of force

  ! judgment for numerical error
  real(PREC), parameter :: INVALID_JUDGE = 1.0e30_PREC
  real(PREC), parameter :: INVALID_VALUE = 1.0e31_PREC
  real(PREC), parameter :: ZERO_JUDGE    = 1.0e-6_PREC
  real(PREC), parameter :: CUTOFF_UNDER_EXP = -50.0_PREC
  real(PREC), parameter :: HIGH_ENERGY = 1.0e10_PREC
  real(PREC), parameter :: HIGH_ENERGY_JUDGE = 1.0e9_PREC

  ! judgement for warning
  ! Output warning if bond angle is larger than WARN_ANGLE.
  real(PREC), parameter :: WARN_ANGLE = F_PI - 0.17e0_PREC  !< Minimum angle [rad]
  ! Output warning if bond length is longer than WARN_BOND.
  real(PREC), parameter :: WARN_BOND        = 5.0e0_PREC !< Maximum bond length [angst.]
  real(PREC), parameter :: WARN_BOND_RNA_SP = 6.0e0_PREC !< Maximum length for RNA S-P bond. [angst.]
  real(PREC), parameter :: WARN_BOND_RNA_SB = 8.0e0_PREC !< Maximum length for RNA S-B bond. [angst.]
  real(PREC), parameter :: WARN_BOND_RNA_PS = 6.0e0_PREC !< Maximum length for RNA P-S bond. [angst.]
  real(PREC), parameter :: WARN_RNA_O3_P = 1.8 !< Maximum atomic distance b/w O3 and P in RNA [angst.]
  real(PREC), parameter :: WARN_RNA_P_O5 = 1.8 !< Maximum atomic distance b/w P and O5 in RNA [angst.]

  ! Mass
  real(PREC), parameter :: MASS_P = 30.973761e0_PREC  !< Phosphorus
  real(PREC), parameter :: MASS_O = 15.9994e0_PREC    !< Oxygen
  real(PREC), parameter :: MASS_C = 12.0107e0_PREC    !< Carbon
  real(PREC), parameter :: MASS_N = 14.0067e0_PREC    !< Nitrogen
  real(PREC), parameter :: MASS_BR= 79.904e0_PREC     !< Boron
  real(PREC), parameter :: MASS_F = 18.9984032e0_PREC !< Fluorine
  real(PREC), parameter :: MASS_S = 32.065e0_PREC     !< Sulfur
  real(PREC), parameter :: MASS_PO4 = MASS_P + MASS_O * 4.0e0_PREC !< Phosphoric acid

  ! Nose-Hoover parameter
  integer,    parameter :: MXCS = 5               !< # of thermal particle
!  real(PREC), parameter :: CSMASS = 100.0e0_PREC  !< Mass of thermal particle
  
  ! O4'+C4'+C3'+C2'+C1'
  real(PREC), parameter :: MASS_RING     = MASS_C * 4 &
                                         + MASS_O * 1

  real(PREC), parameter :: MASS_RIBOSE   = MASS_C * 5 &
                                         + MASS_O * 2

  real(PREC), parameter :: MASS_ADENINE  = MASS_N * 5 &
                                         + MASS_C * 5

  real(PREC), parameter :: MASS_GUANINE  = MASS_N * 5 &
                                         + MASS_C * 5 &
                                         + MASS_O * 1

  real(PREC), parameter :: MASS_URACIL   = MASS_N * 2 &
                                         + MASS_C * 4 &
                                         + MASS_O * 2

  real(PREC), parameter :: MASS_CYTOSINE = MASS_N * 3 &
                                         + MASS_C * 4 &
                                         + MASS_O * 1

  real(PREC), parameter :: MASS_RIN_SUGAR = MASS_P * 1  &
                                          + MASS_C * 10 &
                                          + MASS_O * 9

  real(PREC), parameter :: MASS_SUGAR    = MASS_C * 5 &
                                         + MASS_O * 3

  real(PREC), parameter :: MASS_PHOS     = MASS_P * 1 &
                                         + MASS_O * 2

!  ! Flexible_local_potential
!  real(PREC), parameter :: FBA_MIN_ANG = 1.31e0_PREC
!  real(PREC), parameter :: FBA_MAX_ANG = 2.87e0_PREC
!  real(PREC), parameter :: FBA_MIN_ANG_FORCE = -30.0e0_PREC
!  real(PREC), parameter :: FBA_MAX_ANG_FORCE = 30.0e0_PREC
!  real(PREC), parameter :: FDIH_DEL_MIN_ANG = 2.88e0_PREC ! 2.88 rad = 165 degree
!  real(PREC), parameter :: FDIH_DEL_MAX_ANG = 3.05e0_PREC ! 3.05 rad = 175 degree

endmodule const_physical
