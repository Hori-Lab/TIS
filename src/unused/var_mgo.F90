! var_mgo
!> @brief Modele for defining variables for multiple Go model

! ********************************************************************
! Explanation of arrays
!
! (mgo is omitted in the following descriunition.)
!
! nsystem: number of Multi-Go system. A system has it's own
!          multiple Go-potential, works independently,
!          and interact with other system through non-native interaction.
! nstate(nsystem): number of state of a system. For example,
!                  beta-subunit of F1 has three nucleotide
!                  states: ATP, ADP, and emunity.
! nactnum(nsystem): number of interaction of a system in a state.
! isysmbr(nsystem, nstate, nactnum): members of interaction of a system in a state.
!
! ishadow2real_unit_mgo(nunit): unit number of real unit corresponding
!                               to shadow unit
! iactmat(nunit, nunit): interaction number between iunit and junit
! iact2unit(2, nact): unit numbers involved in a interaction
!
! nact2concalc(nact): number of contact of a interaction,            
!                            and these belong to neighbor list.
! iact2con(nact2con(nact), nact): list of contact number,
!                                  which is in a interaction
! ********************************************************************
module var_mgo

  use const_maxsize

  implicit none

  type input_mgoparameter
    integer :: sz
    integer :: i_multi_mgo
    integer :: nsystem_mgo
    integer :: nstate_mgo(MXSYSTEM_MGO)
    integer :: nstate_max_mgo
    integer :: nactnum_mgo(MXSYSTEM_MGO)
    integer :: nact_mgo
    integer :: isysmbr_mgo(MXSYSTEM_MGO, MXSTATE_MGO, MXACT_MGO)
    integer :: iactmat_mgo(MXUNIT, MXUNIT)
    real(PREC) :: bdemax_mgo
    real(PREC) :: baemax_mgo
    real(PREC) :: dihemax_mgo
    real(PREC) :: enegap(MXSIM, MXSYSTEM_MGO, MXSTATE_MGO)
    real(PREC) :: delta_mgo(MXSYSTEM_MGO, MXSTATE_MGO, MXSTATE_MGO)
  end type input_mgoparameter

  type(input_mgoparameter), save :: inmgo

  integer, save :: iunit2sysmbr_mgo(3, MXUNIT, MXUNIT)
  integer, save :: iact2unit_mgo(2, MXACT_MGO)
  integer, save :: ishadow2real_unit_mgo(MXUNIT)
  integer, save :: ishadow2real_mp_mgo(MXMP)

  integer, allocatable, save :: ibd2sysmbr_mgo(:,:)  !(2, MXBD)
  integer, allocatable, save :: iba2sysmbr_mgo(:,:)  !(2, MXBA)
  integer, allocatable, save :: idih2sysmbr_mgo(:,:) !(2, MXDIH)

  real(PREC), allocatable, save :: enegap_mgo(:,:)  !(MXSYSTEM_MGO, MXSTATE_MGO)
  real(PREC), allocatable, save :: offset_mgo(:,:)  !(MXSYSTEM_MGO, MXSTATE_MGO)
  real(PREC), allocatable, save :: offset_unit(:,:) !(MXUNIT, MXUNIT)
  integer,    allocatable, save :: ncontype_mgo(:)      !(MXCON)
  integer,    allocatable, save :: irefcon_mgo(:)       !(MXCON)
  integer,    allocatable, save :: icon2sysmbr_mgo(:,:) !(2, MXCON)
  !integer, save :: ireal2shadow_con_mgo(MXSTATE_MGO, MXCON)

  real(PREC), allocatable, save :: coef_mgo(:,:)    !(MXSYSTEM_MGO, MXSTATE_MGO)
  real(PREC), allocatable, save :: esystem_mgo(:)   !(MXSYSTEM_MGO)
  real(PREC), allocatable, save :: estate_mgo(:,:)  !(MXSYSTEM_MGO, MXSTATE_MGO)
  real(PREC), allocatable, save :: q_mgo(:,:)       !(MXSYSTEM_MGO, MXSTATE_MGO)
  real(PREC), allocatable, save :: ekai_mgo(:)      !(MXSYSTEM_MGO)

end module var_mgo
