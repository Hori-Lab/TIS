! var_enm
!> @brief Modeul for defining variables for elastic network model

! **************************************************************************
!  variable for elastic network
module var_enm

  use const_maxsize
  implicit none

  type input_enmparameter
     real(PREC) :: cenm !< constant coefficient of the energy function for elastic network model
     real(PREC) :: dfcontact_enm !< cutoff distance to define the native contact for elastic network model
     integer    :: i_enm !< flag variable for elastic network model
     integer    :: sz !< size of the structure
  end type input_enmparameter

  type(input_enmparameter), save :: inenm

end module var_enm
