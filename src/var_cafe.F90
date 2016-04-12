!var_cafe
!> @brief Contains coefficients for analytical tools

module var_cafe

  use const_maxsize
  use const_index

  integer, parameter :: MXNUM = 100000

  type input_reference
     integer :: nmp
     real(PREC) :: xyz(3, MXNUM)
     real(8) :: bfactor(MXNUM)
  end type input_reference

  type input_trajectory
     integer :: ntitle ! the line number of title lines 3 + nunit
     integer :: nblock_size ! block-size 4 + 80*ntitle
     integer :: iat ! the number of atom
     integer :: num 
     integer :: nset ! the number of frames
     integer :: istrt ! starting step number
     integer :: nsavc ! step interval
     integer :: nstep ! the number of steps 
     integer :: nver ! version if CHARMm 24
     integer :: nunit ! the number of unit
     integer :: lunit2mp(MXNUM) ! index atom number of unit

     real(4) :: delta ! time-step
     real(PREC) :: xyz(3, MXNUM) ! coordinate of trajectory
     real(8) :: tempk ! temperature

     character(4) :: hdr ! 'CORD' for coordinate, 'VELD' for velocity
     character(80) :: title ! title
  end type input_trajectory
  

end module var_cafe
