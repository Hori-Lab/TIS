!var_flp
!> @brief Contains parameters and variables which are related &
!>        to flexible local potential

module var_flp
  use const_maxsize
  
  type phi_related_variables
     real(PREC) :: phi                    ! Dihedral angle
     real(PREC) :: vm(3), vn(3)           ! Normal vector
     real(PREC) :: rm, rn                 ! Size of normal vector
     real(PREC) :: vij(3), vkj(3), vkl(3) ! Vector between each particle
     real(PREC) :: vfi(3), vfl(3)
     real(PREC) :: p, q
  end type phi_related_variables
  
  !type(phi_related_variables), save :: phi_related
  
end module var_flp
