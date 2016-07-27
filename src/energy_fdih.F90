!energy_fdih
!> @brief Calculates the energy related to flexible dihedral angle tmer.
!>        Values are added into "energy(E_TYPE%DIHE)" and &
!>        "energy_unit(,,E_TYPE%DIHE)".

subroutine energy_fdih(irep, energy_unit, energy)
  
  use const_maxsize
  use const_index
  use var_struct, only : nfdih, ifdih2mp, fdih_para, iclass_mp, imp2unit, nunit_all, fdih_ener_corr
  use var_setp,   only : inflp
  use var_flp
  use mpiconst

  implicit none

  integer, intent(in) :: irep
  real(PREC), intent(inout) :: energy(E_TYPE%MAX) 
  real(PREC), intent(inout) :: energy_unit(nunit_all, nunit_all, E_TYPE%MAX)
  
  type(phi_related_variables) :: phi_related
  integer :: imp(4)
  integer :: ksta, kend
  integer :: ifdih, i
  real(PREC) :: efull
  integer :: iunit, junit
#ifdef MPI_PAR3
  integer :: klen
#endif

#ifdef MPI_PAR3
  klen=(nfdih-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nfdih)
#else
  ksta = 1
  kend = nfdih
#endif
!$omp do private(i, imp, phi_related, iunit, junit, efull)
  do ifdih = ksta, kend

     do i = 1, 4
        imp(i) = ifdih2mp(i, ifdih)
     end do

     if (iclass_mp( imp(1) ) == CLASS%LIG .AND. &
         iclass_mp( imp(4) ) == CLASS%LIG) cycle

     !call util_dihangle( imp(1), imp(2), imp(3), imp(4), &
     !     th, co_dih, si_dih, xyz_mp_rep(:,:,irep) );
     
     call calc_phi(phi_related, irep, imp)

     !call calc_theta(theta1, irep, imp(1), imp(2), imp(3))
     !call calc_theta(theta2, irep, imp(2), imp(3), imp(4))

     !call calc_u(u1, theta1)
     !call calc_u(u2, theta2)

     efull = fdih_para(1, ifdih)                            &
           + fdih_para(2, ifdih) * cos(    phi_related%phi) &
           + fdih_para(3, ifdih) * sin(    phi_related%phi) &
           + fdih_para(4, ifdih) * cos(2 * phi_related%phi) &
           + fdih_para(5, ifdih) * sin(2 * phi_related%phi) &
           + fdih_para(6, ifdih) * cos(3 * phi_related%phi) &
           + fdih_para(7, ifdih) * sin(3 * phi_related%phi) &
           - fdih_ener_corr(ifdih)

     efull = inflp%k_dih * efull

     !------------------------------------------------------------------------
     ! sum of the energy
     iunit = imp2unit( imp(1) )
     junit = imp2unit( imp(2) )

     energy(E_TYPE%DIHE) = energy(E_TYPE%DIHE) + efull
     energy_unit(iunit, junit, E_TYPE%DIHE) = energy_unit(iunit, junit, E_TYPE%DIHE) + efull
     
  end do
!$omp end do nowait

contains
  ! ----------------------------------------------------------------------
  ! Subroutine for caluculation of phi
  ! ----------------------------------------------------------------------
!  subroutine calc_theta(theta, irep, imp1, imp2, imp3)
!    use const_maxsize
!    use var_struct, only : xyz_mp_rep
!    !use var_flp
!    
!    implicit none
!    
!    ! ----------------------------------------------------------------------
!    ! Arguments
!    real(PREC), intent(inout) :: theta
!    integer, intent(in) :: irep, imp1, imp2, imp3
!    
!    ! ----------------------------------------------------------------------
!    ! Local variables
!    real(PREC) :: v21(3), v32(3)
!    real(PREC) :: c11, c22, c21
!    real(PREC) :: co_theta
!    
!    ! ----------------------------------------------------------------------
!    ! Calculate bond angle
!    
!    v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
!    v32(1:3) = xyz_mp_rep(1:3, imp3, irep) - xyz_mp_rep(1:3, imp2, irep)
!    
!    c11 = v21(1) * v21(1) + v21(2) * v21(2) + v21(3) * v21(3)
!    c22 = v32(1) * v32(1) + v32(2) * v32(2) + v32(3) * v32(3)
!    c21 = v32(1) * v21(1) + v32(2) * v21(2) + v32(3) * v21(3)
!    
!    co_theta = - c21 / sqrt(c11 * c22)
!    
!    if(co_theta > 1.0e0_PREC) then
!       co_theta = 1.0e0_PREC
!    else if(co_theta < -1.0e0_PREC) then
!       co_theta = -1.0e0_PREC
!    end if
!    
!    theta = acos(co_theta)
!    
!  end subroutine calc_theta

!  subroutine calc_u(u, theta)
!    use const_maxsize
!    use const_physical
!
!    implicit none
!
!    ! ----------------------------------------------------------------------
!    ! Arguments
!    real(PREC), intent(inout) :: u
!    real(PREC), intent(in) :: theta
!
!    ! ----------------------------------------------------------------------
!    ! Local variables
!    real(PREC) :: a, b
!    
!    ! ----------------------------------------------------------------------
!    ! Calculate bond angle
!    a = 1.0e0_PREC / (FDIH_DEL_MAX_ANG - FDIH_DEL_MIN_ANG)
!    b = a * FDIH_DEL_MAX_ANG
!
!    if (theta <= FDIH_DEL_MIN_ANG) then
!       u = 1.0e0_PREC
!    else if ( theta > FDIH_DEL_MIN_ANG .and. theta <= FDIH_DEL_MAX_ANG ) then
!       u = -a * theta + b
!       write (*, *) "Warning: u: ", u
!    else if ( theta > FDIH_DEL_MAX_ANG ) then
!       u = 0.0
!       write (*, *) "Warning: u: ", u
!    end if
!       
!  end subroutine calc_u

  subroutine calc_phi(phi_related, irep, imp)
    use const_maxsize
    use var_struct, only : xyz_mp_rep
    use var_flp
    
    implicit none
    
    ! ----------------------------------------------------------------------
    ! Arguments
    type(phi_related_variables), intent(inout) :: phi_related
    integer, intent(in) :: irep, imp(4)
    
    ! ----------------------------------------------------------------------
    ! Local variables
    real(PREC) :: cos_phi
    real(PREC) :: vmvn, vmvm, vnvn
    real(PREC) :: sign
    
    ! ----------------------------------------------------------------------
    ! Calculate dihedral angle
    
    ! Calculate internal vector
    phi_related%vij(1:3) = xyz_mp_rep(1:3, imp(1), irep) - xyz_mp_rep(1:3, imp(2), irep)
    phi_related%vkj(1:3) = xyz_mp_rep(1:3, imp(3), irep) - xyz_mp_rep(1:3, imp(2), irep)
    phi_related%vkl(1:3) = xyz_mp_rep(1:3, imp(3), irep) - xyz_mp_rep(1:3, imp(4), irep)
    
    ! Calculate normal vectors
    call calc_cross_product(phi_related%vm, phi_related%vij, phi_related%vkj)
    call calc_cross_product(phi_related%vn, phi_related%vkj, phi_related%vkl)
    
    ! Calculate dot product of vm and vn
    call calc_dot_product(vmvn, phi_related%vm, phi_related%vn)
    call calc_dot_product(vmvm, phi_related%vm, phi_related%vm)
    call calc_dot_product(vnvn, phi_related%vn, phi_related%vn)
    
    ! Calculate size of normal vectors
    phi_related%rm = vmvm
    phi_related%rn = vnvn
    
    ! Calculate cos(phi)
    cos_phi = vmvn / ( sqrt(phi_related%rm) * sqrt(phi_related%rn) )
    
    if(cos_phi > 1.0e0_PREC) then
       cos_phi = 1.0e0_PREC
    else if(cos_phi < -1.0e0_PREC) then
       cos_phi = -1.0e0_PREC
    end if
    
    phi_related%phi = acos(cos_phi);
    
    call calc_dot_product (sign, phi_related%vij, phi_related%vn)
    
    if (sign < 0.0e0_PREC) then
       phi_related%phi = phi_related%phi * (-1.0e0_PREC)
    end if
    
  end subroutine calc_phi

  ! ----------------------------------------------------------------------
  ! Subroutine for caluculation of cross product
  ! ----------------------------------------------------------------------
  subroutine calc_cross_product(cross_product, vector1, vector2)
    
    use const_maxsize
    
    implicit none
    
    ! ----------------------------------------------------------------------
    ! Arguments
    real(PREC), intent(inout) :: cross_product(3)
    real(PREC), intent(in) :: vector1(3), vector2(3)
    
    cross_product(1) = vector1(2) * vector2(3) - vector1(3) * vector2(2);
    cross_product(2) = vector1(3) * vector2(1) - vector1(1) * vector2(3);
    cross_product(3) = vector1(1) * vector2(2) - vector1(2) * vector2(1);
    
  end subroutine calc_cross_product

  ! ----------------------------------------------------------------------
  ! Subroutine for caluculation of dot product
  ! ----------------------------------------------------------------------
  subroutine calc_dot_product(dot_product, vector1, vector2)
    
    use const_maxsize
    
    implicit none
    
    ! ----------------------------------------------------------------------
    ! Arguments
    real(PREC), intent(inout) :: dot_product
    real(PREC), intent(in) :: vector1(3)
    real(PREC), intent(in) :: vector2(3)
    
    ! ----------------------------------------------------------------------
    ! Local variables
    integer :: i
    
    dot_product = 0.0E0_PREC
    
    do i = 1, 3
       dot_product = dot_product + ( vector1(i) * vector2(i) )
    end do
    
  end subroutine calc_dot_product
  
end subroutine energy_fdih
