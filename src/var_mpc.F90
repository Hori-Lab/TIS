module var_mpc

  use const_maxsize, only : PREC, MXSOLV_MPC,MXMP
  !!  use maxsize, only : PREC, MXSIM, MXUNIT, MXMP, MXBD, MXBA, MXDIH, MXCON, &
  !!       MXSYSTEM_MGO, MXSTATE_MGO, MXACT_MGO
  implicit none
  
  type input_mpcparameter
     integer :: sz
     integer :: i_thermal_mpc !! thermal bath
     integer :: nratio_colli_step  !! the number of ratio_collision_step
     integer :: nratio_vcorrect_step !! step (ratio) number to correct [velocity, temperature]
     
     integer :: n_av_solv_ingrid           !! average solvent number in one grid 
     !! this value must be integer.
     integer :: n_all_solv !! total solvent number in box
     real(PREC) :: rotate_angle_colli !! rotation angle for collision step
     
     integer :: ngrid_x  !! number of grid(cell) x-axis
     integer :: ngrid_y  !! number of grid(cell) y-axis
     integer :: ngrid_z  !! number of grid(cell) z-axis
     
     real(PREC) :: pbox_size_x  !! periodic boundary box size x-axis 
     real(PREC) :: pbox_size_y  !! periodic boundary box size y-axis 
     real(PREC) :: pbox_size_z  !! periodic boundary box size z-axis 
     
     real(PREC) :: pbox_origin_x  !! periodic boundary box origin x-axis 
     real(PREC) :: pbox_origin_y  !! periodic boundary box origin y-axis 
     real(PREC) :: pbox_origin_z  !! periodic boundary box origin z-axis 
     
     real(PREC) :: grid_size_x  !! grid size x-axis 
     real(PREC) :: grid_size_y  !! grid size y-axis 
     real(PREC) :: grid_size_z  !! grid size z-axis 
     
     real(PREC) :: cmass_solv   !! mass of mpc_solvent


     real(PREC) :: nratio_correct_step !! step (ratio) number to correct [velocity, temperature]
     integer :: i_flag_check_mpc, i_flag_check2_mpc !!check flag for mpc

  end type input_mpcparameter
  
  type(input_mpcparameter), save :: inmpc
  
  real(PREC), save :: base_vector_mpc(3,3) !!base_vector for grid
  !! base_vector(1,1:3) !!base vector for grid-1(x) direction  
  !! base_vector(2,1:3) !!base vector for grid-2(y) direction  
  !! base_vector(3,1:3) !!base vector for grid-3(y) direction  
  
  real(PREC), save :: xyz_solv_mpc(3, MXSOLV_MPC) ! coordinate of solvent perticle (MPC)
  real(PREC), save :: velo_solv_mpc(3, MXSOLV_MPC) ! velocity of solvent perticle (MPC)
  
  real(PREC), save :: cmass_solv_mpc(MXSOLV_MPC) ! velocity of solvent perticle (MPC)

  real(PREC), save :: pbox_origin_mpc(3)  !!periodic boudary box coordinate of origin (x,y,z)
  
  real(PREC), save :: pbox_size_mpc(3)  !!periodic boudary box size (x,y,z)
  real(PREC), save :: grid_size_mpc(3)  !!grid size (x,y,z)
!  integer, save :: ngrid_mpc(3)      !!number of grid (x,y,z)
  real(PREC), save :: ngrid_mpc(3)      !!number of grid (x,y,z)
  
  
  integer, save :: isolv2grid_mpc(MXSOLV_MPC) !! mpc solvent particles devied into grid
  integer, save :: imp2grid_mpc(MXMP) !! system particle devied into grid
  
  integer, allocatable, save :: igrid_process(:,:,:)
  integer, allocatable, save :: lgrid_l(:)
  integer, allocatable, save :: igrid_l(:,:,:)

  integer, save :: nis_l, nimp_l
  integer, save :: isolv2grid_mpc_l(2, MXSOLV_MPC) !! mpc solvent particles devied into grid
  integer, save :: imp2grid_mpc_l(2, MXMP) !! system particle devied into grid
  
!  real(PREC), save :: xyz_mp_reduced_mpc(3, MXMP) !! the reduced coordinate of system particle in mpc 
!  real(PREC), save :: xyz_solv_reduced_mpc(3, MXSOLV_MPC) !! the reduced coordinate of solvent particle in mpc 
end module var_mpc
