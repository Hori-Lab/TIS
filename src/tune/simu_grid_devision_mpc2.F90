!!**********************************************
!! input 
!! xyz_solv_mpc(3, MXSOLV_MPC) ! coordinate of solvent perticle (MPC)
!! xyz_mp_rep(3, MXMP)
!!
!! output
!!integer, save :: isolv2grid_mpc(MXSOLV_MPC) !! mpc solvent particles devied into grid
!!integer, save :: imp2grid_mpc(MXMP) !! system particles devied into grid
!!**********************************************
!! solvent_mpc, system_particle devied into grid
subroutine simu_grid_devision_mpc2(irep)

  use const_maxsize, only : PREC, MXSOLV_MPC, MXMP
  use var_setp, only : mts
  use var_struct, only : nmp_real, xyz_mp_rep
  use var_mpc, only: inmpc, xyz_solv_mpc, &
       pbox_origin_mpc, pbox_size_mpc,ngrid_mpc,grid_size_mpc,&
       isolv2grid_mpc, imp2grid_mpc
!  use mtmod
  use mt_stream

  implicit none

  integer, intent(in) :: irep

  ! --------------------------------------------------------------------
  !local variables
  integer :: idimn, is, imp, istream
  integer :: n_solv
  integer :: n_delta(3), i_grid(3)
  real(PREC) :: delta_xyz(3)
  real(PREC) :: vector_grid_shift(3)
  real(PREC) :: rpbox_size(3), rgrid_size(3)
  

  ! -----------------------------------------------------------------
  !! make randam vecotor for grid shift [-a/2: a/2] , where "a" is grid size.
  do idimn = 1, 3
!     vector_grid_shift(idimn) = grid_size_mpc(idimn)*(grnd() - 0.5) &
!          - pbox_origin_mpc(idimn)
     istream = irep
     vector_grid_shift(idimn) = grid_size_mpc(idimn)*(genrand_double1(mts(istream, 0)) - 0.5) &
          - pbox_origin_mpc(idimn)
     rpbox_size(idimn) = 1.0/pbox_size_mpc(idimn)
     rgrid_size(idimn) = 1.0/grid_size_mpc(idimn)
  end do

  !! for mpc_solvent
  !! solvent particle of mpc is devided into grid (cell).
  n_solv = inmpc%n_all_solv
  do is = 1, n_solv
     delta_xyz(1:3) = xyz_solv_mpc(1:3, is) + vector_grid_shift(1:3)
        
     !! periodic boundary boxへ縮約する為のベクトル
     n_delta(1:3) = delta_xyz(1:3) * rpbox_size(1:3)
     delta_xyz(1:3) = delta_xyz(1:3) - n_delta(1:3)*pbox_size_mpc(1:3)
     
     do idimn = 1, 3
        if(delta_xyz(idimn) < 0.0) then
           delta_xyz(idimn) = delta_xyz(idimn) + pbox_size_mpc(idimn)
        end if
     end do
     !! grid_i_x,grid_i_y,grid_i_zをgetする
     i_grid(1:3) = delta_xyz(1:3) * rgrid_size(1:3)

     !! grid numberの確定、格納
     isolv2grid_mpc(is) = ngrid_mpc(1)*ngrid_mpc(2)*i_grid(3) &
          + ngrid_mpc(1)*i_grid(2) + i_grid(1)
  end do

  
  !! for system particle
  !! system particle is devided into grid(cell).
  do imp = 1, nmp_real

     delta_xyz(1:3) = xyz_mp_rep(1:3, imp, irep) + vector_grid_shift(1:3)
     
     !! periodic boundary boxへ縮約する為のベクトル
     n_delta(1:3) = delta_xyz(1:3) * rpbox_size(1:3)
     delta_xyz(1:3) = delta_xyz(1:3) - n_delta(1:3)*pbox_size_mpc(1:3)
     
     do idimn = 1, 3
        if(delta_xyz(idimn) < 0.0d0 ) then
           delta_xyz(idimn) = delta_xyz(idimn) + pbox_size_mpc(idimn)
        end if
     end do

     !! grid_i_x,grid_i_y,grid_i_zをgetする
     i_grid(1:3) = delta_xyz(1:3) * rgrid_size(1:3)

     !! grid numberの確定、格納
     imp2grid_mpc(imp) = ngrid_mpc(1)*ngrid_mpc(2)*i_grid(3) &
          + ngrid_mpc(1)*i_grid(2) + i_grid(1)
  end do
  
end subroutine simu_grid_devision_mpc2
