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
subroutine simu_grid_devision_mpc3(irep)

  use const_maxsize, only : PREC, MXSOLV_MPC, MXMP
  use var_setp, only : irand
  use var_struct, only : nmp_real, xyz_mp_rep
  use var_mpc, only: inmpc, xyz_solv_mpc, &
       pbox_origin_mpc, pbox_size_mpc, ngrid_mpc, grid_size_mpc, &
       igrid_process, lgrid_l, igrid_l, &
       nis_l, nimp_l, isolv2grid_mpc_l, imp2grid_mpc_l

  use mpiconst

  implicit none

  integer, intent(in) :: irep

  ! --------------------------------------------------------------------
  ! function
  real(PREC) :: recipe_rand01

  ! --------------------------------------------------------------------
  !local variables
  integer :: idimn, is, imp, ig, ig_p, ig_l
  integer :: n_solv
  integer :: n_delta(3), i_grid(3)
  real(PREC) :: delta_xyz(3)
  real(PREC) :: vector_grid_shift(3)
  real(PREC) :: rpbox_size(3), rgrid_size(3)
  
  ! -----------------------------------------------------------------
  n_solv = inmpc%n_all_solv

  ! -----------------------------------------------------------------
  !! make randam vecotor for grid shift [-a/2: a/2] , where "a" is grid size.
  do idimn = 1, 3
     vector_grid_shift(idimn) = grid_size_mpc(idimn)*(recipe_rand01(irand) - 0.5) &
          - pbox_origin_mpc(idimn)
     rpbox_size(idimn) = 1.0/pbox_size_mpc(idimn)
     rgrid_size(idimn) = 1.0/grid_size_mpc(idimn)
  end do

  ! -----------------------------------------------------------------
  !! for mpc_solvent
  !! solvent particle of mpc is devided into grid (cell).
  nis_l = 0
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
!     isolv2grid_mpc(is) = ngrid_mpc(1)*ngrid_mpc(2)*i_grid(3) &
!          + ngrid_mpc(1)*i_grid(2) + i_grid(1)
     ig_p = igrid_process(i_grid(1), i_grid(2), i_grid(3))
     if(ig_p == myrank) then
        nis_l = nis_l + 1
        ig_l = igrid_l(i_grid(1), i_grid(2), i_grid(3))
!        isolv2grid_mpc_l(1, nis_l) = ig_l
        ig = ngrid_mpc(1)*ngrid_mpc(2)*i_grid(3) &
             + ngrid_mpc(1)*i_grid(2) + i_grid(1)
        isolv2grid_mpc_l(1, nis_l) = ig
        isolv2grid_mpc_l(2, nis_l) = is
     end if
  end do

  
  !! for system particle
  !! system particle is devided into grid(cell).
  nimp_l = 0
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
!     imp2grid_mpc(imp) = ngrid_mpc(1)*ngrid_mpc(2)*i_grid(3) &
!          + ngrid_mpc(1)*i_grid(2) + i_grid(1)

     ig_p = igrid_process(i_grid(1), i_grid(2), i_grid(3))
     if(ig_p == myrank) then
        nimp_l = nimp_l + 1
        ig_l = igrid_l(i_grid(1), i_grid(2), i_grid(3))
!        imp2grid_mpc_l(1, nimp_l) = ig_l
        ig = ngrid_mpc(1)*ngrid_mpc(2)*i_grid(3) &
             + ngrid_mpc(1)*i_grid(2) + i_grid(1)
        imp2grid_mpc_l(1, nimp_l) = ig
        imp2grid_mpc_l(2, nimp_l) = imp
     end if
  end do

end subroutine simu_grid_devision_mpc3
