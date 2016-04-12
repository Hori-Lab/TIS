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
subroutine simu_grid_devision_mpc(irep)

  use const_maxsize, only : PREC, MXSOLV_MPC, MXMP
  use var_setp, only : irand
  use var_struct, only : nmp_real, xyz_mp_rep
  use var_mpc, only: inmpc, base_vector_mpc, xyz_solv_mpc, &
       pbox_origin_mpc, pbox_size_mpc, ngrid_mpc, grid_size_mpc,&
       isolv2grid_mpc, imp2grid_mpc!,xyz_mp_reduced_mpc,xyz_solv_reduced_mpc
  implicit none

  integer, intent(in) :: irep

  !!real(PREC), save :: xyz_mp_reduced_mpc(3, MXMP) !! the reduced coordinate of system particle  in mpc
  ! --------------------------------------------------------------------
  ! function
  real(PREC) :: recipe_rand01
  ! --------------------------------------------------------------------

  !-----------------
  !local variables
  integer :: idimn, is, imp
  integer :: n_solv
  real(PREC) :: delta_xyz(3),delta_xyz_or(3)
  integer :: n_delta(3),i_grid(3)
  !integer :: n_delta_or(3)
  real(PREC) :: vector_grid_shift(3)
  
  !!real(PREC), save :: base_vector_mpc(3,3) !!base_vector for grid
  !! base_vector(1,1:3) !!base vector for grid-1(x) direction  
  !! base_vector(2,1:3) !!base vector for grid-2(y) direction  
  !! base_vector(3,1:3) !!base vector for grid-3(y) direction  
  !!real(PREC), save :: xyz_solv_mpc(3, MXSOLV_MPC) ! coordinate of solvent perticle (MPC)
  !!real(PREC), save :: pbox_origin_mpc(3)  !!periodic boudary box coordinate of origin (x,y,z)
  !!real(PREC), save :: pbox_size_mpc(3)  !!periodic boudary box size (x,y,z)
  !!real(PREC), save :: grid_size_mpc(3)  !!grid size (x,y,z)
  !!real(PREC), save :: ngrid_mpc(3)      !!number of grid (x,y,z)
  !!integer, save :: isolv2grid_mpc(MXSOLV_MPC) !! mpc solvent particles devied into grid
  !!integer, save :: imp2grid_mpc(MXMP) !! mpc solvent particles devied into grid

  !! make randam vecotor for grid shift [-a/2: a/2] , where "a" is grid size.
  do idimn=1, 3
     !!vector_grid_shift(idimn) = pbox_size_mpc(idimn)/(2.0e0_PREC)*(2.0 * recipe_rand01(irand) - 1.0)
     vector_grid_shift(idimn) = grid_size_mpc(idimn)/(2.0e0_PREC)*(2.0 * recipe_rand01(irand) - 1.0)
     !!vector_grid_shift(idimn) = 0.0e0_PREC
  end do

  !! for mpc_solvent
  !! solvent particle of mpc is devided into grid (cell).
  n_solv=inmpc%n_all_solv
  do is=1, n_solv
     do idimn=1, 3
        delta_xyz_or(idimn)=xyz_solv_mpc(idimn, is)-pbox_origin_mpc(idimn)
        !!delta_xyz(idimn)=xyz_solv_mpc(idimn, is)-pbox_origin_mpc(idimn)+vector_grid_shift(idimn)
        delta_xyz(idimn)=delta_xyz_or(idimn)+vector_grid_shift(idimn)
        
        n_delta(idimn)=delta_xyz(idimn)/pbox_size_mpc(idimn)
        !!n_delta_or(idimn)=delta_xyz_or(idimn)/pbox_size_mpc(idimn)
        
        delta_xyz(idimn)=delta_xyz(idimn)-n_delta(idimn)*pbox_size_mpc(idimn)
        !!delta_xyz_or(idimn)=delta_xyz_or(idimn)-n_delta_or(idimn)*pbox_size_mpc(idimn)
        
        !! periodic boundary boxへ縮約する為のベクトル
        if(delta_xyz(idimn)< 0.0d0 ) then
           delta_xyz(idimn)=delta_xyz(idimn)+pbox_size_mpc(idimn)
        end if
        !! periodic boundary boxへ縮約する為のベクトル
        !!if(delta_xyz_or(idimn)< 0.0d0 ) then
        !!   delta_xyz_or(idimn)=delta_xyz_or(idimn)+pbox_size_mpc(idimn)
        !!end if
        
        !! grid_i_x,grid_i_y,grid_i_zをgetする
        i_grid(idimn)=delta_xyz(idimn)/grid_size_mpc(idimn)
        !! solvent particleのpositionをperiodic boudary boxにreduct (reduce)する
        !!xyz_solv_mpc(idimn, is)=pbox_origin_mpc(idimn)+delta_xyz_or(idimn)
!        xyz_solv_reduced_mpc(idimn, is)=pbox_origin_mpc(idimn)+delta_xyz_or(idimn)
     end do
     !! grid numberの確定、格納
     isolv2grid_mpc(is)= ngrid_mpc(1)*ngrid_mpc(2)*i_grid(3) + ngrid_mpc(1)*i_grid(2) + i_grid(1)
  end do

  
  !! for system particle
  !! system particle is devided into grid(cell).
  do imp=1, nmp_real
     do idimn=1, 3
        delta_xyz_or(idimn) = xyz_mp_rep(idimn, imp, irep) - pbox_origin_mpc(idimn)
        delta_xyz(idimn) = delta_xyz_or(idimn) + vector_grid_shift(idimn)
        !!delta_xyz(idimn)=xyz_mp(idimn, imp)-pbox_origin_mpc(idimn) + vector_grid_shift(idimn)
        
        !!n_delta_or(idimn)=delta_xyz_or(idimn)/pbox_size_mpc(idimn)
        n_delta(idimn)=delta_xyz(idimn)/pbox_size_mpc(idimn)
        
        !!delta_xyz_or(idimn)=delta_xyz_or(idimn)-n_delta_or(idimn)*pbox_size_mpc(idimn)
        delta_xyz(idimn)=delta_xyz(idimn)-n_delta(idimn)*pbox_size_mpc(idimn)
        !! periodic boundary boxへ縮約する為のベクトル


        !!if(delta_xyz_or(idimn)< 0.0d0 ) then
        !!   delta_xyz_or(idimn)=delta_xyz_or(idimn)+pbox_size_mpc(idimn)
        !!end if
        if(delta_xyz(idimn)< 0.0d0 ) then
           delta_xyz(idimn)=delta_xyz(idimn)+pbox_size_mpc(idimn)
        end if
        
        !! grid_i_x,grid_i_y,grid_i_zをgetする
        i_grid(idimn)=delta_xyz(idimn)/grid_size_mpc(idimn) !! i_grid(ii) :[0:ngrid_mpc(ii))
        !proteinは時間発展の為に縮約ベクトルを返す必要はない。が、角運動りょうをzero setする際に必要
!        xyz_mp_reduced_mpc(idimn, imp)=pbox_origin_mpc(idimn)+delta_xyz_or(idimn)
        !! subroutine simu_velo_adjst_settemp_mpc(velo_mp, tempk)で使用する
     end do
     !! grid numberの確定、格納
     imp2grid_mpc(imp)= ngrid_mpc(1)*ngrid_mpc(2)*i_grid(3) + ngrid_mpc(1)*i_grid(2) + i_grid(1)
  end do
  
end subroutine simu_grid_devision_mpc
