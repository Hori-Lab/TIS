!!subroutine setp_para_mpc(velo_mp, tempk)
subroutine setp_para_mpc()


  use const_maxsize
  use const_index, only : SIM
  use const_physical, only : BOLTZC
  use var_setp, only : irand, ifix_mp
  use var_struct, only : cmass_mp
  use var_inp, only : infile, outfile, i_simulate_type
  !!use var_mpc, only : inmpc, base_vector_mpc(3,3),xyz_solv_mpc(3, MXSOLV_MPC),&
  !!     cmass_solv_mpc(MXSOLV_MPC),pbox_origin_mpc(3)

!  use var_mpc, only : inmpc, , base_vector_mpc,xyz_solv_mpc,&
!       cmass_solv_mpc, pbox_origin_mpc, velo_solv_mpc, pbox_size_mpc, grid_size_mpc, ngrid_mpc
  use var_mpc, only : inmpc, pbox_size_mpc, grid_size_mpc, ngrid_mpc, &
                      cmass_solv_mpc, pbox_origin_mpc

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none
  ! -------------------------------------------------------------------
  ! function
  !real(PREC) :: recipe_rand01
  
  ! -----------------------------------------------------------------------
  ! local variables
  integer :: i
  !integer :: is, idimn
  integer :: n_solv
  
  real(PREC) :: pbbox_min_x,pbbox_max_x
  real(PREC) :: pbbox_min_y,pbbox_max_y
  real(PREC) :: pbbox_min_z,pbbox_max_z

  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message


  ! -----------------------------------------------------------------------
!  luninp = ifile(1)
  luninp = infile%inp
  lunout = outfile%data


  !! default variavle
  
  inmpc%i_thermal_mpc=1  !! thermal bath
  inmpc%nratio_colli_step=500  !! the number of ratio_collision_step
  inmpc%nratio_correct_step=500   !! step (ratio) number to correct [velocity, temperature]
  inmpc%nratio_vcorrect_step=500   !! step (ratio) number to correct [velocity, temperature]
  inmpc%n_av_solv_ingrid=4           !! average solvent number in one grid 
  !! this value must be integer.

  inmpc%rotate_angle_colli=90.0e0_PREC !! rotation angle for collision step
  inmpc%ngrid_x=32  !! number of grid(cell) x-axis
  inmpc%ngrid_y=32  !! number of grid(cell) y-axis
  inmpc%ngrid_z=32  !! number of grid(cell) z-axis
     
  inmpc%n_all_solv=inmpc%ngrid_x*inmpc%ngrid_y*inmpc%ngrid_z*inmpc%n_av_solv_ingrid 
  !! total solvent number in box
  
  inmpc%pbox_size_x=200.0e0_PREC  !! periodic boundary box size x-axis 
  inmpc%pbox_size_y=200.0e0_PREC  !! periodic boundary box size y-axis 
  inmpc%pbox_size_z=200.0e0_PREC  !! periodic boundary box size z-axis 
     
  inmpc%grid_size_x=inmpc%pbox_size_x/inmpc%ngrid_x  !! grid size x-axis 
  inmpc%grid_size_y=inmpc%pbox_size_y/inmpc%ngrid_y  !! grid size y-axis 
  inmpc%grid_size_z=inmpc%pbox_size_z/inmpc%ngrid_z  !! grid size z-axis 
  
  inmpc%pbox_origin_x=0.0e0_PREC  !! periodic boundary box origin x-axis 
  inmpc%pbox_origin_y=0.0e0_PREC  !! periodic boundary box origin y-axis 
  inmpc%pbox_origin_z=0.0e0_PREC  !! periodic boundary box origin z-axis 

  inmpc%cmass_solv=10.0e0_PREC   !! mass of mpc_solvent
  !! set of solvent mass
  n_solv=inmpc%n_all_solv
  do i = 1, n_solv
     cmass_solv_mpc(i) = inmpc%cmass_solv
  end do
  
!  write(*,*) "kanada"
  
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  !**********
  if(i_simulate_type == SIM%MPC) then
     call ukoto_uiread2(luninp, lunout, 'mpc_dynamics    ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
  end if
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "mpc_dynamics" field in the input file'
     call util_error(2, error_message)
  end if
  
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     do iequa = 1, nequat

        cvalue = 'pbbox_min_x'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             pbbox_min_x, cvalue)
        
        cvalue = 'pbbox_max_x'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             pbbox_max_x, cvalue)

        cvalue = 'pbbox_min_y'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             pbbox_min_y, cvalue)

        cvalue = 'pbbox_max_y'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             pbbox_max_y, cvalue)

        cvalue = 'pbbox_min_z'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             pbbox_min_z, cvalue)

        cvalue = 'pbbox_max_z'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             pbbox_max_z, cvalue)

        cvalue = 'ngrid_x'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmpc%ngrid_x, cvalue)

        cvalue = 'ngrid_y'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmpc%ngrid_y, cvalue)

        cvalue = 'ngrid_z'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmpc%ngrid_z, cvalue)

        cvalue = 'n_step_collision'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmpc%nratio_colli_step, cvalue)

        cvalue = 'nratio_vcorrect_step'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmpc%nratio_vcorrect_step, cvalue)

        cvalue = 'n_av_solvent'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmpc%n_av_solv_ingrid, cvalue)

        cvalue = 'rmass_solvent'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inmpc%cmass_solv, cvalue)

        cvalue = 'rotate_angle'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inmpc%rotate_angle_colli, cvalue)

        cvalue = 'i_thermal_mpc'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmpc%i_thermal_mpc, cvalue)

        cvalue = 'i_flag_check_mpc'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmpc%i_flag_check_mpc, cvalue)

        cvalue = 'i_flag_check2_mpc'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmpc%i_flag_check2_mpc, cvalue)

     end do
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast(pbbox_min_x, 1, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pbbox_max_x, 1, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pbbox_min_y, 1, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pbbox_max_y, 1, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pbbox_min_z, 1, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pbbox_max_z, 1, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(inmpc, inmpc%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
#endif

  !*********

!  write(*,*) "kanada sssk---1"

  !! base vector for x(1)-axis
  !!base_vector_mpc(1,1)=1.0e0_PREC
  !!base_vector_mpc(1,2)=0.0e0_PREC
  !!base_vector_mpc(1,3)=0.0e0_PREC
  
  !! base vector for y(2)-axis
  !!base_vector_mpc(2,1)=0.0e0_PREC
  !!base_vector_mpc(2,2)=1.0e0_PREC
  !!base_vector_mpc(2,3)=0.0e0_PREC
  
  !! base vector for z(3)-axis
  !!base_vector_mpc(3,1)=0.0e0_PREC
  !!base_vector_mpc(3,2)=0.0e0_PREC
  !!base_vector_mpc(3,3)=1.0e0_PREC
  
  !!*********************************************

  !!write(*,*) "min(x,y,z)",pbbox_min_x, pbbox_min_y, pbbox_min_z
  !!write(*,*) "min(x,y,z)",pbbox_max_x, pbbox_max_y, pbbox_max_z

  !! set of peridoic boudary box origin
  inmpc%pbox_origin_x=pbbox_min_x
  inmpc%pbox_origin_y=pbbox_min_y
  inmpc%pbox_origin_z=pbbox_min_z
  pbox_origin_mpc(1)=inmpc%pbox_origin_x  !! periodic boundary box origin x-axis 
  pbox_origin_mpc(2)=inmpc%pbox_origin_y  !! periodic boundary box origin y-axis 
  pbox_origin_mpc(3)=inmpc%pbox_origin_z  !! periodic boundary box origin z-axis 
  write(*,*) "origin_pbbox(x,y,z)", pbox_origin_mpc(1) , pbox_origin_mpc(2), pbox_origin_mpc(3)

  
  !! number of grid(cell) xyz-axis
  ngrid_mpc(1)=inmpc%ngrid_x !! number of grid(cell) x-axis
  ngrid_mpc(2)=inmpc%ngrid_y !! number of grid(cell) y-axis
  ngrid_mpc(3)=inmpc%ngrid_z !! number of grid(cell) z-axis
  write(*,*) "grid number (x,y,z):", ngrid_mpc(1),ngrid_mpc(2),ngrid_mpc(3)

  !!periodic boudary box size (x,y,z)
  inmpc%pbox_size_x = pbbox_max_x - pbbox_min_x  !!periodic boudary box size (x,y,z)
  inmpc%pbox_size_y = pbbox_max_y - pbbox_min_y  !!periodic boudary box size (x,y,z)
  inmpc%pbox_size_z = pbbox_max_z - pbbox_min_z  !!periodic boudary box size (x,y,z)
  pbox_size_mpc(1)=inmpc%pbox_size_x  !!periodic boudary box size (x,y,z)
  pbox_size_mpc(2)=inmpc%pbox_size_y  !!periodic boudary box size (x,y,z)
  pbox_size_mpc(3)=inmpc%pbox_size_z  !!periodic boudary box size (x,y,z)
  write(*,*) "pbbox size(x,y,z):", pbox_size_mpc(1),pbox_size_mpc(2),pbox_size_mpc(3)
  write(*,*) "pbbox max(x,y,z):", pbox_origin_mpc(1)+pbox_size_mpc(1), pbox_origin_mpc(2)+pbox_size_mpc(2), pbox_origin_mpc(3)+pbox_size_mpc(3)

  !!grid box size (x,y,z)
  inmpc%grid_size_x = pbox_size_mpc(1)/ngrid_mpc(1) 
  inmpc%grid_size_y = pbox_size_mpc(2)/ngrid_mpc(2) 
  inmpc%grid_size_z = pbox_size_mpc(3)/ngrid_mpc(3) 
  grid_size_mpc(1)=inmpc%grid_size_x  !!grid size (x,y,z)
  grid_size_mpc(2)=inmpc%grid_size_y  !!grid size (x,y,z)
  grid_size_mpc(3)=inmpc%grid_size_z  !!grid size (x,y,z)
  write(*,*) "grid box size (x,y,z):", grid_size_mpc(1), grid_size_mpc(2), grid_size_mpc(3)
  
  !!inmpc%nratio_colli_step
  write(*,*) "n_step_collision",inmpc%nratio_colli_step
  
  !!n_av_solvent
  write(*,*) "n_av_solvent in grid",inmpc%n_av_solv_ingrid
  inmpc%n_all_solv=ngrid_mpc(1)*ngrid_mpc(2)*ngrid_mpc(3)*inmpc%n_av_solv_ingrid
  write(*,*) "n_all_solvent in mpc",inmpc%n_all_solv
  
  !! mass of solvent in mpc
  n_solv=inmpc%n_all_solv
  do i = 1, n_solv
     cmass_solv_mpc(i) = inmpc%cmass_solv
  end do
  write(*,*) "mass of solvent particle in mpc",cmass_solv_mpc(5)
  
  write(*,*) "rotation angle in collision step of mpc",inmpc%rotate_angle_colli
  
  !!---------------------------------------------------------------------------
  !! initial position setting for mpc solvent
  !!n_solv=inmpc%n_all_solv
  !!do is = 1, n_solv
  !!do idimn = 1, 3
  !!      xyz_solv_mpc(idimn, is) = pbox_origin_mpc(idimn) + pbox_size_mpc(idimn)*recipe_rand01(irand)
  !!   end do
  !!   write (*,*) "inira",xyz_solv_mpc(1, is),xyz_solv_mpc(2, is),xyz_solv_mpc(3, is)
  !!end do
  !!----------------------------------------------------------------
  !!-----------------------------------------------------------------
  !! initial velocity setting for mpc slovent
  !!n_solv=inmpc%n_all_solv
  !! nmp
  !!do is=1, n_solv
  !!  do idimn = 1, 3
  !!      velo_solv_mpc(idimn, is)= 
  !!   end do
  !!end do
  !!velo_solv_mpc(3, MXSOLV_MPC)=
  !!----------------------------------------------------------------
  

end subroutine setp_para_mpc
