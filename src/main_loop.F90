! main_loop
!> @brief This subroutine is to perform "main loop" for simualation.

! **********************************************************************
subroutine main_loop()
      
  use const_maxsize
  use const_physical
  use const_index
  use var_io,     only : i_go_native_read_style, i_run_mode, &
                          ifile_out_psf, outfile
  use var_setp,    only : insimu, inmisc, inflp, inele
  use var_mgo,     only : inmgo
  use var_fmat,    only : infmat, i_num_sum
  use var_simu,    only : istep_sim, mstep_sim
  use mpiconst

  implicit none
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
#ifdef _DEBUG
  write(*,*) '#### start main_loop'
#endif

  ! -------------------------------------------------------------------
  ! memory allocation
  ! -------------------------------------------------------------------
  call allocate_simu()

  ! -------------------------------------------------------------------
  ! istep_sim: the numbers of potential switch + 1
  ! -------------------------------------------------------------------
  if (i_run_mode == RUN%FMAT) then
     mstep_sim = infmat%n_iter
     if (i_go_native_read_style == NATIVEREAD%INFO) then
        write(error_message, *) 'native-info is not available in fmat mode'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
     if(inmgo%i_multi_mgo >= 1) then
        write(error_message, *) 'multi-go is not available in fmat mode'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  else
     mstep_sim = insimu%n_step_sim
  endif

  do istep_sim = insimu%i_step_sim_init, mstep_sim

     ! -------------------------------------------------------------------
     ! read native information mainly for potential switching
     ! -------------------------------------------------------------------
     if (i_go_native_read_style == NATIVEREAD%INFO) then
        call mloop_nativeinfo(istep_sim)
     endif

     ! -------------------------------------------------------------------
     ! setting up for multiple Go
     ! -------------------------------------------------------------------
     if(inmgo%i_multi_mgo >= 1) then
        call allocate_mgo()
        call mloop_setup_mgo(istep_sim)
     end if

     ! -------------------------------------------------------------------
     ! setting up for flexible local potential
     ! -------------------------------------------------------------------
     if (inflp%i_flp == 1 .or. inmisc%force_flag_local(LINTERACT%L_FLP)) then
        call mloop_flexible_local()
     end if

     ! -------------------------------------------------------------------
     ! delete interaction
     ! -------------------------------------------------------------------
     if (inmisc%i_del_int >= 1) then
        call mloop_del_int()
     end if

     ! -------------------------------------------------------------------
     ! setting up for DTRNA
     ! -------------------------------------------------------------------
     if (inmisc%force_flag(INTERACT%DTRNA)) then
        call mloop_dtrna()
     end if

     ! -------------------------------------------------------------------
     ! setting up for Ewald method
     ! -------------------------------------------------------------------
     if (inmisc%force_flag(INTERACT%ELE)) then
        if (inele%i_function_form == 2) then
           call mloop_ewld()
        endif
     end if

     ! -------------------------------------------------------------------
     ! setting up for Widom method
     ! -------------------------------------------------------------------
     if (i_run_mode == RUN%WIDOM    ) then
        call mloop_widom()
     end if

     ! -------------------------------------------------------------------
     ! writing the native structure
     ! -------------------------------------------------------------------
     if(istep_sim == insimu%i_step_sim_init) then
        call write_seq()
        call write_nativeinfo(outfile%ninfo)
        if (ifile_out_psf == 1) then
           call write_psf()
        endif
        if (i_run_mode == RUN%FMAT) call write_fmat(istep_sim)
     endif 

     ! -------------------------------------------------------------------
     ! molecular dynamics simulation
     ! -------------------------------------------------------------------
     
     call mloop_simulator()

     if (i_run_mode == RUN%FMAT) then
        if (infmat%i_type == FMATTYPE%HOMO) then
           call mloop_fmat_homo(istep_sim)
           call mloop_recalc_coef()
        else if (infmat%i_type == FMATTYPE%HETERO) then
           call mloop_fmat_hetero()
        else
           write(error_message, *) 'undefined inrep%i_type in main_loop'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        call write_fmat(istep_sim)
        i_num_sum = 0 
     else if (i_run_mode == RUN%REPLICA) then
        call write_rep_result()
     endif

     if(inmgo%i_multi_mgo >= 1) then
        call deallocate_mgo()
     endif
  end do

  if (i_run_mode == RUN%FMAT) then
     if (infmat%i_type == FMATTYPE%HETERO) then
        call write_nativeinfo(outfile%fmat)
     endif
  endif

  ! -------------------------------------------------------------------
  ! memory deallocation
  ! -------------------------------------------------------------------
  call deallocate_simu()


#ifdef _DEBUG
  write(*,*) '#### end main_loop'
#endif

end subroutine main_loop
