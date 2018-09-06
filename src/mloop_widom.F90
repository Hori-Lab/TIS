subroutine mloop_widom

   use const_physical
   use const_index
   use var_setp,   only : inwidom, indtrna15, insimu
   use var_struct, only : ntp, iclass_tp, charge_tp, &
                          tp_exv_dt15_eps, tp_exv_dt15_rad
   use var_simu,   only : widom_iw, widom_chp, widom_energy, tempk
   use mpiconst

   integer :: i, itp
   character(CARRAY_MSG_ERROR) :: error_message
   character(CARRAY_MSG_ERROR),parameter :: msg_er_allocate = &
           'failed in memory allocation at widom, PROGRAM STOP'

   ntp = inwidom%n_Mg_add + inwidom%n_Na_add + inwidom%n_K_add + inwidom%n_Cl_add
   if (ntp > MXTP) then 
      error_message = 'Error: ntp > MXTP'
      call util_error(ERROR%STOP_ALL, error_message)
   endif

   iclass_tp(1:ntp) = CLASS%ION

   itp = 0
   do i = 1, inwidom%n_Mg_add
      itp = itp + 1
      charge_tp(itp) = + 2.0
      tp_exv_dt15_eps(itp) = indtrna15%exv_eps( DT15EXV%MG2 )
      tp_exv_dt15_rad(itp) = indtrna15%exv_rad( DT15EXV%MG2 )
   enddo
   do i = 1, inwidom%n_K_add
      itp = itp + 1
      charge_tp(itp) = + 1.0
      tp_exv_dt15_eps(itp) = indtrna15%exv_eps( DT15EXV%K )
      tp_exv_dt15_rad(itp) = indtrna15%exv_rad( DT15EXV%K )
   enddo
   do i = 1, inwidom%n_Na_add
      itp = itp + 1
      charge_tp(itp) = + 1.0
      tp_exv_dt15_eps(itp) = indtrna15%exv_eps( DT15EXV%NA )
      tp_exv_dt15_rad(itp) = indtrna15%exv_rad( DT15EXV%NA )
   enddo
   do i = 1, inwidom%n_Cl_add
      itp = itp + 1
      charge_tp(itp) = - 1.0
      tp_exv_dt15_eps(itp) = indtrna15%exv_eps( DT15EXV%CL )
      tp_exv_dt15_rad(itp) = indtrna15%exv_rad( DT15EXV%CL )
   enddo

   widom_iw = inwidom%iw_init
   if (widom_iw == -1) then
      call read_rst(RSTBLK%WIDOM)
   else if (widom_iw == 0) then
      widom_chp = 0.0e0_PREC
   else
      widom_chp = real(widom_iw,kind=PREC) * exp( -inwidom%chp_init / (BOLTZ_KCAL_MOL * insimu%tempk))
   endif

   if (allocated(widom_energy)) deallocate(widom_energy)
   !allocate( widom_energy(E_TYPE%MAX, 0:nthreads-1), stat=i)
   allocate( widom_energy(E_TYPE%MAX), stat=i)
   if (i/=0) call util_error(ERROR%STOP_ALL, msg_er_allocate)

endsubroutine mloop_widom
