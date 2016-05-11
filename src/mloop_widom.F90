subroutine mloop_widom

   use const_physical
   use const_index
   use var_setp,   only : inwidom, indtrna15
   use var_struct, only : ntp, iclass_tp, charge_tp, &
                          tp_exv_dt15_eps, tp_exv_dt15_rad
   use var_simu,   only : widom_iw, widom_chp

   integer :: i, itp

   ntp = inwidom%n_Mg_add + inwidom%n_Na_add + inwidom%n_K_add + inwidom%n_Cl_add
   if (ntp > MXTP) then
      call util_error(ERROR%STOP_ALL, 'Error: ntp > MXTP')
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

   widom_iw = 0
   widom_chp    = 0.0e0_PREC

endsubroutine mloop_widom
