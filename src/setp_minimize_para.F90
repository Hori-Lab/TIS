! setp_para_emin
!> @brief Reading parameters of "energy_minimization" field in input file

subroutine setp_minimize_para()

   use const_maxsize
   use const_index
   use const_physical
   use var_inp, only : infile, outfile, ifile_out_opt
   use var_setp, only : insimu
   use var_emin, only : inemin
#ifdef MPI_PAR
   use mpiconst
#endif

   implicit none
   logical       :: flg_sd, flg_cg
   integer       :: i
   integer       :: icol, isim
   integer       :: luninp, lunout
   integer       :: iline, nlines, iequa, nequat
   real(PREC)    :: me10 = 1.0e0_PREC
   character(4)  :: kfind
   character(CARRAY_MXCOLM)  :: cwkinp(CARRAY_MXLINE)
   character(CARRAY_MXCOLM)  :: cvalue
   character(CARRAY_MXCOLM)  :: csides(2, CARRAY_MXEQUA)
   character(CARRAY_MXCOLM) :: ctmp02
   character(CARRAY_MSG_ERROR) :: error_message

   ! ---------------------------------------------------------------------
   luninp = infile%inp
   lunout = outfile%data

   inemin%i_method(:) = EMIN_METHOD%VOID
   inemin%i_out = 0   ! default
   inemin%n_out_step = 0
   inemin%eps = INVALID_VALUE  ! will be set to dedault value, 10.0*(machine epsilon)
   inemin%sd_lambda_init = 1.0e-1_PREC  ! default: 0.1
   inemin%sd_rho_accept = 1.2           ! default: 1.2
   inemin%sd_rho_reject = 0.2           ! default: 0.2
   inemin%cg_lambda_init = 1.0e0_PREC  ! default: 1.0 
   inemin%cg_rho = 0.95e0_PREC         ! default: 0.9
   inemin%cg_wolfe_c1 = 1.0e-4_PREC    ! default: 0.0001
   inemin%cg_wolfe_c2 = 1.0e-1_PREC    ! default: 0.1

#ifdef MPI_PAR
   if (myrank == 0) then
#endif
   ! -----------------------------------------------------------------
   ! Find the field
   rewind(luninp)
   call ukoto_uiread2(luninp, lunout, 'energy_minimize ', kfind, &
                      CARRAY_MXLINE, nlines, cwkinp)
   if(kfind /= 'FIND') then
      error_message = 'Error: cannot find "energy_minize" field in the input file'
      call util_error(ERROR%STOP_ALL, error_message)
   end if

   ! -----------------------------------------------------------------
   ! Read the field
   do iline = 1, nlines
      call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
          
      do iequa = 1, nequat
         ctmp02 = csides(1, iequa)
 
         if(ctmp02(1:8) == 'i_method') then
            do icol = 10, CARRAY_MXCOLM
               if(ctmp02(icol:icol) == ')') exit
            end do
            read (ctmp02(10:icol-1), *) isim
            if (isim <= 0 .OR. isim > insimu%n_step_sim) then
               error_message = 'Error: invalid value for i_method in <<<<energy_minimize'
               call util_error(ERROR%STOP_ALL, error_message)
            endif
            cvalue = ctmp02(1:icol)
            call ukoto_ivalue2(lunout, csides(1, iequa), &
                               inemin%i_method(isim), cvalue)
         end if
                            
         cvalue = 'i_out'
         call ukoto_ivalue2(lunout, csides(1, iequa), &
                            inemin%i_out, cvalue)

         cvalue = 'n_out_step'
         call ukoto_ivalue2(lunout, csides(1, iequa), &
                            inemin%n_out_step, cvalue)

         cvalue = 'sd_lambda_init'
         call ukoto_rvalue2(lunout, csides(1, iequa), &
                            inemin%sd_lambda_init, cvalue)

         cvalue = 'sd_rho_accept'
         call ukoto_rvalue2(lunout, csides(1, iequa), &
                            inemin%sd_rho_accept, cvalue)

         cvalue = 'sd_rho_reject'
         call ukoto_rvalue2(lunout, csides(1, iequa), &
                            inemin%sd_rho_reject, cvalue)

         cvalue = 'epsilon'
         call ukoto_rvalue2(lunout, csides(1, iequa), &
                            inemin%eps, cvalue)

         cvalue = 'cg_lambda_init'
         call ukoto_rvalue2(lunout, csides(1, iequa), &
                            inemin%cg_lambda_init, cvalue)

         cvalue = 'cg_rho'
         call ukoto_rvalue2(lunout, csides(1, iequa), &
                            inemin%cg_rho, cvalue)

         cvalue = 'cg_wolfe_c1'
         call ukoto_rvalue2(lunout, csides(1, iequa), &
                            inemin%cg_wolfe_c1, cvalue)

         cvalue = 'cg_wolfe_c2'
         call ukoto_rvalue2(lunout, csides(1, iequa), &
                            inemin%cg_wolfe_c2, cvalue)
      enddo
   enddo
   
   ! -----------------------------------------------------------------
   ! checking input variables
   ! i_out
   if (inemin%i_out == 0) then
      ! Nothing
   else if (inemin%i_out == 1) then
      ! Output to .data file
   else if (inemin%i_out == 2) then
      ! Output to .opt file
      if (ifile_out_opt /= 1) then
         error_message = 'Error: invalid value for i_out in <<<<energy_minimize. '//&
                         '"opt" is needed in the OUTPUT line.'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   else
      error_message = 'Error: invalid value for i_out in <<<<energy_minimize'
      call util_error(ERROR%STOP_ALL, error_message)
   endif
   ! n_out_step
   if (inemin%i_out > 0) then
      if (inemin%n_out_step <= 0) then
         error_message = 'Error: invalid value for n_out_step in <<<<energy_minimize'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endif

   ! i_method
   flg_sd = .false.
   flg_cg = .false.
   do i = 1, insimu%n_step_sim
      if (inemin%i_method(i) == EMIN_METHOD%SD) then
         flg_sd = .true.
      else if (inemin%i_method(i) == EMIN_METHOD%CG) then
         flg_cg = .true.
      else
         error_message = 'Error: invalid value for i_method in <<<<energy_minimize'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   enddo

   ! eps
   me10 = 1.0e1_PREC * epsilon(me10)  ! (machine epsilon) * 10.0
   if (inemin%eps > INVALID_JUDGE .OR. inemin%eps < me10) then
      write(lunout,*) 'epsilon is set to default value: 10*(machine epsilon) =', me10
      inemin%eps = me10
   endif

   if (flg_sd) then
       ! sd_lambda_init
       if (inemin%sd_lambda_init > INVALID_JUDGE) then
          error_message = 'Error: invalid value for sd_lambda_init in <<<<energy_minimize'
          call util_error(ERROR%STOP_ALL, error_message)
       endif
       ! sd_rho_accept
       if (inemin%sd_rho_accept > INVALID_JUDGE) then
          error_message = 'Error: invalid value for sd_rho_accept in <<<<energy_minimize'
          call util_error(ERROR%STOP_ALL, error_message)
       endif
       ! sd_rho_reject
       if (inemin%sd_rho_reject > INVALID_JUDGE) then
          error_message = 'Error: invalid value for sd_rho_reject in <<<<energy_minimize'
          call util_error(ERROR%STOP_ALL, error_message)
       endif
   endif

   if (flg_cg) then
       ! cg_lambda_init
       if (inemin%cg_lambda_init > INVALID_JUDGE) then
          error_message = 'Error: invalid value for cg_lambda_init in <<<<energy_minimize'
          call util_error(ERROR%STOP_ALL, error_message)
       endif
      ! cg_rho
      if (inemin%cg_rho >= 1.0e0_PREC .OR. inemin%cg_rho < 0.0e0_PREC) then
         write(error_message, *) &
         'Error: invalid value for cg_rho in <<<<energy_minimize. ', &
         ' 0.0 < cg_rho < 1.0 should be satisfied.'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
      ! cg_wolfe_c1, cg_wolfe_c2
      if ((inemin%cg_wolfe_c1 <=  0.0e0_PREC) .OR. &
          (inemin%cg_wolfe_c2 <= inemin%cg_wolfe_c1) .OR. &
          (0.5e0_PREC <= inemin%cg_wolfe_c2 ) )then
         write(error_message, *) &
         'Error: invalid value for cg_wolfe_c1/2 in <<<<energy_minimize. ', &
         ' 0 < c1 < c2 < 0.5 shoulde be satisfied.'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endif

#ifdef MPI_PAR
   end if

   call MPI_Bcast (inemin, inemin%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

endsubroutine setp_minimize_para
