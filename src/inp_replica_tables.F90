! inp_replica_tables
!> @brief This subroutine is to make the replica set tables.

subroutine inp_replica_tables()

   use const_maxsize
   use const_index
   use var_io,     only : flg_rst
   use var_replica, only : inrep, rep2lab, lab2rep, lab2val, &
                           lab2step, exchange_step, rep2step,     &
                           flg_rep, rep2val, n_replica_all,  &
                           label_order, &
                           make_type_array, make_exchange_pair_tb
#ifdef MPI_PAR
   use mpiconst
#endif

   implicit none

   ! local
   integer     :: ier
   integer     :: irep, ivar, ido, iseries
   integer     :: label
   integer     :: n_division_pre, n_division, n_continue, n_skip
   real(PREC)  :: rvalue
   integer     :: nstep
   character(CARRAY_MSG_ERROR) :: error_message

   ! ------------------------------------- ! 
#ifdef _DEBUG
   write(*,*) 'inp_replica_tables : START'
#endif

   allocate(label_order(REPTYPE%MAX, n_replica_all), stat=ier)
   if (ier/=0) then
      write(error_message, *) 'defect at allocate(label_order), PROGRAM STOP'
      call util_error(ERROR%STOP_ALL, error_message)
   endif
   label_order(:,:) = 0

   ! ###############################################################
   !  lab2val
   ! ###############################################################
   n_division = 1
   do ivar = 1, REPTYPE%MAX

      if (.not. flg_rep(ivar)) then
         cycle
      endif

      n_division_pre = n_division
      n_division     = n_division * inrep%n_replica(ivar)

      do irep = 1, inrep%n_replica(ivar)
   
         if (inrep%i_style(ivar) == REPVARSTYLE%LINEAR) then
            !rvalue = inrep%lowest(ivar) + (irep-1) * inrep%interval(ivar)
            rvalue = inrep%lowest(ivar) &
            + (irep-1) * (inrep%highest(ivar)-inrep%lowest(ivar))/(real(inrep%n_replica(ivar),PREC)-1.0)

         else if (inrep%i_style(ivar) == REPVARSTYLE%EXPONENTIAL) then
            if (inrep%n_replica(ivar) == 1) then
               rvalue = inrep%lowest(ivar)
            else
               if (irep == 1) then
                  rvalue = inrep%lowest(ivar)
               else if (irep == inrep%n_replica(ivar)) then
                  rvalue = inrep%highest(ivar)
               else
                  rvalue = ((real(irep, PREC) -1.0) ** inrep%exponent_alpha(ivar)) &
                          /((real(inrep%n_replica(ivar), PREC) -1.0) ** inrep%exponent_alpha(ivar))
                  rvalue = inrep%lowest(ivar)                           &
                          * ( (inrep%highest(ivar) / inrep%lowest(ivar)) ** rvalue)
               endif
            endif
         
         else if (inrep%i_style(ivar) == REPVARSTYLE%EXPLICIT) then
            rvalue = inrep%var(irep,ivar)

         else
           write(error_message, *) 'defect at replica_set_table, PROGRAM STOP'
           call util_error(ERROR%STOP_ALL, error_message)
         endif
          
         ! store
         n_continue = n_replica_all / n_division
         n_skip     = n_replica_all / n_division_pre
         do ido = 1, n_division_pre
            do iseries = 1, n_continue

               label = (irep-1)*n_continue + (ido-1)*n_skip + iseries

               lab2val(label, ivar) = rvalue
               label_order(ivar, label) = irep

            enddo
         enddo

      enddo
   enddo

   ! ###############################################################
   !  lab2step
   ! ###############################################################
   if (inrep%n_step_replica > 0) then
      nstep = inrep%n_step_replica
      lab2step(1:n_replica_all) = int(nstep, kind=L_INT)
      exchange_step(1:n_replica_all) = int(nstep, kind=L_INT)
      ! CAUTION: "exchange_step" could be updated in setp_md_info() 
   else
      lab2step(1:n_replica_all) = 1
   endif


   ! ###############################################################
   !  rep2lab, lab2rep
   ! ###############################################################
   do irep = 1, n_replica_all
      rep2lab(irep) = irep
      lab2rep(irep) = irep
   enddo

   if (flg_rst) then
      call read_rst(RSTBLK%REPLICA)
   endif


   ! #####################################################
   !  type_array
   ! #####################################################
   call make_type_array()


   ! ###############################################################
   !  exchagne_pair_tb
   ! ###############################################################
   call make_exchange_pair_tb()


#ifdef _DEBUG
  write(*,*) 'inp_replica_tables : END'
#endif

endsubroutine inp_replica_tables !nh
