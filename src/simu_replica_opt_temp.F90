!>simu_replica_opt_temp
!> @brief Updates the temperature distribution of replicas          &
!>        by feedback-optimization algorithm (Katzgraber et al. 2006)

!
! #############################
!   Feedback-optimized REMD
! #############################
! 
! References:
! #  S. Trebst, M. Troyer and U. H. Hansmann
!    Optimized parallel tempering simulations of proteins
!    J. Chem. Phys. (2006) 124: 174903
! #  H. G. Katzgraber, S. Trebst, D. A. Huse and M. Troyer
!    Feedback-optimized parallel tempering Monte Carlo
!    J. Stat. Mech. (2006) P03018
!
! Aug. 13, 2010, implemented by Naoto Hori
!
subroutine simu_replica_opt_temp(i_current_stage)

   use const_maxsize
   use const_index
   use const_physical
   use var_inp,     only : outfile
   use var_replica, only : n_current_up, n_current_down, &
                           n_replica_all, lab2val, inrep,  &
                           up_or_down, lab2rep
   use mpiconst

   implicit none

   integer, intent(out), optional :: i_current_stage
   integer, save :: iopt_stage = 0
   integer, save :: n_return = 0
   integer    :: i
   integer    :: n_up(REPTYPE%MAX, n_replica_all)
   integer    :: n_down(REPTYPE%MAX, n_replica_all)
   real(PREC) :: target_value
   real(PREC) :: t, dt
   real(PREC) :: f(n_replica_all)
   real(PREC) :: dfodt(n_replica_all)
   real(PREC) :: eta(n_replica_all)
   real(PREC) :: eta_integral
   real(PREC) :: eta_integral_prev
   real(PREC) :: c
   real(PREC) :: delta_T_right
   integer    :: lunout
   integer    :: irep, ilabel
   integer    :: num_add, num_add_i
   logical    :: flg_return
   integer, parameter :: TEMP = REPTYPE%TEMP
   real(PREC) :: lab2val_new(n_replica_all)
   character(CARRAY_MSG_ERROR) :: error_message

   lunout = outfile%data

   n_up(:,:)   = sum(n_current_up,   3) 
   n_down(:,:) = sum(n_current_down, 3)

   ! If any histogram (n_up(T),n_down(T)) has 0 value,
   ! return because the number of steps would be insufficient.
   flg_return = .false.
   if (myrank == 0) then
      write(lunout,*) ''
      write(lunout,*) ''
      write(lunout,*) '#OPT_TEMP histogram'
      write(lunout,*) '#OPT_TEMP ilabel n_up n_down'
      write(lunout,'(i3,1x,i12,1x,i12)') 1, n_up(TEMP,1), n_down(TEMP,1)
   endif

   do ilabel = 2, n_replica_all-1
      if (myrank == 0) then
         write(lunout,'(i3,1x,i12,1x,i12)') ilabel, n_up(TEMP,ilabel), &
                                                    n_down(TEMP,ilabel)
      endif
      if (n_up(TEMP,ilabel)   == 0) flg_return = .true.
      if (n_down(TEMP,ilabel) == 0) flg_return = .true.
   enddo

   if (myrank == 0) then
      write(lunout,'(i3,1x,i12,1x,i12)') n_replica_all, n_up(TEMP,n_replica_all), &
                                                        n_down(TEMP,n_replica_all)
   endif

   if (flg_return) then
      n_return = n_return + 1
      if (myrank == 0) then
         write(lunout,*) ''
         write(lunout,*) 'OPT_TEMP return', n_return
         write(lunout,*) ''
      endif
      ! ########### Return #############
      i_current_stage = iopt_stage
      return
   endif

   ! Double the step interval 
   inrep%n_step_opt_temp = (inrep%n_step_opt_temp * (n_return+1)) * 2
   n_return = 0
   iopt_stage = iopt_stage + 1

   if (myrank == 0) then
      write(lunout,*) ''
      write(lunout,*) ''
      write(lunout,*) '#OPT_TEMP optimization REMD stage',iopt_stage
   
      ! ########################################
      !  f(T)
      !write(lunout,*) '#OPT_TEMP_f'
      !write(lunout,*) '#label T f(T)'
      do ilabel = 1, n_replica_all
         f(ilabel) =  real(n_up(TEMP,ilabel),PREC) &
                    / real((n_up(TEMP,ilabel) + n_down(TEMP,ilabel)),PREC)
         !write(lunout,*) ilabel, lab2val(ilabel,TEMP), f(ilabel)
      enddo
   
      ! ########################################
      !  df(T)/dT
      !write(lunout,*) ''
      !write(lunout,*) '#OPT_TEMP_dfodt'
      !write(lunout,*) '#label T df(T)/dT'
      do ilabel = 1, n_replica_all-1
         delta_T_right= lab2val(ilabel+1, TEMP) - lab2val(ilabel, TEMP)
         dfodt(ilabel) = (f(ilabel+1) - f(ilabel)) / delta_T_right
         !write(lunout,*) ilabel, lab2val(ilabel,TEMP), dfodt(ilabel)
      enddo
   
      ! ########################################
      !   eta(T) (not normalized, here)
      !write(lunout,*) ''
      !write(lunout,*) '#OPT_TEMP_eta'
      !write(lunout,*) '#label T eta(T)'
      do ilabel = 1, n_replica_all-1
         if (dfodt(ilabel) < 0) then
            eta(ilabel) = sqrt( dfodt(ilabel) / (lab2val(ilabel,TEMP)-lab2val(ilabel+1,TEMP)))
         else
            eta(ilabel) = sqrt( dfodt(ilabel) / (lab2val(ilabel+1,TEMP)-lab2val(ilabel,TEMP)))
         endif
         !write(lunout,*) ilabel, lab2val(ilabel,TEMP), eta(ilabel)
      enddo
   
      ! ########################################
      !  Normalizes eta(T)
      c = 0.0_PREC 
      do irep = 1, n_replica_all-1
         c = c + eta(irep) * (lab2val(irep+1,TEMP) - lab2val(irep, TEMP))
      enddo
      c = 1.0_PREC / c
      !write(lunout,*) ''
      !write(lunout,*) '#OPT_TEMP_c ',c
      !write(lunout,*) ''
      !write(lunout,*) '#OPT_TEMP_eta_new'
      !write(lunout,*) '#label T eta(T)'
      do ilabel = 1, n_replica_all-1
         eta(ilabel) = c * eta(ilabel)
      !   write(lunout,*) ilabel, lab2val(ilabel,TEMP), eta(ilabel)
      enddo
   
   
      ! ########################################
      !  Calculates new distribution {T}(new)
      ! ########################################
      write(lunout,*) ''
      write(lunout,*) '#OPT_TEMP new distribution'
      call write_new_variable(lunout, 1, lab2val(1, TEMP))
      lab2val_new(1) = lab2val(1,TEMP)
   
      dt = 1.0e-5_PREC
      t = lab2val(1,TEMP)
      i = 1
      do ilabel = 2, n_replica_all-1
         target_value = 1.0_PREC / real(n_replica_all-1,PREC)
   
         eta_integral = 0.0_PREC
         eta_integral_prev = 0.0_PREC
         num_add = 0
         num_add_i = 0
   
         do while (eta_integral < target_value)
            do while (t >= lab2val(i+1,TEMP))
               if (i < n_replica_all-1) then
                  i = i + 1
                  if (lab2val(i-1,TEMP) > lab2val_new(ilabel-1)) then
                     eta_integral_prev = eta_integral_prev &
                                        + eta(i-1) * (lab2val(i,TEMP)-lab2val(i-1,TEMP))
                  else
                     eta_integral_prev = eta_integral_prev &
                                        + eta(i-1) * (lab2val(i,TEMP)-lab2val_new(ilabel-1))
                  endif
                  num_add_i = 0
               else
                  exit
               endif
            enddo
            
            num_add = num_add + 1
            num_add_i = num_add_i + 1
            eta_integral = eta_integral_prev + eta(i) * (dt * num_add_i) 
            t = lab2val_new(ilabel-1) + dt * num_add
         enddo
   
         lab2val_new(ilabel) = t
         call write_new_variable(lunout, ilabel, t)
      enddo
      call write_new_variable(lunout, n_replica_all, lab2val(n_replica_all, TEMP))
   
      do ilabel = 2, n_replica_all-1
         lab2val(ilabel, TEMP) = lab2val_new(ilabel)
      enddo

   endif

#ifdef MPI_PAR
   call MPI_Bcast(lab2val, MXREPLICA*REPTYPE%MAX, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
#endif

   ! ########################################
   !  Clear histograms
   up_or_down(1:REPTYPE%MAX, 1                ) =  1
   up_or_down(1:REPTYPE%MAX, 2:n_replica_all-1) = -1
   up_or_down(1:REPTYPE%MAX, n_replica_all    ) =  0
   n_current_up  (1:REPTYPE%MAX, 1:n_replica_all, 1:inrep%n_period_prob) = 0
   n_current_down(1:REPTYPE%MAX, 1:n_replica_all, 1:inrep%n_period_prob) = 0

   i_current_stage = iopt_stage

!#########################################################################################
contains

   subroutine write_new_variable(lunout,i,x)
      use const_maxsize
      use const_index
      implicit none
      integer,intent(in)    :: lunout,i
      real(PREC),intent(in) :: x

      if (i<10) then
         write(lunout,'(a8,i1,a6,f12.6)') 'REPLICA(',i,')   = ',x
      elseif (i<100) then
         write(lunout,'(a8,i2,a5,f12.6)') 'REPLICA(',i,')  = ',x
      elseif (i<1000) then
         write(lunout,'(a8,i3,a4,f12.6)') 'REPLICA(',i,') = ',x
      else
        write(error_message,*) 'defect at simu_replica_opt_temp. PROGRAM STOP'
        call util_error(ERROR%STOP_ALL, error_message)
      endif
   endsubroutine write_new_variable

endsubroutine simu_replica_opt_temp
