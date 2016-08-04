!allocate_replica
!> @brief Allocates memory resources for xyz_mp_rep &
!>        and other replica-related variables.

subroutine allocate_replica(xyz_mp_init)

   use const_maxsize
   use const_index
   use const_physical
   use var_io,     only : i_run_mode
   use var_setp,    only : inele
   use var_struct,  only : nmp_all, xyz_mp_rep, pxyz_mp_rep, &
                           xyz_ele_rep, pxyz_ele_rep
   use var_replica, only : hist_combination,    &
                           hist_exchange, hist_attempt,&
                           rate_exchange, &
                           n_replica_all, n_replica_mpi, &
                           inrep, up_or_down, &
                           n_current_up, n_current_down, &
                           n_turnover

   implicit none

   real(PREC), intent(in) :: xyz_mp_init(SDIM, MXMP)

   integer :: irep, imp, ier
   logical :: flg_error
   character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
  write(*,*) '#### start allocate_replica'
#endif

   flg_error = .false.

   ! check
   if (n_replica_all < 1)   flg_error = .true.
   if (allocated(xyz_mp_rep)) flg_error = .true.
   if (allocated(pxyz_mp_rep)) flg_error = .true.
   if (allocated(xyz_ele_rep)) flg_error = .true.
   if (allocated(pxyz_ele_rep)) flg_error = .true.

   if (flg_error) then
      error_message = 'defect at allocate_replica, PROGRAM STOP'
      call util_error(ERROR%STOP_ALL,error_message)
   endif

   !---------------------
   ! allocate xya_mp_rep
   !---------------------
   error_message = 'failed in memory allocation at allocate_replica, PROGRAM STOP'
   allocate( xyz_mp_rep(SDIM, nmp_all, n_replica_mpi),    stat=ier) 
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)

   allocate( pxyz_mp_rep(SDIM, nmp_all, n_replica_mpi),    stat=ier) 
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)

   ! copy initial structure
   do irep = 1, n_replica_mpi
      do imp = 1, nmp_all
         xyz_mp_rep(1:SDIM,imp,irep) = xyz_mp_init(1:SDIM,imp)
      enddo
   enddo

   if(inele%i_calc_method == 1 .or. inele%i_calc_method == 2) then
      !allocate( xyz_ele_rep(SDIM, ncharge, n_replica_mpi),    stat=ier) 
      !if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     ! 
     ! allocate( pxyz_ele_rep(SDIM, ncharge, n_replica_mpi),    stat=ier) 
     ! if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
      error_message = 'allocate_replica: inele%i_calc_method is invalid, PROGRAM STOP'
      call util_error(ERROR%STOP_ALL, error_message)
   end if

   !---------------------
   ! allocation for REM
   !---------------------
   if (i_run_mode == RUN%REPLICA) then
      allocate( hist_combination(n_replica_all, n_replica_all), stat=ier)
      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
      hist_combination(:,:)   = 0
      
      allocate( hist_exchange(n_replica_all, n_replica_all, inrep%n_period_prob), stat=ier)
      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
      hist_exchange(:,:,:)        = 0

      allocate( hist_attempt(n_replica_all, n_replica_all, inrep%n_period_prob), stat=ier)
      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
      hist_attempt(:,:,:) = 0

      allocate( rate_exchange(n_replica_all, n_replica_all, inrep%n_period_prob), stat=ier)
      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
      rate_exchange(:,:,:) = 0.0e0_PREC

      allocate( up_or_down(REPTYPE%MAX, n_replica_all), stat=ier)
      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
      up_or_down(1:REPTYPE%MAX, 1                ) =  1
      up_or_down(1:REPTYPE%MAX, 2:n_replica_all-1) = -1
      up_or_down(1:REPTYPE%MAX, n_replica_all    ) =  0

      allocate( n_current_up(REPTYPE%MAX, n_replica_all, inrep%n_period_prob), stat=ier)
      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
      n_current_up(:,:,:) = 0

      allocate( n_current_down(REPTYPE%MAX, n_replica_all, inrep%n_period_prob), stat=ier)
      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
      n_current_down(:,:,:) = 0

      allocate( n_turnover(REPTYPE%MAX, n_replica_all, inrep%n_period_prob), stat=ier)
      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
      n_turnover(:,:,:) = 0

   endif

#ifdef _DEBUG
  write(*,*) '#### end allocate_replica'
#endif

endsubroutine allocate_replica

!-----------------------------------------------------------------
subroutine deallocate_replica()
   use var_struct, only :  xyz_mp_rep, pxyz_mp_rep, &
                           xyz_ele_rep, pxyz_ele_rep
   use var_replica, only : hist_combination,    &
                           hist_exchange,       &
                           hist_attempt,&
                           n_current_up,       &
                           n_current_down,     &
                           label_order, n_turnover
   implicit none

   if (allocated(xyz_mp_rep) )         deallocate(xyz_mp_rep) 
   if (allocated(pxyz_mp_rep) )        deallocate(pxyz_mp_rep)
   if (allocated(xyz_ele_rep) )         deallocate(xyz_ele_rep) 
   if (allocated(pxyz_ele_rep) )        deallocate(pxyz_ele_rep)
   if (allocated(hist_combination) )   deallocate(hist_combination) 
   if (allocated(hist_exchange) )      deallocate(hist_exchange)
   if (allocated(hist_attempt))        deallocate(hist_attempt)
   if (allocated(n_current_up))        deallocate(n_current_up)
   if (allocated(n_current_down))      deallocate(n_current_down)
   if (allocated(label_order))         deallocate(label_order)
   if (allocated(n_turnover))          deallocate(n_turnover)

endsubroutine deallocate_replica
