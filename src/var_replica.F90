! var_replica
!> @brief Module for defining variables of replica exchange method

module var_replica

  use const_maxsize
  use const_physical
  use const_index

  implicit none

  !----------------------------------------------
  type input_parameter
     integer :: n_step_replica !< # of interval steps of replica exchange
     integer :: n_step_save    !< # of interval steps of saving .rep file
     integer :: npar_rep       !< parameter of parallelization
     integer :: n_replica(REPTYPE%MAX) !< # of replicas
     integer :: i_style(REPTYPE%MAX)   !< how to define replica variables
     integer :: i_loadbalance          !< flag for load-balancing
     integer :: n_step_adjust          !< (not used in current version)
     integer :: n_adjust_interval      !< # of interval steps of load balancing
     integer :: n_pull
     real(PREC) :: lowest(REPTYPE%MAX)   !< lowest value of replica variable
     real(PREC) :: highest(REPTYPE%MAX)  !< highest value of replica variable
     real(PREC) :: exponent_alpha(REPTYPE%MAX)
     !real(PREC) :: interval(REPTYPE%MAX) !< interval value of replica variable
     real(PREC) :: var(MXREPLICA,REPTYPE%MAX)  !< replica variables defined explicitly
     real(PREC) :: window_property(MXREPLICA,WINDTYPE%MAX)  !< window property
     real(PREC) :: window_mp_id(MXREPLICA,2)  !< ids of mass point for window
     integer    :: winz_igrp(MXREPLICA)
     real(PREC) :: winz_z(MXREPLICA)
     real(PREC) :: winz_kxy(MXREPLICA)
     real(PREC) :: winz_kz(MXREPLICA) 
     real(PREC) :: pull_direction(SPACE_DIM, MXPULLING)
     logical :: flg_exchange   !< flag for exchange

     ! for FO-REM
     logical :: flg_opt_temp   !< flag for feedback-optimized REM
     integer(L_INT) :: n_step_opt_temp
     integer :: n_stage_opt_temp

     ! for detail of exchange probability
     integer :: n_period_prob  ! (default = 2)

     integer :: sz    !< own struct size
  end type input_parameter
  type(input_parameter), save :: inrep

  !----------------------------------------------
  logical :: flg_npar_rep
  integer, save :: n_dimension    = -1   !< # of dimensions of replica variable
  integer, save :: n_replica_all  = -1   !< # of replicas
  integer, save :: n_replica_mpi         !< # of replicas per parallelization node
  integer, save :: irep2grep(MXREPLICA)  !< replica(local) => replica(global)
  integer, save :: grep2irep(MXREPLICA)  !< replica(global) => replica(local)
  integer, save :: grep2rank(MXREPLICA)  !< replica(global) => local_rank_rep
  logical, save :: flg_rep(1:REPTYPE%MAX) = .false.
  
  !----------------------------------------------
  ! for load balance
  integer, save :: step_ratio = 20
  integer(L_INT), save :: lab2step(MXREPLICA)  !< # of exchange step for each replica variable (label)
  integer(L_INT), save :: exchange_step(MXREPLICA)!< # of exchange step for each replica

  !----------------------------------------------
  ! Histogram
  integer, save              :: total_attempt = 0     !< Total # of times of exchange attempt
  integer, allocatable, save :: hist_combination(:,:) !< Histogram of existence (label, replica)
  integer, allocatable, save :: hist_exchange(:,:,:)  !< Total # of times of accepted exchange for each (label, period)
  integer, allocatable, save :: hist_attempt(:,:,:) !< Total # of times of exchange attempt for each (label, period)
  real(PREC), allocatable, save :: rate_exchange(:,:,:)  !< (label, label, period)

  !----------------------------------------------
  ! for Feedback-optimized REM
  integer, allocatable, save :: up_or_down(:,:)         !(REPTYPE%MAX, replica)
  integer, allocatable, save :: n_current_up(:,:,:)    !(REPTYPE%MAX, label, period)
  integer, allocatable, save :: n_current_down(:,:,:)  !(REPTEYP%MAX, label, period)
  integer, allocatable, save :: n_turnover(:,:,:)      !(REPTYPE%MAX, repilca, period)

  ! ================================================
  ! Example of parallelization of replica, 
  !    #node=3, #replica=6
  !    -----------------------------------------
  !     Node | replica(global) | replica(local)
  !          |     irep        |    grep
  !    -----------------------------------------
  !      1   |       1         |      1
  !          |       2         |      2
  !      2   |       3         |      1
  !          |       4         |      2
  !      3   |       5         |      1
  !          |       6         |      2
  !    -----------------------------------------
  !     max  | inrep%n_replica | n_replica_mpi
  ! ================================================

  ! Permutation function
  integer,    save :: rep2lab(MXREPLICA)  !< Permutation function (replica => label)
  integer,    save :: lab2rep(MXREPLICA)  !< Permutation function ( label  => replica)
  real(PREC), save :: lab2val(MXREPLICA, REPTYPE%MAX) !< Replica variables

  integer, allocatable, save :: label_order(:,:)  !< (REPTYPE%MAX, n_repica_all) 

  ! ===============================================================
  ! Using these 'Permutation function'
  ! Ex.)
  ! replica => label  (rep2lab)
  ! label => replica  (lab2rep)
  ! --------------------------------
  ! | replica  | 1 | 2 | 3 | 4 | 5 |
  ! | label    | 3 | 1 | 5 | 2 | 4 |
  ! --------------------------------
  ! 
  ! label   => value  (lab2val)
  ! --------------------------------
  ! | label    | 1 | 2 | 3 | 4 | 5 |
  ! | value    |1.0|1.2|1.4|1.6|1.8| (temperature, ionic strength, and so on.)
  ! --------------------------------
  !
  ! As a result, in this situation, 
  ! replica => value  (rep2val (function))
  ! --------------------------------
  ! | replica  | 1 | 2 | 3 | 4 | 5 |
  ! | value    |1.4|1.0|1.8|1.2|1.6| (temperature, ionic strength, and so on.)
  ! --------------------------------
  ! ===============================================================

  !-----------------------------------------------------------
  ! to detect whether current exchange step is odd or even
  ! (These are private variables)
  integer, save, private :: type_array(1:MXREPDIM)
  integer, save, private :: exchange_pair_tb(1:MXREPLICA,1:MXREPDIM*2)
  integer, save, private :: idx_type = 1
  integer, save, private :: idx_pair = 1

! ###########################################################################
contains

  !> rep2val
  !> @brief Permutation function (replica => variable)
  real(PREC) function rep2val (ireplica, vartype)
     implicit none
     integer :: ireplica
     integer :: vartype
     rep2val = lab2val(rep2lab(ireplica), vartype)
  endfunction rep2val


  !> rep2step
  integer function rep2step (ireplica)
     implicit none
     integer :: ireplica
     rep2step = lab2step(rep2lab(ireplica))
  endfunction rep2step


  !> get_pair
  integer function get_pair(i)
     implicit none
     integer, intent(in) :: i
     get_pair = exchange_pair_tb(i, idx_pair)
  endfunction get_pair


  !> set_forward
  !> @brief
  !> Called by simu_replica_exchange 
  !> when each exchange event is finished.
  subroutine set_forward()
     implicit none

     if (idx_type == n_dimension) then
        idx_type = 1
     else
        idx_type = idx_type + 1
     endif

     if (idx_pair == n_dimension*2) then
        idx_pair = 1
     else
        idx_pair = idx_pair + 1
     endif
  endsubroutine set_forward


  !> get_type
  !> @brief
  integer function get_type(i)
     implicit none
     integer, intent(in), optional :: i
     if (present(i)) then
        get_type = type_array(i)
     else
        get_type = type_array(idx_type)
     endif
  endfunction get_type


  !> make_type_array
  !> @brief
  subroutine make_type_array()
     ! called by replica_settable only once
     use const_index
     use const_maxsize
     implicit none
     integer :: idx, ivar

     type_array(:) = 0
     idx = 0
     do ivar = 1, REPTYPE%MAX
        if (flg_rep(ivar)) then
           idx = idx + 1
           if (idx > MXREPDIM) then
              call util_error(ERROR%STOP_ALL, &
              'Error: defect in var_replica::make_type_table')
           endif
           type_array(idx) = ivar
        endif
     enddo
  endsubroutine make_type_array


  !> make_exchange_pair_tb
  !> @brief
  subroutine make_exchange_pair_tb()
     implicit none
     integer :: icounter, ivar, iset
     integer :: idimn, iparity, icycle, ireplica, icontinue
     integer :: ireplica_start
     integer :: n_division_pre, n_division, n_continue
     logical :: flg_post_zero
#ifdef _DEBUG
     integer :: irep
#endif

     exchange_pair_tb(:,:) = 0

     ! ==========================================================
     ! exchange order  (Ex. 3-dimension case)
     !   -----------------------------------------------------
     !    icounter  |  1   |  2   |  3   |  4   |  5   |  6
     !   -----------------------------------------------------
     !    exchange- | dim1 | dim2 | dim3 | dim1 | dim2 | dim3 
     !      pair    | odd  | odd  | odd  | even | even | even
     !   -----------------------------------------------------
     ! ==========================================================

     icounter       = 0
     do iparity = 1, 2   ! odd or even
                         ! 1 odd  : exchange 1-2, 3-4, 5-6 ,.....
                         ! 2 even : exchange 2-3, 4-5, 6-7 ,.....

        n_division     = 1
        do idimn = 1, n_dimension 

           icounter = icounter + 1
           ivar = type_array(idimn)
           n_division_pre = n_division
           n_division     = n_division * inrep%n_replica(ivar) 
           n_continue     = n_replica_all / n_division

           iset = 0
           do icycle = 1, n_division_pre

              flg_post_zero = .false.
              if (iparity == 1) then
                 if (mod(inrep%n_replica(ivar),2) == 1) then 
                    ! #replica is odd
                    flg_post_zero = .true.
                 endif
                 ireplica_start = 1
              else
                 ! zeroing 1st replica
                 do icontinue = 1, n_continue
                    iset = iset + 1
                    exchange_pair_tb(iset, icounter) = 0
                 enddo
                 if (mod(inrep%n_replica(ivar),2) == 0) then
                    ! #replica is even
                    flg_post_zero = .true.
                  endif
                 ireplica_start = 2
              endif

              do ireplica = ireplica_start, (inrep%n_replica(ivar)-1), 2

                 do icontinue = 1, n_continue*2
                    iset = iset + 1

                    if (icontinue <= n_continue) then
                       exchange_pair_tb(iset, icounter) = iset + n_continue
                    else
                       exchange_pair_tb(iset, icounter) = iset - n_continue
                    endif
                 enddo

              enddo ! ireplica

              ! zeroing last replica
              if (flg_post_zero) then
                 do icontinue = 1, n_continue
                    iset = iset + 1
                    exchange_pair_tb(iset, icounter) = 0
                 enddo
              endif

           enddo ! icycle
        enddo ! iparity
     enddo ! idimn

#ifdef _DEBUG
     do icounter = 1, MXREPDIM*2
        write(*,*) '#################'
        write(*,*) 'icounter = ',icounter
        do irep = 1, n_replica_all
           write(*,*) 'exchange_pair_tb(',irep,',',icounter,')=', exchange_pair_tb(irep, icounter)
        enddo
     enddo
#endif

  endsubroutine make_exchange_pair_tb

endmodule var_replica
