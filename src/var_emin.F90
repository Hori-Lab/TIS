! var_emin
!> @brief Modeul for defining variables for energy minimization

! **************************************************************************
!  variable for energy minimiazation
module var_emin

   use const_maxsize
   implicit none
 
   type input_emin
      integer :: i_method(MXSIM) ! = 1 : Steepest Descent (with Arbitrary Step Approach)
                                 ! = 2 : Conjugate Gradient
                                 !  (with the Polak-Riviere method and the line-search approach)
      integer :: i_out ! Whether to output time-series of minimization parameter (optional)
                       ! = 0 : nothing (default) 
                       ! = 1 : output to .data
                       ! = 2 : output to .opt ("opt" must be specified in the OUTPUT line above)
      integer :: n_out_step ! How often information is saved in output files
 
      real(PREC) :: eps  ! Criterion of convergence (common to both methods)
 
      !### FOR THE STEEPEST DESCENT METHOD ###
      real(PREC) :: sd_lambda_init ! The initial value of the maximum displacement
      real(PREC) :: sd_rho_accept  ! Change rates of lambda in each step
      real(PREC) :: sd_rho_reject
      
      !### FOR THE CONJUGATE GRADIENT METHOD ###
      real(PREC) :: cg_lambda_init ! The initial value of the maximum displacement
      real(PREC) :: cg_rho ! Change rate of lambda when rejected
      real(PREC) :: cg_wolfe_c1 ! Coefficients of strong Wolfe condition
      real(PREC) :: cg_wolfe_c2 ! (0.0 < c1 < c2 < 0.5 should be satisfied)
 
      integer    :: sz
   endtype input_emin
 
   type(input_emin), save :: inemin

   integer, save :: lunout
   real(PREC), save :: lambda
   real(PREC), save :: norm_max  ! Maximun norm
   real(PREC), allocatable, save :: pvec(:,:)

! ###########################################################################
contains

   logical function func_check_norm()
#ifdef MPI_PAR
      use mpiconst
#endif
      implicit none

      func_check_norm = .FALSE.
      if (norm_max <= inemin%eps) then
#ifdef MPI_PAR
         if (local_rank_mpi == 0) then
#endif
         write(*,'(a)') '### norm reached under the epsilon value. minimization ends.'
         write(*,'(a8,g22.15,a7,g22.15)') '## norm=',norm_max,'lambda=',lambda
         if (inemin%i_out > 0) then
            write(lunout, '(a)') '### norm reached under the epsilon value. minimization ends.'
            write(lunout, '(a8,g22.15,a7,g22.15)') '## norm=',norm_max,'lambda=',lambda
         endif
#ifdef MPI_PAR
         endif
#endif
         func_check_norm = .TRUE.
      endif
   endfunction func_check_norm


   logical function func_check_lambda()
#ifdef MPI_PAR
      use mpiconst
#endif
      implicit none

      func_check_lambda = .FALSE.
      if (lambda <= inemin%eps) then
#ifdef MPI_PAR
         if (local_rank_mpi == 0) then
#endif
         write(*,'(a)') '### lambda reached under the epsilon value. minimization ends.'
         write(*,'(a8,g22.15,a7,g22.15)') '## norm=',norm_max,'lambda=',lambda
         if (inemin%i_out > 0) then
            write(lunout, '(a)') '### lambda reached under the epsilon value. minimization ends.'
            write(lunout, '(a8,g22.15,a7,g22.15)') '## norm=',norm_max,'lambda=',lambda
         endif
#ifdef MPI_PAR
         endif
#endif
         func_check_lambda = .TRUE.
      endif
   endfunction func_check_lambda


end module var_emin
