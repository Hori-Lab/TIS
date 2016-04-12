! simu_searchingtf
!> @brief This subroutine is to search the folding transition temperature (T_f) automatically.

subroutine simu_searchingtf(istep, ntstep, qscore, tempk)
  
  use const_maxsize
  use const_index
  use var_inp, only : outfile
  use var_setp, only : insear

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer(L_INT), intent(in)    :: istep, ntstep
  real(PREC), intent(inout) :: qscore, tempk

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: lunout
  integer, save :: n_non_native, n_native, n_phasetrans, istate
  real(PREC), save :: tneigh_lower, tneigh_upper
  real(PREC) :: wide
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
  lunout = outfile%data

  if(istep == 0) then
     tneigh_lower = insear%tempk_lower
     tneigh_upper = insear%tempk_upper
     tempk = (tneigh_lower + tneigh_upper) / 2.0

     n_non_native = 0
     n_native = 0
     n_phasetrans = 0

  else if(istep == ntstep) then

#ifdef MPI_PAR
  if (myrank == 0) then
#endif
     write (lunout, "(72('*'))")
     write (lunout, "(a6, 4a8)") &
             'tf_out', 'tempk', 'n_state', 'd_state', 'p_trans'
     write (lunout, "(a6, f8.3, 3i8)") &
          'tf_out', tempk, n_native, n_non_native, n_phasetrans
     write (lunout, "(72('*'))")
     write (lunout, *)''
     write (lunout, *)''
#ifdef MPI_PAR
  end if
#endif
     
     if(n_native >= n_non_native) then
        wide = abs(tneigh_upper - tempk)
        tneigh_lower = tempk
     else
        wide = abs(tneigh_lower - tempk)
        tneigh_upper = tempk
     end if

     tempk = (tneigh_lower + tneigh_upper) / 2.0

     if(wide < 1.0e0_PREC) then
        !  This program stop

!#ifdef MPI_PAR
!  call MPI_Finalize(ierr)
!#endif

        error_message = 'PROGRAM STOP'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     n_native = 0
     n_non_native = 0
     n_phasetrans = 0

  else
     if(istep == 1) then
        if(qscore >= 0.5) then
           istate = 1
           n_native = n_native + 1
        else
           istate = -1
           n_non_native = n_non_native + 1
        end if
     else
        if(qscore >= 0.5) then
           if(istate == -1) then
              n_phasetrans = n_phasetrans + 1
           end if
           istate = 1
           n_native = n_native + 1
        else
           if(istate == 1) then
              n_phasetrans = n_phasetrans + 1
           end if
           istate = -1
           n_non_native = n_non_native + 1
        end if
     end if
     
  end if
      
end subroutine simu_searchingtf
