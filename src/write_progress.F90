subroutine write_progress(istep, ntstep )

   use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT 
   use const_maxsize
   use var_setp, only : insimu

   implicit none
   
   integer(L_INT), intent(in) :: istep
   integer(L_INT), intent(in) :: ntstep

   real(PREC), parameter :: sec_in_day = 60*60*24.0
   logical, save :: flg_first = .true.
   integer(L_INT) :: tt, t_rate
   integer(L_INT), save :: step_pre, step_ini
   real(PREC), save :: fn_step_progress
   real(PREC), save :: t_pre, t_ini
   real(PREC) :: mstep_per_day_pre, mstep_per_day_ini
   real(PREC) :: days_from_pre, days_from_ini, remain_hours
   real(PREC) :: t

   if (flg_first) then
      write(OUTPUT_UNIT,'(a,i10,a)') '# PROGRESS shwon every ', insimu%n_step_progress, ' steps'
      write(OUTPUT_UNIT,'(a)') '#1: step'
      write(OUTPUT_UNIT,'(a)') '#2: million steps / day (latest cycle)'
      write(OUTPUT_UNIT,'(a)') '#3: million steps / day (averaged)'
      write(OUTPUT_UNIT,'(a)') '#4: remaining time / hours'

      fn_step_progress = real(insimu%n_step_progress, kind=PREC)

      call system_clock(tt,t_rate)
      t_pre = tt / real(t_rate, kind=PREC)
      t_ini = t_pre
      step_ini = istep

      flg_first = .false.
      return
   endif

   call system_clock(tt,t_rate)
   t = tt / real(t_rate, kind=PREC)

   days_from_ini = (t - t_ini) / sec_in_day
   days_from_pre = (t - t_pre) / sec_in_day

   mstep_per_day_ini = real(istep - step_ini, kind=PREC) / fn_step_progress / days_from_ini
   mstep_per_day_pre = real(istep - step_pre, kind=PREC) / fn_step_progress / days_from_pre

   remain_hours = real(ntstep - istep, kind=PREC) / fn_step_progress / mstep_per_day_ini * 24.0

   write(OUTPUT_UNIT,'(i12,1x,g13.6,1x,g13.6,1x,f6.1)') istep, mstep_per_day_pre, mstep_per_day_ini, remain_hours
   flush(OUTPUT_UNIT)

   t_pre = t
   step_pre = istep

endsubroutine write_progress
