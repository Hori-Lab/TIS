subroutine write_progress(istep, ntstep )

   use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT 
   use const_maxsize
   use var_setp, only : insimu

   implicit none
   
   integer(L_INT), intent(in) :: istep
   integer(L_INT), intent(in) :: ntstep

   real(PREC), parameter :: sec_in_day = 60*60*24.0
   logical, save :: flg_first = .true.
   integer(L_INT) :: clock
   integer(L_INT), save :: step_pre, step_ini, t_rate
   real(PREC), parameter :: fn_million = 1000000.0_PREC
   !real(PREC), save :: t_pre, t_ini
   integer(L_INT), save :: clock_pre, clock_ini
   real(PREC), save :: rt_rate
   real(PREC) :: mstep_per_day_pre, mstep_per_day_ini
   real(PREC) :: days_from_pre, days_from_ini, remain_hours

   if (flg_first) then
      write(OUTPUT_UNIT,'(a,i10,a)') '# PROGRESS shwon every ', insimu%n_step_progress, ' steps'
      write(OUTPUT_UNIT,'(a)') '#1: step'
      write(OUTPUT_UNIT,'(a)') '#2: million steps / day (latest cycle)'
      write(OUTPUT_UNIT,'(a)') '#3: million steps / day (averaged)'
      write(OUTPUT_UNIT,'(a)') '#4: remaining time / hours'

      call system_clock(clock,t_rate)
      clock_pre = clock
      clock_ini = clock
      step_ini = istep
      step_pre = istep
      rt_rate = 1.0_PREC / real(t_rate, kind=PREC)

      flg_first = .false.
      return
   endif

   call system_clock(clock)

   days_from_ini = rt_rate * (clock - clock_ini) / sec_in_day
   days_from_pre = rt_rate * (clock - clock_pre) / sec_in_day

   mstep_per_day_ini = real(istep - step_ini, kind=PREC) / fn_million / days_from_ini
   mstep_per_day_pre = real(istep - step_pre, kind=PREC) / fn_million / days_from_pre

   remain_hours = real(ntstep - istep, kind=PREC) / fn_million / mstep_per_day_ini * 24.0

   write(OUTPUT_UNIT,'(i12,1x,g13.6,1x,g13.6,1x,f6.1)') istep, mstep_per_day_pre, mstep_per_day_ini, remain_hours
   flush(OUTPUT_UNIT)

   clock_pre = clock
   step_pre = istep

endsubroutine write_progress
