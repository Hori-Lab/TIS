subroutine write_progress(istep, ntstep )

   use const_physical
   use var_setp, only : insimu

   implicit none
   
   integer(L_INT), intent(in) :: istep
   integer(L_INT), intent(in) :: ntstep

   real(PREC), parameter :: sec_in_day = 60*60*24.0
   logical, save :: flg_first = .true.
   integer :: tt, t_rate
   integer(L_INT), save :: step_pre, step_ini
   real(PREC), save :: t_pre, t_ini
   real(PREC) :: mstep_per_day_pre, mstep_per_day_ini
   real(PREC) :: days_from_pre, days_from_ini, remain_hours
   real(PREC) :: t

   if (flg_first) then
      write(6,'(a,i10,a)') '# PROGRESS shwon every ', insimu%n_step_progress, ' steps'
      write(6,'(a)') '#1: step'
      write(6,'(a)') '#2: million steps / day (latest cycle)'
      write(6,'(a)') '#3: million steps / day (averaged)'
      write(6,'(a)') '#4: remaining time / hours'

      call system_clock(tt,t_rate)
      t_pre = tt / dble(t_rate)
      t_ini = t_pre
      step_ini = istep

      flg_first = .false.
      return
   endif

   call system_clock(tt,t_rate)
   t = tt / dble(t_rate)

   days_from_ini = (t - t_ini) / sec_in_day
   days_from_pre = (t - t_pre) / sec_in_day

   mstep_per_day_ini = real(istep - step_ini, kind=PREC) / 1000000.0 / days_from_ini
   mstep_per_day_pre = real(istep - step_pre, kind=PREC) / 1000000.0 / days_from_pre

   remain_hours = real(ntstep - istep, kind=PREC) / 1000000.0 / mstep_per_day_ini * 24.0

   write(6,'(i12,1x,g12.6,1x,g12.6,1x,f6.1)') istep, mstep_per_day_pre, mstep_per_day_ini, remain_hours

   t_pre = t
   step_pre = istep

endsubroutine write_progress
