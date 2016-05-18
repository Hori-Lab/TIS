subroutine widom

   use const_physical
   use const_index
   use var_io,    only : outfile
   use var_setp,  only : inwidom, mts, inele, inperi
   use var_struct,only : ntp, xyz_tp
   use var_simu,  only : istep, tempk, widom_iw, widom_chp, widom_flg_exv_inf
   use mt_stream

   implicit none

   integer, parameter :: irep = 1
   integer :: iw
   integer :: tp1
   real(PREC) :: random(SDIM)
   real(PREC) :: energy_test(E_TYPE%MAX), test_total, chp
   real(PREC) :: kT

interface
   subroutine energy_exv_dt15_tp(irep, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: energy(:)
   endsubroutine energy_exv_dt15_tp
   subroutine energy_ele_coulomb_tp(irep, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: energy(:)
   endsubroutine energy_ele_coulomb_tp
   subroutine energy_ele_coulomb_ewld_tp(irep, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: energy(:)
   endsubroutine energy_ele_coulomb_ewld_tp
endinterface

   kT = BOLTZC * tempk
   
   do iw = 1, inwidom%n_trial

      widom_iw = widom_iw + 1

      ! Generate coordinates
      do tp1 = 1, ntp
         random(1) = genrand_double1(mts(0,0))
         random(2) = genrand_double1(mts(0,0))
         random(3) = genrand_double1(mts(0,0))
         xyz_tp(1:3, tp1) = (random(1:3) - 0.5e0_PREC) * inperi%psize(1:3)
      enddo
   
      ! Calculate energy
      energy_test(:) = 0.0e0_PREC

      widom_flg_exv_inf = .False.
      call energy_exv_dt15_tp(irep, energy_test)
      if (widom_flg_exv_inf) then
         write(outfile%chp,'(i15,1x, f10.4,1x, i10,1x, i5,1x,a9,1x,a9,1x,a9)') &
                  widom_iw, -kT*log(widom_chp/real(widom_iw)), istep, iw, 'inf', 'inf', 'NaN'
         cycle
      endif

      if (inele%i_function_form == 1) then ! Coulomb potential
         call energy_ele_coulomb_tp(irep, energy_test)
      elseif (inele%i_function_form == 2) then ! Coulomb potential
         call energy_ele_coulomb_ewld_tp(irep, energy_test)
      else
         call util_error(ERROR%STOP_ALL, 'Error in widom.F90')
      endif

      test_total = energy_test(E_TYPE%EXV_DT15) + energy_test(E_TYPE%ELE)

      ! Calculate chemical potential
      chp = exp( - test_total / kT)
      widom_chp =  widom_chp + chp

      ! Output
      write(outfile%chp,'(i15,1x, f10.4,1x, i10,1x, i5,1x,g10.4,1x,g10.4,1x,g10.4)') &
                  widom_iw, -kT*log(widom_chp/real(widom_iw)), istep, iw, &
                  test_total, energy_test(E_TYPE%EXV_DT15), energy_test(E_TYPE%ELE)
   enddo

endsubroutine widom
