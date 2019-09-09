subroutine widom

   use const_physical
   use const_index
   use var_io,    only : outfile
   use var_setp,  only : inwidom, mts, inele, inperi
   use var_struct,only : ntp, xyz_tp
   use var_simu,  only : istep, tempk, widom_iw, widom_chp, widom_count_exv_inf, widom_energy
   use mt_stream
   use mpiconst

   implicit none

   integer, parameter :: irep = 1
   integer :: iw
   integer :: tp1
   real(PREC) :: random(SDIM)
   real(PREC) :: test_total, chp
   real(PREC) :: kT
   character(CARRAY_MSG_ERROR) :: error_message

interface
   subroutine energy_exv_dt15_tp(irep, energy)
      use const_maxsize
      use const_index
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(inout) :: energy(E_TYPE%MAX)
   endsubroutine energy_exv_dt15_tp
   subroutine energy_ele_coulomb_tp(irep, energy)
      use const_maxsize
      use const_index
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: energy(E_TYPE%MAX)
   endsubroutine energy_ele_coulomb_tp
   subroutine energy_ele_coulomb_ewld_tp(irep, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: energy(:)
   endsubroutine energy_ele_coulomb_ewld_tp
endinterface

   kT = BOLTZ_KCAL_MOL * tempk

   do iw = 1, inwidom%n_trial

      widom_iw = widom_iw + 1

      ! Generate coordinates
      do tp1 = 1, ntp
         random(1) = genrand_double1(mts(0,0))
         random(2) = genrand_double1(mts(0,0))
         random(3) = genrand_double1(mts(0,0))
         xyz_tp(1:3, tp1) = (random(1:3) - 0.5e0_PREC) * inperi%psize(1:3)
      enddo

      widom_count_exv_inf = 0
   
      ! Call energy routines one by one
      widom_energy(:) = 0.0e0_PREC

      call energy_exv_dt15_tp(irep, widom_energy(:))
      if (widom_count_exv_inf > 0) then
         write(outfile%chp(irep),'(i15,1x, f10.4,1x, i10,1x, i5,1x,a9,1x,a9,1x,a9)') &
                  widom_iw, -kT*log(widom_chp/real(widom_iw,kind=PREC)), istep, iw, 'inf', 'inf', 'NaN'
         cycle
      endif

      if (inele%i_function_form == 1) then ! Coulomb potential
         call energy_ele_coulomb_tp(irep, widom_energy(:))
      elseif (inele%i_function_form == 2) then ! Coulomb potential
         call energy_ele_coulomb_ewld_tp(irep, widom_energy(:))
      else
         error_message = 'Error in widom.F90'
         call util_error(ERROR%STOP_ALL, error_message)
      endif

      test_total = widom_energy(E_TYPE%EXV_DT15) + widom_energy(E_TYPE%ELE)

      ! Calculate chemical potential
      chp = exp( - test_total / kT)
      widom_chp =  widom_chp + chp

      ! Output
      write(outfile%chp(irep),'(i15,1x, f10.4,1x, i10,1x, i5,1x,g11.4,1x,g11.4,1x,g11.4)') &
                  widom_iw, -kT*log(widom_chp/real(widom_iw,kind=PREC)), istep, iw, &
                  test_total, widom_energy(E_TYPE%EXV_DT15), widom_energy(E_TYPE%ELE)
   enddo

endsubroutine widom
