subroutine widom

   use const_physical
   use const_index
   use var_inp,    only : inperi, outfile
   use var_setp,   only : inwidom, mts
   use var_struct, only : ntp, xyz_tp, iclass_tp, charge_tp, &
                          tp_exv_dt15_eps, tp_exv_dt15_rad
   use var_simu,   only : istep, tempk, &
                          widom_iw, widom_energy, widom_chp
   use mt_stream

   implicit none

   integer, parameter :: irep = 1
   integer :: iw
   integer :: tp1, tp2
   real(PREC) :: random(SDIM)
   real(PREC) :: energy_test(E_TYPE%MAX), test_total, chp, delta 
   real(PREC) :: kT

   kT = BOLTZC * tempk
   
   do iw = 1, inwidom%n_trial

      ! Generate coordinates
      do tp1 = 1, ntp
         random(1) = genrand_double1(mts(0,0))
         random(2) = genrand_double1(mts(0,0))
         random(3) = genrand_double1(mts(0,0))
         xyz_tp(1:3, tp1) = (random(1:3) - 0.5e0_PREC) * inperi%psize(1:3)
      enddo
   
      ! Calculate energy
      energy_test(:) = 0.0e0_PREC
      call energy_exv_dt15_tp(irep, energy_test)
      call energy_ele_coulomb_tp(irep, energy_test)
      test_total = energy_test(E_TYPE%EXV_DT15) + energy_test(E_TYPE%ELE)

      ! Calculate chemical potential
      !delta = test_total - e_total
      !chp = exp( - delta / kT)
      widom_iw = widom_iw + 1
      chp = exp( - test_total / kT)
      widom_chp =  widom_chp + chp

      ! Output
      write(outfile%chp,'(i10,1x, f10.4,1x, i10,1x, i10,1x,f9.4,1x,f9.4,1x,f9.4)') &
                  widom_iw, -kT*log(widom_chp/real(widom_iw)), istep, iw, &
                  test_total, energy_test(E_TYPE%EXV_DT15), energy_test(E_TYPE%ELE)
   enddo

endsubroutine widom
