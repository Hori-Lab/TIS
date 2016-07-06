! simu_replica_exchange
!> @brief Main subroutine for replica exchange methods

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

#ifdef _DUMP_REP_EX
#define DUMP
#else
#undef DUMP
#endif

subroutine simu_replica_exchange(velo_mp, replica_energy, tempk)

   use const_maxsize
   use const_physical
   use const_index
   use var_setp,    only : insimu, mts
   use var_replica, only : n_replica_all,                 &
                           rep2lab, lab2rep, lab2val,     &
                           total_attempt, hist_combination, &
                           hist_attempt, hist_exchange, &
                           rate_exchange, &
                           flg_rep, get_pair,  &
                           set_forward, inrep, &
                           n_current_up, n_current_down, &
                           up_or_down, label_order, &
                           n_turnover
   use mt_stream
   use mpiconst
   use time

   implicit none

   !---------------------------------------------------------------------------
   real(PREC), intent(inout) :: velo_mp(:,:,:)      ! (SDIM, mp, replica)
   real(PREC), intent(in)    :: replica_energy(:,:) ! (2, replica)
   real(PREC), intent(in)    :: tempk

   !---------------------------------------------------------------------------
   ! local
   real(PREC) :: delta  ! energy difference
   real(8)    :: random0(n_replica_all) ! for Metropolis
   integer    :: nrand

   ! exchanged targets are represented by "i" and "j"
   integer    :: i
   integer    :: ivar
   integer    :: rep_i, rep_j
   integer    :: l_i, l_j
   real(PREC) :: temp_m, temp_n
   real(PREC) :: pot_i_m, pot_i_n
   real(PREC) :: pot_j_m, pot_j_n

   integer, save :: iperiod = 0
   integer, save :: nexchange_period
   integer :: ntotal_exchange

   ! ###########################################
   !  Initialize  (only once)
   ! ###########################################
   if (iperiod == 0) then
      iperiod = 1
      if (inrep%flg_opt_temp) then
         ntotal_exchange = int(inrep%n_step_opt_temp / inrep%n_step_replica) ! 500 = 1000 / 10
         nexchange_period = int(ntotal_exchange / inrep%n_period_prob)    ! 50 = 100 / 2
      else
         ntotal_exchange = int(insimu%n_tstep_all / inrep%n_step_replica) ! 1000 = 10000 / 10
         nexchange_period = int(ntotal_exchange / inrep%n_period_prob)    ! 333 = 1000 / 3
      endif
      !write(*,*) 'ntotal_exchange=',ntotal_exchange
      !write(*,*) 'nexchange_period=',nexchange_period
   endif


   ! ###########################################
   !  Generate random numbers
   ! ###########################################
   TIME_S( tm_random)
   nrand = 0
   do i = 1, n_replica_all, 2
      nrand = nrand + 1
      !random0(nrand) = grnd()
      random0(nrand) = genrand_double1(mts(0, 0))
   enddo
   TIME_E( tm_random)
#ifdef MPI_REP
   call MPI_BCAST(random0, nrand, MPI_REAL8, 0, mpi_comm_rep, ierr)
#endif


   ! ###########################################
   !  Main operation (judgement and exchange)
   ! ###########################################
   !  l_i   : label i
   !  l_j   : label j
   !  rep_i : replica i  (a replica which has label i)
   !  rep_j : replica j  (a replica which has label j)
   ! 
   ! PRE-EXCHANGE
   !           label = temperature & potential
   !   rep_i :  l_i  =   temp_m    &  pot_i_m
   !   rep_j :  l_j  =   temp_n    &  pot_j_n
   !
   ! CANDIDATE
   !           label = temperature & potential
   !   rep_i :  l_j  =   temp_n    &  pot_i_n
   !   rep_j :  l_i  =   temp_m    &  pot_j_m
   !
#ifdef DUMP
   write(6,*) ''
#endif
   nrand = 0
   do l_i = 1, n_replica_all

      l_j = get_pair(l_i)

      if (l_j < l_i) then
         cycle
      endif

      nrand = nrand + 1

      rep_i     = lab2rep( l_i)
      rep_j     = lab2rep( l_j)
#ifdef DUMP
      write(6,*) 'exchange label: ', l_i, l_j
      write(6,*) 'exchange replica: ', rep_i, rep_j
#endif

      if (flg_rep(REPTYPE%TEMP)) then
         temp_m = lab2val( l_i ,REPTYPE%TEMP) 
         temp_n = lab2val( l_j ,REPTYPE%TEMP)
         ! Note that, temp_m < temp_n
         !   i.e.      T(rep_i) < T(rep_j)
      else
         temp_m = tempk
         temp_n = tempk
      endif
#ifdef DUMP
      write(6,*) 'temp_m, temp_n = ', temp_m, temp_n
#endif

      pot_i_m = replica_energy(1, rep_i)
      pot_i_n = replica_energy(2, rep_i)
      pot_j_n = replica_energy(1, rep_j)
      pot_j_m = replica_energy(2, rep_j)

      ! Y. Sugita et al. J Chem phys (2000)  eq.[14]
      delta = ( (pot_j_m - pot_i_m) / temp_m    &
               -(pot_j_n - pot_i_n) / temp_n )  &
             / BOLTZC
#ifdef DUMP
      write(6,*) 'pot_i_m, pot_i_n = ',pot_i_m,pot_i_n
      write(6,*) 'pot_j_m, pot_j_n = ',pot_j_m,pot_j_n
      write(6,*) 'delta = ',delta
#endif

      ! Judgment
      if (delta < 0.0e0_PREC) then
         call exchange_permutation()
#ifdef DUMP
         write(6,*) 'Exchange (delta < 0.0)'
#endif
         hist_exchange(l_i, l_j, iperiod) = hist_exchange(l_i, l_j, iperiod) + 1     
         if (flg_rep(REPTYPE%TEMP)) then
            call scale_velo(temp_m, temp_n)
         endif

      else
         if (random0(nrand) <= exp(-delta)) then
            call exchange_permutation()
            hist_exchange(l_i, l_j, iperiod) = hist_exchange(l_i, l_j, iperiod) + 1     
            if (flg_rep(REPTYPE%TEMP)) then
               call scale_velo(temp_m, temp_n)
            endif
#ifdef DUMP
         write(6,*) 'Exchange (random(',random0(nrand),') <= delta'
#endif
         else
#ifdef DUMP
         write(6,*) 'NOT Exchange (random(',random0(nrand),') > delta'
#endif
         endif

         rate_exchange(l_i, l_j, iperiod) = rate_exchange(l_i, l_j, iperiod) - delta

      endif
      hist_attempt(l_i, l_j, iperiod) = hist_attempt(l_i, l_j, iperiod) + 1

   enddo ! l_i 


   ! ###########################################
   !  Aftertreatment
   ! ###########################################
   call set_forward()

   ! Update histogram
   total_attempt = total_attempt + 1
   do l_i = 1, n_replica_all
      rep_i = lab2rep(l_i)
      hist_combination(l_i, rep_i)   &
                    = hist_combination(l_i, rep_i) + 1
   enddo

   ! turn ON(up)/OFF(down) flags
   do l_i = 1, n_replica_all
      do ivar = 1, REPTYPE%MAX
         if (.not. flg_rep(ivar)) then
            cycle
         endif

         if (label_order(ivar, l_i) == 1) then
            if (up_or_down(ivar, lab2rep(l_i)) == 0) then
               ! increment
               n_turnover(ivar, lab2rep(l_i), iperiod) = n_turnover(ivar, lab2rep(l_i), iperiod) + 1
            endif
            up_or_down(ivar, lab2rep(l_i)) = 1  ! up
         else if (label_order(ivar, l_i) == inrep%n_replica(ivar)) then
            up_or_down(ivar, lab2rep(l_i)) = 0  ! down
         endif
      enddo
   enddo
   
   ! check flag state
   do l_i = 1, n_replica_all
      do ivar = 1, REPTYPE%MAX
         if (.not. flg_rep(ivar)) then
            cycle
         endif

         if (up_or_down(ivar, lab2rep(l_i)) == 1) then
            n_current_up(ivar, l_i, iperiod) = n_current_up(ivar, l_i, iperiod) + 1
         else if (up_or_down(ivar, lab2rep(l_i)) == 0) then
            n_current_down(ivar, l_i, iperiod) = n_current_down(ivar, l_i, iperiod) + 1
         endif
      enddo
   enddo

   ! next period
   if (mod(total_attempt, nexchange_period) == 0) then
      if (iperiod < inrep%n_period_prob) then
         iperiod = iperiod + 1
      endif
   endif

!---------------------------------------------------------------------------
contains

   !---------------------------------!
   !-- exchange permutation arrays --!
   !---------------------------------!
   subroutine exchange_permutation()
      rep2lab(rep_i) = l_j
      rep2lab(rep_j) = l_i
   
      lab2rep(l_i) = rep_j
      lab2rep(l_j) = rep_i
   end subroutine exchange_permutation

   !---------------------------------!
   !-- scale velocity ---------------!
   !---------------------------------!
   subroutine scale_velo(temperature_i, temperature_j)
      use var_replica, only : irep2grep, n_replica_mpi
      implicit none
      real(PREC), intent(in) :: temperature_i, temperature_j
      real(PREC) :: scale_factor
      integer    :: irep, grep

      scale_factor = sqrt( temperature_j / temperature_i )
!      write(6,*) '##scale velo'
!      write(6,*) 'Ti = ',temperature_i
!      write(6,*) 'Tj = ',temperature_j
!      write(6,*) 'Tj / Ti = ', temperature_j / temperature_i 
!      write(6,*) 'SQRT(Tj/Ti) = ', sqrt( temperature_j / temperature_i )

      do irep = 1, n_replica_mpi
         grep = irep2grep(irep)

         if (rep_i == grep) then
!            write(6,*) 'scale_velo i pre : ',velo_mp(1,1,irep)
            velo_mp(:,:, irep) = scale_factor * velo_mp(:,:, irep)
!            write(6,*) 'scale_velo i post: ',velo_mp(1,1,irep)
         endif
         if (rep_j == grep) then
!            write(6,*) 'scale_velo j pre : ',velo_mp(1,1,irep)
            velo_mp(:,:, irep) = velo_mp(:,:, irep) / scale_factor
!            write(6,*) 'scale_velo j post: ',velo_mp(1,1,irep)
         endif
      enddo
   end subroutine

end subroutine simu_replica_exchange
#undef DUMP
