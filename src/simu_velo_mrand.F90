! simu_velo_mrand
!> @brief This subroutine is to generate the initial velocities &
!>        for all of mass-points under specified temperature.

! ***********************************************************************
! Maxwell distribution
! ***********************************************************************
subroutine simu_velo_mrand(tempk_in)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : fix_mp
  use var_struct,  only : nmp_real, cmass_mp
  use var_replica, only : flg_rep, rep2val, irep2grep, &
                          n_replica_all, n_replica_mpi
  use var_simu, only : velo_mp
  implicit none

  ! --------------------------------------------------------------------
  real(PREC), intent(in) :: tempk_in

  ! --------------------------------------------------------------------
  ! function
  real(PREC) :: rfunc_boxmuller

  ! --------------------------------------------------------------------
  ! local variables
  integer :: irep, grep, imp, idimn, istream
  real(PREC) :: coef
  real(PREC) :: tempk
  real(PREC) :: r_boxmuller(SDIM, nmp_real, n_replica_all)

  ! --------------------------------------------------------------------
#ifdef _DEBUG
  write(*,*) 'start simu_velo_mrand : START'
#endif

  tempk = tempk_in

  do irep = 1, n_replica_mpi
     istream = irep
     do imp= 1, nmp_real
        do idimn = 1, SDIM
           r_boxmuller(idimn, imp, irep) = rfunc_boxmuller(istream, 0)
        enddo
     enddo
  enddo

  !  do irep = 1, n_replica_all
  do irep = 1, n_replica_mpi

     grep = irep2grep(irep)
     if (flg_rep(REPTYPE%TEMP)) then
        tempk = rep2val(grep, REPTYPE%TEMP)
     endif 

     do imp = 1, nmp_real
        if(fix_mp(imp)) then
           velo_mp(1:3, imp, irep) = 0.0
        else
           coef = sqrt(tempk * BOLTZ_KCAL_MOL / cmass_mp(imp))
           do idimn = 1, 3
              velo_mp(idimn, imp, irep) = coef * r_boxmuller(idimn, imp, irep)
           end do
        end if
     end do
  enddo

#ifdef _DEBUG
  write(6,*) 'simu_velo_mrand : END'
#endif

end subroutine simu_velo_mrand
