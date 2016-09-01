!force_sumup
!> @brief Subroutine for force calculation

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! ************************************************************************
! subroutine for the force
! ************************************************************************
subroutine force_sumup(force_mp, &  ! [ o]
                       irep      &  ! [i ]
                       )
            
  use const_maxsize
  use const_physical
  use const_index
  use var_struct, only : nmp_all  !, nunit_all
  use mpiconst

  implicit none

  real(PREC), intent(out) :: force_mp(SDIM, nmp_all)
  integer,    intent(in)  :: irep

  integer    :: tn, n
  real(PREC) :: force_mp_l(SDIM, nmp_all, 0:nthreads-1)

!$omp parallel private(tn)
  tn = 0
!$  tn = omp_get_thread_num()

  force_mp_l(1:SDIM,1:nmp_all  ,tn) = 0.0_PREC

  call force_bond  (irep, force_mp_l(1,1,tn))

  call force_bangle(irep, force_mp_l(1,1,tn))

        call force_dtrna_stack_nlocal(irep, force_mp_l(1,1,tn))
        call force_dtrna_stack(irep, force_mp_l(1,1,tn))
        call force_dtrna_hbond15(irep, force_mp_l(1,1,tn))

        call force_exv_dt15 (irep, force_mp_l(1,1,tn))

        call force_ele_coulomb(irep, force_mp_l(1,1,tn))

!$omp end parallel

!  do n = 1, nthreads-1
!    force_mp_l(1:SDIM,1:nmp_all,0) = &
!    force_mp_l(1:SDIM,1:nmp_all,0) + &
!    force_mp_l(1:SDIM,1:nmp_all,n)
!  end do
!
!  force_mp(1:SDIM,1:nmp_all) = force_mp_l(1:SDIM,1:nmp_all,0) 

  force_mp(1:SDIM, 1:nmp_all) = force_mp_l(1:SDIM,1:nmp_all,0)
  do n = 1, nthreads-1
    force_mp(1:SDIM,1:nmp_all) = &
    force_mp(1:SDIM,1:nmp_all) + force_mp_l(1:SDIM,1:nmp_all,n)
  end do

end subroutine force_sumup
