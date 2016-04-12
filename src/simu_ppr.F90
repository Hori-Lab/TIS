! simu_ppr
! @brief   Control the Push-Pull-Release (PPR) protocol.
subroutine simu_ppr()

  use const_maxsize
  use var_setp, only : inmisc
  use var_simu, only : flg_ppr_release
  use var_struct, only : xyz_mp_rep, grp

  implicit none

  integer :: ip, icom, n
  integer :: igrp1, igrp2, i1, i2, imp1, imp2
  real(PREC) :: d(3), dist2
  real(PREC) :: Rt

  integer, save :: isgn(MXBRIDGE)
  integer, save :: istep(MXBRIDGE)
  integer, save :: ioffset(MXBRIDGE)
  real(PREC), save :: slope(MXBRIDGE)
  real(PREC), save :: Rref(MXBRIDGE)
  logical, save :: flg_first = .true.

  integer, parameter :: IREP=1  !! Replica is not supported in PPR.

  !           Push-Pull cycle
  !      Rt
  !      |
  ! Rmax | _ _ _ _ _ _ _ _ _ _ _ _ 
  !      |\P        l/\P        l/
  !      | \u      l/  \u      l/
  !      |  \s    u/    \s    u/
  ! Rcut |_ _\h_ P/_ _ _ \h_ P/_ _
  !      |    \  /        \  /
  ! Rmin | _ _ \/_ _ _ _ _ \/_ _ _
  !      |                        ----> step
  !       <--cycle--> <--cycle-->
  !
  !   Push: Rt =   Rmax   +   (-1)   *   slope   * (  step    -    0      ))
  !   Pull: Rt =   Rmin   +   (+1)   *   slope   * (  step    -  cycle/2  ))
  !   Code: Rt = Rref(ip) + isgn(ip) * slope(ip) * (istep(ip) - ioffset(ip))
  !
  !  * Slope = (Rmax - Rmin) / (cycle / 2)
  !  * step is reset to zero in each cycle.
  !  * When Rt > Rcut, the reference of PPR potential is equal to Rt.
  !    On the other hand, when Rt <= Rcut, it is equal to Rcut. (flattened)
  !
  !  * inmisc%ibrid_ppr_opt
  !    = 1: The reference is not flattened at Rcut, and thus always equal to Rt.

  if (flg_first) then
     n = inmisc%nbrid_ppr
     slope(1:n) = 2.0*(inmisc%brid_ppr_rmax(1:n) - inmisc%brid_ppr_rmin(1:n)) &
                 / real(inmisc%ibrid_ppr_cycl(1:n))
     istep(1:n) = 0
     ioffset(1:n) = 0
     isgn(1:n) = -1
     Rref(1:n) = inmisc%brid_ppr_rmax(1:n)
     flg_first = .false.
  endif

  do ip = 1, inmisc%nbrid_ppr

     if (isgn(ip) == -1 .and. istep(ip) >= inmisc%ibrid_ppr_cycl(ip) / 2) then
        isgn(ip) = 1
        Rref(ip) = inmisc%brid_ppr_rmin(ip)
        ioffset(ip) = istep(ip)
   
     else if (istep(ip) == inmisc%ibrid_ppr_cycl(ip)) then
        isgn(ip) = -1
        Rref(ip) = inmisc%brid_ppr_rmax(ip)
        istep(ip) = 0
        ioffset(ip) = 0

     endif
     
     Rt = Rref(ip) + isgn(ip) * slope(ip) *  (istep(ip) - ioffset(ip))
     icom = inmisc%ibrid_ppr_com(ip)
     inmisc%brid_com_dist(icom) = Rt

     flg_ppr_release(icom) = .false.
     if (Rt <= inmisc%brid_ppr_rcut(ip)) then
     
        if (inmisc%ibrid_ppr_opt(ip) /= 1) then
           inmisc%brid_com_dist(icom) = inmisc%brid_ppr_rcut(ip)
        endif

        igrp1 = inmisc%ibrid_ppr_gid_r(1,ip)
        igrp2 = inmisc%ibrid_ppr_gid_r(2,ip)

#ifndef _OPENMP
        loop_dist: do i1 = 1, grp%nmp(igrp1)
           imp1 = grp%implist(i1, igrp1)
           do i2 = 1, grp%nmp(igrp2)
              imp2 = grp%implist(i2, igrp2)
              d(:) = xyz_mp_rep(:,imp1,IREP) - xyz_mp_rep(:,imp2,IREP)
              dist2 = dot_product(d,d)
              if (dist2 < inmisc%brid_ppr_rzero2(ip)) then
                 flg_ppr_release(icom) = .true.
                 exit loop_dist
              endif
           enddo
        enddo loop_dist
#else
        !! NOTE: When using openmp, "exit" can not be used in parallelized loop.
        do i1 = 1, grp%nmp(igrp1)
           imp1 = grp%implist(i1, igrp1)
!$omp parallel do private(i2,imp2,d,dist2)
           do i2 = 1, grp%nmp(igrp2)
              imp2 = grp%implist(i2, igrp2)
              d(:) = xyz_mp_rep(:,imp1,IREP) - xyz_mp_rep(:,imp2,IREP)
              dist2 = dot_product(d,d)
              if (dist2 < inmisc%brid_ppr_rzero2(ip)) then
                 flg_ppr_release(icom) = .true.
              endif
           enddo
!$omp end parallel do
           if (flg_ppr_release(icom)) exit
        enddo
#endif

     endif
  enddo

  istep(1:ip) = istep(1:ip) + 1

endsubroutine simu_ppr
