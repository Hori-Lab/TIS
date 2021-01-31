! force_ele
!> @brief This subroutine calculates the force of electrostatic interaction.

subroutine force_ele_DH(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inele, inperi, inpmf
  use var_struct, only : pxyz_mp_rep, lele, iele2mp, coef_ele, nmp_all
  use var_replica,only : irep2grep
  use var_simu,   only : istep, pmfdh_force
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SDIM, nmp_all)

  integer :: imp1, imp2
  integer :: ksta, kend
  integer :: grep, iele1
  integer :: imir, ipmf, ibin
  real(PREC) :: dist1, dist2, rdist1
  real(PREC) :: dvdw_dr, rcdist, cutoff2
  real(PREC) :: ek, ek_simu
  real(PREC) :: v21(3), for(3)
#ifdef MPI_PAR
  integer :: klen
#endif
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
#ifdef _DEBUG
#ifdef MPI_PAR
  call MPI_barrier(MPI_COMM_local,ierr)
#endif
  write(*,*) '#### start force_ele_DH'
  call flush(6)
#endif

  ! for speed up
  grep = irep2grep(irep)
  cutoff2 = (inele%cutoff_ele * inele%cdist(grep))**2
  rcdist = 1.0e0_PREC / inele%cdist(grep)

  ek_simu = inele%diele
  ek = inele%diele_water

#ifdef _DEBUG
  write(*,*) 'lele(irep), ',lele(irep)
#endif
#ifdef MPI_PAR
#ifdef SHARE_NEIGH
  klen=(lele(irep)-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,lele(irep))
#else
  ksta = 1
  kend = lele(irep)
#endif
#ifdef MPI_DEBUG
  print *,"exv2_6       = ", kend-ksta+1
#endif

#else
  ksta = 1
  kend = lele(irep)
#endif

!$omp do private(imp1,imp2,v21,dist2,dist1,rdist1,dvdw_dr,for,imir, ipmf, ibin)
  do iele1=ksta, kend
     imp1 = iele2mp(1, iele1, irep)
     imp2 = iele2mp(2, iele1, irep)
     imir = iele2mp(3, iele1, irep)
     ipmf = iele2mp(4, iele1, irep)

     v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imir)

     dist2 = dot_product(v21,v21)
     if(dist2 > cutoff2) cycle

     ! -----------------------------------------------------------------
     dist1 = sqrt(dist2)
     rdist1 = 1.0 / dist1

     ! PMF + DH (semiexplicit model)
     if (ipmf > 0 .and. dist1 <= inpmf%Rmax(ipmf)) then 

        if (dist1 < inpmf%Rmin(PMFTYPE%MG_P)) then
           write(error_message,*) 'force_ele_DH(r < Rmin)', istep, irep, imp1, imp2, dist1
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        !dvdw_dr = rdist1 * force_pmfdh(dist1, grep)

        ibin = floor( (dist1 - inpmf%Rmin(PMFTYPE%MG_P)) * inpmf%Rbininv(PMFTYPE%MG_P)) + 1

        dvdw_dr = rdist1 * pmfdh_force(ibin, grep, PMFTYPE%MG_P)

     ! DH
     else
        dvdw_dr = coef_ele(iele1, irep) * rdist1 * rdist1 &
                 * (rdist1 + rcdist) * exp(-dist1 * rcdist)

     endif
     
     if(dvdw_dr > DE_MAX) then
        write(error_message,*) 'force_ele_DH > DE_MAX', istep, irep, imp1, imp2, dist1, dvdw_dr, DE_MAX
        call util_error(ERROR%WARN_ALL, error_message)
        dvdw_dr = DE_MAX
     end if
        
     for(1:3) = dvdw_dr * v21(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
  end do
!$omp end do nowait
#ifdef _DEBUG
  write(*,*) '#### end force_ele_DH'
  call flush(6)
#ifdef MPI_PAR
  call MPI_barrier(MPI_COMM_local,ierr)
#endif
#endif

contains
   
!   real(PREC) function force_pmfdh(r, grep)
!      use var_simu, only : pmfdh_force
!      implicit none
!      real(PREC),intent(in) :: r
!      integer, intent(in) :: grep
!      integer :: ibin
!
!      ! Rmin = 2.6
!      ! Rbin = 0.025
!      !
!      ! r = 2.624 ==> bin 1
!      ! r = 2.625 ==> bin 2
!      ! r = 2.626 ==> bin 2
!      ! r = 2.649 ==> bin 2
!      ! r = 2.650 ==> bin 3
!
!      ! If r = 2.624
!      ! (r - Rmin) = 0.024
!      ! (r - Rmin) * Rbininv = 0.96  ==floor==>  0
!      ! floor((r - Rmin) * Rbininv) + 1 = 1
!
!      ! If r = 2.625
!      ! (r - Rmin) = 0.025
!      ! (r - Rmin) * Rbininv = 1.0  ==floor==>  1
!      ! floor((r - Rmin) * Rbininv) + 1 = 2
!
!      ! If r = 2.626
!      ! (r - Rmin) = 0.026
!      ! (r - Rmin) * Rbininv = 1.04  ==floor==>  1
!      ! floor((r - Rmin) * Rbininv) + 1 = 2
!
!      ! If r = 2.649
!      ! (r - Rmin) = 0.049
!      ! (r - Rmin) * Rbininv = 1.96 ==floor==>  1
!      ! floor((r - Rmin) * Rbininv) + 1 = 2
!
!      ! If r = 2.650
!      ! (r - Rmin) = 0.050
!      ! (r - Rmin) * Rbininv = 2.0 ==floor==>  2
!      ! floor((r - Rmin) * Rbininv) + 1 = 3
!
!      ibin = floor( (r - Rmin) * Rbininv) + 1
!
!      force_pmfdh = pmfdh_force(ibin, grep, PMFTYPE%MG_P)
!
!   endfunction force_pmfdh

end subroutine force_ele_DH
