subroutine simu_force_ele(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : outfile
  use var_setp,   only : inele, inmisc
  use var_struct, only : xyz_mp_rep, lpnl, ipnl2mp, &
                         lele, iele2mp, coef_ele,   &
                         nmp_all
  use var_replica,only : irep2grep
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

#ifdef MPI_PAR
  integer :: klen, ksta, kend
#endif

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

  ! --------------------------------------------------------------------
  ! local variables                                                             
  integer :: imp1, imp2, imp3, imp4
  integer :: grep
  integer :: iele1
  real(PREC) :: dist1, dist2, rdist1
  real(PREC) :: dvdw_dr, coef, rcdist, cutoff2
  real(PREC) :: v21(3), for(3)
  real(PREC) :: v211,v212,v213, for1, for2, for3
  real(PREC) :: v2111,v2121,v2131
  real(PREC) :: v2112,v2122,v2132
  real(PREC) :: force_mp11,force_mp21,force_mp31
  real(PREC) :: force_mp12,force_mp22,force_mp32
  real(PREC) :: dist21, dist22
  integer    :: ielem

#ifdef _DEBUG
  write(*,*) '#### start simu_force_ele'
#endif

  if (inmisc%force_flag(INTERACT%ELE)) then
     ! for speed up
     grep = irep2grep(irep)
     cutoff2 = (inele%cutoff_ele * inele%cdist(grep))**2
     rcdist = 1.0e0_PREC / inele%cdist(grep)
#ifdef _DEBUG
     write(*,*) 'electrostatic'
     write(*,*) 'lele(irep), ',lele(irep)
#endif

#ifdef MPI_PAR
     klen=(lele(irep)-1+npar_mpi)/npar_mpi
     ksta=1+klen*local_rank_mpi
     kend=min(ksta+klen-1,lele(irep))

#ifdef SHARE_NEIGH
     ielem = mod(kend-ksta+1, 2)
     do iele1=ksta, ksta+ielem-1
#else
     ielem = mod(lele(irep), 2)
     do iele1=1, ielem
#endif

#else
     ielem = mod(lele(irep), 2)
     do iele1 = 1, ielem
#endif
        imp1 = iele2mp(1, iele1, irep)
        imp2 = iele2mp(2, iele1, irep)
   
        v211 = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
        v212 = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
        v213 = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)
        dist2 = v211*v211 + v212*v212 + v213*v213

        if(dist2 > cutoff2) cycle
   
        force_mp11 = force_mp(1, imp1)
        force_mp21 = force_mp(2, imp1)
        force_mp31 = force_mp(3, imp1)

        ! -----------------------------------------------------------------
        dist1 = sqrt(dist2)
        dvdw_dr = coef_ele(iele1, irep) / dist1 / dist1 * &
             (1.0e0_PREC / dist1 + rcdist) * &
             exp(-dist1 * rcdist)

        force_mp12 = force_mp(1, imp2)
        force_mp22 = force_mp(2, imp2)
        force_mp32 = force_mp(3, imp2)

        if(dvdw_dr > DE_MAX) then
   !        write (*, *) "electrostatic", imp1, imp2, dvdw_dr
           dvdw_dr = DE_MAX
        end if
        !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
        
        for1 = dvdw_dr * v211
        for2 = dvdw_dr * v212
        for3 = dvdw_dr * v213
        force_mp(1, imp1) = force_mp11 - for1
        force_mp(2, imp1) = force_mp21 - for2
        force_mp(3, imp1) = force_mp31 - for3

        force_mp(1, imp2) = force_mp12 + for1
        force_mp(2, imp2) = force_mp22 + for2
        force_mp(3, imp2) = force_mp32 + for3
     end do

#ifdef MPI_PAR
!$omp do private(imp1,imp2,imp3,imp4, &
!$omp            v2111,v2121,v2131, &
!$omp            v2112,v2122,v2132, &
!$omp            dist1,dist21,dist22,rdist1, &
!$omp&           dvdw_dr,for1,for2,for3)

#ifdef SHARE_NEIGH
     do iele1=ksta+ielem, kend, 2
#else
     do iele1=ielem+1, lele(irep), 2
#endif

#else
     do iele1 =ielem+1, lele(irep), 2
#endif
        imp1 = iele2mp(1, iele1, irep)
        imp2 = iele2mp(2, iele1, irep)
   
        v2111  = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
        v2121  = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
        v2131  = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)
        dist21 = v2111*v2111 + v2121*v2121 + v2131*v2131

        imp3 = iele2mp(1, iele1+1, irep)
        imp4 = iele2mp(2, iele1+1, irep)

        v2112  = xyz_mp_rep(1, imp4, irep) - xyz_mp_rep(1, imp3, irep)
        v2122  = xyz_mp_rep(2, imp4, irep) - xyz_mp_rep(2, imp3, irep)
        v2132  = xyz_mp_rep(3, imp4, irep) - xyz_mp_rep(3, imp3, irep)
        dist22 = v2112*v2112 + v2122*v2122 + v2132*v2132

        if(dist21 > cutoff2) goto 200

        ! -----------------------------------------------------------------
        dist1 = sqrt(dist21)
        dvdw_dr = coef_ele(iele1, irep) / dist1 / dist1 * &
             (1.0e0_PREC / dist1 + rcdist) * &
             exp(-dist1 * rcdist)

        if(dvdw_dr > DE_MAX) then
   !        write (*, *) "electrostatic", imp1, imp2, dvdw_dr
           dvdw_dr = DE_MAX
        end if
        !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
   
        force_mp11 = force_mp(1, imp1)
        force_mp21 = force_mp(2, imp1)
        force_mp31 = force_mp(3, imp1)
        
        force_mp12 = force_mp(1, imp2)
        force_mp22 = force_mp(2, imp2)
        force_mp32 = force_mp(3, imp2)

        for1 = dvdw_dr * v2111
        for2 = dvdw_dr * v2121
        for3 = dvdw_dr * v2131

        force_mp(1, imp1) = force_mp11 - for1
        force_mp(2, imp1) = force_mp21 - for2
        force_mp(3, imp1) = force_mp31 - for3

        force_mp(1, imp2) = force_mp12 + for1
        force_mp(2, imp2) = force_mp22 + for2
        force_mp(3, imp2) = force_mp32 + for3

 200    continue 

        if(dist22 > cutoff2) cycle
   
        ! -----------------------------------------------------------------
        dist1 = sqrt(dist22)
        dvdw_dr = coef_ele(iele1+1, irep) / dist1 / dist1 * &
             (1.0e0_PREC / dist1 + rcdist) * &
             exp(-dist1 * rcdist)

        if(dvdw_dr > DE_MAX) then
   !        write (*, *) "electrostatic", imp1, imp2, dvdw_dr
           dvdw_dr = DE_MAX
        end if
        !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
        
        force_mp11 = force_mp(1, imp3)
        force_mp21 = force_mp(2, imp3)
        force_mp31 = force_mp(3, imp3)

        force_mp12 = force_mp(1, imp4)
        force_mp22 = force_mp(2, imp4)
        force_mp32 = force_mp(3, imp4)

        for1 = dvdw_dr * v2112
        for2 = dvdw_dr * v2122
        for3 = dvdw_dr * v2132

        force_mp(1, imp3) = force_mp11 - for1
        force_mp(2, imp3) = force_mp21 - for2
        force_mp(3, imp3) = force_mp31 - for3

        force_mp(1, imp4) = force_mp12 + for1
        force_mp(2, imp4) = force_mp22 + for2
        force_mp(3, imp4) = force_mp32 + for3
     end do
#ifdef MPI_PAR
!$omp end do nowait
#endif
  endif

end subroutine simu_force_ele
