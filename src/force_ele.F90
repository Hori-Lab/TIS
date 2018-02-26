! force_ele
!> @brief This subroutine calculates the force of electrostatic interaction.

subroutine force_ele(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inele, inperi
  use var_struct, only : pxyz_mp_rep, lele, iele2mp, coef_ele, nmp_all
  use var_replica,only : irep2grep
  use var_simu, only : istep
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SDIM, nmp_all)

  integer :: imp1, imp2
  integer :: ksta, kend
  integer :: grep, iele1
  integer :: imirror
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
  write(*,*) '#### start force_ele'
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

!$omp do private(imp1,imp2,v21,dist2,dist1,rdist1,dvdw_dr,for,imirror)
  do iele1=ksta, kend
     imp1 = iele2mp(1, iele1, irep)
     imp2 = iele2mp(2, iele1, irep)
     imirror = iele2mp(3, iele1, irep)

     !if(inperi%i_periodic == 0) then
     !   v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     !else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     !end if

     ! v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = dot_product(v21,v21)
     if(dist2 > cutoff2) cycle
   
     ! -----------------------------------------------------------------
     dist1 = sqrt(dist2)
     rdist1 = 1.0 / dist1
     
     dvdw_dr = coef_ele(iele1, irep) * rdist1 * rdist1 * &
          (rdist1 + rcdist) * &
          exp(-dist1 * rcdist)
     
     if(dvdw_dr > DE_MAX) then
        write(error_message,*) 'force_ele > DE_MAX', istep, imp1, imp2, dist1, dvdw_dr, DE_MAX
        call util_error(ERROR%WARN_ALL, error_message)
        dvdw_dr = DE_MAX
     end if
     ! if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
        
     for(1:3) = dvdw_dr * v21(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
  end do
!$omp end do nowait
#ifdef _DEBUG
  write(*,*) '#### end force_ele'
  call flush(6)
#ifdef MPI_PAR
  call MPI_barrier(MPI_COMM_local,ierr)
#endif
#endif

end subroutine force_ele
