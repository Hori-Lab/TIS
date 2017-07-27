! neighbor_list_hb
!> @brief This subroutine is to make neighborlist especially for the hydrogen-bonding interaction

subroutine neighbor_list_hb(irep)
  
  use const_maxsize
  use const_index
  use var_setp,    only : inmisc
  use var_struct,  only : ndtrna_hb, idtrna_hb2mp, nhbneigh, ineigh2hb, nhbsite, & 
                          idtrna_hb2hbsite, list_hb_at_hbsite, num_hb_at_hbsite, &
                          xyz_mp_rep, nmp_all, dtrna_hb_neigh_dist2
  use mpiconst

  implicit none

  integer, intent(in) :: irep
  ! -------------------------------------------------------------------
  ! local variables
  integer :: i, ihb, ineigh, ihbsite
  integer :: ksta, kend
  integer :: ineigh2hb_l(1:MXMPHBNEIGHBOR*nmp_all/2)
#ifdef MPI_PAR
  integer :: klen
  integer :: nhbneigh_l
  integer :: nhbneigh_lall(0:npar_mpi-1)
  integer :: recvcounts(0:npar_mpi-1), displs(0:npar_mpi-1)
#endif
  real(PREC) :: v12(3), d1212
!  real(PREC) :: dist_cut_sq

  ! -------------------------------------------------------------------
  if (inmisc%i_dtrna_model /= 2015) then
     !nhbneigh(:) = 0   ! this is not allocated
     return
  end if

  ! --------------------------------------------------------------------
!  if (inmisc%i_neigh_dynamic == 1) then
!     dist_cut_sq = (dtrna_hb_longest + indtrna15%hb_cutoff_dist + inpara%neigh_margin) ** 2
!  else
!     dist_cut_sq = 20.0e0_PREC ** 2
!  endif

  ineigh = 0
  !ineigh2hb(:, irep) = 0
  ineigh2hb_l(:) = 0

#ifdef MPI_PAR2
  klen=(ndtrna_hb-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,ndtrna_hb)
#else
  ksta = 1
  kend = ndtrna_hb
#endif

  do ihb=ksta,kend

     !===== calc vectors =====
     v12 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
          -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)

     d1212 = dot_product(v12,v12)

     !if (d1212 < dist_cut_sq) then
     if (d1212 < dtrna_hb_neigh_dist2(ihb)) then
        ineigh = ineigh + 1
        ineigh2hb_l(ineigh) = ihb
     endif
  enddo

#ifdef MPI_PAR2
  nhbneigh_l = ineigh

  call mpi_allgather( nhbneigh_l,    1, MPI_INTEGER, &
                      nhbneigh_lall, 1, MPI_INTEGER, &
                      mpi_comm_local, ierr)

  nhbneigh(irep) = sum( nhbneigh_lall(0:npar_mpi-1) )
  
  displs(0) = 0
  recvcounts(0) = nhbneigh_lall(0)

  do i=1, npar_mpi-1
     displs(i) = displs(i-1) + nhbneigh_lall(i-1)
     recvcounts(i) = nhbneigh_lall(i)
  enddo

  call mpi_allgatherv( ineigh2hb_l,       nhbneigh_l, MPI_INTEGER,&
                       ineigh2hb(1,irep), recvcounts,  displs,  MPI_INTEGER, &
                       mpi_comm_local, ierr)
#else
  nhbneigh(irep) = ineigh
  ineigh2hb(:,irep) = ineigh2hb_l(:)
#endif

     
  num_hb_at_hbsite(1:nhbsite, irep) = 0
  list_hb_at_hbsite(:, :, irep) = 0

!#ifdef MPI_PAR2
!  klen=(ndtrna_hb-1+npar_mpi)/npar_mpi
!  ksta=1+klen*local_rank_mpi
!  kend=min(ksta+klen-1,ndtrna_hb)
!#else
  ksta = 1
  kend = nhbneigh(irep)
!#endif

  do ineigh=ksta, kend
     ihb = ineigh2hb(ineigh, irep)

     do i=1,3
        ihbsite = idtrna_hb2hbsite(i, 1, ihb)
        if (ihbsite > 0) then
           num_hb_at_hbsite(ihbsite, irep) = num_hb_at_hbsite(ihbsite, irep) + 1
           list_hb_at_hbsite( num_hb_at_hbsite(ihbsite, irep), ihbsite, irep) = ihb
        endif

        ihbsite = idtrna_hb2hbsite(i, 2, ihb)
        if (ihbsite > 0) then
           num_hb_at_hbsite(ihbsite, irep) = num_hb_at_hbsite(ihbsite, irep) + 1
           list_hb_at_hbsite( num_hb_at_hbsite(ihbsite, irep), ihbsite, irep) = ihb
        endif
     enddo
  enddo

end subroutine neighbor_list_hb
