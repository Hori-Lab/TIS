! neighbor_list_hb
!> @brief This subroutine is to make neighborlist especially for the hydrogen-bonding interaction

subroutine neighbor_list_bbr(irep)
  
  use if_util
  use const_physical
  use const_maxsize
  use const_index
  use var_setp,    only : inmisc, inperi, indtrna, inpara
  use var_struct,  only : nbbr_bd, lbbr_bd, nbbr, ibbr2mp, &
                          ibd2mp, xyz_mp_rep, ibd2type
  use mpiconst

  implicit none

  integer, intent(in) :: irep
  ! -------------------------------------------------------------------
  ! local variables
  integer :: i, ibd, imp1, imp2
  integer :: j, jbd, jmp1, jmp2
  integer :: ibbr
  integer :: ksta, kend
  integer :: ibbr2mp_l(4, MXBBR)
#ifdef MPI_PAR
  integer :: klen
  integer :: nbbr_l
  integer :: nbbr_lall(0:npar_mpi-1)
  integer :: recvcounts(0:npar_mpi-1), displs(0:npar_mpi-1)
#endif
  real(PREC) :: vij(SDIM), vi(SDIM), vj(SDIM)
  real(PREC) :: d2, dist_cut_sq
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  ! --------------------------------------------------------------------
  if (inmisc%i_neigh_dynamic == 1) then
     dist_cut_sq = (indtrna%bbr_cutoff + inpara%neigh_margin) ** 2
  else
     dist_cut_sq = (indtrna%bbr_cutoff + 5.0) ** 2
     write(error_message,*) 'WARNING: in neighbor_list_BBR, examine this cutoff carefully', indtrna%bbr_cutoff + 5.0
     call util_error(ERROR%WARN_ALL, error_message)
  endif

  ibbr = 0
  ibbr2mp_l(:,:) = 0

#ifdef MPI_PAR2
  klen=(nbbr_bd-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nbbr_bd)
#else
  ksta = 1
  kend = nbbr_bd
#endif

  do i = ksta, kend

     ibd = lbbr_bd(i)
     imp1 = ibd2mp(1,ibd)
     imp2 = ibd2mp(2,ibd)

     do j = i+1, nbbr_bd

        jbd = lbbr_bd(j)

        ! Consider only P-S vs P-S, and P-S and S-P
        ! Excluded volume (DT15) works well for S-P vs S-P
        if (ibd2type(ibd) == BDTYPE%RNA_SP .and. ibd2type(jbd) == BDTYPE%RNA_SP) then
           cycle
        endif

        jmp1 = ibd2mp(1,jbd)
        jmp2 = ibd2mp(2,jbd)

        if (imp1 == jmp1 .or. imp1 == jmp2 .or. imp2 == jmp1 .or. imp2 == jmp2) then
           cycle
        endif

        if(inperi%i_periodic == 0) then
           !   0.5 * (imp1 + imp2) - 0.5 * (jmp1 + jmp2)
           ! = 0.5 * {(imp1 + imp2) - (jmp1 + jmp2)}
           vi(1:3) = xyz_mp_rep(1:3, imp1, irep) + xyz_mp_rep(1:3,imp2,irep)
           vj(1:3) = xyz_mp_rep(1:3, jmp1, irep) + xyz_mp_rep(1:3,jmp2,irep)
           vij(1:3) = 0.5 * (vi(1:3) - vj(1:3))
        else
           !   0.5 * (imp1 + imp2) - 0.5 * (jmp1 + jmp2)
           !  (Because we need to use util_pbneighbor)
           vi(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
           call util_pbneighbor(vi)
           vi(1:3) = xyz_mp_rep(1:3, imp1, irep) + 0.5 * vi(1:3)
           call util_pbneighbor(vi)
   
           vj(1:3) = xyz_mp_rep(1:3, jmp2, irep) - xyz_mp_rep(1:3, jmp1, irep)
           call util_pbneighbor(vj)
           vj(1:3) = xyz_mp_rep(1:3, jmp1, irep) + 0.5 * vj(1:3)
           call util_pbneighbor(vj)
   
           vij(1:3) = vi(1:3) - vj(1:3)
           call util_pbneighbor(vij)
        end if

        d2 = dot_product(vij,vij)
      
        if (d2 <= dist_cut_sq) then
           ibbr = ibbr + 1
           if (ibbr > MXBBR) then
              write(error_message,*) 'Error: ibbr > MXBBR in neighbor_list_BBR.'
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           ibbr2mp_l(1,ibbr) = imp1
           ibbr2mp_l(2,ibbr) = imp2
           ibbr2mp_l(3,ibbr) = jmp1
           ibbr2mp_l(4,ibbr) = jmp2
           !write(*,*) 'neighbor_list_BBR: ',ibbr, imp1,imp2,jmp1,jmp2
        endif

     enddo
  enddo

#ifdef MPI_PAR2
  nbbr_l = ibbr

  call mpi_allgather( nbbr_l,    1, MPI_INTEGER, &
                      nbbr_lall, 1, MPI_INTEGER, &
                      mpi_comm_local, ierr)

  nbbr(irep) = sum( nbbr_lall(0:npar_mpi-1) )
  
  displs(0) = 0
  recvcounts(0) = 4*nbbr_lall(0)

  do i=1, npar_mpi-1
     displs(i) = displs(i-1) + 4*nbbr_lall(i-1)
     recvcounts(i) = 4*nbbr_lall(i)
  enddo

  call mpi_allgatherv( ibbr2mp_l,       4*nbbr_l, MPI_INTEGER,&
                       ibbr2mp(1,1,irep), recvcounts,  displs,  MPI_INTEGER, &
                       mpi_comm_local, ierr)
#else
  nbbr(irep) = ibbr
  ibbr2mp(:,:,irep) = ibbr2mp_l(:,:)
#endif

end subroutine neighbor_list_bbr
