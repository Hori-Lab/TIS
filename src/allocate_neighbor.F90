! allocate_neighbor
!> @brief Allocate/Deallocate arrays of neighborling list

subroutine allocate_neighbor()

  use const_maxsize
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inmisc, inele
  use var_struct,  only : lpnl, ipnl2mp, lele, iele2mp, &
                          coef_ele, &
                          nhpneigh, lhp2neigh, ineigh2hp, cutoff_dmin_hp, &
                          cutoff_dmax_hp, ncharge, nmp_all, nhp, &
                          nhbneigh, ineigh2hb
  use var_replica, only : n_replica_mpi

#ifdef MPI_PAR2
  use mpiconst
#endif

  implicit none

  ! -------------------------------------------------------------------
  ! local variables
  integer :: ier
  integer :: ncharge_mpi
  integer :: n_index, n_index_tail
  logical :: flg_error
  character(CARRAY_MSG_ERROR) :: error_message
  
  ! -------------------------------------------------------------------
  ncharge_mpi = ncharge
#ifdef MPI_PAR2
  ncharge_mpi = ncharge_l
#endif
  
  n_index = 2 + inperi%n_mirror_index
  n_index_tail = 1 + inperi%n_mirror_index
  
  ! check
  flg_error = .false.
  if (allocated(lpnl))           flg_error = .true.
  if (allocated(ipnl2mp))        flg_error = .true.
  if (allocated(lele))           flg_error = .true.
  if (allocated(iele2mp))        flg_error = .true.
  if (allocated(coef_ele))       flg_error = .true.
  if (allocated(nhpneigh))       flg_error = .true.
  if (allocated(lhp2neigh))   flg_error = .true.
  if (allocated(ineigh2hp))   flg_error = .true.
  if (allocated(cutoff_dmax_hp)) flg_error = .true.
  if (allocated(cutoff_dmin_hp)) flg_error = .true.
  if (allocated(nhbneigh))       flg_error = .true.
  if (allocated(ineigh2hb))      flg_error = .true.
  
  if (flg_error) then
     error_message = 'defect at allocate_neighbor, PROGRAM STOP'
     call util_error(ERROR%STOP_ALL,error_message)
  endif
  
  !-----------
  ! allocate
  !-----------
  error_message = 'failed in memory allocation at allocate_neighbor, PROGRAM STOP'
  
  allocate( lpnl(2, E_TYPE%MAX, n_replica_mpi),             stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  lpnl(:,:,:) = 0

#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_PNL
  allocate( ipnl2mp(n_index, MXMPNEIGHBOR*nmp_all, n_replica_mpi), stat=ier)
#else
  allocate( ipnl2mp(n_index, MXMPNEIGHBOR*nmp_all/npar_mpi+1, n_replica_mpi), stat=ier)
#endif

#else
  allocate( ipnl2mp(n_index, MXMPNEIGHBOR*nmp_all, n_replica_mpi), stat=ier)
#endif

  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  ipnl2mp(:,:,:) = 0
  
  ! for electrostatic interaction
  if (inmisc%force_flag(INTERACT%ELE)) then
     !      if(inele%i_calc_method == 0) then
     allocate( lele(n_replica_mpi),                           stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     lele(:) = 0
     
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH
     allocate( iele2mp(n_index, MXMPELE*ncharge, n_replica_mpi),    stat=ier)
#else
     allocate( iele2mp(n_index, MXMPELE*ncharge/npar_mpi+1, n_replica_mpi),    stat=ier)
#endif

#else
     allocate( iele2mp(n_index, MXMPELE*ncharge, n_replica_mpi),    stat=ier)
#endif

     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     iele2mp(:,:,:) = 0
     
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH
     allocate( coef_ele(MXMPELE*ncharge, n_replica_mpi),      stat=ier)
#else
     allocate( coef_ele(MXMPELE*ncharge/npar_mpi+1, n_replica_mpi),      stat=ier)
#endif

#else
     allocate( coef_ele(MXMPELE*ncharge, n_replica_mpi),      stat=ier)
#endif

     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     coef_ele(:,:) = 0.0e0_PREC

     if(inele%i_calc_method /= 0) then
        call util_error(ERROR%STOP_ALL, 'i_calc_method /= 0 in allocate_neighbor')
     endif
  endif
  
  ! for hydrophobic interaction
  if (inmisc%force_flag(INTERACT%HP)) then
     allocate( nhpneigh(n_replica_mpi),              stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message) 
     nhpneigh(:) = 0
     
     allocate( lhp2neigh(2, nhp, n_replica_mpi),   stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message) 
     lhp2neigh(:,:,:) = 0

#ifdef MPI_PAR2
     
#ifdef SHARE_NEIGH_HP
     allocate( ineigh2hp(MXMPHP*nhp, n_replica_mpi),   stat=ier)
#else
     allocate( ineigh2hp(MXMPHP*nhp/npar_mpi + 1, n_replica_mpi),   stat=ier)
#endif

#else
     allocate( ineigh2hp(MXMPHP*nhp, n_replica_mpi),   stat=ier)
#endif

     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message) 
     ineigh2hp(:,:) = 0
     
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_HP
     allocate( cutoff_dmin_hp(MXMPHP*nhp, n_replica_mpi),    stat=ier)
#else
     allocate( cutoff_dmin_hp(MXMPHP*nhp/npar_mpi + 1, n_replica_mpi),    stat=ier)
#endif

#else
     allocate( cutoff_dmin_hp(MXMPHP*nhp, n_replica_mpi),    stat=ier)
#endif

     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message) 
     cutoff_dmin_hp(:,:) = 0
     
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_HP
     allocate( cutoff_dmax_hp(MXMPHP*nhp, n_replica_mpi),    stat=ier)
#else
     allocate( cutoff_dmax_hp(MXMPHP*nhp/npar_mpi + 1, n_replica_mpi),    stat=ier)
#endif

#else
     allocate( cutoff_dmax_hp(MXMPHP*nhp, n_replica_mpi),    stat=ier)
#endif

     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message) 
     cutoff_dmax_hp(:,:) = 0
  end if

  ! for hydrogen-bonding interaction in DTRNA2015
  if (inmisc%i_dtrna_model == 2015) then
     allocate( nhbneigh(n_replica_mpi), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)

     allocate( ineigh2hb(MXMPHBNEIGHBOR*nmp_all/2 , n_replica_mpi), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  endif

end subroutine allocate_neighbor

!######################################################################################

subroutine deallocate_neighbor

   use var_struct,  only : lpnl, ipnl2mp, lele, iele2mp, coef_ele, &
                           nhpneigh, ineigh2hp, lhp2neigh, cutoff_dmax_hp, cutoff_dmin_hp, &
                           nhbneigh, ineigh2hb

   if (allocated(lpnl))          deallocate(lpnl)
   if (allocated(ipnl2mp))       deallocate(ipnl2mp)
   if (allocated(lele))          deallocate(lele)
   if (allocated(iele2mp))       deallocate(iele2mp)
   if (allocated(coef_ele))      deallocate(coef_ele)
   if (allocated(nhpneigh))       deallocate(nhpneigh)
   if (allocated(ineigh2hp))  deallocate(ineigh2hp)
   if (allocated(lhp2neigh))  deallocate(lhp2neigh)
   if (allocated(cutoff_dmax_hp))deallocate(cutoff_dmax_hp)
   if (allocated(cutoff_dmin_hp))deallocate(cutoff_dmin_hp)
   if (allocated(nhbneigh))      deallocate(nhbneigh)
   if (allocated(ineigh2hb))     deallocate(ineigh2hb)
   
endsubroutine deallocate_neighbor
