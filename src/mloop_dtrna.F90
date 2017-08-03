subroutine mloop_dtrna()

  use const_maxsize
  use const_index
  use var_io,      only : outfile, flg_file_out
  use var_setp,    only : inmisc, indtrna15, inpara
  use var_struct,  only : ndtrna_hb, ndtrna_st, ndtrna_tst, &
                          idtrna_st2mp, idtrna_tst2st, idtrna_tst2side, idtrna_tst2mp, &
                          dtrna_hb_nat, dtrna_hb_neigh_dist2, nhbsite
  use var_replica, only : n_replica_mpi
  use var_simu,    only : hb_energy, hb_status, flg_hb_energy, st_status, hbsite_excess, &
                          ene_st, ene_tst

  implicit none
  integer :: ier = 0
  integer :: lunout
  integer :: ist, itst, imp1, imp2, imp1_st, imp2_st, ihb
  logical :: flg_found
  character(CARRAY_MSG_ERROR) :: error_message
! **********************************************************************

  imp2 = 0  ! to supress compiler message
  lunout   = outfile%data
  error_message = 'failed in memory allocation at mloop_dtrna, PROGRAM STOP'

  if (allocated(st_status)) deallocate(st_status)
  if (allocated(ene_st )) deallocate(ene_st)
  if (allocated(ene_tst)) deallocate(ene_tst)

  allocate( st_status(1:ndtrna_st, 1:n_replica_mpi), stat=ier)
  if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
  st_status(:,:) = .True.

  ! ene_st and ene_tst
  if (flg_file_out%st  .or. flg_file_out%stall .or. &
      flg_file_out%tst .or. flg_file_out%tstall ) then

     allocate( ene_st(ndtrna_st, n_replica_mpi), stat=ier)
     if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
     ene_st(:,:) = 0.0e0_PREC

     allocate( ene_tst(ndtrna_tst, n_replica_mpi), stat=ier)
     if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
     ene_tst(:,:) = 0.0e0_PREC
  endif


  if (inmisc%i_dtrna_model == 2015) then
     
     if (allocated(hb_energy)) deallocate(hb_energy)
     if (allocated(hb_status)) deallocate(hb_status)
   
     allocate( hb_energy(1:ndtrna_hb, 1:n_replica_mpi), stat=ier)
     if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
     hb_energy(:,:) = 0.0e0_PREC
     flg_hb_energy = .False.
   
     allocate( hb_status(1:ndtrna_hb, 1:n_replica_mpi), stat=ier)
     if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
     hb_status(:,:) = .False.

     allocate( hbsite_excess(1:nhbsite), stat=ier)
     if (ier /= 0) call util_error(ERROR%STOP_ALL,error_message)
     hbsite_excess(:) = 0

     ! Assign idtrna_tst2st
     do itst = 1, ndtrna_tst
   
        imp1 = idtrna_tst2mp(1,itst)
        if (idtrna_tst2side(1,itst) == 1) then
           imp2 = imp1 - 3
        else if (idtrna_tst2side(1,itst) == 2) then
           imp2 = imp1 + 3
        else
           write(error_message, *) "invalid idtrna_tst2side in mloop_dtrna"
           call util_error(ERROR%STOP_ALL, error_message)
        endif
   
        flg_found = .False.
        do ist = 1, ndtrna_st
           imp1_st = idtrna_st2mp(1,ist)
           imp2_st = idtrna_st2mp(2,ist)
           if ((imp1 == imp1_st .and. imp2 == imp2_st) .or. &
               (imp1 == imp2_st .and. imp2 == imp1_st)) then
              idtrna_tst2st(1,itst) = ist
              flg_found = .True.
              exit
           endif
        enddo
   
        if (.not. flg_found) then
           write(error_message, *) "WARNING: corresponding stack was not found in mloop_dtrna, itst=", itst
           call util_error(ERROR%WARN_ALL, error_message)
        endif
   
   
        imp1 = idtrna_tst2mp(2,itst)
        if (idtrna_tst2side(2,itst) == 1) then
           imp2 = imp1 - 3
        else if (idtrna_tst2side(2,itst) == 2) then
           imp2 = imp1 + 3
        else
           write(error_message, *) "invalid idtrna_tst2side in mloop_dtrna"
           call util_error(ERROR%STOP_ALL, error_message)
        endif
   
        flg_found = .False.
        do ist = 1, ndtrna_st
           imp1_st = idtrna_st2mp(1,ist)
           imp2_st = idtrna_st2mp(2,ist)
           if ((imp1 == imp1_st .and. imp2 == imp2_st) .or.&
               (imp1 == imp2_st .and. imp2 == imp1_st)) then
              idtrna_tst2st(2,itst) = ist
              flg_found = .True.
              exit
           endif
        enddo
   
        if (.not. flg_found) then
           write(error_message, *) "WARNING: corresponding stack was not found in mloop_dtrna, itst=", itst
           call util_error(ERROR%WARN_ALL, error_message)
        endif
   
     enddo

     ! calc square of distance for neighbor list
     if (inmisc%i_neigh_dynamic == 1) then
        do ihb = 1, ndtrna_hb
           dtrna_hb_neigh_dist2(ihb) = (dtrna_hb_nat(1,ihb) + indtrna15%hb_cutoff_dist + inpara%neigh_margin) ** 2
        enddo
     else
        dtrna_hb_neigh_dist2(:) = 20.0e0_PREC ** 2
     endif

  endif

end subroutine mloop_dtrna
