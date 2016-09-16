!setp_dtrna15
!> @brief
!>       

subroutine setp_dtrna15()

  use const_index
  use const_maxsize
  use var_struct, only : nmp_all, cmp2seq, cmp2atom, iclass_mp,&
                         nhbsite, imp2hbsite, nvalence_hbsite,&
                         list_hb_at_hbsite, num_hb_at_hbsite, &
                         exv_radius_mp, exv_epsilon_mp
  use var_replica, only : n_replica_mpi
  use var_setp,    only : indtrna15
  implicit none
  
  ! -------------------------------------------------------------
  integer :: imp, ihbsite
  integer :: ier
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------

  allocate( imp2hbsite(2,nmp_all) )

  imp2hbsite(1:2, 1:nmp_all) = 0

  ihbsite = 0

  do imp = 1, nmp_all
     
     if (iclass_mp(imp) /= CLASS%RNA) then
        cycle
     endif

     if (cmp2atom(imp) == ' S  ') then
        imp2hbsite(1,imp) = ihbsite + 1
        imp2hbsite(2,imp) = ihbsite + 2
        ihbsite = ihbsite + 2
     else if (cmp2atom(imp) == ' P  ') then
        imp2hbsite(1,imp) = ihbsite + 1
        imp2hbsite(2,imp) = ihbsite + 2
        ihbsite = ihbsite + 2
     else if (cmp2atom(imp) == ' Ab ') then
        imp2hbsite(1,imp) = ihbsite + 1
        imp2hbsite(2,imp) = ihbsite + 4
        ihbsite = ihbsite + 4
     else if (cmp2atom(imp) == ' Cb ') then
        imp2hbsite(1,imp) = ihbsite + 1
        imp2hbsite(2,imp) = ihbsite + 3
        ihbsite = ihbsite + 3
     else if (cmp2atom(imp) == ' Gb ') then
        imp2hbsite(1,imp) = ihbsite + 1
        imp2hbsite(2,imp) = ihbsite + 5
        ihbsite = ihbsite + 5
     else if (cmp2atom(imp) == ' Ub ') then
        imp2hbsite(1,imp) = ihbsite + 1
        imp2hbsite(2,imp) = ihbsite + 3
        ihbsite = ihbsite + 3
     endif

  enddo

  nhbsite = ihbsite

  allocate( nvalence_hbsite(1:nhbsite) , stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  nvalence_hbsite(:) = 0
  allocate( list_hb_at_hbsite(MXMPHBNEIGHBOR, 1:nhbsite, n_replica_mpi) , stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  list_hb_at_hbsite(:,:,:) = 0
  allocate( num_hb_at_hbsite(1:nhbsite, n_replica_mpi) , stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  num_hb_at_hbsite(:,:) = 0

  ihbsite = 0
  do imp = 1, nmp_all
     
     if (iclass_mp(imp) == CLASS%ION) then

        if (cmp2seq(imp) == 'Mg ') then
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%MG2))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%MG2)
        else if (cmp2seq(imp) == 'K  ') then
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%K))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%K)
        else if (cmp2seq(imp) == 'Na ') then
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%NA))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%NA)
        else if (cmp2seq(imp) == 'Cl ') then
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%CL))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%CL)
        else if (cmp2seq(imp) == 'Ca2') then
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%CA2))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%CA2)
        endif

     else if (iclass_mp(imp) == CLASS%RNA) then

        if (cmp2atom(imp) == ' S  ') then
           nvalence_hbsite(ihbsite+1) = 2
           nvalence_hbsite(ihbsite+2) = 1
           ihbsite = ihbsite + 2
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%S))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%S)
        else if (cmp2atom(imp) == ' P  ') then
           nvalence_hbsite(ihbsite+1) = 1
           nvalence_hbsite(ihbsite+2) = 1
           ihbsite = ihbsite + 2
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%P))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%P)
        else if (cmp2atom(imp) == ' Ab ') then
           nvalence_hbsite(ihbsite+1) = 1
           nvalence_hbsite(ihbsite+2) = 1
           nvalence_hbsite(ihbsite+3) = 2
           nvalence_hbsite(ihbsite+4) = 1
           ihbsite = ihbsite + 4
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%A))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%A)
        else if (cmp2atom(imp) == ' Cb ') then
           nvalence_hbsite(ihbsite+1) = 2
           nvalence_hbsite(ihbsite+2) = 1
           nvalence_hbsite(ihbsite+3) = 2
           ihbsite = ihbsite + 3
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%C))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%C)
        else if (cmp2atom(imp) == ' Gb ') then
           nvalence_hbsite(ihbsite+1) = 1
           nvalence_hbsite(ihbsite+2) = 2
           nvalence_hbsite(ihbsite+3) = 1
           nvalence_hbsite(ihbsite+4) = 2
           nvalence_hbsite(ihbsite+5) = 1
           ihbsite = ihbsite + 5
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%G))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%G)
        else if (cmp2atom(imp) == ' Ub ') then
           nvalence_hbsite(ihbsite+1) = 1
           nvalence_hbsite(ihbsite+2) = 1
           nvalence_hbsite(ihbsite+3) = 1
           ihbsite = ihbsite + 3
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%U))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%U)
        endif

     else if (iclass_mp(imp) == CLASS%LIG) then

        if (cmp2atom(imp) == ' X1 ') then
           exv_epsilon_mp(imp) = sqrt(indtrna15%exv_eps(DT15EXV%X1))
           exv_radius_mp(imp)  = indtrna15%exv_rad(DT15EXV%X1)
        endif

     else
        cycle
     endif

  enddo

endsubroutine setp_dtrna15
