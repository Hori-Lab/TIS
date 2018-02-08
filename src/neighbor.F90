! neighbor
!> @brief This subroutine is to make neighborlist for the whole of non-local interaction.

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine neighbor(irep)
  
  use const_maxsize
  use const_index
  use var_setp,   only : inmisc, inperi
  use var_struct, only : lexv, nmp_all, nhbneigh, lele, nbbr
  use var_io,     only : flg_file_out, outfile
  use var_replica,only : irep2grep

  use time
  use mpiconst
  implicit none

  integer, intent(in) :: irep

  integer :: grep
  integer :: lmp2neigh((nmp_l+nthreads-1)/nthreads  ,0:nthreads-1)
  integer :: ineigh2mp(MXMPNEIGHBOR*nmp_all/nthreads,0:nthreads-1)

  ! -------------------------------------------------------------------
#ifdef _DEBUG
  write(6,*) '####### start neighbor'
#endif

  if(inperi%i_periodic == 1) then
     call util_periodic(irep)
  end if

!  if (inmisc%force_flag(INTERACT%ENM)) then
!     lexv(1, 1:E_TYPE%MAX,:) = 1
!     lexv(2, 1:E_TYPE%MAX,:) = 0
!     return
!  end if

  TIME_S( tm_neighbor_exv )
  if (inmisc%force_flag(INTERACT%EXV_DT15) .OR. inmisc%force_flag(INTERACT%EXV_WCA) .OR.&
      inmisc%force_flag(INTERACT%EXV12) .OR. inmisc%force_flag(INTERACT%EXV6) .OR.&
      inmisc%force_flag(INTERACT%EXV_GAUSS) .OR. &
      inmisc%force_flag(INTERACT%GO) .OR. inmisc%force_flag(INTERACT%LJ) ) then
      !inmisc%force_flag(INTERACT%PAIR_RNA) .OR.& inmisc%force_flag(INTERACT%STACK_RNA) .OR.&
      !inmisc%force_flag(INTERACT%MORSE) .OR. inmisc%force_flag(INTERACT%SASA) .OR.&
      !inmisc%force_flag(INTERACT%AICG1) .OR. inmisc%force_flag(INTERACT%AICG2)) then
     ! make neighborlist
     call neighbor_list(irep, ineigh2mp, lmp2neigh)

     ! assign iconcal2con and iexv2mp
     call neighbor_assign(irep, ineigh2mp, lmp2neigh)

  else
     lexv(1, 1:E_TYPE%MAX,:) = 1
     lexv(2, 1:E_TYPE%MAX,:) = 0
  endif
  TIME_E( tm_neighbor_exv )
  
  TIME_S( tm_neighbor_ele )
  if (inmisc%force_flag(INTERACT%ELE)) then
     ! make neighborlist for electrostatic interaction
     call neighbor_list_ele(irep)
     !if(inele%i_calc_method == 0) then
     !   call neighbor_list_ele(irep)
     !else if(inele%i_calc_method == 1) then
     !   call neighbor_list_ele2(irep)
     !end if
  endif
  TIME_E( tm_neighbor_ele )

!  TIME_S( tm_neighbor_hp )
!  if (inmisc%force_flag(INTERACT%HP)) then
!     ! make neighborlist for hydrophobic interaction
!     call neighbor_list_hp(irep)
!  endif
!  TIME_E( tm_neighbor_hp )

  if (inmisc%i_dtrna_model == 2015) then
     TIME_S( tm_neighbor_hb )
     call neighbor_list_hb(irep)
     TIME_E( tm_neighbor_hb )
  endif

  if (inmisc%i_BBR == 1) then
     call neighbor_list_BBR(irep)
  endif

#ifdef MPI_PAR
  if (local_rank_mpi == 0) then
#endif
  if (flg_file_out%neigh) then
     grep = irep2grep(irep)
     if (inmisc%force_flag(INTERACT%EXV12)) then
        write(outfile%neigh(grep),'(1xi6)',advance='no') lexv(2,E_TYPE%EXV12,irep) - lexv(1,E_TYPE%EXV12,irep) + 1
     endif
     if (inmisc%force_flag(INTERACT%EXV6)) then
        write(outfile%neigh(grep),'(1xi6)',advance='no') lexv(2,E_TYPE%EXV6,irep) - lexv(1,E_TYPE%EXV6,irep) + 1
     endif
     if (inmisc%force_flag(INTERACT%EXV_WCA)) then
        write(outfile%neigh(grep),'(1xi6)',advance='no') lexv(2,E_TYPE%EXV_WCA,irep) - lexv(1,E_TYPE%EXV_WCA,irep) + 1
     endif
     if (inmisc%force_flag(INTERACT%EXV_DT15)) then
        write(outfile%neigh(grep),'(1xi6)',advance='no') lexv(2,E_TYPE%EXV_DT15,irep) - lexv(1,E_TYPE%EXV_DT15,irep) + 1
     endif
     if (inmisc%force_flag(INTERACT%EXV_GAUSS)) then
        write(outfile%neigh(grep),'(1xi6)',advance='no') lexv(2,E_TYPE%EXV_GAUSS,irep) - lexv(1,E_TYPE%EXV_GAUSS,irep) + 1
     endif
     if (inmisc%i_dtrna_model == 2015) then
        write(outfile%neigh(grep),'(1xi6)',advance='no') nhbneigh(irep)
     endif
     if (inmisc%force_flag(INTERACT%ELE)) then
        write(outfile%neigh(grep),'(1xi6)',advance='no') lele(irep)
     endif
     if (inmisc%i_BBR == 1) then
        write(outfile%neigh(grep),'(1xi6)',advance='no') nbbr(irep)
     endif
     write(outfile%neigh(grep), *) ''
     flush(outfile%neigh(grep))
  endif
#ifdef MPI_PAR
  endif
#endif

#ifdef _DEBUG
  write(6,*) '####### end neighbor'
#endif
end subroutine neighbor
