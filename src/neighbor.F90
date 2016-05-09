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
  use var_inp,     only : inperi
  use var_setp,    only : inmisc, inele
  use var_struct,  only : lexv, nmp_all

  use time
  use mpiconst
  implicit none

  integer, intent(in) :: irep
  ! -------------------------------------------------------------------
  ! local variables
  integer :: lmp2neigh((nmp_l+nthreads-1)/nthreads  ,0:nthreads-1)
  integer :: ineigh2mp(MXMPNEIGHBOR*nmp_all/nthreads,0:nthreads-1)

  ! -------------------------------------------------------------------

  if(inperi%i_periodic == 1) then
     call util_periodic(irep)
  end if

  if (inmisc%force_flag(INTERACT%ENM)) then
     lexv(1, 1:E_TYPE%MAX,:) = 1
     lexv(2, 1:E_TYPE%MAX,:) = 0
     return
  end if

  TIME_S( tm_neighbor_exv )
  ! make neighborlist
  call neighbor_list(irep, ineigh2mp, lmp2neigh)

  ! assign iconcal2con and iexv2mp
  call neighbor_assign(irep, ineigh2mp, lmp2neigh)
  TIME_E( tm_neighbor_exv )
  
  TIME_S( tm_neighbor_ele )
  if (inmisc%force_flag(INTERACT%ELE)) then
     ! make neighborlist for electrostatic interaction
     if(inele%i_calc_method == 0) then
        call neighbor_list_ele(irep)
!     if(inele%i_calc_method == 0) then
     else if(inele%i_calc_method == 1) then
        call neighbor_list_ele2(irep)
     end if
  endif
  TIME_E( tm_neighbor_ele )

  TIME_S( tm_neighbor_hp )
  if (inmisc%force_flag(INTERACT%HP)) then
     ! make neighborlist for hydrophobic interaction
     call neighbor_list_hp(irep)
  endif
  TIME_E( tm_neighbor_hp )

  if (inmisc%i_dtrna_model == 2015) then
     call neighbor_list_hb(irep)
  endif

end subroutine neighbor
