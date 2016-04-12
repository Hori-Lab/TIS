subroutine dump_var_struct_neighbor(lunout)

  use const_maxsize
  use var_struct

  integer :: i,j,k
  integer :: ni,nj,nk

  write(lunout,*) ''
  write(lunout,*) '### var_struct neighbor list'

  ! electrostatic DNA
  write(lunout,*) '# neighbor list of electrostatic'
  if (allocated(lele)) then
     ni = ubound(lele,1)
     do i = 1, ni
        write(lunout,*) 'lele(',i,'),',lele(i)
     enddo
  else
     write(lunout,*) 'lele is not allocated.'
  endif
  if (allocated(iele2mp)) then
     nj = ubound(iele2mp,2)
     nk = ubound(iele2mp,3)
     do k = 1, nk
        do j = 1, nj
           write(lunout,*) 'iele2mp(:,',j,',',k,'),',iele2mp(:,j,k)
        enddo
     enddo
  else 
     write(lunout,*) 'iele2mp is not allocated.'
  endif
  if (allocated(coef_ele)) then
     ni = ubound(coef_ele,1)
     nj = ubound(coef_ele,2)
     do j = 1, nj
        do i = 1, ni
           write(lunout,*) 'coef_ele(',i,',',j,'),',coef_ele(i,j)
        enddo
     enddo
  else
     write(lunout,*) 'coef_ele is not allocated.'
  endif

  ! solvation DNA
  write(lunout,*) '# neighbor list of solvation DNA'
  if (allocated(lsolv)) then
     ni = ubound(lsolv,1)
     do i = 1, ni
        write(lunout,*) 'lele(',i,'),',lsolv(i)
     enddo
  else
     write(lunout,*) 'lsolv is not allocated.'
  endif
  if (allocated(isolv2mp)) then
     nj = ubound(isolv2mp,2)
     nk = ubound(isolv2mp,3)
     do k = 1, nk
        do j = 1, nj
           write(lunout,*) 'isolv2mp(:,',j,',',k,'),',isolv2mp(:,j,k)
        enddo
     enddo
  else
     write(lunout,*) 'isolv2mp is not allocated.'
  endif

  ! ----------------------------------------------------------------
  ! neighbor_list
  write(lunout,*) '# neighbor list'
  if (allocated(lpnl)) then
     nj = ubound(lpnl,2)
     nk = ubound(lpnl,3)
     do k = 1, nk
        do j = 1, nj
           write(lunout,*) 'lpnl(:,',j,',',k,'),',lpnl(:,j,k)
        enddo
     enddo
  else 
     write(lunout,*) 'lpnl is not allocated.'
  endif
  if (allocated(ipnl2mp)) then
     nj = ubound(ipnl2mp,2)
     nk = ubound(ipnl2mp,3)
     do k = 1, nk
        do j = 1, nj
           write(lunout,*) 'ipnl2mp(:,',j,',',k,'),',ipnl2mp(:,j,k)
        enddo
     enddo
  else 
     write(lunout,*) 'ipnl2mp is not allocated.'
  endif
  
  ! ----------------------------------------------------------------
  ! int lipid (Noguchi)
  write(lunout,*) '# lipid (Noguchi)'
  if (allocated(itail2mp)) then
     ni = ubound(itail2mp,1)
     nj = ubound(itail2mp,2)
     do j = 1, nj
        do i = 1, ni
           write(lunout,*) 'itail2mp(',i,',',j,'),',itail2mp(i,j)
        enddo
     enddo
  else
     write(lunout,*) 'itail2mp is not allocated.'
  endif
  if (allocated(ltail2mp)) then
     nj = ubound(ltail2mp,2)
     nk = ubound(ltail2mp,3)
     do k = 1, nk
        do j = 1, nj
           write(lunout,*) 'ltail2mp(:,',j,',',k,'),',ltail2mp(:,j,k)
        enddo
     enddo
  else 
     write(lunout,*) 'ltail2mp is not allocated.'
  endif

  ! -------------------------------------------------------------------
  ! hydrophobic
  write(lunout,*) '# hydrophobic'
  if (allocated(nhpneigh)) then
     ni = ubound(nhpneigh,1)
     do i = 1, ni
        write(lunout,*) 'nhpneigh(',i,'),',nhpneigh(i)
     enddo
  else
     write(lunout,*) 'nhpneigh is not allocated.'
  endif
  if (allocated(lhp2neigh)) then
     nj = ubound(lhp2neigh,2)
     nk = ubound(lhp2neigh,3)
     do k = 1, nk
        do j = 1, nj
           write(lunout,*) 'lhp2neigh(:,',j,',',k,'),',lhp2neigh(:,j,k)
        enddo
     enddo
  else 
     write(lunout,*) 'lhp2neigh is not allocated.'
  endif
  if (allocated(ineigh2hp)) then
     nj = ubound(ineigh2hp,1)
     nk = ubound(ineigh2hp,2)
     do k = 1, nk
        do j = 1, nj
           write(lunout,*) 'ineigh2hp(',j,',',k,'),',ineigh2hp(j,k)
        enddo
     enddo
  else 
     write(lunout,*) 'ineigh2hp is not allocated.'
  endif
  if (allocated(cutoff_dmin_hp)) then
     ni = ubound(cutoff_dmin_hp,1)
     nj = ubound(cutoff_dmin_hp,2)
     do j = 1, nj
        do i = 1, ni
           write(lunout,*) 'cutoff_dmin_hp(',i,',',j,'),',cutoff_dmin_hp(i,j)
        enddo
     enddo
  else
     write(lunout,*) 'cutoff_dmin_hp is not allocated.'
  endif
  if (allocated(cutoff_dmax_hp)) then
     ni = ubound(cutoff_dmax_hp,1)
     nj = ubound(cutoff_dmax_hp,2)
     do j = 1, nj
        do i = 1, ni
           write(lunout,*) 'cutoff_dmax_hp(',i,',',j,'),',cutoff_dmax_hp(i,j)
        enddo
     enddo
  else
     write(lunout,*) 'cutoff_dmax_hp is not allocated.'
  endif
endsubroutine dump_var_struct_neighbor
