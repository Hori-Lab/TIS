subroutine dump_var_struct(lunout)

  use const_maxsize
  use var_struct

  integer :: i,j,k
  integer :: ni,nj,nk
  integer, intent(in) :: lunout

  write(lunout,'(a)') ''
  write(lunout,'(a)') '### var_struct'
  write(lunout,*) 'nunit_real,',nunit_real
  write(lunout,*) 'nunit_all,', nunit_all
  write(lunout,*) 'nmp_real,',  nmp_real
  write(lunout,*) 'nmp_all,',   nmp_all
  write(lunout,*) 'nres,',       nres
  do i = 1, MXUNIT
     write(lunout,*) 'lunit2mp(:,',i,'),', lunit2mp(1, i),lunit2mp(2, i)
  enddo
  do i = 1, MXMP
     write(lunout,*) 'ires_mp(',i,'),', ires_mp(i)
  enddo
  do i = 1, MXUNIT
     write(lunout,*) 'iclass_unit(',i,'),', iclass_unit(i)
  enddo
  do i = 1, MXMP
     write(lunout,*) 'iclass_mp(',i,'),', iclass_mp(i)
  enddo

! xyz_mp_rep
  write(lunout,*) '# xyz_mp_rep'
  if (allocated(xyz_mp_rep)) then
     nj = ubound(xyz_mp_rep,2)
     nk = ubound(xyz_mp_rep,3)
     do k = 1, nk
        do j = 1, nj
           write(lunout,*) 'xyz_mp_rep(:,',j,',',k,'),'
           write(lunout,*) xyz_mp_rep(:,j,k)
        enddo
     enddo
  else
     write(lunout,*) '#xyz_mp_rep(:,:,:) is not allocated.'
  endif

  do i = 1, MXMP
     write(lunout,*) 'cmass_mp(',i,'),',cmass_mp(i)
  enddo
  do i = 1, MXMP
     write(lunout,*) 'fric_mp(',i,'),',fric_mp(i)
  enddo
  do i = 1, MXMP
     write(lunout,*) 'cmp2seq(',i,'),',cmp2seq(i)
  enddo
  do i = 1, MXMP
     write(lunout,*) 'cmp2atom(',i,'),',cmp2atom(i)
  enddo
  do i = 1, MXMP
     write(lunout,*) 'imp2type(',i,'),',imp2type(i)
  enddo
  do i = 1, MXMP
     write(lunout,*) 'imp2unit(',i,'),',imp2unit(i)
  enddo

! secondary structure
  do i = 1, MXMP
     write(lunout,*) 'istype_mp(',i,'),',istype_mp(i)
  enddo

  ! ----------------------------------------------------------------
  ! protein structure
  ! bond length
  write(lunout,*) ''
  write(lunout,*) '# bond length'
  write(lunout,*) 'nbd,',nbd
  do i = 1, MXMPBD*nmp_all
     write(lunout,*) 'ibd2mp(:,',i,'),',ibd2mp(1, i),ibd2mp(2,i)
  enddo
  do i = 1, MXMPBD*nmp_all
     write(lunout,*) 'bd_nat(',i,'),',bd_nat(i)
  enddo
  do i = 1, MXMPBD*nmp_all
     write(lunout,*) 'factor_bd(',i,'),',factor_bd(i)
  enddo
  do i = 1, MXMPBD*nmp_all
     write(lunout,*) 'coef_bd(:,',i,'),',coef_bd(1, i),coef_bd(2,i)
  enddo
  do i = 1, MXMPBD*nmp_all
     write(lunout,*) 'correct_bd_mgo(',i,'),',correct_bd_mgo(i)
  enddo

  ! bond angle
  write(lunout,*) ''
  write(lunout,*) '# bond angle'
  write(lunout,*) 'nba,',nba
  do i = 1, MXMPBA*nmp_all
     write(lunout,*) 'iba2mp(:,',i,'),',iba2mp(1, i),iba2mp(2,i),iba2mp(3,i)
  enddo
  do i = 1, MXMPBA*nmp_all
     write(lunout,*) 'ba_nat(',i,'),',ba_nat(i)
  enddo
  do i = 1, MXMPBA*nmp_all
     write(lunout,*) 'factor_ba(',i,'),',factor_ba(i)
  enddo
  do i = 1, MXMPBA*nmp_all
     write(lunout,*) 'coef_ba(:,',i,'),',coef_ba(1, i),coef_ba(2,i)
  enddo
  do i = 1, MXMPBA*nmp_all
     write(lunout,*) 'correct_ba_mgo(',i,'),',correct_ba_mgo(i)
  enddo

  ! dihedral angle
  write(lunout,*) ''
  write(lunout,*) '# dihedral angle'
  write(lunout,*) 'ndih,',ndih
  do i = 1, MXMPDIH*nmp_all
     write(lunout,*) 'idih2mp(:,',i,'),',idih2mp(1, i),idih2mp(2,i),idih2mp(3,i),idih2mp(4,i)
  enddo
  do i = 1, MXMPDIH*nmp_all
     write(lunout,*) 'dih_nat(',i,'),',dih_nat(i)
  enddo
  do i = 1, MXMPDIH*nmp_all
     write(lunout,*) 'factor_dih(',i,'),',factor_dih(i)
  enddo
  do i = 1, MXMPDIH*nmp_all
     write(lunout,*) 'coef_dih(:,',i,'),',coef_dih(1, i),coef_dih(2,i)
  enddo
  do i = 1, MXMPDIH*nmp_all
     write(lunout,*) 'dih_sin_nat(',i,'),',dih_sin_nat(i)
  enddo
  do i = 1, MXMPDIH*nmp_all
     write(lunout,*) 'dih_cos_nat(',i,'),',dih_cos_nat(i)
  enddo
  do i = 1, MXMPDIH*nmp_all
     write(lunout,*) 'correct_dih_mgo(',i,'),',correct_dih_mgo(i)
  enddo

  ! go
  write(lunout,*) ''
  write(lunout,*) '# go'
  write(lunout,*) 'ncon,',ncon
  do i = 1, MXCON
     write(lunout,*) 'icon2mp(:,',i,'),',icon2mp(1,i),icon2mp(2,i)
  enddo
  do i = 1, MXMP
     write(lunout,*) 'lmp2con(',i,'),',lmp2con(i)
  enddo
  do i = 1, MXCON
     write(lunout,*) 'icon2unit(:,',i,'),',icon2unit(1,i),icon2unit(2,i)
  enddo
  do i = 1, MXCON
     write(lunout,*) 'icon_dummy_mgo(',i,'),',icon_dummy_mgo(i)
  enddo
  do i = 1, MXCON
     write(lunout,*) 'go_nat(',i,'),',go_nat(i)
  enddo
  do i = 1, MXCON
     write(lunout,*) 'go_nat2(',i,'),',go_nat2(i)
  enddo
  do i = 1, MXCON
     write(lunout,*) 'factor_go(',i,'),',factor_go(i)
  enddo
  do i = 1, MXCON
     write(lunout,*) 'coef_go(',i,'),',coef_go(i)
  enddo

  ! qscore
  write(lunout,*) '# qsocre'
  do i = 1,MXUNIT
     do j = 1,MXUNIT
        write(lunout,*) 'ncon_unit(',i,',',j,'),',ncon_unit(i,j)
     enddo
  enddo

  ! electrostatic
  write(lunout,*) '# electrostatic'
  write(lunout,*) 'ncharge,',ncharge
  do i = 1, MXCHARGE
     write(lunout,*) 'icharge2mp(',i,'),',icharge2mp(i)
  enddo
  do i = 1, MXCHARGE
     nj = ubound(coef_charge,2)
     do j = 1, nj
        write(lunout,*) 'coef_charge(',i,',',j,'),',coef_charge(i,j)
     enddo
  enddo
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
  

  ! hydrophobic and binding site 
  write(lunout,*) '# hydrophobic'
  write(lunout,*) 'nhp,',nhp
  do i = 1, MXHP
     write(lunout,*) 'ihp2mp(',i,'),',ihp2mp(i)
  enddo
  do i = 1, MXUNIT
     write(lunout,*) 'lunit2hp(1:2,',i,'),',lunit2hp(1:2,i)
  enddo
  do i = 1, MXHP
     write(lunout,*) 'ncoor_hp(',i,'),',ncoor_hp(i)
  enddo
  do i = 1, MXHP
     write(lunout,*) 'ncoor_max_hp(',i,'),',ncoor_max_hp(i)
  enddo
  do i = 1, MXHP
     write(lunout,*) 'coef_aa_hp(',i,'),',coef_aa_hp(i)
  enddo
  if (allocated(nhpneigh)) then
     ni = ubound(nhpneigh,1)
     do i = 1, ni
        write(lunout,*) 'nhpneigh(',i,'),',nhpneigh(i)
     enddo
  else
     write(lunout,*) 'nhpneigh is not allocated.'
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

endsubroutine dump_var_struct
