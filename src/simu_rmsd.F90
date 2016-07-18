! simu_rmsd
!> @brief Calculate the root mean squared deviation (RMSD)

subroutine simu_rmsd(rmsd_unit, rmsd)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : fix_mp
  use var_struct,  only : nunit_all, nmp_real, lunit2mp, &
                          xyz_mp_rep, xyz_ref_mp, iclass_unit
  use var_replica, only : n_replica_mpi

  implicit none
  
  real(PREC), intent(out) :: rmsd(:)               ! (replica)
  real(PREC), intent(out) :: rmsd_unit(:,:)        ! (unit, replica)

  integer :: imp, iunit, irep, num, i1, i2
  integer :: list1(nmp_real), list2(nmp_real)
  real(PREC) :: xyz_init_mp(3, nmp_real), xyz_in_mp(3, nmp_real), xyz_out_mp(3, nmp_real)
  logical :: flg_all

  ! -------------------------------------------------------------------

  flg_all = .True.

  do imp = 1, nmp_real
     list1(imp) = imp
     list2(imp) = imp
  end do

  do irep = 1, n_replica_mpi

     do iunit = 1, nunit_all
        if (iclass_unit(iunit) == CLASS%ION) then
            rmsd_unit(iunit,irep) = 0.0
            flg_all = .False.
        else
           i1 = lunit2mp(1, iunit)
           i2 = lunit2mp(2, iunit)

           if (all(fix_mp(i1:i2))) then
              rmsd_unit(iunit,irep) = 0.0
              flg_all = .False.
           else
              num = i2 - i1 + 1
              xyz_init_mp(1:3, 1:num) = xyz_ref_mp(1:3, i1:i2)
              xyz_in_mp(1:3, 1:num) = xyz_mp_rep(1:3, i1:i2, irep)
              call util_bestfit(xyz_init_mp, num, xyz_in_mp, num, num, xyz_out_mp, &
                   list1, list2, rmsd_unit(iunit, irep))
           endif
        endif
     end do

     if (flg_all) then
        i1 = 1
        i2 = nmp_real
        num = i2 - i1 + 1
        xyz_init_mp(1:3, i1:i2) = xyz_ref_mp(1:3, i1:i2)
        xyz_in_mp(1:3, i1:i2) = xyz_mp_rep(1:3, i1:i2, irep)
        call util_bestfit(xyz_init_mp, num, xyz_in_mp, num, num, xyz_out_mp, &
             list1, list2, rmsd(irep))
     else
        rmsd(irep) = 0.0
     endif

  end do

end subroutine simu_rmsd
