! simu_rmsd
!> @brief Calculate the root mean squared deviation (RMSD)

subroutine simu_rmsd(rmsd_unit, rmsd)

  use const_maxsize
  use const_physical
  use var_struct,  only : nunit_all, nmp_real, lunit2mp, &
                          xyz_mp_rep, xyz_ref_mp
  use var_replica, only : n_replica_mpi

  implicit none
  
  ! -------------------------------------------------------------------
  real(PREC), intent(out) :: rmsd(:)               ! (replica)
  real(PREC), intent(out) :: rmsd_unit(:,:)        ! (unit, replica)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: imp, iunit, irep, num, i1, i2
  integer :: list1(nmp_real), list2(nmp_real)
  real(PREC) :: xyz_init_mp(3, nmp_real), xyz_in_mp(3, nmp_real), xyz_out_mp(3, nmp_real)

  ! -------------------------------------------------------------------
  do imp = 1, nmp_real
     list1(imp) = imp
     list2(imp) = imp
  end do

  do irep = 1, n_replica_mpi

     i1 = 1
     i2 = nmp_real
     num = i2 - i1 + 1
     xyz_init_mp(1:3, i1:i2) = xyz_ref_mp(1:3, i1:i2)
     xyz_in_mp(1:3, i1:i2) = xyz_mp_rep(1:3, i1:i2, irep)
     call util_bestfit(xyz_init_mp, num, xyz_in_mp, num, num, xyz_out_mp, &
          list1, list2, rmsd(irep))

     do iunit = 1, nunit_all
        i1 = lunit2mp(1, iunit)
        i2 = lunit2mp(2, iunit)
        num = i2 - i1 + 1
        xyz_init_mp(1:3, 1:num) = xyz_ref_mp(1:3, i1:i2)
        xyz_in_mp(1:3, 1:num) = xyz_mp_rep(1:3, i1:i2, irep)
        call util_bestfit(xyz_init_mp, num, xyz_in_mp, num, num, xyz_out_mp, &
             list1, list2, rmsd_unit(iunit, irep))
     end do

  end do

end subroutine simu_rmsd
