! simu_energy_box
!> @brief Calculates the boxing energy, when ``i_in_box=1" in ``<<<<md_information" block.

! ****************************************************************
subroutine simu_energy_box(irep, pnle_unit, pnlet)

  use const_maxsize
  use const_index
  use var_setp,    only : inmisc
  use var_struct,  only : nmp_real, xyz_mp_rep, imp2unit
  use var_replica, only : n_replica_mpi

  implicit none
  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: pnlet(:)         ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: pnle_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)

  ! ------------------------------------------------------------
  ! local variables
  integer :: imp, idimn, iunit, junit
  real(PREC) :: boxsigma, thre(2), boxsize(3)
  real(PREC) :: dbox, coef, kbox, efull
  real(PREC) :: rdbox, rdbox2, rdbox4, rdbox8, rdbox12
  real(PREC) :: crdbox, crdbox2, crdbox4, crdbox12, emin
  
  ! ------------------------------------------------------------
  boxsize(1) = inmisc%xbox * 0.5e0_PREC
  boxsize(2) = inmisc%ybox * 0.5e0_PREC
  boxsize(3) = inmisc%zbox * 0.5e0_PREC

  kbox = 10.0e0_PREC
  boxsigma = inmisc%boxsigma

  coef = kbox
  thre(1) = 3.0e0_PREC * boxsigma
 ! thre(2) = 0.5e0_PREC * boxsigma
  thre(2) = 0.8e0_PREC * boxsigma
  crdbox = boxsigma / thre(2)
  crdbox2 = crdbox * crdbox
  crdbox4 = crdbox2 * crdbox2
  crdbox12 = coef * crdbox4 * crdbox4 * crdbox4
 !calculating bottom energy to subtract 
  crdbox = boxsigma / thre(1)
  crdbox2 = crdbox * crdbox
  crdbox4 = crdbox2 * crdbox2
  emin = coef * crdbox4 * crdbox4 * crdbox4
  
  do imp = 1, nmp_real
     do idimn = 1, 3
  
        if(xyz_mp_rep(idimn, imp, irep) > 0.0e0_PREC) then
           dbox = boxsize(idimn) - xyz_mp_rep(idimn, imp, irep)
        else
           dbox = xyz_mp_rep(idimn, imp, irep) + boxsize(idimn)
        end if

        if(dbox < thre(1)) then
           if(dbox < thre(2)) then
              !efull = linear - bottom
              efull = crdbox12 * (1.0e0_PREC + 12.0e0_PREC * (thre(2) - dbox) / thre(2)) - emin
           else
              !efull = 1/r**12 - bottom
              rdbox = boxsigma / dbox
              rdbox2 = rdbox * rdbox
              rdbox4 = rdbox2 * rdbox2 
              rdbox8 = rdbox4 * rdbox4
              rdbox12 = rdbox4 * rdbox8
              efull = coef*rdbox12 - emin
           end if

           pnlet(E_TYPE%BOX) = pnlet(E_TYPE%BOX) + efull
           iunit = imp2unit(imp)
           junit = iunit
           pnle_unit(iunit, junit, E_TYPE%BOX) = pnle_unit(iunit, junit, E_TYPE%BOX) + efull
        end if

     end do
  end do

end subroutine simu_energy_box
