! force_box
!> @brief Calculate force of box option

! ****************************************************************
! ****************************************************************
subroutine force_box(irep, force_mp)

  use const_maxsize
  use var_setp,   only : inmisc
  use var_struct, only : nmp_real, xyz_mp_rep, nmp_all
  implicit none

  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  ! ------------------------------------------------------------
  ! local variables
  integer :: imp, idimn
  real(PREC) :: boxsigma, thre(2), boxsize(3)
  real(PREC) :: fbox, dbox, boxsign, coef, kbox
  real(PREC) :: rdbox, rdbox2, rdbox4, rdbox8, rdbox13
  real(PREC) :: crdbox, crdbox2, crdbox4, crdbox13
  
  ! ------------------------------------------------------------
  boxsize(1) = inmisc%xbox * 0.5e0_PREC
  boxsize(2) = inmisc%ybox * 0.5e0_PREC
  boxsize(3) = inmisc%zbox * 0.5e0_PREC

  kbox = 10.0e0_PREC
  boxsigma = inmisc%boxsigma

  coef = 12.0e0_PREC * kbox / boxsigma
  thre(1) = 3.0e0_PREC * boxsigma
!  thre(2) = 0.5e0_PREC * boxsigma
  thre(2) = 0.8e0_PREC * boxsigma
  crdbox = boxsigma / thre(2)
  crdbox2 = crdbox * crdbox
  crdbox4 = crdbox2 * crdbox2
  crdbox13 =  coef * crdbox * crdbox4 * crdbox4 * crdbox4

  ! ------------------------------------------------------------
  do imp = 1, nmp_real
     do idimn = 1, 3
        if(xyz_mp_rep(idimn, imp, irep) > 0.0e0_PREC) then
           dbox = boxsize(idimn) - xyz_mp_rep(idimn, imp, irep)
           boxsign = 1.0e0_PREC
        else
           dbox = xyz_mp_rep(idimn, imp, irep) + boxsize(idimn)
           boxsign = -1.0e0_PREC
        end if

        if(dbox < thre(1)) then
           if(dbox < thre(2)) then
              fbox = - crdbox13
           else
              rdbox = boxsigma / dbox
              rdbox2 = rdbox * rdbox
              rdbox4 = rdbox2 * rdbox2 
              rdbox8 = rdbox4 * rdbox4
              rdbox13 = rdbox * rdbox4 * rdbox8
              fbox = - coef*rdbox13
           end if

           force_mp(idimn, imp) = force_mp(idimn, imp) + boxsign * fbox
        end if
     end do
  end do
      
  return

end subroutine force_box
