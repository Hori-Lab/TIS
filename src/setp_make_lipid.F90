! setp_make_lipid
!> @brief  Constructs the initial structure for lipid, and &
!>         stores it into the argument array.

! ****************************************************************
subroutine setp_make_lipid(xyz_lip)

  use const_maxsize
  use const_index
  use var_setp,   only : inlip
  use var_struct, only : nmp_real, nmp_all, imp2unit, lunit2mp,  &
                         nres, ires_mp, cmp2seq, cmp2atom, iclass_mp
  implicit none
  
  ! -------------------------------------------------------------
  real(PREC), intent(inout) :: xyz_lip(3, MXMP)

  ! -------------------------------------------------------------
  integer :: imp, jmp, kmp, lmp, iunit, ires
  integer :: lcore, lint
  integer :: ilayer
  real(PREC) :: sigma
  real(PREC) :: xyz_frame(3)

  ! -------------------------------------------------------------
  sigma = inlip%sigma_lipid
  lcore = inlip%num_lip_core
  lint = lcore + inlip%num_lip_int

  ! -------------------------------------------------------------
  lmp  = 0
  ires = 0

  ! -------------------------------------------------------------
  iunit = 1
  lunit2mp(1, iunit) = lmp + 1
  do ilayer = 1, inlip%nlayer_lipid

     xyz_frame(1) = -0.5 * ((inlip%nmp_transverse_lipid - 1) * sigma)
     xyz_frame(3) = inlip%z_coord_lipid(ilayer) * sigma
     do imp = 1, inlip%nmp_transverse_lipid
        xyz_frame(2) = -0.5 * ((inlip%nmp_longitudinal_lipid - 1) * sigma)
        do jmp = 1, inlip%nmp_longitudinal_lipid
           ires = ires + 1
           do kmp = 1, inlip%num_lip_total
              lmp = lmp + 1
              xyz_lip(1, lmp) = xyz_frame(1)
              xyz_lip(2, lmp) = xyz_frame(2)
              xyz_lip(3, lmp) = xyz_frame(3) - kmp * sigma

              if(mod(ilayer,2) == 0)then
                 xyz_lip(3, lmp) = xyz_frame(3) + kmp * sigma
              else if(mod(ilayer,2) == 1)then
                 xyz_lip(3, lmp) = xyz_frame(3) - kmp * sigma
              end if

              imp2unit(lmp) = iunit
              ires_mp(lmp) = ires
              if(kmp <= lcore) then
                 cmp2seq(lmp) = 'COR'
                 cmp2atom(lmp) = ' C  '
              else if(kmp <= lint) then
                 cmp2seq(lmp) = 'INT'
                 cmp2atom(lmp) = ' O  '
              else
                 cmp2seq(lmp) = 'TAI'
                 cmp2atom(lmp) = ' N  '
              end if
              iclass_mp(lmp) = CLASS%LIP
           end do
           xyz_frame(2) = xyz_frame(2) + inlip%grid_size_lipid*sigma
        end do
        xyz_frame(1) = xyz_frame(1) + inlip%grid_size_lipid*sigma
     end do

  end do
  lunit2mp(2, iunit) = lmp

  ! -------------------------------------------------------------
  nmp_real = lmp
  nmp_all  = lmp
  nres = ires

end subroutine setp_make_lipid
