! energy_rest1d
!> @brief Calculate energy of 1D-restraint option

! ************************************************************************
subroutine energy_rest1d(irep, e_exv_unit, e_exv)

  use const_maxsize
  use const_index
  use const_physical
  use var_setp, only : inmisc
  use var_struct, only : xyz_mp_rep, imp2unit
  implicit none

  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)
  real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: irest, irest_center, imp, iunit, junit
  real(PREC) :: norm, dot, efull
  real(PREC) :: xyz_sa(SPACE_DIM), xyz_sb(SPACE_DIM)
  real(PREC) :: v_ab(SPACE_DIM), v_ai(SPACE_DIM)
  real(PREC) :: xyz_com(SPACE_DIM)
  
  ! ----------------------------------------------------------------------
  do irest = 1, inmisc%nrest1d
   
     imp = inmisc%irest1d2mp(irest)
     xyz_sa(:) = xyz_mp_rep(:,inmisc%irest1d_mp_sa(irest),irep)
     xyz_sb(:) = xyz_mp_rep(:,inmisc%irest1d_mp_sb(irest),irep)

     v_ab = xyz_sb - xyz_sa
     v_ai = xyz_mp_rep(:,imp,irep) - xyz_sa

     norm = sqrt((xyz_sb(1) - xyz_sa(1)) **2  &
                +(xyz_sb(2) - xyz_sa(2)) **2  &
                +(xyz_sb(3) - xyz_sa(3)) **2)

     inmisc%rest1d_s(irest) = dot_product(v_ab,v_ai) / norm

     efull = inmisc%coef_rest1d(irest) &
           * (inmisc%rest1d_s(irest) - inmisc%rest1d_s0(irest)) ** 2

     e_exv(E_TYPE%REST1D) = e_exv(E_TYPE%REST1D) + efull
     iunit = imp2unit(imp)
     junit = iunit
     e_exv_unit(iunit, junit, E_TYPE%REST1D) =  &
          e_exv_unit(iunit, junit, E_TYPE%REST1D) + efull

  end do

  do irest_center = 1, inmisc%nrest1d_center

     if (inmisc%rest1d_center_init_flag(irest_center) == 1) then
        xyz_sa(:) = xyz_mp_rep(:,inmisc%irest1d_center_mp_sa(irest_center),irep)
        xyz_sb(:) = xyz_mp_rep(:,inmisc%irest1d_center_mp_sb(irest_center),irep)
        norm = sqrt((xyz_sb(1) - xyz_sa(1)) **2 &
                   +(xyz_sb(2) - xyz_sa(2)) **2 &
                   +(xyz_sb(3) - xyz_sa(3)) **2)
        v_ab = xyz_sb - xyz_sa
        v_ab(:) = v_ab(:) / norm
        inmisc%rest1d_center_v(irest_center, :) = v_ab(:)
        inmisc%rest1d_center_origin(irest_center, :) = xyz_sa(:)
        inmisc%rest1d_center_init_flag(irest_center) = 0
     end if

     xyz_com(:) = 0.0
     do imp = 1, inmisc%nrest1d_center_mp(irest_center)
        xyz_com(1) = xyz_com(1) + xyz_mp_rep(1, inmisc%irest1d_center2mp(irest_center, imp), irep) / real( inmisc%nrest1d_center_mp(irest_center) )
        xyz_com(2) = xyz_com(2) + xyz_mp_rep(2, inmisc%irest1d_center2mp(irest_center, imp), irep) / real( inmisc%nrest1d_center_mp(irest_center) )
        xyz_com(3) = xyz_com(3) + xyz_mp_rep(3, inmisc%irest1d_center2mp(irest_center, imp), irep) / real( inmisc%nrest1d_center_mp(irest_center) )
     end do

     v_ai = xyz_com(:) - inmisc%rest1d_center_origin(irest_center, :)

     dot = inmisc%rest1d_center_v(irest_center, 1) * v_ai(1) &
         + inmisc%rest1d_center_v(irest_center, 2) * v_ai(2) &
         + inmisc%rest1d_center_v(irest_center, 3) * v_ai(3) 

     efull = inmisc%coef_rest1d_center(irest_center) &
             * ( dot - inmisc%rest1d_center_s0(irest_center) ) ** 2

     e_exv(E_TYPE%REST1D) = e_exv(E_TYPE%REST1D) + efull
     !iunit = imp2unit(imp)
     !junit = iunit
     !e_exv_unit(iunit, junit, E_TYPE%REST1D) =  &
     !e_exv_unit(iunit, junit, E_TYPE%REST1D) + efull
     
  end do

end subroutine energy_rest1d
