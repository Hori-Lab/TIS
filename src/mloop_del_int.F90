! mloop_del_int
!> @brief This subroutine is deleting interaction

! *****************************************************************************
subroutine mloop_del_int()

  use const_maxsize
  !use var_io, only : ius2unit
  use var_setp, only : inmisc
  use var_struct, only : nba, iba2mp, factor_ba, coef_ba, &
                         ndih, idih2mp, factor_dih, coef_dih, &
                         ncon, icon2mp, coef_go
!  use var_mgo, only : inmgo, ishadow2real_mp_mgo
  use mpiconst

  implicit none

  integer :: idel, ini, las, ini2, las2
  integer :: iba, idih, icon !,irna_bp, imorse
  integer :: imp1, imp2, imp3, imp4
!  integer :: imp_real1, imp_real2
  integer :: imp
!  integer :: imp_shadow, iunit_real, iunit_shadow, instate
  integer :: idel_mp(MXMP)

  ! -----------------------------------------------------------------

  idel_mp(:) = 0

  if(inmisc%ndel_lgo > 0) then
     
     do idel = 1, inmisc%ndel_lgo
        ini = inmisc%idel_lgo(1, idel)
        las = inmisc%idel_lgo(2, idel)
        do imp = ini, las
           idel_mp(imp) = 1

!           if(inmgo%i_multi_mgo >= 1) then
!              iunit_real = imp2unit(imp)
!              do instate = 2, MXSTATE_MGO
!                 iunit_shadow = ius2unit(iunit_real, instate)
!                 if(iunit_shadow <= 1) then
!                    exit
!                 else
!                    imp_shadow = imp + lunit2mp(1, iunit_shadow) - lunit2mp(1, iunit_real)
!
!                    idel_mp(imp_shadow) = 1
!                 end if
!              end do
!           end if

        end do
     end do


     ! del bond angle interaction
     do iba = 1, nba
        imp1 = iba2mp(1, iba)
        imp2 = iba2mp(2, iba)
        imp3 = iba2mp(3, iba)

        if(idel_mp(imp1) == 1 .or. idel_mp(imp2) == 1 .or. idel_mp(imp3) == 1) then
           factor_ba(iba) = 0.0e0_PREC
           coef_ba(1, iba) = 0.0e0_PREC
           coef_ba(2, iba) = 0.0e0_PREC
        end if
     end do

     ! del dihedral angle interaction
     do idih = 1, ndih
        imp1 = idih2mp(1, idih)
        imp2 = idih2mp(2, idih)
        imp3 = idih2mp(3, idih)
        imp4 = idih2mp(4, idih)

        if(idel_mp(imp1) == 1 .or. idel_mp(imp2) == 1 .or. idel_mp(imp3) == 1 .or. idel_mp(imp4) == 1) then
           factor_dih(idih) = 0.0e0_PREC
           coef_dih(1, idih) = 0.0e0_PREC
           coef_dih(2, idih) = 0.0e0_PREC
        end if
     end do
  end if


  ! ----------------------------------------------------------------------
  ! del nonlocal Go interaction

  if(inmisc%ndel_go > 0) then
     do idel = 1, inmisc%ndel_go
        ini = inmisc%idel_go(1, idel)
        las = inmisc%idel_go(2, idel)
        ini2 = inmisc%idel_go(3, idel)
        las2 = inmisc%idel_go(4, idel)

        do icon = 1, ncon
           imp1 = icon2mp(1, icon)
           imp2 = icon2mp(2, icon)
           if((imp1 >= ini  .and. imp1 <= las .and. &
                imp2 >= ini2 .and. imp2 <= las2) .or. &
                (imp1 >= ini2 .and. imp1 <= las2 .and. &
                imp2 >= ini  .and. imp2 <= las)) then
              coef_go(icon) = 0.0
           end if

!           if(inmgo%i_multi_mgo >= 1) then
!              imp_real1 = ishadow2real_mp_mgo(imp1)
!              imp_real2 = ishadow2real_mp_mgo(imp2)
!              if((imp_real1 >= ini  .and. imp_real1 <= las .and. &
!                   imp_real2 >= ini2 .and. imp_real2 <= las2) .or. &
!                   (imp_real1 >= ini2 .and. imp_real1 <= las2 .and. &
!                   imp_real2 >= ini  .and. imp_real2 <= las)) then
!                 coef_go(icon) = 0.0
!              end if
!           end if
        end do

!        do irna_bp = 1, nrna_bp
!           imp1 = irna_bp2mp(1, irna_bp)
!           imp2 = irna_bp2mp(2, irna_bp)
!           if((imp1 >= ini  .and. imp1 <= las .and. &
!                imp2 >= ini2 .and. imp2 <= las2) .or. &
!                (imp1 >= ini2 .and. imp1 <= las2 .and. &
!                imp2 >= ini  .and. imp2 <= las)) then
!              coef_rna_bp(irna_bp) = 0.0
!              coef_rna_bp_a(irna_bp) = 0.0
!              coef_rna_bp_fD(irna_bp) = 0.0
!           end if
!        end do

!        do imorse = 1, nmorse
!           imp1 = imorse2mp(1, imorse)
!           imp2 = imorse2mp(2, imorse)
!           if((imp1 >= ini  .and. imp1 <= las .and. &
!                imp2 >= ini2 .and. imp2 <= las2) .or. &
!                (imp1 >= ini2 .and. imp1 <= las2 .and. &
!                imp2 >= ini  .and. imp2 <= las)) then
!              coef_morse_a(imorse) = 0.0
!              coef_morse_fD(imorse) = 0.0
!           end if
!        end do

     end do
  end if

end subroutine mloop_del_int
