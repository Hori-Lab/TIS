! energy_ele_coulomb
!> @brief Calculate the energy of electrostatic interaction 

subroutine energy_ele_coulomb_brute(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_replica,only : irep2grep
  use var_setp,   only : inele, inperi
  use var_struct, only : icharge2mp, coef_charge, imp2unit, pxyz_mp_rep, ncharge
  use var_simu,   only : ewld_f_n, ewld_h
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: icharge, jcharge, imp, jmp, iunit, junit, ih, grep
  real(PREC) :: coef, coef_i, ene
  real(PREC) :: v21(SDIM), pv21(SDIM)

  grep = irep2grep(irep)

!  do icharge = 1, ncharge
!     imp = icharge2mp(icharge)
!     iunit = imp2unit(imp)
!     coef_i = coef_charge(icharge,grep) * inele%coef(grep)
!
!     do jcharge = 1, ncharge
!        jmp = icharge2mp(jcharge)
!        junit = imp2unit(jmp)
!        coef = coef_i * coef_charge(jcharge, grep)
!
!        v21(1:SDIM) = pxyz_mp_rep(1:SDIM, jmp, irep) - pxyz_mp_rep(1:SDIM, imp, irep)
!
!        ene = 0.0e0_PREC
!        ! Self
!        if (jmp > imp) then
!           ene = ene + coef / sqrt(dot_product(v21,v21))
!        endif
!
!        ! Periodic
!        do ih = 1, ewld_f_n
!           pv21(1:SDIM) = v21(:) + inperi%psize(:) * ewld_h(:,ih)
!           ene = ene + 0.5 * coef / sqrt(dot_product(pv21,pv21))
!        enddo
!
!        energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene
!        energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + ene
!     enddo
!
!  end do
!
!  ene = 0.0e0_PREC
!  do ih = 1, ewld_f_n
!     do icharge = 1, ncharge
!        imp = icharge2mp(icharge)
!        !iunit = imp2unit(imp)
!        coef_i = coef_charge(icharge,grep) * inele%coef(grep)
!
!        do jcharge = 1, ncharge
!           jmp = icharge2mp(jcharge)
!           !junit = imp2unit(jmp)
!           coef = coef_i * coef_charge(jcharge, grep)
!
!           v21(1:SDIM) = pxyz_mp_rep(1:SDIM, jmp, irep) - pxyz_mp_rep(1:SDIM, imp, irep)
!
!           ! Self
!           if (ih == 1 .and. jmp /= imp) then
!              ene = ene + coef / sqrt(dot_product(v21,v21))
!           endif
!
!           ! Periodic
!           pv21(1:SDIM) = v21(:) + inperi%psize(:) * ewld_h(:,ih)
!           ene = ene + coef / sqrt(dot_product(pv21,pv21))
!        enddo
!     enddo
!  end do
!  energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + 0.5 * ene
!  !energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + 0.5 * ene

  !
  ! This would be fastest among three
  !
!!$omp do private(imp,jmp,iunit,junit,coef_i,coef,v21,pv21,ene,iunit,junit)
  do icharge = 1, ncharge
     imp = icharge2mp(icharge)
     iunit = imp2unit(imp)
     coef_i = coef_charge(icharge,grep) * inele%coef(grep)

     do jcharge = icharge, ncharge
        jmp = icharge2mp(jcharge)
        junit = imp2unit(jmp)
        coef = coef_i * coef_charge(jcharge, grep)

        v21(1:SDIM) = pxyz_mp_rep(1:SDIM, jmp, irep) - pxyz_mp_rep(1:SDIM, imp, irep)

        ! Self
        ene = 0.0
        if (jcharge > icharge) then
           ene = coef / sqrt(dot_product(v21,v21))
        endif

        ! Periodic
        do ih = 1, ewld_f_n
           pv21(1:SDIM) = v21(:) + inperi%psize(:) * ewld_h(:,ih)
           ene = ene + 0.5 * coef / sqrt(dot_product(pv21,pv21))

           if (jcharge > icharge) then
              pv21(1:SDIM) = -v21(:) + inperi%psize(:) * ewld_h(:,ih)
              ene = ene + 0.5 * coef / sqrt(dot_product(pv21,pv21))
           endif
        enddo

        energy(E_TYPE%ELE) = energy(E_TYPE%ELE) + ene
        energy_unit(iunit, junit, E_TYPE%ELE) = energy_unit(iunit, junit, E_TYPE%ELE) + ene
     enddo
  end do
!!!!$omp end do nowait

end subroutine energy_ele_coulomb_brute
