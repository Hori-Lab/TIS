!simu_mc_implig
!> @brief Changes the state of implicit_ligand(istate_implig) by MC (Monte Calro).

subroutine simu_mc_implig(irep, istep, tempk)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : mts
  use var_implig, only : inimplig, Etbind_implig, istate_implig, p_implig
  use mt_stream
  implicit none

!-------------------------------------------------  
  integer,    intent(in) :: irep
  integer(L_INT), intent(in) :: istep
  real(PREC), intent(in) :: tempk

!-------------------------------------------------
  integer :: ist, istream
  integer :: rdm_offset
  real(PREC) :: p
! -----------------------------------------------
  
  rdm_offset = 0

  do ist=1,inimplig%nsite_implig
     !! generate [0:1] uniform random number: grnd().
     ! p = grnd()
     istream = irep
     p = genrand_double1(mts(istream, 0))

     !! for binding
     if(istate_implig(ist, irep)== IMPLIGBOUND_STATE%UN_BOUND .and. &
          mod(istep,inimplig%istep_implig)==0)then
 
        ! for debug
        ! write(*,*)'istep=',istep,'p=',p,'p_implig=',p_implig(ist)
        if(p<p_implig(ist))then
           if(Etbind_implig(ist, irep) < 0.0_PREC)then
              !istate_implig(ist, irep) = 1
              istate_implig(ist, irep) = IMPLIGBOUND_STATE%BOUND
           else
              if(p<exp(-Etbind_implig(ist, irep)/(BOLTZ_KCAL_MOL*tempk)))then
                 !istate_implig(ist, irep) = 1
                 istate_implig(ist, irep) = IMPLIGBOUND_STATE%BOUND
              endif
           endif
        endif

     !! for unbinding
     elseif(istate_implig(ist, irep)== IMPLIGBOUND_STATE%BOUND .and. &
          mod(istep,inimplig%istep_un_implig)==0)then

        if(Etbind_implig(ist, irep)>0.0_PREC)then
           !istate_implig(ist, irep) = 0
           istate_implig(ist, irep) =  IMPLIGBOUND_STATE%UN_BOUND
        else
           if(p<exp(Etbind_implig(ist, irep)/(BOLTZ_KCAL_MOL*tempk)))then
              !istate_implig(ist, irep) = 0
              istate_implig(ist, irep) = IMPLIGBOUND_STATE%UN_BOUND
           endif
        endif

     endif
     !     if(istate_implig(ist, irep)==0)then
     !        Etbind_implig(ist, irep)=0.0e0_PREC
     !     endif
  enddo
  
  ! for debug
  !  do ist=1,inimplig%nsite_implig
  !     write(*,*)ist,'istate_implig=',istate_implig(ist, irep),'Etbind_implig=',Etbind_implig(ist, irep)
  !  enddo
  !  return
end subroutine simu_mc_implig
