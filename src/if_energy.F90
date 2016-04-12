!if_energy
!> @brief Interface module for subroutines which called by simu_energy.

module if_energy
interface

   subroutine simu_energy_implig(irep, pnle_unit, pnlet, iflag_for_mc)   
     use const_maxsize
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: pnle_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)
     real(PREC), intent(inout) :: pnlet(:)         ! (E_TYPE%MAX)
     integer,    intent(in)    :: iflag_for_mc 
   endsubroutine simu_energy_implig

   subroutine simu_energy_anchor(irep, pnle_unit, pnlet)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: pnle_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: pnlet(:)         ! (E_TYPE%MAX, replica)
   endsubroutine simu_energy_anchor

   subroutine simu_energy_rest1d(irep, pnle_unit, pnlet)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: pnle_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: pnlet(:)         ! (E_TYPE%MAX, replica)
   endsubroutine simu_energy_rest1d
   
   subroutine simu_energy_bangle(irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
      real(PREC), intent(inout) :: pnlet(:)
   endsubroutine simu_energy_bangle
   
   subroutine simu_energy_bond  (irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_bond
   
   subroutine simu_energy_box(irep, pnle_unit, pnlet)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: pnlet(:)
     real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_box
   
   subroutine simu_energy_bridge(irep, pnle_unit, pnlet)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: pnle_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: pnlet(:)         ! (E_TYPE%MAX, replica)
   endsubroutine simu_energy_bridge
   
   subroutine simu_energy_dih   (irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_dih
  
   subroutine simu_energy_dih_harmonic (irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_dih_harmonic

   subroutine simu_energy_rna_stack(irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_rna_stack

   subroutine simu_energy_dtrna_stack(irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_dtrna_stack

   subroutine simu_energy_dtrna_hbond13(irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_dtrna_hbond13

   subroutine simu_energy_dtrna_hbond15(irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_dtrna_hbond15

   subroutine simu_energy_dtrna_stack_nlocal(irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_dtrna_stack_nlocal

   subroutine simu_energy_exv_wca (irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: pnlet(:)
      real(PREC), intent(out) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_exv_wca

   subroutine simu_energy_exv_dt15 (irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: pnlet(:)
      real(PREC), intent(out) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_exv_dt15
   
   subroutine simu_energy_enm(irep, now_con, pnlet, pnle_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_con(:,:)
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_enm
   
   subroutine simu_energy_mgo(pnle_unit, pnlet)
      use const_maxsize
      implicit none
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_mgo
   
   subroutine simu_energy_nlocal_go(irep, now_con, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_con(:,:)
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_nlocal_go
   
   subroutine simu_energy_nlocal_morse(irep, now_morse, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_morse(:,:)
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_nlocal_morse
   
   subroutine simu_energy_nlocal_rna_bp(irep, now_rna_bp, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_rna_bp(:,:)
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_nlocal_rna_bp
   
   subroutine simu_energy_nlocal_mgo(irep, now_con, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      integer,    intent(out) :: now_con(:,:)
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_nlocal_mgo
   
   subroutine simu_energy_orderpara(irep, now_allcon)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      integer,    intent(in)  :: now_allcon(:,:)
   endsubroutine simu_energy_orderpara
   
   subroutine simu_energy_pnl (irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: pnlet(:)
      real(PREC), intent(out) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_pnl
   
   subroutine simu_energy_pnl_restype (irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: pnlet(:)
      real(PREC), intent(out) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_pnl_restype
   
   subroutine simu_energy_pnl2(irep, pnlet, pnle_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: pnlet(:)
      real(PREC), intent(out)   :: pnle_unit(:,:,:)
   endsubroutine simu_energy_pnl2
   
   subroutine simu_energy_pnl3 (irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
      real(PREC), intent(inout) :: pnlet(:)
   endsubroutine simu_energy_pnl3

   subroutine simu_energy_ion (irep, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
      real(PREC), intent(inout) :: pnlet(:)
   endsubroutine simu_energy_ion

   subroutine simu_energy_ele(irep, pnlet, pnle_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: pnlet(:)
      real(PREC), intent(out)   :: pnle_unit(:,:,:)
   endsubroutine simu_energy_ele
   
   subroutine simu_energy_ele2(irep, pnlet, pnle_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: pnlet(:)
      real(PREC), intent(out)   :: pnle_unit(:,:,:)
   endsubroutine simu_energy_ele2
   
   subroutine simu_energy_ele3(irep, pnlet, pnle_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: pnlet(:)
      real(PREC), intent(out)   :: pnle_unit(:,:,:)
   endsubroutine simu_energy_ele3

   subroutine simu_energy_ele_coulomb(irep, pnlet, pnle_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: pnlet(:)
      real(PREC), intent(out)   :: pnle_unit(:,:,:)
   endsubroutine simu_energy_ele_coulomb
   
   subroutine simu_energy_hp (irep, pnlet, pnle_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
      real(PREC), intent(inout) :: pnlet(:)
   endsubroutine simu_energy_hp

!sasa
   subroutine simu_energy_sasa (irep, pnlet)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
!      real(PREC), intent(inout) :: pnle_unit(:,:,:)
      real(PREC), intent(inout) :: pnlet(:)
   endsubroutine simu_energy_sasa
   
   subroutine simu_energy_pulling(irep, pnle_unit, pnlet)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: pnlet(:)
     real(PREC), intent(inout) :: pnle_unit(:,:,:)
   end subroutine simu_energy_pulling
   
   subroutine simu_energy_velo(velo_mp, pnle_unit, pnlet)
      use const_maxsize
      implicit none
      real(PREC), intent(in)    :: velo_mp(:,:)
      real(PREC), intent(inout) :: pnlet(:)
      real(PREC), intent(inout) :: pnle_unit(:,:,:)
   endsubroutine simu_energy_velo
endinterface
endmodule if_energy
