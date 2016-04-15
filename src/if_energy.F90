!if_energy
!> @brief Interface module for subroutines which called by simu_energy.

module if_energy
interface

   subroutine simu_energy_implig(irep, e_exv_unit, e_exv, iflag_for_mc)   
     use const_maxsize
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)
     real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX)
     integer,    intent(in)    :: iflag_for_mc 
   endsubroutine simu_energy_implig

   subroutine simu_energy_anchor(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX, replica)
   endsubroutine simu_energy_anchor

   subroutine simu_energy_rest1d(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX, replica)
   endsubroutine simu_energy_rest1d
   
   subroutine simu_energy_bangle(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine simu_energy_bangle
   
   subroutine simu_energy_bond  (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_bond
   
   subroutine simu_energy_box(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv(:)
     real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_box
   
   subroutine simu_energy_bridge(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX, replica)
   endsubroutine simu_energy_bridge
   
   subroutine simu_energy_dih   (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_dih
  
   subroutine simu_energy_dih_harmonic (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_dih_harmonic

   subroutine simu_energy_rna_stack(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_rna_stack

   subroutine simu_energy_dtrna_stack(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_dtrna_stack

   subroutine simu_energy_dtrna_hbond13(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_dtrna_hbond13

   subroutine simu_energy_dtrna_hbond15(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_dtrna_hbond15

   subroutine simu_energy_dtrna_stack_nlocal(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_dtrna_stack_nlocal

   subroutine simu_energy_exv_wca (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: e_exv(:)
      real(PREC), intent(out) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_exv_wca

   subroutine simu_energy_exv_dt15 (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: e_exv(:)
      real(PREC), intent(out) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_exv_dt15
   
   subroutine simu_energy_enm(irep, now_con, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_con(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_enm
   
   subroutine simu_energy_mgo(e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_mgo
   
   subroutine simu_energy_nlocal_go(irep, now_con, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_con(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_nlocal_go
   
   subroutine simu_energy_nlocal_morse(irep, now_morse, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_morse(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_nlocal_morse
   
   subroutine simu_energy_nlocal_rna_bp(irep, now_rna_bp, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_rna_bp(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_nlocal_rna_bp
   
   subroutine simu_energy_nlocal_mgo(irep, now_con, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      integer,    intent(out) :: now_con(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_nlocal_mgo
   
   subroutine simu_energy_orderpara(irep, now_allcon)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      integer,    intent(in)  :: now_allcon(:,:)
   endsubroutine simu_energy_orderpara
   
   subroutine simu_energy_exv (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: e_exv(:)
      real(PREC), intent(out) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_exv
   
   subroutine simu_energy_exv_restype (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: e_exv(:)
      real(PREC), intent(out) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_exv_restype
   
   subroutine simu_energy_exv2(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_exv2
   
   subroutine simu_energy_exv3 (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine simu_energy_exv3

   subroutine simu_energy_ion (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine simu_energy_ion

   subroutine simu_energy_ele(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_ele
   
   subroutine simu_energy_ele2(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_ele2
   
   subroutine simu_energy_ele3(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_ele3

   subroutine simu_energy_ele_coulomb(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_ele_coulomb
   
   subroutine simu_energy_hp (irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine simu_energy_hp

!sasa
   subroutine simu_energy_sasa (irep, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
!      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine simu_energy_sasa
   
   subroutine simu_energy_pulling(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv(:)
     real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   end subroutine simu_energy_pulling
   
   subroutine simu_energy_velo(velo_mp, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      real(PREC), intent(in)    :: velo_mp(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine simu_energy_velo
endinterface
endmodule if_energy
