!if_energy
!> @brief Interface module for subroutines which called by simu_energy.

module if_energy
interface

   subroutine energy_implig(irep, e_exv_unit, e_exv, iflag_for_mc)   
     use const_maxsize
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)
     real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX)
     integer,    intent(in)    :: iflag_for_mc 
   endsubroutine energy_implig

   subroutine energy_anchor(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX, replica)
   endsubroutine energy_anchor

   subroutine energy_rest1d(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX, replica)
   endsubroutine energy_rest1d
   
   subroutine energy_bangle(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine energy_bangle
   
   subroutine energy_bond  (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_bond
   
   subroutine energy_box(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv(:)
     real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_box
   
   subroutine energy_bridge(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX, replica)
   endsubroutine energy_bridge
   
   subroutine energy_dih   (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_dih
  
   subroutine energy_dih_harmonic (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_dih_harmonic

   subroutine energy_rna_stack(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_rna_stack

   subroutine energy_dtrna_stack(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_dtrna_stack

   subroutine energy_dtrna_hbond13(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_dtrna_hbond13

   subroutine energy_dtrna_hbond15(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_dtrna_hbond15

   subroutine energy_dtrna_stack_nlocal(irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_dtrna_stack_nlocal

   subroutine energy_exv_wca (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: e_exv(:)
      real(PREC), intent(out) :: e_exv_unit(:,:,:)
   endsubroutine energy_exv_wca

   subroutine energy_exv_dt15 (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: e_exv(:)
      real(PREC), intent(out) :: e_exv_unit(:,:,:)
   endsubroutine energy_exv_dt15
   
   subroutine energy_enm(irep, now_con, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_con(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_enm
   
   subroutine energy_mgo(e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_mgo
   
   subroutine energy_nlocal_go(irep, now_con, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_con(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_nlocal_go
   
   subroutine energy_nlocal_morse(irep, now_morse, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_morse(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_nlocal_morse
   
   subroutine energy_nlocal_rna_bp(irep, now_rna_bp, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_rna_bp(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_nlocal_rna_bp
   
   subroutine energy_nlocal_mgo(irep, now_con, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      integer,    intent(out) :: now_con(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_nlocal_mgo
   
   subroutine energy_orderpara(irep, now_allcon)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      integer,    intent(in)  :: now_allcon(:,:)
   endsubroutine energy_orderpara
   
   subroutine energy_exv (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: e_exv(:)
      real(PREC), intent(out) :: e_exv_unit(:,:,:)
   endsubroutine energy_exv
   
   subroutine energy_exv_restype (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: e_exv(:)
      real(PREC), intent(out) :: e_exv_unit(:,:,:)
   endsubroutine energy_exv_restype
   
   subroutine energy_exv2(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine energy_exv2
   
   subroutine energy_exv3 (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine energy_exv3

   subroutine energy_ion (irep, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine energy_ion

   subroutine energy_ele(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine energy_ele
   
   subroutine energy_ele2(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine energy_ele2
   
   subroutine energy_ele3(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine energy_ele3

   subroutine energy_ele_coulomb(irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: e_exv(:)
      real(PREC), intent(out)   :: e_exv_unit(:,:,:)
   endsubroutine energy_ele_coulomb
   
   subroutine energy_hp (irep, e_exv, e_exv_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine energy_hp

!sasa
   subroutine energy_sasa (irep, e_exv)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
!      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
      real(PREC), intent(inout) :: e_exv(:)
   endsubroutine energy_sasa
   
   subroutine energy_pulling(irep, e_exv_unit, e_exv)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: e_exv(:)
     real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   end subroutine energy_pulling
   
   subroutine energy_velo(velo_mp, e_exv_unit, e_exv)
      use const_maxsize
      implicit none
      real(PREC), intent(in)    :: velo_mp(:,:)
      real(PREC), intent(inout) :: e_exv(:)
      real(PREC), intent(inout) :: e_exv_unit(:,:,:)
   endsubroutine energy_velo
endinterface
endmodule if_energy
