!if_energy
!> @brief Interface module for subroutines which called by simu_energy.

module if_energy
interface

!   subroutine energy_implig(irep, energy_unit, energy, iflag_for_mc)   
!     use const_maxsize
!     integer,    intent(in)    :: irep
!     real(PREC), intent(inout) :: energy_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)
!     real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX)
!     integer,    intent(in)    :: iflag_for_mc 
!   endsubroutine energy_implig

   subroutine energy_anchor(irep, energy_unit, energy)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: energy_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX, replica)
   endsubroutine energy_anchor

   subroutine energy_rest1d(irep, energy_unit, energy)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: energy_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX, replica)
   endsubroutine energy_rest1d
   
   subroutine energy_bangle(irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy_unit(:,:,:)
      real(PREC), intent(inout) :: energy(:)
   endsubroutine energy_bangle
   
   subroutine energy_bond  (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_bond
   
   subroutine energy_fene  (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_fene

   subroutine energy_rouse  (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_rouse
   
!   subroutine energy_box(irep, energy_unit, energy)
!     use const_maxsize
!     implicit none
!     integer,    intent(in)    :: irep
!     real(PREC), intent(inout) :: energy(:)
!     real(PREC), intent(inout) :: energy_unit(:,:,:)
!   endsubroutine energy_box
   
   subroutine energy_bridge(irep, energy_unit, energy)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: energy_unit(:,:,:) ! (unit, unit, E_TYPE%MAX, replica)
     real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX, replica)
   endsubroutine energy_bridge
   
!   subroutine energy_dih   (irep, energy_unit, energy)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!      real(PREC), intent(inout) :: energy(:)
!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!   endsubroutine energy_dih
  
!   subroutine energy_dih_harmonic (irep, energy_unit, energy)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!      real(PREC), intent(inout) :: energy(:)
!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!   endsubroutine energy_dih_harmonic

!   subroutine energy_rna_stack(irep, energy_unit, energy)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!      real(PREC), intent(inout) :: energy(:)
!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!   endsubroutine energy_rna_stack

   subroutine energy_dtrna_stack(irep, energy_unit, energy, ene_st)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
      real(PREC), intent(inout) :: ene_st(:)
   endsubroutine energy_dtrna_stack

   subroutine energy_dtrna_hbond13(irep, energy_unit, energy, ene_hb)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
      real(PREC), intent(out)   ::  ene_hb(:)
   endsubroutine energy_dtrna_hbond13

   subroutine energy_dtrna_hbond15(irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_dtrna_hbond15

   subroutine energy_dtrna_stack_nlocal(irep, energy_unit, energy, ene_tst, st_status)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
      real(PREC), intent(out)   :: ene_tst(:)
      logical,    intent(out)   :: st_status(:)
   endsubroutine energy_dtrna_stack_nlocal

   subroutine energy_exv_wca (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: energy(:)
      real(PREC), intent(out) :: energy_unit(:,:,:)
   endsubroutine energy_exv_wca

   subroutine energy_exv_dt15 (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: energy(:)
      real(PREC), intent(out) :: energy_unit(:,:,:)
   endsubroutine energy_exv_dt15
   
   subroutine energy_enm(irep, now_con, energy, energy_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_con(:,:)
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_enm
   
!   subroutine energy_mgo(energy_unit, energy)
!      use const_maxsize
!      implicit none
!      real(PREC), intent(inout) :: energy(:)
!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!   endsubroutine energy_mgo
   
!   subroutine energy_nlocal_go(irep, now_con, energy_unit, energy)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!      integer,    intent(out)   :: now_con(:,:)
!      real(PREC), intent(inout) :: energy(:)
!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!   endsubroutine energy_nlocal_go
   
   subroutine energy_LJ(irep, now_LJ, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_LJ(:,:)
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_LJ
   
   subroutine energy_LJ_1210(irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_LJ_1210

   subroutine energy_HPS(irep, now_HPS, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      integer,    intent(out)   :: now_HPS(:,:)
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_HPS

   !subroutine energy_wca(irep, now_wca, energy_unit, energy)
   subroutine energy_wca(irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      !integer,    intent(out)   :: now_wca(:,:)
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_wca
   
   subroutine energy_con_gauss(irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_con_gauss
   
!   subroutine energy_nlocal_morse(irep, now_morse, energy_unit, energy)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!      integer,    intent(out)   :: now_morse(:,:)
!      real(PREC), intent(inout) :: energy(:)
!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!   endsubroutine energy_nlocal_morse
   
!   subroutine energy_nlocal_rna_bp(irep, now_rna_bp, energy_unit, energy)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!      integer,    intent(out)   :: now_rna_bp(:,:)
!      real(PREC), intent(inout) :: energy(:)
!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!   endsubroutine energy_nlocal_rna_bp
   
!   subroutine energy_nlocal_mgo(irep, now_con, energy_unit, energy)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)  :: irep
!      integer,    intent(out) :: now_con(:,:)
!      real(PREC), intent(inout) :: energy(:)
!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!   endsubroutine energy_nlocal_mgo
   
   subroutine energy_orderpara(irep, now_allcon)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      integer,    intent(in)  :: now_allcon(:,:)
   endsubroutine energy_orderpara
   
   subroutine energy_exv_rep12 (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: energy(:)
      real(PREC), intent(out) :: energy_unit(:,:,:)
   endsubroutine energy_exv_rep12
   
   subroutine energy_exv_rep6 (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: energy(:)
      real(PREC), intent(out) :: energy_unit(:,:,:)
   endsubroutine energy_exv_rep6
   
   subroutine energy_exv_restype (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: energy(:)
      real(PREC), intent(out) :: energy_unit(:,:,:)
   endsubroutine energy_exv_restype
   
   subroutine energy_exv_gauss (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: energy(:)
      real(PREC), intent(out) :: energy_unit(:,:,:)
   endsubroutine energy_exv_gauss

   subroutine energy_ele_DH(irep, energy, energy_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: energy(:)
      real(PREC), intent(out)   :: energy_unit(:,:,:)
   endsubroutine energy_ele_DH
   
!   subroutine energy_ele2(irep, energy, energy_unit)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!      real(PREC), intent(out)   :: energy(:)
!      real(PREC), intent(out)   :: energy_unit(:,:,:)
!   endsubroutine energy_ele2
!   
!   subroutine energy_ele3(irep, energy, energy_unit)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!      real(PREC), intent(out)   :: energy(:)
!      real(PREC), intent(out)   :: energy_unit(:,:,:)
!   endsubroutine energy_ele3

   subroutine energy_ele_coulomb(irep, energy, energy_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: energy(:)
      real(PREC), intent(out)   :: energy_unit(:,:,:)
   endsubroutine energy_ele_coulomb

   subroutine energy_ele_coulomb_ewld(irep, energy, energy_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: energy(:)
      real(PREC), intent(out)   :: energy_unit(:,:,:)
   endsubroutine energy_ele_coulomb_ewld

   subroutine energy_ele_coulomb_brute(irep, energy, energy_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)    :: irep
      real(PREC), intent(out)   :: energy(:)
      real(PREC), intent(out)   :: energy_unit(:,:,:)
   endsubroutine energy_ele_coulomb_brute
   
!   subroutine energy_hp (irep, energy, energy_unit)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!      real(PREC), intent(inout) :: energy(:)
!   endsubroutine energy_hp

!   subroutine energy_sasa (irep, energy)
!      use const_maxsize
!      implicit none
!      integer,    intent(in)    :: irep
!!      real(PREC), intent(inout) :: energy_unit(:,:,:)
!!      real(PREC), intent(inout) :: energy(:)
!   endsubroutine energy_sasa
   
   subroutine energy_pulling(irep, energy_unit, energy)
     use const_maxsize
     implicit none
     integer,    intent(in)    :: irep
     real(PREC), intent(inout) :: energy(:)
     real(PREC), intent(inout) :: energy_unit(:,:,:)
   end subroutine energy_pulling
   
   subroutine energy_kinetic(velo_mp, energy_unit, energy)
      use const_maxsize
      implicit none
      real(PREC), intent(in)    :: velo_mp(:,:)
      real(PREC), intent(inout) :: energy(:)
      real(PREC), intent(inout) :: energy_unit(:,:,:)
   endsubroutine energy_kinetic

   subroutine energy_BBR (irep, energy_unit, energy)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(out) :: energy(:)
      real(PREC), intent(out) :: energy_unit(:,:,:)
   endsubroutine energy_BBR
endinterface
endmodule if_energy
