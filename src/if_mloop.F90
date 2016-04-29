! if_mloop
!> @brief Interface module for subroutines called by mloop_simulator

module if_mloop
interface

   subroutine simu_boxmuller(r_boxmuller,nsize)
      use const_maxsize
      implicit none
      real(PREC), intent(inout) :: r_boxmuller(:)
      integer,    intent(in)    :: nsize
   endsubroutine
 
   subroutine energy_allrep(energy_unit, energy, &
                                 velo_mp, replica_energy, flg_replica, tempk)
      use const_maxsize
      implicit none
      real(PREC), intent(out) :: energy_unit(:,:,:,:)  ! (MXUNIT, MXUNIT, E_TYPE%MAX, replica)
      real(PREC), intent(out) :: energy(:,:)          ! (E_TYPE%MAX,replica)
      real(PREC), intent(in)  :: velo_mp(:,:,:)      ! (3, MXMP, replica)
      real(PREC), intent(out) :: replica_energy(:,:) ! (2, replica)
      real(PREC), intent(in)  :: tempk
      logical, intent(in)  :: flg_replica
   endsubroutine energy_allrep

   subroutine energy_sumup(irep, velo_mp, energy, energy_unit)
      use const_maxsize
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(in)  :: velo_mp(:,:)      ! (3, nmp_real)
      real(PREC), intent(out) :: energy(:)          ! (E_TYPE%MAX)
      real(PREC), intent(out) :: energy_unit(:,:,:)  ! (nunit_all, nunit_all, E_TYPE%MAX)
   endsubroutine energy_sumup

   subroutine write_traject_file(ibefore_time, istep, tempk, velo_mp)
      use const_maxsize
      implicit none
      integer(L_INT), intent(in) :: ibefore_time
      integer(L_INT), intent(in) :: istep
      real(PREC), intent(in) :: tempk
      real(PREC), intent(in) :: velo_mp(:,:,:)
   endsubroutine write_traject_file

   subroutine write_record_file(istep, velo_mp)
      use const_maxsize
      implicit none
      integer(L_INT), intent(in) :: istep
      real(PREC), intent(in) :: velo_mp(:,:,:)
   endsubroutine write_record_file

   subroutine simu_radiusg(rg_unit, rg)
      use const_maxsize
      implicit none
      real(PREC), intent(out) :: rg(:)               ! (replica)
      real(PREC), intent(out) :: rg_unit(:,:)        ! (unit, replica)
   endsubroutine simu_radiusg

   subroutine simu_rmsd(rmsd_unit, rmsd)
      use const_maxsize
      implicit none
      real(PREC), intent(out) :: rmsd(:)             ! (replica)
      real(PREC), intent(out) :: rmsd_unit(:,:)      ! (unit, replica)
   endsubroutine simu_rmsd

   subroutine simu_replica_exchange(velo_mp, replica_energy, tempk)
      use const_maxsize
      implicit none
      real(PREC), intent(inout)  :: velo_mp(:,:,:)      ! (SDIM, MXMP, replica)
      real(PREC), intent(in)     :: replica_energy(:,:) ! (2, replica)
      real(PREC), intent(in)     :: tempk
   endsubroutine simu_replica_exchange

   subroutine simu_velo_adjst(velo_mp, irep)
      use const_maxsize
      implicit none
      real(PREC), intent(inout) :: velo_mp(:,:,:)
      integer,    intent(in)    :: irep
   endsubroutine simu_velo_adjst

!   subroutine simu_velo_mrand(velo_mp, tempk_in)
!      use const_maxsize
!      implicit none
!      real(PREC), intent(out) :: velo_mp(:,:,:)
!      real(PREC), intent(in)  :: tempk_in
!   endsubroutine simu_velo_mrand

   subroutine simu_velo_settemp(velo_mp, irep, tempk)
      use const_maxsize
      implicit none
      real(PREC), intent(in)    :: tempk
      real(PREC), intent(inout) :: velo_mp(:,:,:)
      integer,    intent(in)    :: irep
   endsubroutine simu_velo_settemp

   subroutine simu_velo_nosehoover(velo_mp, irep, tempk, velo_yojou)
      use const_maxsize
      implicit none
      real(PREC), intent(in) :: tempk
      integer,    intent(in) :: irep
      real(PREC), intent(out) :: velo_yojou
      real(PREC), intent(inout) :: velo_mp(:,:,:)
   endsubroutine simu_velo_nosehoover

   subroutine simu_mc_implig(irep, istep, tempk)
      use const_maxsize
      use var_implig, only : inimplig
      use var_replica, only : n_replica_all
      implicit none
      integer,    intent(in) :: irep
      integer(L_INT), intent(in) :: istep
      real(PREC), intent(in) :: tempk
   endsubroutine simu_mc_implig

   subroutine simu_replica_opt_temp(i_current_stage)
      implicit none
      integer, intent(out), optional :: i_current_stage
   endsubroutine simu_replica_opt_temp

endinterface 
endmodule if_mloop
