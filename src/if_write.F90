! if_write
!> @brief Interface module for output subroutine

module if_write
interface

   !subroutine write_tseries(ibefore_time, ntstep, istep, &
   subroutine write_tseries(ibefore_time, istep, &
                            rg_unit, rg, rmsd_unit, rmsd, &
                            pnle_unit, pnlet, temp_in, &
                            flg_header)
      use const_maxsize
      integer(L_INT), intent(in) :: ibefore_time
   !   integer,    intent(in) :: ntstep
      integer(L_INT), intent(in) :: istep
      real(PREC), intent(in) :: rg_unit(:,:)        ! (unit, replica)
      real(PREC), intent(in) :: rg(:)               ! (replica)
      real(PREC), intent(in) :: rmsd_unit(:,:)      ! (unit, replica)
      real(PREC), intent(in) :: rmsd(:)             ! (replica)
      real(PREC), intent(in) :: pnle_unit(:,:,:,:)  ! (unit, unit, E_TYPE%MAX, replica)
      real(PREC), intent(in) :: pnlet(:,:)          ! (E_TYPE%MAX, replica)
      real(PREC), intent(in) :: temp_in
      logical, intent(in), optional :: flg_header
   endsubroutine write_tseries

   subroutine write_xyz_dcd(i_coor_velo, ibefore_time, istep, ntstep, tempk, velo_mp)
      use const_maxsize
      implicit none
      integer,    intent(in) :: i_coor_velo
      integer(L_INT), intent(in) :: ibefore_time
      integer(L_INT), intent(in) :: istep
      integer(L_INT), intent(in) :: ntstep
      real(PREC), intent(in) :: tempk
      real(PREC), intent(in) :: velo_mp(:,:,:)
   endsubroutine write_xyz_dcd

   subroutine write_crd(i_record_file, istep, velo_mp)
      use const_maxsize
      implicit none
      integer, intent(in) :: i_record_file
      integer(L_INT), intent(in) :: istep
!      real(PREC), intent(in), optional :: velo_mp(:,:,:)
      real(PREC), intent(in) :: velo_mp(:,:,:)
   endsubroutine write_crd


   subroutine read_crd(i_record_file, velo_mp)
     use const_maxsize
     implicit none
      integer, intent(in) :: i_record_file
!     real(PREC), intent(out), optional :: velo_mp(:,:,:)
     real(PREC), intent(out) :: velo_mp(:,:,:)
   endsubroutine read_crd

endinterface
endmodule if_write
