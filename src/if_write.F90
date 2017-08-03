! if_write
!> @brief Interface module for output subroutine

module if_write
interface

   !subroutine write_tseries(ibefore_time, ntstep, istep, &
   subroutine write_tseries(ibefore_time, istep, &
                            rg_unit, rg, rmsd_unit, rmsd, &
                            energy_unit, energy, temp_in, &
                            flg_header)
      use const_maxsize
      integer(L_INT), intent(in) :: ibefore_time
   !   integer,    intent(in) :: ntstep
      integer(L_INT), intent(in) :: istep
      real(PREC), intent(in) :: rg_unit(:,:)        ! (unit, replica)
      real(PREC), intent(in) :: rg(:)               ! (replica)
      real(PREC), intent(in) :: rmsd_unit(:,:)      ! (unit, replica)
      real(PREC), intent(in) :: rmsd(:)             ! (replica)
      real(PREC), intent(in) :: energy_unit(:,:,:,:)  ! (unit, unit, E_TYPE%MAX, replica)
      real(PREC), intent(in) :: energy(:,:)          ! (E_TYPE%MAX, replica)
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

   subroutine write_stack(ene_st, ene_tst, tempk)
      use const_maxsize
      implicit none
      real(PREC), intent(in) :: ene_st(:,:)
      real(PREC), intent(in) :: ene_tst(:,:)
      real(PREC), intent(in) :: tempk
   endsubroutine write_stack

endinterface
endmodule if_write
