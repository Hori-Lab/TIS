! if_neighbor
!> @brief This is the interface for the neighborlist.

module if_neighbor
interface

   subroutine simu_neighbor_pre(xyz_mp, ineigh_unit)
      use const_maxsize
      implicit none
      real(PREC), intent(in)  :: xyz_mp(:,:)
      integer,    intent(out) :: ineigh_unit(MXUNIT, MXUNIT)
   endsubroutine simu_neighbor_pre

   subroutine simu_neighbor_sort(irep, npnl, ipnl2mp_in, ipnl2mp_out, npnl_lall)
     use const_maxsize
     use mpiconst
     use var_struct, only : nmp_all
     integer, intent(in)  :: irep
     integer, intent(in)  :: npnl
     integer, intent(in)  :: ipnl2mp_in (3,MXMPNEIGHBOR*nmp_all)
     integer, intent(in) ,optional :: npnl_lall(0:npar_mpi-1)
     integer, intent(out),optional :: ipnl2mp_out(3,MXMPNEIGHBOR*nmp_all)
   endsubroutine simu_neighbor_sort

endinterface
endmodule if_neighbor
