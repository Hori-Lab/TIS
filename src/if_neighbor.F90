! if_neighbor
!> @brief This is the interface for the neighborlist.

module if_neighbor
interface

   subroutine neighbor_pre(xyz_mp, ineigh_unit)
      use const_maxsize
      implicit none
      real(PREC), intent(in)  :: xyz_mp(:,:)
      integer,    intent(out) :: ineigh_unit(MXUNIT, MXUNIT)
   endsubroutine neighbor_pre

   subroutine neighbor_sort(irep, nexv, iexv2mp_in, iexv2mp_out, nexv_lall)
     use const_maxsize
     use mpiconst
     use var_struct, only : nmp_all
     integer, intent(in)  :: irep
     integer, intent(in)  :: nexv
     integer, intent(in)  :: iexv2mp_in (3,MXMPNEIGHBOR*nmp_all)
     integer, intent(in) ,optional :: nexv_lall(0:npar_mpi-1)
     integer, intent(out),optional :: iexv2mp_out(3,MXMPNEIGHBOR*nmp_all)
   endsubroutine neighbor_sort

endinterface
endmodule if_neighbor
