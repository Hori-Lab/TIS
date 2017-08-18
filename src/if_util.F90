! if_write
!> @brief Interface module for utility subroutines

module if_util
interface

   subroutine util_pbneighbor(vx, imirror)
      use const_maxsize
      real(PREC), intent(inout) :: vx(3)
      integer, intent(out), optional :: imirror
   endsubroutine util_pbneighbor

endinterface
endmodule if_util
