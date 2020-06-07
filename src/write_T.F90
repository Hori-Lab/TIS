subroutine write_T(tempk, energy)
 
   use const_maxsize
   use const_physical
   use const_index
   use var_io,      only : outfile
   use var_struct,  only : nmp_real
   use var_replica, only : n_replica_mpi, irep2grep

   implicit none

   real(PREC), intent(in) :: tempk
   real(PREC), intent(in) :: energy(:,:)

   integer :: irep, grep
   real(PREC) :: t

   do irep = 1, n_replica_mpi

      grep = irep2grep(irep)

      t = 2.0 / (3.0 * nmp_real * BOLTZ_KCAL_MOL) * energy(E_TYPE%VELO, irep)

      write(outfile%T(grep), '(f7.3,1x,f7.3)') tempk, t
   enddo

endsubroutine write_T
