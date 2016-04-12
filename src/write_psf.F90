!write_psf
!> @breif Writes '.psf' file which is used in VMD to specify  
!>        molecular topology.

subroutine write_psf()

#include "format.F90"

   use const_index
   use var_struct, only : nmp_real, imp2unit, iclass_unit,      &
                          cmp2seq, ires_mp, cmp2atom, imp2type, &
                          cmass_mp, nbd, ibd2mp, lmp2charge, coef_charge
   use var_inp,    only : outfile
   use var_setp,   only : inmisc
#ifdef MPI_PAR
   use mpiconst
#endif
   
   implicit none
   integer :: imp, ibd, icounter
   integer :: lunout

#ifdef MPI_PAR
   if (myrank == 0) then
#endif
   lunout = outfile%psf

!!!!! atom
!         1         2         3         4         5         6         7
!1234567890123456789012345678901234567890123456789012345678901234567890
!       6 !NATOM
!       1 S1   1    STEA S    CTL2  -0.180000       12.0110           0
!       2 S1   1    STEA B    HAL2   0.090000        1.0080           0
!       3 S1   1    STEA P    HAL2   0.090000        1.0080           0
!       4 S1   1    STEA S    CTL2  -0.180000       12.0110           0
!       5 S1   1    STEA B    HAL2   0.090000        1.0080           0
!       6 S1   1    STEA P    HAL2   0.090000        1.0080           0

! See also http://www.ks.uiuc.edu/Research/namd/wiki/index.cgi?charmm2namd
! '   %5d %-4s %-4s %-4s %-4s %4s  %9.6f      %8.4f           %d'

   ! header
   write(lunout,'(i8)', ADVANCE='NO') nmp_real
   write(lunout,'(a)') ' !NATOM'
   ! data
   do imp = 1, nmp_real
      write(lunout,'(i8)', ADVANCE='NO') imp
      write(lunout,'(2xa3)', ADVANCE='NO') cmp2seq(imp)
      write(lunout,'(1xi4)', ADVANCE='NO') mod(ires_mp(imp), 10000)
      write(lunout,'(2xa3)', ADVANCE='NO') cmp2seq(imp)
      write(lunout,'(1xa4)', ADVANCE='NO') cmp2atom(imp)
      if (iclass_unit(imp2unit(imp)) == CLASS%RNA) then
         if (imp2type(imp) == MPTYPE%RNA_PHOS) then
            write(lunout,'(1xa4)', ADVANCE='NO') _STR_PSF_ATOM_RNA_P_
         else if (imp2type(imp) == MPTYPE%RNA_BASE) then
            write(lunout,'(1xa4)', ADVANCE='NO') _STR_PSF_ATOM_RNA_B_
         else if (imp2type(imp) == MPTYPE%RNA_SUGAR) then
            write(lunout,'(1xa4)', ADVANCE='NO') _STR_PSF_ATOM_RNA_S_
         endif
      else 
         write(lunout,'(1xa4)', ADVANCE='NO') _STR_PSF_ATOM_PRO_
      endif
      if (lmp2charge(imp) == 0) then
         write(lunout, '(2xf9.6)', ADVANCE='NO') 0.0
      else
         write(lunout, '(2xf9.6)', ADVANCE='NO') coef_charge( lmp2charge(imp) , 1)
      endif
      write(lunout, '(6xf8.4)', ADVANCE='NO') cmass_mp(imp)
      write(lunout, '(11xa1)') '0'
   enddo

!!!!! bond
!         1         2         3         4         5         6         7
!1234567890123456789012345678901234567890123456789012345678901234567890
!       5 !NBOND: bonds
!       1       2       1       3       3       4       4       5
!       4       6 
   ! header
   write(lunout, '(i8)', ADVANCE='NO') nbd
   write(lunout, '(a)') ' !NBOND: bonds'
   ! data
   icounter = 0
   do ibd = 1, nbd
      icounter = icounter + 1
      if (icounter < 4) then
         write(lunout, '(i8,i8)', ADVANCE='NO') ibd2mp(1,ibd), ibd2mp(2,ibd)
      else
         write(lunout, '(i8,i8)') ibd2mp(1,ibd), ibd2mp(2,ibd)
         icounter = 0
      endif
   enddo
   if (icounter /= 0) then
      write(lunout,*) ''
   endif

#ifdef MPI_PAR
  endif
#endif

end subroutine write_psf
