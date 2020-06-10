subroutine write_T()
 
   use const_maxsize
   use const_physical
   use const_index
   use var_io,      only : outfile
   use var_struct,  only : nmp_real, xyz_mp_rep
   use var_replica, only : n_replica_mpi, irep2grep
   use var_setp,    only : inmisc, fix_mp
   use var_simu,    only : tempk, energy, force_mp

   implicit none

   integer :: irep, grep, nmp_fix, imp
   real(PREC) :: t, tc, tc_com
   real(PREC) :: com(3)
   logical, save :: flg_init = .True.

   if (flg_init) then
      do irep = 1, n_replica_mpi
         grep = irep2grep(irep)
         write(outfile%T(grep), '(a1,a6,1x,a8,1x,a12,1x,a12)') '#', 'Tset', 'Tkin', 'Tconf', 'Tconf*'
      enddo
      flg_init = .False.
   endif

   nmp_fix = 0
   if (inmisc%i_fix == 1) then
      nmp_fix = count(fix_mp)
   endif

   do irep = 1, n_replica_mpi
      
      ! Kinetic temperature
      t = 2.0 / (3.0 * (nmp_real - nmp_fix) * BOLTZ_KCAL_MOL) * energy(E_TYPE%VELO, irep)

      ! Configurational temperature (simpler form, not applicable for periodic boundary)
      com(:) = 0.0
      do imp = 1, nmp_real
         com(:) = com(:) + xyz_mp_rep(:,imp,irep)
      enddo
      com(:) = com(:) / real(nmp_real,kind=PREC)

      !! It seems that it is converged quicker when the translation (center of mass movement)
      !! is removed (tc_com).
      tc = 0.0
      tc_com = 0.0
      do imp = 1, nmp_real
         if(fix_mp(imp)) cycle
         tc = tc + dot_product(xyz_mp_rep(:,imp,irep), -force_mp(:,imp,irep))
         tc_com = tc_com + dot_product((xyz_mp_rep(:,imp,irep) - com(:)), -force_mp(:,imp,irep))
      enddo

      tc = tc / (3.0 * (nmp_real - nmp_fix) * BOLTZ_KCAL_MOL)
      tc_com = tc_com / (3.0 * (nmp_real - nmp_fix - 1) * BOLTZ_KCAL_MOL)

      grep = irep2grep(irep)
      write(outfile%T(grep), '(f7.3,1x,f8.3,1x,f12.3,1x,f12.3)') tempk, t, tc, tc_com

   enddo

endsubroutine write_T
