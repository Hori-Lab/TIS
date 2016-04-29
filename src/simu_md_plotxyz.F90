! simu_md_plotxyz
!> @brief Make a random initial configuration

! ***********************************************************************
! input into xyz_mp the random structure or a structure
! ***********************************************************************
subroutine simu_md_plotxyz()

  use const_maxsize
  use const_physical
  use const_index
  use var_setp, only : mts
  use var_struct,  only : nmp_real, xyz_mp_rep, imp2type, iclass_mp
  use var_replica, only : n_replica_mpi
  use mt_stream
  implicit none

  ! -------------------------------------------------------------------
  integer    :: imp, irep, istream
  real(PREC) :: radxy, radyz
  real(PREC) :: pi2
  real(PREC) :: center_P(SDIM)
  real(PREC) :: grnd_array(2, nmp_real, n_replica_mpi)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  pi2 = F_PI * 70.0e0_PREC / 180.0e0_PREC

  ! -------------------------
  ! prepare random numbers 
  ! -------------------------
  do irep = 1, n_replica_mpi
     istream = irep
     do imp = 1, nmp_real-1
        grnd_array(1, imp, irep) = genrand_double1(mts(istream, 0))
        grnd_array(2, imp, irep) = genrand_double1(mts(istream, 0))
     enddo
  enddo

  if (iclass_mp(1) == CLASS%RNA) then
     select case (imp2type(1))
     case (MPTYPE%RNA_PHOS)
        xyz_mp_rep(1:SDIM, 1, 1:n_replica_mpi) = 0.0e0_PREC

        do irep = 1, n_replica_mpi
           do imp = 2, nmp_real, 3
              ! S
              radxy = (pi2) * grnd_array(1, imp-1, irep)
              radyz = (pi2) * grnd_array(2, imp-1, irep)
              xyz_mp_rep(1, imp, irep) = 5.0e0_PREC * sin(radyz) * cos(radxy) + xyz_mp_rep(1, imp - 1, irep)
              xyz_mp_rep(2, imp, irep) = 5.0e0_PREC * sin(radyz) * sin(radxy) + xyz_mp_rep(2, imp - 1, irep)
              xyz_mp_rep(3, imp, irep) = 5.0e0_PREC * cos(radyz)              + xyz_mp_rep(3, imp - 1, irep)

              ! P
              if (imp+2 > nmp_real) cycle
              radxy = (pi2) * grnd_array(1, imp, irep)
              radyz = (pi2) * grnd_array(2, imp, irep)
              xyz_mp_rep(1, imp+2, irep) = 5.0e0_PREC * sin(radyz) * cos(radxy) + xyz_mp_rep(1, imp, irep)
              xyz_mp_rep(2, imp+2, irep) = 5.0e0_PREC * sin(radyz) * sin(radxy) + xyz_mp_rep(2, imp, irep)
              xyz_mp_rep(3, imp+2, irep) = 5.0e0_PREC * cos(radyz)              + xyz_mp_rep(3, imp, irep)
           end do

           do imp = 3, nmp_real, 3
              if (imp+1 > nmp_real) then
                 ! B (last)
                 radxy = (pi2) * grnd_array(1, imp-1, irep)
                 radyz = (pi2) * grnd_array(2, imp-1, irep)
                 xyz_mp_rep(1, imp, irep) = 5.0e0_PREC * sin(radyz) * cos(radxy) + xyz_mp_rep(1, imp-1, irep)
                 xyz_mp_rep(2, imp, irep) = 5.0e0_PREC * sin(radyz) * sin(radxy) + xyz_mp_rep(2, imp-1, irep)
                 xyz_mp_rep(3, imp, irep) = 5.0e0_PREC * cos(radyz)              + xyz_mp_rep(3, imp-1, irep)
                 cycle
              endif
              ! B
              center_P(1:3) = (xyz_mp_rep(1:3, imp-2, irep) + xyz_mp_rep(1:3, imp+1, irep)) / 2.0e0_PREC
              xyz_mp_rep(1:3, imp, irep) = -3.0e0_PREC *  center_P(1:3) + 4.0e0_PREC * xyz_mp_rep(1:3, imp-1, irep)
           enddo
        enddo

     case (MPTYPE%RNA_SUGAR)
        ! B (imp=2)
        xyz_mp_rep(1:SDIM, 2, 1:n_replica_mpi) = 0.0e0_PREC

        do irep = 1, n_replica_mpi
           ! S (imp=1)
           radxy = (pi2) * grnd_array(1, 1, irep)
           radyz = (pi2) * grnd_array(2, 1, irep)
           xyz_mp_rep(1, 1, irep) = 5.0e0_PREC * sin(radyz) * cos(radxy) + xyz_mp_rep(1, 2, irep)
           xyz_mp_rep(2, 1, irep) = 5.0e0_PREC * sin(radyz) * sin(radxy) + xyz_mp_rep(2, 2, irep)
           xyz_mp_rep(3, 1, irep) = 5.0e0_PREC * cos(radyz)              + xyz_mp_rep(3, 2, irep)
           do imp = 3, nmp_real, 3
              ! P
              radxy = (pi2) * grnd_array(1, imp-1, irep)
              radyz = (pi2) * grnd_array(2, imp-1, irep)
              xyz_mp_rep(1, imp, irep) = 5.0e0_PREC * sin(radyz) * cos(radxy) + xyz_mp_rep(1, imp - 2, irep)
              xyz_mp_rep(2, imp, irep) = 5.0e0_PREC * sin(radyz) * sin(radxy) + xyz_mp_rep(2, imp - 2, irep)
              xyz_mp_rep(3, imp, irep) = 5.0e0_PREC * cos(radyz)              + xyz_mp_rep(3, imp - 2, irep)

              ! S
              if (imp+1 > nmp_real) cycle
              radxy = (pi2) * grnd_array(1, imp, irep)
              radyz = (pi2) * grnd_array(2, imp, irep)
              xyz_mp_rep(1, imp+1, irep) = 5.0e0_PREC * sin(radyz) * cos(radxy) + xyz_mp_rep(1, imp, irep)
              xyz_mp_rep(2, imp+1, irep) = 5.0e0_PREC * sin(radyz) * sin(radxy) + xyz_mp_rep(2, imp, irep)
              xyz_mp_rep(3, imp+1, irep) = 5.0e0_PREC * cos(radyz)              + xyz_mp_rep(3, imp, irep)
           end do

           do imp = 5, nmp_real, 3
              ! B
              if (imp+1 > nmp_real) then
                 ! B (last)
                 radxy = (pi2) * grnd_array(1, imp-1, irep)
                 radyz = (pi2) * grnd_array(2, imp-1, irep)
                 xyz_mp_rep(1, imp, irep) = 5.0e0_PREC * sin(radyz) * cos(radxy) + xyz_mp_rep(1, imp-1, irep)
                 xyz_mp_rep(2, imp, irep) = 5.0e0_PREC * sin(radyz) * sin(radxy) + xyz_mp_rep(2, imp-1, irep)
                 xyz_mp_rep(3, imp, irep) = 5.0e0_PREC * cos(radyz)              + xyz_mp_rep(3, imp-1, irep)
                 cycle
              endif
              center_P(1:3) = (xyz_mp_rep(1:3, imp-2, irep) + xyz_mp_rep(1:3, imp+1, irep)) / 2.0e0_PREC
              xyz_mp_rep(1:3, imp, irep) = -3.0e0_PREC *  center_P(1:3) + 4.0e0_PREC * xyz_mp_rep(1:3, imp-1, irep)
           enddo
        enddo

     case default
        error_message = 'Error: first mp should be P(phosphate) or S(sugar)'
        call util_error(ERROR%STOP_ALL, error_message)
     endselect

  else
     !xyz_mp_rep(1:SDIM, 1, 1:n_replica_all) = 0.0e0_PREC
     xyz_mp_rep(1:SDIM, 1, 1:n_replica_mpi) = 0.0e0_PREC

     ! -------------------------------------------------------------------
     ! plot the random structure
     !do irep = 1, n_replica_all
     do irep = 1, n_replica_mpi
        do imp = 2, nmp_real
           radxy = (pi2) * grnd_array(1, imp-1, irep)
           radyz = (pi2) * grnd_array(2, imp-1, irep)
           xyz_mp_rep(1, imp, irep) = 3.8e0_PREC * sin(radyz) * cos(radxy) + xyz_mp_rep(1, imp - 1, irep)
           xyz_mp_rep(2, imp, irep) = 3.8e0_PREC * sin(radyz) * sin(radxy) + xyz_mp_rep(2, imp - 1, irep)
           xyz_mp_rep(3, imp, irep) = 3.8e0_PREC * cos(radyz) + xyz_mp_rep(3, imp - 1, irep)
        end do
     enddo
  endif

end subroutine simu_md_plotxyz
