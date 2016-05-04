subroutine simu_initial()

  use if_write
  use const_maxsize
  use const_physical
  use const_index
  use var_inp, only : i_initial_state, i_initial_velo, &
                      i_seq_read_style, i_simulate_type, flg_rst
  use var_setp, only : insimu, inmisc
  use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep, xyz_ref_mp
  use var_replica, only : n_replica_mpi
  use var_simu, only : istep_sim, velo_mp, tempk
  use mt_stream
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  integer :: irep

  ! -----------------------------------------------------------------
#ifdef _DEBUG
  write(*,*) '#### start simu_initial'
#endif
  
  pxyz_mp_rep(1:SDIM,1:nmp_all,1:n_replica_mpi) = 0.0

  ! -----------------------------------------------------------------
  ! initial structure
  ! -----------------------------------------------------------------
  if(istep_sim == 1) then
     ! Random configuration
     if(i_initial_state == INISTAT%RANDOM) then
        call simu_md_plotxyz()
        
     ! Native configuration
     else if(i_initial_state == INISTAT%NATIVE) then
        call simu_copyxyz_replica()

     ! Configuration given in the input
     else if(i_initial_state == INISTAT%INPUT .or. &
             i_initial_state == INISTAT%CG    .or. &
             i_initial_state == INISTAT%CARD  ) then
        call simu_initial_xyz()
        call simu_copyxyz_replica()
        
     else if(i_initial_state == INISTAT%RST) then
        call read_rst(RSTBLK%XYZ)
        
     end if

     if(insimu%i_com_zeroing_ini == 1 .or. insimu%i_com_zeroing_ini == 2) then
        call simu_xyz_adjst()
     end if

     call simu_copyxyz(0)

     if(i_seq_read_style /= SEQREAD%PDB) then
        xyz_ref_mp(1:SDIM, 1:nmp_all) = xyz_mp_rep(1:SDIM, 1:nmp_all, 1)
     end if
  end if
  
  if (inmisc%class_flag(CLASS%ION)) then
     do irep = 1, n_replica_mpi
        call simu_initial_ion(irep)
     end do
  end if

  if (inmisc%i_reset_struct == 1 .and. .not. flg_rst) then
     call simu_copyxyz_ref()
  endif
  
  pxyz_mp_rep(1:SDIM,1:nmp_all,1:n_replica_mpi) = xyz_mp_rep(1:SDIM,1:nmp_all,1:n_replica_mpi)

  ! -----------------------------------------------------------------
  ! initial velocity
  ! -----------------------------------------------------------------
  if(istep_sim == 1 .OR. inmisc%i_reset_struct == 1) then
     ! Random velocities 
     if(i_initial_velo == INIVELO%MAXWELL) then
        call simu_velo_mrand(tempk)

     ! Given in .velo file
     else if(i_initial_velo == INIVELO%CARD) then
        call read_crd(RECORD_FILE%VELO, velo_mp)

     else if(i_initial_velo == INIVELO%RST) then
        call read_rst(RSTBLK%VELO)

     endif
  end if

  ! -----------------------------------------------------------------
  ! stacking and H-bonding in DTRNA15
  ! -----------------------------------------------------------------
  if (flg_rst) then
     if (inmisc%class_flag(CLASS%RNA) .and. inmisc%i_dtrna_model==2015) then
        call read_rst(RSTBLK%DTRNA15)
     endif
  endif

#ifdef _DEBUG
!  do irep=1,n_replica_mpi
!     do imp=1, nmp_real
!        write(6,'(i5,1pd12.5)'),imp,xyz_mp_rep(1,imp,irep)
!     enddo
!  enddo
#endif

#ifdef _DEBUG
  write(*,*) '#### end simu_initial'
#endif

end subroutine simu_initial
