subroutine simu_initial()

  use if_write
  use const_maxsize
  use const_physical
  use const_index
  use var_io, only : i_initial_state, i_initial_velo, i_seq_read_style, flg_rst
  use var_setp, only : insimu, inmisc
  use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep, xyz_ref_mp, nmp_all
  use var_replica, only : n_replica_mpi, irep2grep
  use var_simu, only : istep_sim, tempk
  use mt_stream
  use mpiconst

  implicit none

  integer :: irep
#ifdef _DEBUG
  integer :: imp, grep
#endif

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
#ifdef _DEBUG
        write(*,*) '#####simu_initial: xyz_mp_rep after read_rst'
        do irep = 1, n_replica_mpi
           grep = irep2grep(irep)
           do imp = 1, nmp_all
              write(6,'(2i5,3g14.7)') grep, imp, xyz_mp_rep(1:3,imp,irep)
           enddo
        enddo
#endif
        
     end if

     if(insimu%i_com_zeroing_ini == 1 .or. insimu%i_com_zeroing_ini == 2) then
        call simu_xyz_adjst()
     end if

!     call simu_copyxyz(0)

     if(i_seq_read_style /= SEQREAD%PDB) then
        xyz_ref_mp(1:SDIM, 1:nmp_all) = xyz_mp_rep(1:SDIM, 1:nmp_all, 1)
     end if
  end if
  

  if (inmisc%class_flag(CLASS%ION)) then
     if (i_initial_state == INISTAT%RST) then
        ! Do not generate if restarted
        continue

     else
        do irep = 1, n_replica_mpi
           call simu_initial_ion(irep)
        end do
     endif
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

     else if(i_initial_velo == INIVELO%RST) then
        call read_rst(RSTBLK%VELO)

     endif
  end if

  ! -----------------------------------------------------------------
  ! stacking and H-bonding in DTRNA15
  ! -----------------------------------------------------------------
  if (flg_rst) then
     if (inmisc%class_flag(CLASS%RNA) .and. &
         (inmisc%i_dtrna_model==2015 .or. inmisc%i_dtrna_model==2019)) then
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
