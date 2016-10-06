! inp_energy_func
!> @brief Subroutine for reading "energy_function"

! **********************************************************************
! This routine is for reading the itype_nlocal from input file.
! itype_nlocal : types of the force and the energy of Go
! **********************************************************************
subroutine inp_energy_func()

  use const_maxsize
  use const_index
  use const_physical
  use var_io,    only : infile, ius2unit, outfile !, i_aicg
  use var_setp,   only : inmisc !, inflp
  use var_struct, only : nunit_all
!  use var_mgo,    only : inmgo
  use mpiconst

  implicit none

  integer :: k
  integer :: i1 = 1
  integer :: i2 = 1
  integer :: isw, itype, icol, jcol, i_ninfo
  integer :: iu, ju, iunit, junit, i_local_type
!  integer :: isys, istat, iactnum, iact, isw2
  integer :: inunit(2), jnunit(2), instate, instate2
  integer :: inum
  integer :: luninp, lunout
  integer :: iline, nlines
  integer :: iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: char00
  character(12) :: char12
  character(CARRAY_MXCOLM)  :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! ----------------------------------------------------------------------
  luninp  = infile%inp
  lunout  = outfile%data
  
  ! ----------------------------------------------------------------------
  ! default setting
  inmisc%flag_local_unit(1:MXUNIT, 1:MXUNIT, 1:LINTERACT%MAX) = .FALSE.
  inmisc%flag_nlocal_unit(1:MXUNIT, 1:MXUNIT, 1:INTERACT%MAX) = .FALSE.
  inmisc%flag_nlocal_unit(1:MXUNIT, 1:MXUNIT, INTERACT%NOTHING) = .TRUE.

!  i_aicg = 1  ! aicg
!  inflp%i_flp = 0 ! flexible_local_potential

  inmisc%i_use_atom_protein = 0
  inmisc%i_residuenergy_radii = 0
  inmisc%i_output_energy_style = 0
  inmisc%i_triple_angle_term = 1   ! default
  inmisc%flg_coef_from_ninfo = .true.  ! default
  inmisc%i_temp_independent = 0   ! default
  inmisc%i_dtrna_model = 2013  ! default

  ! ----------------------------------------------------------------------
  ! read energy function
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'energy_function ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "energy_function" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     char00 = cwkinp(iline)

     ! LOCAL
     
     if(char00(1:5) == 'LOCAL') then
        i_local_type = 0
        do icol = 7, CARRAY_MXCOLM
           if(char00(icol:icol) == '/') exit
           if(icol == CARRAY_MXCOLM) then
              i_local_type = 1
           end if
        end do
        do jcol = 7, CARRAY_MXCOLM
           if(char00(jcol:jcol) == ')') exit
        end do

        if(i_local_type == 1) then
           read (char00(7:jcol-1), '(a)') char12
           call util_unitstate(char12, inunit, instate)
           read (char00(7:jcol-1), '(a)') char12
           call util_unitstate(char12, jnunit, instate2)
        else
           read (char00(7:icol-1), '(a)') char12
           call util_unitstate(char12, inunit, instate)
           read (char00(icol+1:jcol-1), '(a)') char12
           call util_unitstate(char12, jnunit, instate2)
        end if

        isw = 0
        i1 = 0
        i2 = 0
        do k = jcol + 1, CARRAY_MXCOLM
           if (cwkinp(iline)(k:k) /= ' ' .and. cwkinp(iline)(k:k) /= ',') then
              if(isw == 0) then
                 i1 = k
                 isw = 1
              end if
              i2 = k

           else if (isw == 1) then

              call lchar2itype(char00, i1, i2, itype)
              if(itype /= 0) then
                 
                 do iu = inunit(1), inunit(2)
                    
                    
                    if (iu > MXUNIT) then
                       error_message = 'Error: MXUNIT is too small.'
                       call util_error(ERROR%STOP_ALL, error_message)
                    endif
                    do ju = jnunit(1), jnunit(2)
                       if (iu > MXUNIT) then
                          error_message = 'Error: MXUNIT is too small.'
                          call util_error(ERROR%STOP_ALL, error_message)
                       endif
                       iunit = ius2unit(iu, instate)
                       junit = ius2unit(ju, instate2)
                       inmisc%flag_local_unit(iunit, junit, itype) = .TRUE.
                    end do
                 end do
              end if

              isw = 0
           end if
        end do
     end if

     ! NLOCAL
     if(char00(1:6) == 'NLOCAL') then
        do icol = 8, CARRAY_MXCOLM
           if(char00(icol:icol) == '/') exit
        end do
        do jcol = 8, CARRAY_MXCOLM
           if(char00(jcol:jcol) == ')') exit
        end do

        read (char00(8:icol-1), '(a)') char12
        call util_unitstate(char12, inunit, instate)
        read (char00(icol+1:jcol-1), '(a)') char12
        call util_unitstate(char12, jnunit, instate2)

        isw = 0
        i1 = 0
        i2 = 0
        do k = jcol + 1, CARRAY_MXCOLM
           if (cwkinp(iline)(k:k) /= ' ' .and. cwkinp(iline)(k:k) /= ',') then
              if(isw == 0) then
                 i1 = k
                 isw = 1
              end if
              i2 = k

           else if (isw == 1) then

              call char2itype(char00, i1, i2, itype)
              
              if(itype /= 0) then
                 do iu = inunit(1), inunit(2)
                    if (iu > MXUNIT) then
                       error_message = 'Error: MXUNIT is too small.'
                       call util_error(ERROR%STOP_ALL, error_message)
                    endif
                    do ju = jnunit(1), jnunit(2)
                       if (iu > MXUNIT) then
                          error_message = 'Error: MXUNIT is too small.'
                          call util_error(ERROR%STOP_ALL, error_message)
                       endif
                       iunit = ius2unit(iu, instate)
                       junit = ius2unit(ju, instate2)
                       inmisc%flag_nlocal_unit(iunit, junit, INTERACT%NOTHING) = .FALSE.
                       inmisc%flag_nlocal_unit(iunit, junit, itype) = .TRUE.
                    end do
                 end do
              end if

              isw = 0
           end if
        end do
     end if
  end do

  ! check nonlocal energy function
  call check_nlocal()

  ! check and modify local energy function
  call check_local()

#ifdef MPI_PAR
  end if
  call MPI_Bcast(inmisc%flag_local_unit,MXUNIT*MXUNIT*LINTERACT%MAX,MPI_LOGICAL,0 ,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc%flag_nlocal_unit,MXUNIT*MXUNIT*INTERACT%MAX,MPI_LOGICAL,0 ,MPI_COMM_WORLD,ierr)
#endif


  ! ----------------------------------------------------------------------
  ! for check the force type
  inmisc%force_flag(1:INTERACT%MAX) = .false.
  do iunit = 1, nunit_all
     do junit = iunit, nunit_all
        do inum = 1, INTERACT%MAX
           if(inmisc%flag_nlocal_unit(iunit, junit, inum)) then
              inmisc%force_flag(inum) = .true.
           end if
        end do
     end do
  end do

  ! ----------------------------------------------------------------------
  ! for check the local force type   !AICG
  inmisc%force_flag_local(1:LINTERACT%MAX) = .false.
  do iunit = 1, nunit_all
     do inum = 1, LINTERACT%MAX
        if(inmisc%flag_local_unit(iunit, iunit, inum)) then
           inmisc%force_flag_local(inum) = .true.
        end if
     end do
  end do

  ! ----------------------------------------------------------------------
  ! for enm
!  if(inmisc%force_flag(INTERACT%ENM)) then
!     do inum = 1, INTERACT%MAX
!        if(inmisc%force_flag(inum) .and. inum /= INTERACT%ENM .and. inum /= INTERACT%NOTHING) then
!           error_message = 'Error: elastic network model cannot be mixed'
!           call util_error(ERROR%STOP_ALL, error_message)
!        end if
!     end do
!  end if


!  ! ----------------------------------------------------------------------
!  ! for Multi Go
!  if(inmgo%i_multi_mgo == 1) then
!
!#ifdef MPI_PAR
!     if (myrank == 0) then
!#endif
!
!     inmgo%nsystem_mgo = 0
!     inmgo%isysmbr_mgo(1:MXSYSTEM_MGO, 1:MXSTATE_MGO, 1:MXACT_MGO) = 0
!     inmgo%nstate_mgo(1:MXSYSTEM_MGO) = 0
!     inmgo%nactnum_mgo(1:MXSYSTEM_MGO) = 0
!     inmgo%iactmat_mgo(1:MXUNIT, 1:MXUNIT) = 0
!
!     iact = 0
!     do iline = 1, nlines
!        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
!        char00 = cwkinp(iline)
!
!        if(char00(1:14) == 'MULTIGO_SYSTEM') then
!           do icol = 16, CARRAY_MXCOLM
!              if(char00(icol:icol) == ')') exit
!           end do
!           read (char00(16:icol-1), '(a)') char12
!           call util_unitstate(char12, inunit, instate)
!           isys = inunit(1)
!           istat = instate
!           if(inmgo%nsystem_mgo < isys) inmgo%nsystem_mgo = isys
!           if(inmgo%nstate_mgo(isys) < istat) inmgo%nstate_mgo(isys) = istat
!
!           iactnum = 0
!           isw = 0
!           i1 = 0
!           i2 = 0
!           do k = icol + 1, CARRAY_MXCOLM
!              if (cwkinp(iline)(k:k) /= ' ' .and. cwkinp(iline)(k:k) /= ',') then
!                 if(isw == 0) then
!                    i1 = k
!                    isw = 1
!                 end if
!                 i2 = k
!
!              else if (isw == 1) then
!                 isw2 = 0
!                 do jcol = i1, i2
!                    if(char00(jcol:jcol) == '/') isw2 = jcol
!                 end do
!
!                 read (char00(i1:isw2-1), '(a)') char12
!                 call util_unitstate(char12, inunit, instate)
!                 read (char00(isw2+1:i2), '(a)') char12
!                 call util_unitstate(char12, jnunit, instate2)
!
!                 do iu = inunit(1), inunit(2)
!                    do ju = jnunit(1), jnunit(2)
!                       iact = iact + 1
!                       iactnum = iactnum + 1
!                       iunit = ius2unit(iu, instate)
!                       junit = ius2unit(ju, instate2)
!                       inmgo%isysmbr_mgo(isys, istat, iactnum) = iact
!                       inmgo%iactmat_mgo(iunit, junit) = iact
!                    end do
!                 end do
!              
!                 i1 = 0
!                 isw = 0
!              end if
!           end do
!
!           inmgo%nactnum_mgo(isys) = iactnum
!        end if
!     end do
!
!     inmgo%nact_mgo = iact
!
!#ifdef MPI_PAR
!     end if
!     call MPI_Bcast (inmgo, inmgo%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
!#endif
!
!  end if


!  ! ----------------------------------------------------------------------
!  ! flexible local potentail
!#ifdef MPI_PAR
!  if (myrank == 0) then
!#endif
!
!  do iline = 1, nlines
!     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
!     
!     do iequa = 1, nequat
!        ! flexible local potential
!        cvalue = 'i_flp'
!        call ukoto_ivalue2(lunout, csides(1, iequa), &
!             inflp%i_flp, cvalue)
!        
!     end do
!  end do
!
!#ifdef MPI_PAR
!  end if
!  call MPI_Bcast(inflp%i_flp,     1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!#endif

  ! ----------------------------------------------------------------------
  ! atom selection, energy output style, ....
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)

     do iequa = 1, nequat
        
        cvalue = 'i_use_atom_protein'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmisc%i_use_atom_protein, cvalue)

        cvalue = 'i_residuenergy_radii'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmisc%i_residuenergy_radii, cvalue)

        cvalue = 'i_output_energy_style'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmisc%i_output_energy_style, cvalue)

        cvalue = 'i_triple_angle_term'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmisc%i_triple_angle_term, cvalue)

        cvalue = 'i_coef_from_ninfo'
        if(csides(1, iequa) == cvalue) then
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                i_ninfo, cvalue)
           !if (i_ninfo == 0) then
           !   inmisc%flg_coef_from_ninfo = .false.
           !else if (i_ninfo == 1) then
           if (i_ninfo == 1) then
              inmisc%flg_coef_from_ninfo = .true.
           else
              error_message = 'Error: i_coef_from_ninfo must be 1'
              call util_error(ERROR%STOP_ALL, error_message)
           end if
        end if

        cvalue = 'i_temp_independent'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmisc%i_temp_independent, cvalue)

        cvalue = 'i_dtrna_model'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inmisc%i_dtrna_model, cvalue)

     end do
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc%flg_coef_from_ninfo,   1,MPI_LOGICAL,0 ,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc%i_use_atom_protein,    1,MPI_INTEGER,0 ,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc%i_residuenergy_radii,   1,MPI_INTEGER,0 ,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc%i_output_energy_style, 1,MPI_INTEGER,0 ,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc%i_triple_angle_term,   1,MPI_INTEGER,0 ,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc%i_temp_independent,    1,MPI_INTEGER,0 ,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc%i_dtrna_model,         1,MPI_INTEGER,0 ,MPI_COMM_WORLD,ierr)
  
#endif

!  ! ----------------------------------------------------------------------
!  ! ----------------------------------------------------------------------
!  ! read i_aicg  !AICG
!#ifdef MPI_PAR
!  if (myrank == 0) then
!#endif
!
!  if (inmisc%force_flag(INTERACT%AICG1) .OR. &
!      inmisc%force_flag(INTERACT%AICG2) .OR. &
!      inmisc%force_flag_local(LINTERACT%L_AICG1) .OR. &
!      inmisc%force_flag_local(LINTERACT%L_AICG2) .OR. &
!      inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
!
!     rewind(luninp)
!     call ukoto_uiread2(luninp, lunout, 'aicg            ', kfind, &
!          CARRAY_MXLINE, nlines, cwkinp)
!     if(kfind /= 'FIND') then 
!        error_message = 'Error: cannot find "aicg " field in the input file'
!        call util_error(ERROR%STOP_ALL, error_message)
!     end if
!     
!     do iline = 1, nlines
!        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
!        
!        do iequa = 1, nequat
!           ! aicg
!           cvalue = 'i_aicg'
!           call ukoto_ivalue2(lunout, csides(1, iequa), &
!                i_aicg, cvalue)
!        end do
!     end do
!     
!  end if
!
!  if (inmisc%force_flag_local(LINTERACT%L_AICG2)) then
!     error_message = 'Warning: L_AICG2 interaction may cause mirror image problem. ' // &
!                     'You should use L_AICG2_PLUS interaction instead of L_AICG2 interaction.'
!     call util_error(ERROR%WARN_ALL, error_message)
!  end if
!
!#ifdef MPI_PAR
!  end if
!  call MPI_Bcast(i_aicg,          1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!#endif


  ! ----------------------------------------------------------------------
  ! write in lunout 

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  ! -----------------------------------------------------------------
  ! output energy function
!  write (lunout, '(72(1H*))')
!  write (lunout, '(a)') &
!       '**** interaction type and coefficient between units ****'
!  write (lunout, '(a)') '<type>'
!  write (lunout, '(a)') '--------------------------------------------'
!  write (lunout, '(a4, a2, 200i4)') 'unit','|', (i , i = 1, nunit_all)
!  write (lunout, '(a)') '--------------------------------------------'
!  do iunit = 1, nunit_all
!     write (lunout, '(i4, a2, 200i4)') &
!          iunit, '|', (itype_nlocal(junit, iunit), junit = 1, iunit)
!  end do
!
!
!  write (lunout, '(a)') '--------------------------------------------'
!  write (lunout, '(i2, a)') 1, ' : no interation'
!  write (lunout, '(i2, a)') INTERACT%GO, ' : 12-10 Go potential'
!  write (lunout, '(i2, a)') INTERACT%EXV, ' : (c/r)**12 repulsion'
!  write (lunout, '(i2, a)') INTERACT%ELE, ' : electrostatic interaction'
!  write (lunout, '(i2, a)') INTERACT%ENM, ' : elastic network model'
!  write (lunout, '(i2, a)') INTERACT%HP, ' : hydrophobic interaction'
!  write (lunout, '(i2, a)') INTERACT%SASA, ' : sasa interaction'  !sasa
!  write (lunout, '(a)') ''


!  ! -----------------------------------------------------------------
!  ! AICG
!  write (lunout, *) 'i_aicg = ', i_aicg
!  if(i_aicg == 1) then
!     write (lunout, *) 'aicg parameters are given by cafemol'
!  else if(i_aicg == 2) then
!     write (lunout, *) 'aicg parameters are given by users'
!  else
!     error_message = 'Error: invalid value about i_aicg'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if
  ! -----------------------------------------------------------------
  ! using atoms of protein, inmisc%i_use_atom_protein
  if(inmisc%i_use_atom_protein == 0) then
     write (lunout, *) 'using Ca'
  else if(inmisc%i_use_atom_protein == 1) then
     write (lunout, *) 'using Cb'
  else if(inmisc%i_use_atom_protein == 2) then
     write (lunout, *) 'using the center of mass '
  else
     error_message = 'Error: invalid value for inmisc%i_use_atom_protein'
     call util_error(ERROR%STOP_ALL, error_message)
  end if


  ! -----------------------------------------------------------------
  ! parameter set for excluded volume interaction: inmisc%i_residuenergy_radii
  if(inmisc%i_residuenergy_radii == 0) then
     write (lunout, *) 'using constant radii for excluded volume interactions'
  else if(inmisc%i_residuenergy_radii == 1) then
     write (lunout, *) 'using residue type dependent radii for excluded volume interactions'
  else
     error_message = 'Error: invalid value for inmisc%i_residuenergy_radii'
     call util_error(ERROR%STOP_ALL, error_message)
  end if


  ! -----------------------------------------------------------------
  if(inmisc%i_output_energy_style == 0) then
     write (lunout, *) 'output inter energy are summed up intra energy half and half'
  else if(inmisc%i_output_energy_style == 1) then
     write (lunout, *) 'output intra and inter energy separately'
  else
     error_message = 'Error: invalid value for inmisc%i_output_energy_style'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! -----------------------------------------------------------------
  if(inmisc%i_triple_angle_term == 0) then
     write (lunout, *) 'no triple-angle term in dihedral angle function'
  else if (inmisc%i_triple_angle_term == 1 .or. inmisc%i_triple_angle_term == 2) then
     continue !default
  else
     error_message = 'Error: invalid value for inmisc%i_triple_angle_term'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  ! -----------------------------------------------------------------
  if (inmisc%flg_coef_from_ninfo) then
     write (lunout, *) 'i_coef_from_ninfo = 1; "coef" values in ninfo file will be used'    
  endif

  ! -----------------------------------------------------------------
  if(inmisc%i_temp_independent == 0) then
     continue !default
  else if(inmisc%i_temp_independent == 1) then
     write (lunout, *) 'i_temp_independent = 1: Calculate temperature-independent energy'
  else if(inmisc%i_temp_independent == 2) then
     write (lunout, *) 'i_temp_independent = 2: Calculate temperature-independent energy'
  else
     error_message = 'Error: invalid value for inmisc%i_temp_independent'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  ! -----------------------------------------------------------------
  if(inmisc%i_dtrna_model == 2013) then
     write (lunout, *) 'i_dtrna_model = 2013'
  else if(inmisc%i_dtrna_model == 2015) then
     write (lunout, *) 'i_dtrna_model = 2015'
  else
     error_message = 'Error: invalid value for inmisc%i_dtrna_model'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

#ifdef MPI_PAR
  end if
#endif


contains


  subroutine lchar2itype(char00, i1, i2, itype)

    implicit none

    character(CARRAY_MXCOLM) :: char00
    integer, intent(in) :: i1, i2
    integer, intent(out) :: itype

    if(char00(i1:i2) == 'NOTHING') then
       itype = LINTERACT%NOTHING

    else if(char00(i1:i2) == 'L_GO') then
       itype = LINTERACT%L_GO

    else if(char00(i1:i2) == 'L_FENE') then
       itype = LINTERACT%L_FENE

!    else if(char00(i1:i2) == 'L_AICG1') then
!       itype = LINTERACT%L_AICG1

!    else if(char00(i1:i2) == 'L_AICG2') then
!       itype = LINTERACT%L_AICG2

!    else if(char00(i1:i2) == 'L_AICG2_PLUS') then
!       itype = LINTERACT%L_AICG2_PLUS

!    else if(char00(i1:i2) == 'L_FLP') then
!       itype = LINTERACT%L_FLP

    else if(char00(i1:i2) == 'L_BOND') then
       itype = LINTERACT%L_BOND

    else if(char00(i1:i2) == 'L_ROUSE') then
       itype = LINTERACT%L_ROUSE

!    else if(char00(i1:i2) == 'L_ENM') then
!       itype = LINTERACT%L_ENM

    else if(char00(i1:i2) == 'L_RIGID_LIG') then
       itype = LINTERACT%L_RIGID_LIG

    else if(char00(i1:i2) == 'L_DTRNA') then
       itype = LINTERACT%L_DTRNA

    else
       itype = 0

    end if

  end subroutine lchar2itype


  subroutine char2itype(char00, i1, i2, itype)

    implicit none

    character(CARRAY_MXCOLM) :: char00
    integer, intent(in) :: i1, i2
    integer, intent(out) :: itype

    if(char00(i1:i2) == 'NOTHING') then
       itype = INTERACT%NOTHING

    else if(char00(i1:i2) == 'GO') then
       itype = INTERACT%GO

    else if(char00(i1:i2) == 'LJ') then
       itype = INTERACT%LJ

    else if(char00(i1:i2) == 'EXV12') then
       itype = INTERACT%EXV12

    else if(char00(i1:i2) == 'EXV6') then
       itype = INTERACT%EXV6

    else if(char00(i1:i2) == 'ELE') then
       itype = INTERACT%ELE

!    else if(char00(i1:i2) == 'ENM') then
!       itype = INTERACT%ENM

!    else if(char00(i1:i2) == 'HP') then
!       itype = INTERACT%HP

!    else if(char00(i1:i2) == 'MORSE') then
!       itype = INTERACT%MORSE

!    else if(char00(i1:i2) == 'PAIR_RNA') then
!       itype = INTERACT%PAIR_RNA
!
!    else if(char00(i1:i2) == 'AICG1') then
!       itype = INTERACT%AICG1

!    else if(char00(i1:i2) == 'AICG2') then
!       itype = INTERACT%AICG2

!    else if(char00(i1:i2) == 'SASA') then
!       itype = INTERACT%SASA

    else if(char00(i1:i2) == 'DTRNA') then
       itype = INTERACT%DTRNA

    else if(char00(i1:i2) == 'EXV_WCA') then
       itype = INTERACT%EXV_WCA

    else if(char00(i1:i2) == 'EXV_DT15') then
       itype = INTERACT%EXV_DT15

    else
       itype = 0

    end if

  end subroutine char2itype


  subroutine check_nlocal()

    use const_index
    use var_struct, only : iclass_unit
    implicit none

    integer :: n_go
    integer :: iforce

    do iunit = 1, nunit_all
       do junit = 1, nunit_all

          n_go = 0
          do iforce = 1, LINTERACT%MAX
             if(.not. inmisc%flag_nlocal_unit(iunit, junit, iforce)) then
                cycle
             end if
             
             !if(iforce == INTERACT%GO .or. iforce == INTERACT%LJ .or. &
             !   iforce == INTERACT%AICG1 .or. iforce == INTERACT%AICG2) then
             if(iforce == INTERACT%GO .or. iforce == INTERACT%LJ) then
                n_go = n_go + 1
             endif
                
!             if(iforce == INTERACT%ENM) then
!                if(iclass_unit(iunit) /= CLASS%PRO) then
!                   error_message = 'Error: ENM interaction is only applicable for protein'
!                   call util_error(ERROR%STOP_ALL, error_message)
!                end if
                
!             else if(iforce == INTERACT%AICG1 .or. iforce == INTERACT%AICG2) then
!                if(iclass_unit(iunit) /= CLASS%PRO) then
!                   error_message = 'Error: AICG interaction is only applicable for protein'
!                   call util_error(ERROR%STOP_ALL, error_message)
!                end if
                
!             else if(iforce == INTERACT%PAIR_RNA) then
!                if(iclass_unit(iunit) /= CLASS%RNA) then
!                   error_message = 'Warning: PAIR_RNA interaction is only applicable for RNA'
!                   call util_error(ERROR%WARN_ALL, error_message)
!                end if
                
             !else if(iforce == INTERACT%DTRNA) then
             if(iforce == INTERACT%DTRNA) then
                if(iclass_unit(iunit) /= CLASS%RNA) then
                   error_message = 'Error: DTRNA interactions are only applicable for RNA'
                   call util_error(ERROR%STOP_ALL, error_message)
                end if
             end if
          

          end do

          if(n_go >= 2) then
             error_message = 'Error: should not specify more than one Go-type interaction for each unit-unit interaction'
             call util_error(ERROR%STOP_ALL, error_message)
          end if

       end do
    end do

  end subroutine check_nlocal


  subroutine check_local()

    use const_index
    use var_struct, only : iclass_unit
    implicit none

    integer :: n_local_input
    integer :: iforce

    do iunit = 1, nunit_all

       n_local_input = 0
       do iforce = 1, LINTERACT%MAX
          if(inmisc%flag_local_unit(iunit, iunit, iforce)) then
             n_local_input = n_local_input + 1
          end if
       end do

       ! default setting of local interaction from nonlocal_interaction
       if(n_local_input == 0) then

!          if(iclass_unit(iunit) == CLASS%PRO) then
!             if(inmisc%flag_nlocal_unit(iunit, iunit, INTERACT%GO)) then
!                inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_GO) = .TRUE.
!             else if(inmisc%flag_nlocal_unit(iunit, iunit, INTERACT%AICG1)) then
!                inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG1) = .TRUE.
!             else if(inmisc%flag_nlocal_unit(iunit, iunit, INTERACT%AICG2)) then
!                inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) = .TRUE.
!                inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS) = .TRUE.
!             else if(inmisc%flag_nlocal_unit(iunit, iunit, INTERACT%ENM)) then
!                inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_ENM) = .TRUE.
!             else
!!                inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_BOND) = .TRUE.
!                inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_GO) = .TRUE.
!             end if
!
!          else if(iclass_unit(iunit) == CLASS%RNA) then
!             if(inmisc%flag_nlocal_unit(iunit, iunit, INTERACT%GO)) then
!                inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_GO) = .TRUE.
!             else
!                inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_BOND) = .TRUE.
!             end if
!             
!          else if(iclass_unit(iunit) == CLASS%LIG) then
!             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_RIGID_LIG) = .TRUE.
!
!          else if(iclass_unit(iunit) == CLASS%ION) then
!             inmisc%flag_local_unit(iunit, iunit, LINTERACT%NOTHING) = .TRUE.
!
!          end if

          error_message = 'Error: should specify one local interaction for each unit'
          call util_error(ERROR%STOP_ALL, error_message)

       else if(n_local_input == 1) then

          if(inmisc%flag_local_unit(iunit,iunit,LINTERACT%L_DTRNA)) then
             if(iclass_unit(iunit) /= CLASS%RNA) then
                error_message = 'Error: L_DTRNA interactions is only applicable for RNA'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          endif

       else if(n_local_input >= 2) then

          error_message = 'Error: should not specify more than one local interaction for each unit'
          call util_error(ERROR%STOP_ALL, error_message)
          
       end if

!       ! add flp to local AICG2
!       if(iclass_unit(iunit) == CLASS%PRO) then  !AICG2
!          if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .OR. &
!             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS)) then !AICG2
!             inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_FLP) = .TRUE.  !AICG2
!          end if !AICG2
!       end if !AICG2

    end do

  end subroutine check_local

end subroutine inp_energy_func
