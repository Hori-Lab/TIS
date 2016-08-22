! mloop_flexible_local
!> @brief This subroutine is to read and set the parameters for flexible local 
!>        potential

! *****************************************************************************
subroutine mloop_flexible_local()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : inflp
  use var_struct, only: nres, nfba, nfdih, nmp_all, &
       fba_para_x, fba_para_y, fba_para_y2, fdih_para, &
       ifba2mp, ifdih2mp, &
       fdih_ener_corr, fba_ener_corr, &
       fba_min_th, fba_max_th, fba_min_th_ener, fba_max_th_ener, &
       ifba2ba, ifdih2dih !flp_mgo

#ifdef MPI_PAR
  use mpiconst
  use var_struct, only : factor_ba, factor_dih, coef_ba, coef_dih, &
                         factor_aicg13, coef_aicg13_gauss, &
                         factor_aicg14, coef_aicg14_gauss, coef_dih_gauss
#endif

  implicit none

  !----------------------------------------------------------------------------
  ! Local variables
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  integer :: ier
  integer :: i_ang1, i_ang2, i_dih
  integer :: iba, idih, ip
  integer :: iflp_mp(MXMP), idel_mp(MXMP)

  character(5) :: header
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  real(PREC), allocatable, save :: fba_mult_para_y1(:,:)  !(10, nres)
  real(PREC), allocatable, save :: fba_mult_para_y2(:,:)  !(10, nres)
  integer,    allocatable, save :: fba_mult_resID1(:)     !(nres)
  integer,    allocatable, save :: fba_mult_resID2(:)     !(nres)
  real(PREC), allocatable, save :: fdih_mult_para(:,:) !(7, nres)
  integer,    allocatable, save :: fdih_mult_resID(:)  !(nres)

  !----------------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  ! Allocate ifba2mp and ifdih2mp
  allocate( ifba2mp(3, MXMPBA*nmp_all), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  ifba2mp(:,:) = 0

  allocate( ifdih2mp(4, MXMPDIH*nmp_all), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  ifdih2mp(:,:) = 0

  ! Allocate ifba2ba and ifdih2dih  !flp_mgo
  allocate( ifba2ba(MXMPBA*nmp_all), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  ifba2ba(:) = 0

  allocate( ifdih2dih(MXMPDIH*nmp_all), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  ifdih2dih(:) = 0
  
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     iflp_mp(:) = 0
     idel_mp(:) = 0
     ! Set ifba2mp and ifdih2mp (Note that nfba and nfdih also change)
     call set_flp_ba_dih(iflp_mp, idel_mp) ! Set the region flexible local potential will be used
     call set_ifba2mp(iflp_mp, idel_mp)
     call set_ifdih2mp(iflp_mp, idel_mp)

#ifdef MPI_PAR
  end if   

  ! Transfer changed variables to all node
  call MPI_Bcast(iflp_mp, MXMP, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(idel_mp, MXMP, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ifba2mp, 3*MXMPBA*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ifdih2mp, 4*MXMPDIH*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nfba, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nfdih, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(factor_ba, MXMPBA*nmp_all, PREC_MPI,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(factor_dih, MXMPDIH*nmp_all, PREC_MPI,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(factor_aicg13, MXMPBA*nmp_all, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(factor_aicg14, MXMPDIH*nmp_all,PREC_MPI, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_ba, 2*MXMPBA*nmp_all, PREC_MPI,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_dih, 2*MXMPDIH*nmp_all, PREC_MPI,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_aicg13_gauss, MXMPBA*nmp_all, PREC_MPI,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_aicg14_gauss, MXMPDIH*nmp_all, PREC_MPI,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_dih_gauss, MXMPDIH*nmp_all, PREC_MPI,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ifba2ba, MXMPBA*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   !flp_mgo
  call MPI_Bcast(ifdih2dih, MXMPDIH*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)  !flp_mgo

#endif

  ! Allocate arrays which contain parameters for flexible local potential
  allocate( fba_para_x(10, nfba), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fba_para_x(:,:) = 0
  allocate( fba_para_y(10, nfba), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fba_para_y(:,:) = 0
  allocate( fba_para_y2(10, nfba), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fba_para_y2(:,:) = 0
  allocate( fba_ener_corr(nfba), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fba_ener_corr(:) = 0
  allocate( fba_max_th(nfba), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fba_max_th(:) = 0
  allocate( fba_min_th(nfba), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fba_min_th(:) = 0
  allocate( fba_max_th_ener(nfba), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fba_max_th_ener(:) = 0
  allocate( fba_min_th_ener(nfba), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fba_min_th_ener(:) = 0
  allocate( fdih_para(7, nfdih), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fdih_para(:,:) = 0
  allocate( fdih_ener_corr(nfdih), stat=ier)
  if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
  fdih_ener_corr(:) = 0
  
#ifdef MPI_PAR
  if (myrank == 0) then
#endif
  
     call set_para_fba()
     call set_para_fdih()
     
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'flexible_local  ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)

     if (kfind /= 'FIND') then
        error_message = 'Error: cannot find "flexible_local" field in setp_flexible_local'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     error_message = 'failed in memory allocation at allocate_nativestruct, PROGRAM STOP'

     allocate( fba_mult_resID1(nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fba_mult_resID1(:) = 0
     allocate( fba_mult_resID2(nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fba_mult_resID2(:) = 0
     allocate( fdih_mult_resID(nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fdih_mult_resID(:) = 0
     allocate( fba_mult_para_y1(10, nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fba_mult_para_y1(:,:) = 0
     allocate( fba_mult_para_y2(10, nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fba_mult_para_y2(:,:) = 0
     allocate( fdih_mult_para(7, nres), stat=ier)
     if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
     fdih_mult_para(:,:) = 0

     i_ang1 = 1
     i_ang2 = 1
     i_dih  = 1

     ! Default value setting
     inflp%k_ang = 1.0e0_PREC
     inflp%k_dih = 1.0E0_PREC
     
     inflp%boundary = 0
     inflp%boundary_lower_limit = 1.20e0_PREC
     inflp%boundary_upper_limit = 2.80e0_PREC
     inflp%boundary_max_force   = 2.00e2_PREC
     
     do iline = 1, nlines
        ctmp00 = cwkinp(iline)
        if (ctmp00(1:5) == 'ANGY1') then
           read (ctmp00, *) header, fba_mult_resID1(i_ang1), &
                fba_mult_para_y1(1,  i_ang1), &
                fba_mult_para_y1(2,  i_ang1), &
                fba_mult_para_y1(3,  i_ang1), &
                fba_mult_para_y1(4,  i_ang1), &
                fba_mult_para_y1(5,  i_ang1), &
                fba_mult_para_y1(6,  i_ang1), &
                fba_mult_para_y1(7,  i_ang1), &
                fba_mult_para_y1(8,  i_ang1), &
                fba_mult_para_y1(9,  i_ang1), &
                fba_mult_para_y1(10, i_ang1)
           i_ang1 = i_ang1 + 1
        else if (ctmp00(1:5) == 'ANGY2') then
           read (ctmp00, *) header, fba_mult_resID2(i_ang2), &
                fba_mult_para_y2(1,  i_ang2), &
                fba_mult_para_y2(2,  i_ang2), &
                fba_mult_para_y2(3,  i_ang2), &
                fba_mult_para_y2(4,  i_ang2), &
                fba_mult_para_y2(5,  i_ang2), &
                fba_mult_para_y2(6,  i_ang2), &
                fba_mult_para_y2(7,  i_ang2), &
                fba_mult_para_y2(8,  i_ang2), &
                fba_mult_para_y2(9,  i_ang2), &
                fba_mult_para_y2(10, i_ang2)
           i_ang2 = i_ang2 + 1
        else if (ctmp00(1:4) == 'DIH') then
           read (ctmp00, *) header, fdih_mult_resID(i_dih), &
                fdih_mult_para(1,  i_dih), &
                fdih_mult_para(2,  i_dih), &
                fdih_mult_para(3,  i_dih), &
                fdih_mult_para(4,  i_dih), &
                fdih_mult_para(5,  i_dih), &
                fdih_mult_para(6,  i_dih), &
                fdih_mult_para(7,  i_dih)
           i_dih = i_dih + 1
        else 
           call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
           
           do iequa = 1, nequat
              cvalue = 'k_ang'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inflp%k_ang, cvalue)

              cvalue = 'k_dih'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inflp%k_dih, cvalue)

              cvalue = 'boundary'
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   inflp%boundary, cvalue)

              cvalue = 'boundary_lover_limit'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inflp%boundary_lower_limit, cvalue)

              cvalue = 'boundary_upper_limit'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inflp%boundary_upper_limit, cvalue)

              cvalue = 'boundary_cutoff'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inflp%boundary_cutoff, cvalue)

              cvalue = 'boundary_max_force'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inflp%boundary_max_force, cvalue)
           end do
        end if
     end do

     do iba = 1, nres
        if (fba_mult_resID1(iba) /= 0) then
           do ip = 1, 10
              fba_para_y(ip, fba_mult_resID1(iba)) = fba_mult_para_y1(ip, iba)
           end do
        end if
     end do
    
     do iba = 1, nres
        if (fba_mult_resID2(iba) /= 0) then
           do ip = 1, 10
              fba_para_y2(ip, fba_mult_resID2(iba)) = fba_mult_para_y2(ip, iba)
           end do
        end if
     end do

     do idih = 1, nres
        if (fdih_mult_resID(idih) /= 0) then
           do ip = 1, 7
              fdih_para(ip, fdih_mult_resID(idih)) = fdih_mult_para(ip, idih)
           end do
        end if
     end do

     ! Set the values of correction term array for setting the energy minimum to 0
     ! For dihedral angle
     call set_fdih_ener_corr()
     ! For bond angle
     call set_fba_ener_corr()

!     do idih = 1, nfba
!        write(*, *) fba_min_th(idih), fba_max_th(idih), fba_min_th_ener(idih), fba_max_th_ener(idih)
!     end do

#ifdef MPI_PAR
  end if

  !  call MPI_Bcast (inflp%k_ang, 1, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
  !  call MPI_Bcast (inflp%k_dih, 1, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(inflp, inflp%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_para_x, 10*nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_para_y, 10*nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_para_y2, 10*nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_ener_corr, nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_max_th, nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_min_th, nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_max_th_ener, nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fba_min_th_ener, nfba, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fdih_para, 7*nfdih, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fdih_ener_corr, nfdih, PREC_MPI, 0, MPI_COMM_WORLD, ierr)

#endif

end subroutine mloop_flexible_local

subroutine set_flp_ba_dih(iflp_mp, idel_mp)

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile, ius2unit
  use var_setp, only : inflp, inmisc
  use var_struct, only : nmp_real, nunit_all, lunit2mp, imp2unit
  use var_mgo, only : inmgo
  
  implicit none

  integer, intent(inout)  :: iflp_mp(MXMP), idel_mp(MXMP)

  !----------------------------------------------------------------------------
  ! Local variables
  integer :: icol
  integer :: luninp, lunout
  integer :: iline, nlines, nequat
  integer :: inunit(2), instate
  integer :: imp, iunit
  integer :: imp_shadow, iunit_real, iunit_shadow
  integer :: i_del_lgo, i_add_flp, i_del_flp, num_char

  character(4) :: kfind
  character(12) :: char12
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message
  !----------------------------------------------------------------------------
  
  luninp = infile%inp
  lunout = outfile%data
!  inflp%nflp = 0

  do iunit = 1, nunit_all
     if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_FLP)) then
        do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
           iflp_mp(imp) = 1
!           idel_mp(imp) = 1
        end do
     end if
  end do


  rewind(luninp)
  ! Find out flexible_local field
  call ukoto_uiread2(luninp, lunout, 'flexible_local  ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "flexible_local" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! Read file
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     ctmp00 = cwkinp(iline)

     i_del_lgo = 0
     i_add_flp = 0
     i_del_flp = 0
     num_char = 0
     if(ctmp00(1:15) == 'DEL_LGO_ADD_FLP') then
        i_del_lgo = 1
        i_add_flp = 1
        i_del_flp = 0
        num_char = 15
     else if(ctmp00(1:7) == 'ADD_FLP') then
        i_del_lgo = 0
        i_add_flp = 1
        i_del_flp = 0
        num_char = 7
     else if(ctmp00(1:7) == 'DEL_FLP') then
        i_del_lgo = 0
        i_add_flp = 0
        i_del_flp = 1
        num_char = 7
     end if

     if(i_del_lgo == 1 .or. i_add_flp == 1 .or. i_del_flp == 1) then
        do icol = num_char + 2, CARRAY_MXCOLM
           if (ctmp00(icol:icol) == ')') exit
        end do
        read (ctmp00(num_char + 2:icol-1), *) char12
        write (lunout, '(4a)') '---reading flp: ', ctmp00(1:num_char), ' ', trim(char12)
        call util_unitstate(char12, inunit, instate)


        if(inunit(2) > nmp_real) then
           error_message = 'Error: flp should not specify shadow chain'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        if(inflp%i_flp /= 1) then
           do imp = inunit(1), inunit(2)
              iunit = imp2unit(imp)
              if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_FLP)) then
                 
              else
                 error_message = 'Error: should be use i_flp option'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end do
        end if

        do imp = inunit(1), inunit(2)
           if(i_del_lgo == 1) then
              idel_mp(imp) = 1
           end if

           if(i_add_flp == 1) then
              iflp_mp(imp) = 1
           else if(i_del_flp == 1) then
              iflp_mp(imp) = 0
           end if

           if(inmgo%i_multi_mgo >= 1) then
              iunit_real = imp2unit(imp)
              do instate = 2, MXSTATE_MGO
                 iunit_shadow = ius2unit(iunit_real, instate)
                 if(iunit_shadow <= 1) then
                    exit
                 else
                    imp_shadow = imp + lunit2mp(1, iunit_shadow) - lunit2mp(1, iunit_real)

                    if(i_del_lgo == 1) then
                       idel_mp(imp_shadow) = 1
                    end if
                    
                    if(i_add_flp == 1) then
                       iflp_mp(imp_shadow) = 1
                    else if(i_del_flp == 1) then
                       iflp_mp(imp_shadow) = 0
                    end if
                 end if
              end do
           end if

        end do
     end if
     
  end do

end subroutine set_flp_ba_dih

subroutine set_ifba2mp(iflp_mp, idel_mp)

  use const_maxsize
  use var_struct, only: nba, nfba, ifba2mp, iba2mp, factor_ba, coef_ba, &
                        factor_aicg13, coef_aicg13_gauss, ifba2ba !flp_mgo
  
  implicit none

  integer, intent(inout)  :: iflp_mp(MXMP), idel_mp(MXMP)

  !----------------------------------------------------------------------------
  ! Local variables
  integer :: iba
  integer :: imp1, imp2, imp3
  !----------------------------------------------------------------------------
  nfba = 0
 
  do iba = 1, nba
     imp1 = iba2mp(1, iba)
     imp2 = iba2mp(2, iba)
     imp3 = iba2mp(3, iba)
     if(idel_mp(imp1) == 1 .or. idel_mp(imp2) == 1 .or. idel_mp(imp3) == 1) then
        factor_aicg13(iba) = 0.0e0_PREC
        coef_aicg13_gauss(iba) = 0.0e0_PREC

        factor_ba(iba) = 0.0e0_PREC
        coef_ba(1, iba) = 0.0e0_PREC
        coef_ba(2, iba) = 0.0e0_PREC
     end if
     if(iflp_mp(imp1) == 1 .or. iflp_mp(imp2) == 1 .or. iflp_mp(imp3) == 1) then
        nfba = nfba + 1
        ifba2mp(1, nfba) = imp1
        ifba2mp(2, nfba) = imp2
        ifba2mp(3, nfba) = imp3
        ifba2ba(nfba) = iba   !flp_mgo
     end if
  end do

end subroutine set_ifba2mp

subroutine set_ifdih2mp(iflp_mp, idel_mp)

  use const_maxsize
  use var_struct, only: ndih, nfdih, idih2mp, ifdih2mp, factor_dih, coef_dih, &
                        factor_aicg14, coef_aicg14_gauss, coef_dih_gauss, ifdih2dih !flp_mgo
  
  implicit none

  integer, intent(inout)  :: iflp_mp(MXMP), idel_mp(MXMP)

  !----------------------------------------------------------------------------
  ! Local variables
  integer :: idih
  integer :: imp1, imp2, imp3, imp4
  !----------------------------------------------------------------------------
  nfdih = 0
  
  do idih = 1, ndih
     imp1 = idih2mp(1, idih)
     imp2 = idih2mp(2, idih)
     imp3 = idih2mp(3, idih)
     imp4 = idih2mp(4, idih)
     if(idel_mp(imp1) == 1 .or. idel_mp(imp2) == 1 .or. idel_mp(imp3) == 1 .or. idel_mp(imp4) == 1) then
        factor_aicg14(idih) = 0.0e0_PREC
        coef_aicg14_gauss(idih) = 0.0e0_PREC
        coef_dih_gauss(idih) = 0.0e0_PREC
        factor_dih(idih) = 0.0e0_PREC
        coef_dih(1, idih) = 0.0e0_PREC
        coef_dih(2, idih) = 0.0e0_PREC
     end if
     if(iflp_mp(imp1) == 1 .or. iflp_mp(imp2) == 1 .or. iflp_mp(imp3) == 1 .or. iflp_mp(imp4) == 1) then
        nfdih = nfdih + 1
        ifdih2mp(1, nfdih) = imp1
        ifdih2mp(2, nfdih) = imp2
        ifdih2mp(3, nfdih) = imp3
        ifdih2mp(4, nfdih) = imp4
        ifdih2dih(nfdih) = idih     !flp_mgo
     end if
  end do

end subroutine set_ifdih2mp

subroutine set_para_fba()
  use const_maxsize
  use const_index
  use var_setp, only: inflp
  use var_struct, only: nfba, ifba2mp, &
                        fba_para_x, fba_para_y, fba_para_y2, &
                        cmp2seq
  
  implicit none
  !----------------------------------------------------------------------------
  ! Local variables
  !integer :: ier
  integer :: iba, id, ip
  !character(CARRAY_MSG_ERROR) :: error_message
  !----------------------------------------------------------------------------
  ! Function
  integer :: ifunc_seq2id
  
  do iba = 1, nfba
     id = ifunc_seq2id( cmp2seq( ifba2mp(2, iba) ) )
     do ip = 1, 10
        fba_para_x(ip, iba) = inflp%ang_para_x(ip)
        fba_para_y(ip, iba)  = inflp%ang_para_y(id, ip)
        fba_para_y2(ip, iba) = inflp%ang_para_y2(id, ip)
     end do
  end do

end subroutine set_para_fba

subroutine set_para_fdih()
  use const_maxsize
  use const_index
  use var_setp, only: inflp
  use var_struct, only: nfdih, fdih_para, ifdih2mp, cmp2seq
    
  implicit none
  !----------------------------------------------------------------------------
  ! Local variables
  !integer :: ier
  integer :: idih, ip
  integer :: id1, id2
  !character(CARRAY_MSG_ERROR) :: error_message
  !----------------------------------------------------------------------------
  ! Function
  integer :: ifunc_seq2id
  
  do idih = 1, nfdih
     id1 = ifunc_seq2id( cmp2seq( ifdih2mp(2, idih) ) )
     id2 = ifunc_seq2id( cmp2seq( ifdih2mp(3, idih) ) )

     do ip = 1, 7
        fdih_para(ip, idih) = inflp%dih_para(id1, id2, ip)
     end do
  end do

end subroutine set_para_fdih

subroutine set_fdih_ener_corr()

  use const_maxsize
  use const_physical
  use var_struct, only: nfdih, fdih_para, fdih_ener_corr

  implicit none

  !----------------------------------------------------------------------------
  ! Local variables
  integer :: ifdih
  real(PREC) :: th, dth
  real(PREC) :: energy, min_energy
  !----------------------------------------------------------------------------
  ! initialize local variables
  dth = 1e-4_PREC
  !----------------------------------------------------------------------------

  do ifdih = 1, nfdih
     th = -F_PI
     min_energy = calc_energy(ifdih, th)
     
     do
        th = th + dth
        if (th > F_PI) exit
        energy = calc_energy(ifdih, th)
        min_energy = min(min_energy, energy)
     end do
     
     fdih_ener_corr(ifdih) =  min_energy
     
  end do
contains
  function calc_energy(ifdih, th)

    !----------------------------------------------------------------------------
    real(PREC) :: calc_energy
    integer, intent(in) :: ifdih
    real(PREC), intent(in) :: th
    
    !----------------------------------------------------------------------------
    ! local variables
    real(PREC) :: energy

    !----------------------------------------------------------------------------
    
    energy = fdih_para(1, ifdih)               &
           + fdih_para(2, ifdih) * cos(    th) &
           + fdih_para(3, ifdih) * sin(    th) &
           + fdih_para(4, ifdih) * cos(2 * th) &
           + fdih_para(5, ifdih) * sin(2 * th) &
           + fdih_para(6, ifdih) * cos(3 * th) &
           + fdih_para(7, ifdih) * sin(3 * th)           
    
    calc_energy = energy
  end function
  
end subroutine set_fdih_ener_corr

subroutine set_fba_ener_corr()

  use const_maxsize
  use const_physical
  use var_struct, only: nfba,            &
                        fba_ener_corr,   &
                        fba_min_th,      &
                        fba_max_th,      &
                        fba_min_th_ener, &   
                        fba_max_th_ener

  implicit none

  !----------------------------------------------------------------------------
  ! Local variables
  real(PREC) :: center
  integer    :: ifba
  real(PREC) :: th, dth
  real(PREC) :: energy, min_energy
  real(PREC) :: min_th, max_th
  real(PREC) :: min_th_ener, max_th_ener
  real(PREC) :: force

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! initialize local variables
  dth = 1e-4_PREC
  !----------------------------------------------------------------------------

  center = ( FBA_MAX_ANG - FBA_MIN_ANG ) / 2.0e0_PREC

  do ifba = 1, nfba
     th     = FBA_MIN_ANG
     min_th = FBA_MIN_ANG
     max_th = FBA_MIN_ANG
     call splint (th, ifba, min_energy)
     call dsplint(th, ifba, force)
     min_th_ener = min_energy
     
     do
        th = th + dth
        if (th > FBA_MAX_ANG) exit
        call splint (th, ifba, energy)
        min_energy = min(min_energy, energy)
        call dsplint(th, ifba, force)
        if (force < FBA_MIN_ANG_FORCE) then
           min_th = th
           min_th_ener = energy
        end if
        if (th > center .and. max_th == FBA_MIN_ANG .and. force > FBA_MAX_ANG_FORCE) then
           max_th = th
           max_th_ener = energy
        end if
     end do

     fba_min_th(ifba)      = min_th
     fba_max_th(ifba)      = max_th
     fba_min_th_ener(ifba) = min_th_ener
     fba_max_th_ener(ifba) = max_th_ener
     fba_ener_corr(ifba)   = min_energy

  end do

end subroutine set_fba_ener_corr
