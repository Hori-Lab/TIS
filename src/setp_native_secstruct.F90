! setp_native_secstruct
!> @brief Calculate the secondary structure of protein based on &
!>        the given PDB structure

subroutine setp_native_secstruct()

  use const_maxsize
  use const_physical
  use const_index
  use var_io,    only : outfile, ifile_pdb, num_file
  use var_struct, only : nmp_all, imp2unit, istype_mp
#ifdef MPI_PAR
  use mpiconst
#endif
  
  implicit none

#ifdef MPI_PAR
  integer :: iflag2
#endif

  ! -----------------------------------------------------------------
  ! intent(out) :: istype_mp

  ! -----------------------------------------------------------------
  ! local variables
  integer :: i, iclass
  integer :: input_status, lunout, lunpdb, npdb
  integer :: iunit, junit, imp, impread, imunitemp, jmp
  integer :: ihbond, jhbond, lhbond
  integer :: iatom, ico, ico_2, idev, iflag, inh, inh_2, jfile

  ! to be determined
  ! isec shows whether the inputed pdbfile include N,H,O data
  integer :: nhbond
  integer :: lhbond2mp(2, MXHBOND), ihbond_stype(MXHBOND)
  real(PREC) :: hbonde(MXHBOND)

  real(PREC) :: e, f, rch, rcn, roh, ron, q1, q2
  real(PREC) :: x, y, z, tempfactor, occupancy
  real(PREC) :: xyz_h(SDIM, MXMP)
  real(PREC) :: xyz_n(SDIM, MXMP)
  real(PREC) :: xyz_ca(SDIM, MXMP)
  real(PREC) :: xyz_c(SDIM, MXMP)
  real(PREC) :: xyz_o(SDIM, MXMP)
  character(72) :: char72
  character(1) :: multistruct
  character(4) :: nameatom
  character(6) :: nameid
  character(3) :: nameofmp
  character(2) :: chainid
  character(CARRAY_MSG_ERROR) :: error_message

  ! -----------------------------------------------------------------
  lunout = outfile%data

  ! -----------------------------------------------------------------
  ! for bench mark simulation
  ! skip this subroutine
  do imp = 1, nmp_all
     istype_mp(imp) = 1
  end do
  return
  ! -----------------------------------------------------------------

  npdb = num_file%pdb
  do imp = 1, nmp_all
     istype_mp(imp) = -1
  end do

  imp = 0
  imunitemp = -999
  iflag = 3 

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     do jfile = 1, npdb
        iclass = ifile_pdb(1, jfile)
        if (iclass /= CLASS%PRO) cycle
        lunpdb = ifile_pdb(2, jfile)

        rewind(lunpdb)
   
        do
           read (lunpdb, '(a72)', iostat = input_status) char72  

#ifdef MPI_PAR
           iflag2 = 0
#endif

           if(input_status < 0) then
              exit
           else if(input_status > 0) then
              error_message = 'Error: input error in setp_native_secstruct'
              call util_error(ERROR%STOP_ALL, error_message)
           end if
   
           if(char72(1:4) == 'ATOM') then
              read (char72, '(a6, i5, 1x, a4, a1, a3, a2, i4, 4x, 3f8.3, 2f6.2)') &
                   nameid, iatom, nameatom, multistruct, &
                   nameofmp, chainid, &
                   impread, x, y, z, occupancy, tempfactor
   
              if(multistruct /= ' ' .and. multistruct /= 'A') cycle
              
              if(nameatom == ' N  ' .or. nameatom == ' O  ' .or. &
                   nameatom == ' C  ' .or. nameatom == ' CA ') then
                 
                 if(imunitemp /= impread) then
                    if(iflag /= 3) then
                       write (lunout, *) 'cannot determine secondary structure'
                       write (lunout, *) 'There are few data.'
#ifdef MPI_PAR
                       iflag2 = -1
#endif
                       return
                    endif
                    iflag = 0 
                    imunitemp = impread
                    imp = imp + 1
                 end if
                 
                 if(nameatom == ' N  ') then
                    iflag = iflag + 1 
                    xyz_n(1, imp) = x
                    xyz_n(2, imp) = y
                    xyz_n(3, imp) = z
                 else if(nameatom == ' O  ') then
                    iflag = iflag + 1 
                    xyz_o(1, imp) = x
                    xyz_o(2, imp) = y
                    xyz_o(3, imp) = z
                 else if(nameatom == ' C  ') then
                    iflag = iflag + 1 
                    xyz_c(1, imp) = x
                    xyz_c(2, imp) = y
                    xyz_c(3, imp) = z
                 else if(nameatom == ' CA ') then
                    xyz_ca(1, imp) = x
                    xyz_ca(2, imp) = y
                    xyz_ca(3, imp) = z
                 end if
                 
              end if
           end if
        end do
     end do

#ifdef MPI_PAR
  endif

  call MPI_Bcast(iflag2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( iflag2 .eq. -1 ) then
     return
  endif

  call MPI_Bcast(xyz_n, 3*MXMP,PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(xyz_o, 3*MXMP,PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(xyz_c, 3*MXMP,PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(xyz_ca,3*MXMP,PREC_MPI,0,MPI_COMM_WORLD,ierr)
#endif

  call protonate(nmp_all, xyz_n, xyz_ca, xyz_c, xyz_o, xyz_h)

  ! -----------------------------------------------------------------
  ! find the hbond
  q1 = 0.20e0_PREC
  q2 = 0.42e0_PREC
  f = 332.0e0_PREC
  
  ihbond = 0
  
  do imp = 1, nmp_all - 2
     do jmp = imp + 2, nmp_all
            
        iunit = imp2unit(imp)
        junit = imp2unit(jmp)

        if(iunit == junit) then

           ! ---------------------------------------------------------
           ! between C=O(imp) and N-H(jmp)
           rch = sqrt((xyz_c(1, imp) - xyz_h(1, jmp))**2 + &
                      (xyz_c(2, imp) - xyz_h(2, jmp))**2 + &
                      (xyz_c(3, imp) - xyz_h(3, jmp))**2)
           rcn = sqrt((xyz_c(1, imp) - xyz_n(1, jmp))**2 + &
                      (xyz_c(2, imp) - xyz_n(2, jmp))**2 + &
                      (xyz_c(3, imp) - xyz_n(3, jmp))**2)
           roh = sqrt((xyz_o(1, imp) - xyz_h(1, jmp))**2 + &
                      (xyz_o(2, imp) - xyz_h(2, jmp))**2 + &
                      (xyz_o(3, imp) - xyz_h(3, jmp))**2)
           ron = sqrt((xyz_o(1, imp) - xyz_n(1, jmp))**2 + &
                      (xyz_o(2, imp) - xyz_n(2, jmp))**2 + &
                      (xyz_o(3, imp) - xyz_n(3, jmp))**2)
           
           e = q1 * q2 * (1.0 / ron + 1.0 / rch - 1.0 / roh - 1.0 / rcn) * f
               
           if(e < - 0.5) then
              ihbond = ihbond + 1
              hbonde(ihbond) = e
              lhbond2mp(1, ihbond) = imp
              lhbond2mp(2, ihbond) = jmp       
           end if

           ! ---------------------------------------------------------
           ! between C=O(jmp) and N-H(imp)
           rch = sqrt((xyz_c(1, jmp) - xyz_h(1, imp))**2 + &
                      (xyz_c(2, jmp) - xyz_h(2, imp))**2 + &
                      (xyz_c(3, jmp) - xyz_h(3, imp))**2)
           rcn = sqrt((xyz_c(1, jmp) - xyz_n(1, imp))**2 + &
                      (xyz_c(2, jmp) - xyz_n(2, imp))**2 + &
                      (xyz_c(3, jmp) - xyz_n(3, imp))**2)
           roh = sqrt((xyz_o(1, jmp) - xyz_h(1, imp))**2 + &
                      (xyz_o(2, jmp) - xyz_h(2, imp))**2 + &
                      (xyz_o(3, jmp) - xyz_h(3, imp))**2)
           ron = sqrt((xyz_o(1, jmp) - xyz_n(1, imp))**2 + &
                      (xyz_o(2, jmp) - xyz_n(2, imp))**2 + &
                      (xyz_o(3, jmp) - xyz_n(3, imp))**2)
               
               
           e = q1 * q2 * (1.0e0_PREC / ron + 1.0e0_PREC / rch - 1.0 / roh - 1.0e0_PREC / rcn) * f
           if(e < - 0.5e0_PREC) then
                  
              ihbond = ihbond + 1
              hbonde(ihbond) = e
                  
              lhbond2mp(1, ihbond) = jmp
              lhbond2mp(2, ihbond) = imp
           end if
               
        end if

     end do
  end do
      
#ifdef MPI_PAR
! call MPI_Bcast(ihbond,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

  nhbond = ihbond

  ! -------------------------------------------------------------------
  do ihbond = 1, nhbond
     ihbond_stype(ihbond) = 0
  end do
  do imp = 1, nmp_all
     istype_mp(imp)= 0
  end do

  ! -------------------------------------------------------------------
  ! define the alpha helix
  ihbond = 1
  do
     if(ihbond >= nhbond) then
        exit
     else
        ico = lhbond2mp(1, ihbond)
        inh = lhbond2mp(2, ihbond)
        idev = inh - ico
     
        if(idev == 4) then
           
           lhbond = ihbond
           do jhbond = ihbond + 1, nhbond
              ico_2 = lhbond2mp(1, jhbond)
              inh_2 = lhbond2mp(2, jhbond)
              idev = inh_2 - ico_2
           
              if(idev == 4 .and. ico_2 == ico + 1 &
                   .and. inh_2 == inh + 1) then
                 ihbond_stype(ihbond) = 1              
                 ihbond_stype(jhbond) = 1
                 istype_mp(ico) = 1
                 istype_mp(inh) = 1
                 istype_mp(ico_2) = 1
                 istype_mp(inh_2) = 1
                 ico = ico_2
                 inh = inh_2
                 lhbond = jhbond
                 
                 ! -----------------------------------------------------------
                 !start of alpha helix is not helix
                 i = lhbond2mp(1, ihbond)
                 istype_mp(i) = 0
              end if
           end do
        
           inh_2 = lhbond2mp(2, lhbond)
           istype_mp(inh_2) = 0
           
           ihbond = lhbond + 1
           
        else
           ihbond = ihbond + 1
        end if
     end if
  end do

  ! -------------------------------------------------------------------
  ! define beta sheet 
  do ihbond = 1, nhbond - 1
     ico = lhbond2mp(1, ihbond)
     inh = lhbond2mp(2, ihbond)
     
     do jhbond = ihbond + 1, nhbond
        ico_2 = lhbond2mp(1, jhbond)
        inh_2 = lhbond2mp(2, jhbond)

        ! --------------------------------------------------------------
        ! anti parallel
        if(ico == inh_2 .and. inh == ico_2) then
           ihbond_stype(ihbond) = 2
           ihbond_stype(jhbond) = 2
           istype_mp(ico) = 2
           istype_mp(inh) = 2               
        end if
            
        if(inh_2-2 == ico .and. ico_2 + 2 == inh) then
           ihbond_stype(ihbond) = 2
           ihbond_stype(jhbond) = 2
           istype_mp(ico) = 2
           istype_mp(inh) = 2
           istype_mp(ico + 1) = 2
           istype_mp(inh - 1) = 2
           istype_mp(ico + 2) = 2
           istype_mp(inh - 2) = 2
        end if

        ! --------------------------------------------------------------
        ! parallel
        if(inh == ico_2 .and. ico + 2 == inh_2) then
           ihbond_stype(ihbond) = 2
           ihbond_stype(jhbond) = 2
           istype_mp(ico) = 3  
           istype_mp(inh) = 3
           istype_mp(ico_2) = 3
           istype_mp(inh_2) = 3
           istype_mp(ico + 1) = 3
        end if
            
     end do
  end do
         
  ! ------------------------------------------------------------------- 
  ! if the edge of sheet is parallel,the aminacid is not sheet 

  do jmp = 1, nmp_all - 1  
     if(istype_mp(jmp) == 3 .and. istype_mp(jmp + 1) == 0) then
        istype_mp(jmp) = 0
     elseif(istype_mp(jmp) == 0 .and. istype_mp(jmp + 1) == 3) then
        istype_mp(jmp+1) = 0
     end if
  end do

  ! ------------------------------------------------------------------- 
  ! define the sheet 2            
  do jmp = 1, nmp_all
     if(istype_mp(jmp) == 3) then
        istype_mp(jmp) = 2
     end if
  end do
      
contains

  ! ***********************************************************************
  ! The xyz coordinate of H atom from N on the main chain 
  ! are caluculated 
  ! ***********************************************************************
  subroutine  protonate(nmp_all, xyz_n, xyz_ca, xyz_c, xyz_o, xyz_h)
    
    use var_io, only : outfile
    implicit none

    ! -------------------------------------------------------------------
    integer, intent(in) :: nmp_all
    real(PREC), intent(in) :: xyz_n(3, MXMP), xyz_ca(3, MXMP)
    real(PREC), intent(in) :: xyz_c(3, MXMP), xyz_o(3, MXMP)
    ! to be determined
    real(PREC), intent(out) :: xyz_h(3, MXMP)

    ! local variables
    ! -------------------------------------------------------------------
    integer :: lunout
    integer :: imp
    real(PREC) :: bond, bondangle, dihed
    real(PREC) :: coorda(3), coordb(3), coordc(3), coordd(3)
    
    ! -------------------------------------------------------------------
    lunout = outfile%data

    bond = 1.0e0_PREC

    bondangle = 118.2e0_PREC * F_PI / 180.e0_PREC
    dihed = 0.0e0_PREC
    
    coordc(1:3) = xyz_n(1:3, 1)
    coordb(1:3) = xyz_ca(1:3, 1)
    coorda(1:3) = xyz_c(1:3, 1)
      
    call util_uintra2xyz(lunout, bond, bondangle, dihed, 0, &
         coorda, coordb, coordc, coordd)
      
    xyz_h(1:3, 1) = coordd(1:3)
    
    bond = 1.0e0_PREC
    bondangle = 119.5e0_PREC * F_PI / 180.0e0_PREC
    dihed = F_PI
    
    do imp = 2, nmp_all
       coordc(1:3) = xyz_n(1:3, imp)
       coordb(1:3) = xyz_c(1:3, imp - 1)
       coorda(1:3) = xyz_ca(1:3, imp)
       
       call util_uintra2xyz(lunout, bond, bondangle, dihed, 0, &
            coorda, coordb, coordc, coordd)
       
       xyz_h(1:3, imp) = coordd(1:3)
       
    end do

  end subroutine protonate
  
end subroutine setp_native_secstruct
