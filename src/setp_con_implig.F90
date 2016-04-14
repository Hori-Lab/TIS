!setp_con_implig
!> @brief  Fixes the implicit-ligand mediated contact pair.

subroutine setp_con_implig(xyz_mp_init, &
                           iatomnum, xyz) 

  use const_maxsize
  use const_physical
  use var_implig, only : inimplig, naa_site_implig, p_implig,&
                         ncon_implig, icon2mp_implig, vdwrad_implig, inimplig_bindsite
  use var_struct, only : lunit2mp, icon2mp, ncon, imp2unit
  use var_inp,    only : outfile


#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none
  !-----------------------------------------------------------------
  real(PREC), intent(in) :: xyz_mp_init(SPACE_DIM, MXMP)
  integer,    intent(in) :: iatomnum(MXMP)
  real(PREC), intent(in) :: xyz(SPACE_DIM, MXATOM_MP, MXMP)

  
  !-----------------------------------------------------------------
  ! local variables
  integer :: i, ist, ich, iaa, iaa_ini, iaa_las, jaa, inum, jnum, nconsum
  integer :: icon, lunout
  integer :: isep_contact
  integer :: n_total_ncon
  real(PREC) :: dfcontact2_implig, dist2
  integer :: imp1, imp2
  integer :: iunit1, iunit2
  integer :: imp1un, imp2un
  !!  real(PREC) :: dfcontact2_go, dfcontact2_implig, dist2
  !!  integer, intent(in) :: iatomnum(MXMP)
  !!  real(PREC), intent(in) :: xyz(3, MXMP, 25)
  ! -----------------------------------------------------------
  lunout = outfile%data

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  
     dfcontact2_implig = inimplig%dfcontact_implig**2
     !!dfcontact_implig = 10.0e0_PREC

     isep_contact = inimplig%isep_contact_implig
     !!isep_contact = inpara%n_sep_contact
     !!isep_contact = 4
     


!!     react_rate_implig(i) = 4.0e0_PREC*inimplig%diff_const_implig(i)*inimplig%react_rad_implig(i) ! k_D
!!     react_rate_implig(i): reaction rate for implicit-ligand(i)     

  
  !------------------------------------------------------------------   
  ! calculate p_implig(i)
  ! p_implig(i): binding probability for implicit-ligand(i) per ``istep_implig" MD time steps
  ! bind_rate_implig(i): the ligand-binding rate constant [k_b]  for implicit-ligand(i)   
  do i=1,inimplig%nsite_implig
     p_implig(i) = inimplig%bind_rate_implig(i)*inimplig%istep_implig

!!   p_implig(i) = react_rate_implig(i)*inimplig%c_implig(i)*inimplig%istep_implig
     ! for debug
     !     write(*,*)'react_rate_implig=',react_rate_implig(i)
     !     write(*,*)'p_implig=',p_implig(i)
     !     write(*,*)'diff_const_implig=',inimplig%diff_const_implig(i)
     !     write(*,*)'react_rad_implig=',inimplig%react_rad_implig(i)
     !     write(*,*)'c_implig=',inimplig%c_implig(i)
     !     write(*,*)'istep_implig=',inimplig%istep_implig
  enddo


  
  ! fix the list of binding residues
  do ist=1,inimplig%nsite_implig
     iaa_ini = 0
     iaa_las = 0
     do ich=1,inimplig_bindsite%nunit_implig(ist)
        iaa_ini = iaa_las + 1
        iaa_las = iaa_ini + inimplig_bindsite%naa_implig(ich,ist) - 1
        do iaa=iaa_ini,iaa_las
           inimplig_bindsite%list_aa_implig(iaa,ist) =&
               inimplig_bindsite%list_aa_implig(iaa,ist) + &
               lunit2mp(1,inimplig_bindsite%idunit_implig(ich,ist)) - 1
        enddo
     enddo
  enddo
  ! for debug
  !  do ist=1,inimplig%nsite_implig
  !     write(*,*) ist
  !     write(*,*) (inimplig_bindsite%list_aa_implig(iaa,ist),iaa=1,naa_site_implig(ist))
  !  enddo
  


  ! -----------------------------------------------------------
  ! define the impligand-mediated contacts
  nconsum = 0
  do ist=1,inimplig%nsite_implig
     ncon_implig(ist) = 0
     do iaa=1,naa_site_implig(ist)-1
        jaa_loop:do jaa=iaa+1,naa_site_implig(ist)
           if(inimplig_bindsite%list_aa_implig(jaa,ist)- &
                inimplig_bindsite%list_aa_implig(iaa,ist) >= isep_contact)then
              do icon=1,ncon !! excluding original contacts (native contact)
                 if(icon2mp(1,icon)==inimplig_bindsite%list_aa_implig(iaa,ist) .and. &
                      icon2mp(2,icon)==inimplig_bindsite%list_aa_implig(jaa,ist))then
                    cycle jaa_loop
                 endif
              enddo
              do inum=1,iatomnum(inimplig_bindsite%list_aa_implig(iaa,ist))
                 if (xyz(1,inum,inimplig_bindsite%list_aa_implig(iaa,ist)) > INVALID_JUDGE) cycle
                 do jnum=1,iatomnum(inimplig_bindsite%list_aa_implig(jaa,ist))
                    if (xyz(1,jnum,inimplig_bindsite%list_aa_implig(jaa,ist)) > INVALID_JUDGE) cycle
                    dist2 = (xyz(1,jnum,inimplig_bindsite%list_aa_implig(jaa,ist)) - &
                         xyz(1,inum,inimplig_bindsite%list_aa_implig(iaa,ist)))**2 + &
                         (xyz(2,jnum,inimplig_bindsite%list_aa_implig(jaa,ist)) - &
                         xyz(2,inum,inimplig_bindsite%list_aa_implig(iaa,ist)))**2 + &
                         (xyz(3,jnum,inimplig_bindsite%list_aa_implig(jaa,ist)) - &
                         xyz(3,inum,inimplig_bindsite%list_aa_implig(iaa,ist)))**2
                    if(dist2<=dfcontact2_implig)then
                       ncon_implig(ist) = ncon_implig(ist) + 1
                       nconsum = nconsum + 1
                       icon2mp_implig(1,nconsum) = inimplig_bindsite%list_aa_implig(iaa,ist)
                       icon2mp_implig(2,nconsum) = inimplig_bindsite%list_aa_implig(jaa,ist)
                       vdwrad_implig(nconsum) = &
                         sqrt( ( xyz_mp_init(1,inimplig_bindsite%list_aa_implig(jaa,ist))       &
                                -xyz_mp_init(1,inimplig_bindsite%list_aa_implig(iaa,ist)) )**2  &
                              +( xyz_mp_init(2,inimplig_bindsite%list_aa_implig(jaa,ist))       &
                                -xyz_mp_init(2,inimplig_bindsite%list_aa_implig(iaa,ist)) )**2  &
                              +( xyz_mp_init(3,inimplig_bindsite%list_aa_implig(jaa,ist))       &
                                -xyz_mp_init(3,inimplig_bindsite%list_aa_implig(iaa,ist)) )**2 )
                       cycle jaa_loop
                    endif
                 enddo
              enddo
           endif
        enddo jaa_loop
     enddo
  enddo

  
  !------------
  ! output
  !--------------------------------------------------------------------
  ! write the ligand-mediated contact information to data_file
  !-------------------------------------------------------------------
  write (lunout, '(a)') '************************************************************************'
  write (lunout, '(a)') '<<<< ligand-mediated contacts for implicit ligand model'
  ! write (lunout, '(a)') '<<<< native contact '

  ! calculate total number of ligand-mediated contact
  n_total_ncon = 0
  do ist=1,inimplig%nsite_implig
     n_total_ncon = n_total_ncon + ncon_implig(ist)
  end do
  write (lunout, '(a, i6)') '** total_mediated_contact = ', n_total_ncon
  write (lunout, '(a, f10.2, a)') '** definition_of_ligand-mediated_contact (constant) = ', sqrt(dfcontact2_implig), ' A'
  write (lunout, '(a)') ''
  
  nconsum = 0
  do ist=1,inimplig%nsite_implig
     write (lunout, '(a, i3)') '** ligand-mediated contact for ligand_', ist
     write (lunout, '(a, i6)') '** total_mediated_contact = ', ncon_implig(ist)
     write (lunout, '(a)', ADVANCE='NO') '**        icon  ligand  iunit1-iunit2 '
     write (lunout, '(a)', ADVANCE='NO') 'imp1 - imp2 imp1un-imp2un'
     write (lunout, '(a)') '    mediated_contact_length'
     do i=1,ncon_implig(ist)
        nconsum = nconsum + 1
        imp1 = icon2mp_implig(1,nconsum)
        imp2 = icon2mp_implig(2,nconsum)
        iunit1 = imp2unit(imp1)
        iunit2 = imp2unit(imp2)
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit2) + 1
        write (lunout, "(a7, 8(1xi6), 1(f12.4))") &
             'contact', nconsum, ist,iunit1, iunit2, imp1, imp2, &
             imp1un, imp2un, &
             vdwrad_implig(nconsum)
     enddo
     write (lunout, '(a)') ''
  end do
  write (lunout, '(a4)') '>>>>'
  write (lunout, '(a)') '************************************************************************'
  
  ! for debug
  ! ----------------------------------------------------------------------
  ! output to data_file (ligand-mediated contact information)
  ! ----------------------------------------------------------------------
  !  write (lunout, '(a)') '****************************************************'
  !  write (lunout, '(a)') '<<<< ligand-mediated contacts for implicit ligand model'
  !  ! calculate total number of ligand-mediated contact
  !  n_total_ncon = 0
  !  do ist=1,inimplig%nsite_implig
  !     n_total_ncon = n_total_ncon + ncon_implig(ist)
  !  end do
  !
  !  nconsum = 0
  !  do ist=1,inimplig%nsite_implig
  !     write (lunout, '(a, i6, a, i6)') '# of BINDING SITE ', ist, ' CONTACTS=', ncon_implig(ist)
  !     do i=1,ncon_implig(ist)
  !        nconsum = nconsum + 1
  !        write (lunout, '(a, i6, a, i6, a, i6, a, i6)') 'BINDING SITE ', ist, ' CONTACT', i, ':', &
  !            icon2mp_implig(1,nconsum), '-', icon2mp_implig(2,nconsum)
  !     enddo
  !  enddo
  !  write (lunout, '(a)') '>>>>'
  !  write (lunout, '(a)') '****************************************************'

#ifdef MPI_PAR
  endif
  call MPI_Bcast(inimplig_bindsite, inimplig_bindsite%sz, MPI_BYTE,   0,MPI_COMM_WORLD,ierr)  
  call MPI_Bcast(ncon_implig,       MXSITE_IMPLIG,        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(icon2mp_implig,    2*MXCON_IMPLIG,       MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!call MPI_Bcast(react_rate_implig, MXSITE_IMPLIG,        PREC_MPI,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(p_implig,          MXSITE_IMPLIG,        PREC_MPI,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(vdwrad_implig,     MXCON_IMPLIG,         PREC_MPI,   0,MPI_COMM_WORLD,ierr)
#endif

!  return
end subroutine setp_con_implig
