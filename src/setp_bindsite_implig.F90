! setp_bindsite_implig
!> @brief Read binding-site (residues) from input-file &
!>       (<<<< binding_site field) for implicit ligand simulation

subroutine setp_bindsite_implig()

  use const_maxsize
  use const_index
  use var_io,    only : infile, ius2unit, outfile
  use var_struct, only : nunit_all


  use var_implig, only : inimplig, naa_site_implig, inimplig_bindsite

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: isite
  integer :: iligsite_tmp(2)
  
  integer :: list_aa_implig_prepare(MXAA_IMPLIG,MXUNIT,MXSITE_IMPLIG)
  integer :: naa_implig_prepare(MXUNIT,MXSITE_IMPLIG)
  integer :: iaa, iaa_tmp
  integer :: iunit_implig
  integer :: k, ii, i_start,i_end
  integer :: i1 = 1
  integer :: i2 = 1

  integer :: isw
  integer :: isw_end, jj_start
  integer :: iu, iunit
  integer :: inunit(2), instate
  integer :: inu_tmp(2),ins_tmp
  integer :: iresi_tmp(2),iresi_state_tmp
  integer :: iresin

  integer :: luninp, lunout, ioutput
  integer :: iline, nlines, kprint

  character(4) :: kfind
  character(12) :: char12
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MSG_ERROR) :: error_message

  ! ----------------------------------------------------------------------
  luninp  = infile%inp
  lunout  = outfile%data
  ioutput = 0
! kprint  = 0
  kprint  = 1

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  
  ! ----------------------------------------------------------------------
  ! read data from imput_file (<<<< binding_site section) [IMPLIGSITE]
  ! ----------------------------------------------------------------------
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'binding_site    ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "<<<< binding_site" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     
     
     !! initialization of local variables{ naa_implig_prepare(iunit,isite), list_aa_implig_prepare(iaa, iunit, isite)}
     do isite = 1, MXSITE_IMPLIG
        do iunit = 1, nunit_all
           naa_implig_prepare(iunit,isite) = 0 
           !naa_implig_prepare(iunit,isite): max number of implig-binding site at (iunit,isite)
           do iaa = 1, MXAA_IMPLIG
              list_aa_implig_prepare(iaa, iunit, isite) = 0
              !list_aa_implig_prepare(iaa, iunit, isite): implig-binding site number for (iaa, iunit,isite) [iaa=1..naa_implig_prepare(iunit,isite)]
           enddo
        enddo
     enddo
     

     do iline = 1, nlines
        ctmp00 = cwkinp(iline)
        
        if(ctmp00(1:10) == 'IMPLIGSITE') then
           
           write (lunout, '(2a)') '---reading implicit ligand site: ', trim(ctmp00)
           !! read implicit ligand site number
           isw = 0
           isw_end = 0
           i1 = 0
           i2 = 0
           jj_start = 10
           do k = jj_start+1, CARRAY_MXCOLM
              if(isw_end == 0) then
                 if (cwkinp(iline)(k:k) /= ' ' .and. cwkinp(iline)(k:k) /= ',') then
                    if(isw == 0) then
                       i1 = k
                       isw = 1
                    end if
                    i2 = k
                 else if (isw == 1) then
                    read (ctmp00(i1:i2), *) char12
                    call util_unitstate(char12, inu_tmp, ins_tmp)
                    if(ins_tmp /= 0)then
                       error_message = 'Error: implicit ligand site number should be integer.'
                       call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    if(inu_tmp(2) > inimplig%nsite_implig) then
                       error_message = 'Error: implicit ligand site number should be less than nsite_implig.'
                       call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    iligsite_tmp(1)=inu_tmp(1) !! start_number for ligind site
                    iligsite_tmp(2)=inu_tmp(2) !! end_number   for ligind site

                    jj_start = i2
                    isw_end = 1
                 end if
              end if
           end do
           
           
           !! read unit&state for each implicit ligand site
           isw = 0
           isw_end = 0
           i1 = 0
           i2 = 0
           do k = jj_start+1, CARRAY_MXCOLM
              if(isw_end == 0) then
                 if (cwkinp(iline)(k:k) /= ' ' .and. cwkinp(iline)(k:k) /= ',') then
                    if(isw == 0) then
                       i1 = k
                       isw = 1
                    end if
                    i2 = k
                 else if (isw == 1) then
                    read (ctmp00(i1:i2), *) char12
                    call util_unitstate(char12, inu_tmp, ins_tmp)
                    inunit(1) = inu_tmp(1) !! start unit number for ligind binding
                    inunit(2) = inu_tmp(2) !! end unit number for ligind binding
                    instate   = ins_tmp    !! state of unit for ligind binding

                    jj_start = i2
                    isw_end = 1
                 end if
              end if
           end do
           

           !! read the sign ("u") for specification of intra-unit mass-point-id 
           isw = 0
           isw_end = 0
           i1 = 0
           i2 = 0
           do k = jj_start+1, CARRAY_MXCOLM
              if(isw_end == 0) then
                 if (cwkinp(iline)(k:k) /= ' ' .and. cwkinp(iline)(k:k) /= ',') then
                    if(isw == 0) then
                       i1 = k
                       isw = 1
                    end if
                    i2 = k
                 else if (isw == 1) then
                    if(ctmp00(i1:i2) == 'u') then
                       !! specification is OK
                    else
                       error_message = 'Error: the specification "u" which stands'// &
                       ' for intra-unit mass-point-id (residues-id) should be necessary.'
                       call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    
                    jj_start = i2
                    isw_end = 1
                 end if
              end if
           end do
           


           !! read the implicit-ligand residues for each binding-site
           isw = 0
           i1 = 0
           i2 = 0
           
           do k = jj_start + 1, CARRAY_MXCOLM
              if (cwkinp(iline)(k:k) /= ' ' .and. cwkinp(iline)(k:k) /= ',') then
                 if(isw == 0) then
                    i1 = k
                    isw = 1
                 end if
                 i2 = k
              else if (isw == 1) then
                 !read (ctmp00(i1:i2), *) itype
                 read (ctmp00(i1:i2), *) char12
                 call util_unitstate(char12, iresi_tmp, iresi_state_tmp)
                 if(iresi_state_tmp /= 0)then
                    error_message = 'Error: residue number should be integer.'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif
                 
                 do isite = iligsite_tmp(1), iligsite_tmp(2)
                    if (isite > MXSITE_IMPLIG) then
                       error_message = 'Error: MXSITE_IMPLIG is too small.'
                       call util_error(ERROR%STOP_ALL, error_message)
                    endif
                    do iu = inunit(1), inunit(2)
                       iunit = ius2unit(iu, instate)
                       if (iunit < 1 .or. iunit > nunit_all)  then
                          error_message = 'Error: unit&state specification for implicit-ligand is incorrect.'
                          call util_error(ERROR%STOP_ALL, error_message)
                       endif
                       do iresin = iresi_tmp(1), iresi_tmp(2)
                          naa_implig_prepare(iunit, isite) = naa_implig_prepare(iunit, isite) + 1
                          list_aa_implig_prepare(naa_implig_prepare(iunit, isite), iunit, isite) = iresin
                       enddo
                    enddo
                 enddo
                 isw = 0
              end if
           end do
        end if
     end do
     
     
     !check the order of list_aa_implig_prepare(iaa_tmp, iunit, isite)
     do isite = 1, inimplig%nsite_implig
        do iunit = 1, nunit_all
           do iaa_tmp = 1, naa_implig_prepare(iunit, isite)-1
              if(list_aa_implig_prepare(iaa_tmp, iunit, isite) >= list_aa_implig_prepare(iaa_tmp + 1, iunit, isite)) then
                 error_message = &
                      'Error: at "binding_site" field, residue number for each ligand site should be arranged in ascending orders.'
                 call util_error(ERROR%STOP_ALL, error_message)
              endif
           end do
        enddo
     end do
  

  
!**********************************************************************
! initialization for following global-parameters 
!                    inimplig_bindsite%nunit_implig(isite), 
!                    inimplig_bindsite%idunit_implig(iunit_implig, isite), 
!                    inimplig_bindsite%list_aa_implig(iaa ,isite)
!                    inimplig_bindsite%naa_implig(iunit_implig, isite)
!                    naa_site_implig(isite)
 do isite = 1, MXSITE_IMPLIG
     inimplig_bindsite%nunit_implig(isite)=0
     naa_site_implig(isite) = 0
     do iunit_implig = 1, MXUNIT_IMPLIG
        inimplig_bindsite%idunit_implig(iunit_implig, isite) = 0
        inimplig_bindsite%naa_implig(iunit_implig, isite) = 0
     enddo
     do iaa =1, MXAA_IMPLIG
        inimplig_bindsite%list_aa_implig(iaa ,isite) = 0 
     enddo
  enddo


!***********************************************************************
! data conversion 
! <<from>> 
!   naa_implig_prepare(iunit, ligand_isite) !! number of ligand-binding residue for (iunit&state, ligand_isite)
!   list_aa_implig_prepare(i, iunit, ligand_isite) !! ligand-binding residue number for (i, iunit&state, ligand_isite) [i=1, naa_implig_prepare(iunit, ligand_isite)] 
! <<to>> 
!   inimplig_bindsite%nunit_implig(ligand_isite) !! number of unit&state (chain) for ligand_isite
!   naa_site_implig(ligand_isite) !! number of ligand-binding residue for ligand_isite
!   inimplig_bindsite%idunit_implig(iunit_implig, ligand_isite) !! unit&state (chain) number for (i, ligand_isite)   [i=1, inimplig_bindsite%nunit_implig(ligand_isite)] 
!   inimplig_bindsite%naa_implig(iunit_implig, isite)    !! number of ligand-binding residue for (i, ligand_isite)   [i=1, inimplig_bindsite%nunit_implig(ligand_isite)] 
!   inimplig_bindsite%list_aa_implig(iaa ,ligand_isite) !! ligand-binding residue number for (i, ligand_isite) [i=1, naa_site_implig(ligand_isite)]
!*************************************************************************
  do isite = 1, inimplig%nsite_implig
     iaa = 0
     do iunit = 1, nunit_all
        if(naa_implig_prepare(iunit, isite) > 0) then
           inimplig_bindsite%nunit_implig(isite) = inimplig_bindsite%nunit_implig(isite) + 1
           inimplig_bindsite%idunit_implig(inimplig_bindsite%nunit_implig(isite), isite) = iunit
           inimplig_bindsite%naa_implig(inimplig_bindsite%nunit_implig(isite), isite) = naa_implig_prepare(iunit, isite)
           do iaa_tmp = 1, naa_implig_prepare(iunit, isite)
              iaa = iaa + 1
              inimplig_bindsite%list_aa_implig(iaa, isite) = list_aa_implig_prepare(iaa_tmp, iunit, isite)
              naa_site_implig(isite) = iaa
           enddo
        endif
     enddo
  enddo



  !------------------------------------
  !output to data_file
  !------------------------------------
  do isite=1,inimplig%nsite_implig
     write(lunout, '(a, i3, a ,i4)') 'ligand_',isite, '    total number of binding_site = ',naa_site_implig(isite)

     iaa = 0
     do iunit=1, inimplig_bindsite%nunit_implig(isite)
        write(lunout,'(a, i3, a, i3)') '<< ligand_',isite,',  unit&state =', &
             inimplig_bindsite%idunit_implig(iunit, isite)
        
        write(lunout,'(a, i4)') 'number of binding_site =',inimplig_bindsite%naa_implig(iunit,isite)
        write(lunout,'(a)') 'list of binding_site:'
        i_start = iaa + 1
        i_end   = iaa + inimplig_bindsite%naa_implig(iunit,isite)
        write(lunout, *) (inimplig_bindsite%list_aa_implig(ii,isite),ii=i_start,i_end)
        iaa = iaa + inimplig_bindsite%naa_implig(iunit,isite)
        write(lunout,'(a)') '>>'

        ! for debug
        !write(lunout,'(a, i3, a, i2, a, i4)') 'number of binding residues for [unit&state=',&
        !    inimplig_bindsite%idunit_implig(iunit,isite),', ligand(',isite,') ] =',inimplig_bindsite%naa_implig(iunit,isite)
        !        write(lunout,'(a, i3, a, i2, a, i4)')'number of binding residues for [unit&state=',&
        !             inimplig_bindsite%idunit_implig(iunit,isite),', ligand(',isite,') ] =',inimplig_bindsite%naa_implig(iunit,isite)
        !write(lunout,'(a, i3, a, i3)') '>> ligand_',isite,',  unit&state =', &
        !     inimplig_bindsite%idunit_implig(iunit, isite)

     enddo
     write(lunout,'(a)') ''
  enddo


#ifdef MPI_PAR
   end if
   call MPI_Bcast (inimplig_bindsite, inimplig_bindsite%sz, MPI_BYTE,   0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast (naa_site_implig,   MXSITE_IMPLIG,        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

!  return
end subroutine setp_bindsite_implig
