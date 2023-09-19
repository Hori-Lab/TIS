! mloop_neigh_dist
!> @brief calculate rneighbordist2_unit when i_neigh_dynamic == 1

! ************************************************************************
subroutine mloop_neigh_dist()

  use const_maxsize
  use const_index
  use var_io,     only : outfile
  use var_setp,   only : insimu, inmisc, indtrna13, indtrna15, &! inrna, &
                         inpro, inligand, inexv
  use var_struct, only : nunit_all, iclass_unit, exv_radius_mp, lunit2mp, imp2type
                          
  use mpiconst

  implicit none

  integer :: i, iunit, junit, imp, jmp, iclass, jclass
  integer :: lunout
  real(PREC) :: cut, max_cut
  character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
  write(*,*) '##### start mloop_neigh_dist'
#endif
  ! ----------------------------------------------------------------------

  do iunit = 1, nunit_all
     iclass = iclass_unit(iunit)

     do junit = 1, nunit_all
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%NOTHING)) then
           cycle
        endif

        jclass = iclass_unit(junit)
        max_cut = 0.0

        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV_WCA)) then
           if (indtrna13%exv_dist > max_cut) then
              max_cut = indtrna13%exv_dist
           endif
        endif

        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV_DT15)) then
           if (iclass == CLASS%RNA .AND. jclass == CLASS%RNA) then
              cut  = indtrna15%exv_dist
              if (cut > max_cut) then
                 max_cut = cut
              endif
              cut  = indtrna15%exv_dist_PS
              if (cut > max_cut) then
                 max_cut = cut
              endif
           else
              do imp = lunit2mp(1,iunit), lunit2mp(2, iunit)
                 do jmp = lunit2mp(1,junit), lunit2mp(2, junit)
                    cut  = exv_radius_mp(imp)  + exv_radius_mp(jmp)
                    if (cut > max_cut) then
                       max_cut = cut
                    endif
                 enddo
              enddo
           endif
        endif

        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV_GAUSS)) then
           if (inpro%cutoff_exv_gauss > max_cut) then
              max_cut = inpro%cutoff_exv_gauss
           endif
        endif

!        if (inmisc%i_residuenergy_radii == 0) then
!           if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV12)) then
!              !if (iclass == CLASS%RNA .AND. jclass == CLASS%RNA) then
!              !   cut = inrna%cutoff_exvol * inrna%cdist_rep12
!              !else if ((iclass == CLASS%RNA .AND. jclass == CLASS%PRO) .OR. &
!              !         (iclass == CLASS%PRO .AND. jclass == CLASS%RNA)) then
!              !   cut = (inrna%cutoff_exvol*inrna%cdist_rep12 + inpro%cutoff_exvol*inpro%cdist_rep12) / 2.0
!              !else if (iclass == CLASS%LIG .OR. jclass == CLASS%LIG) then
!              if (iclass == CLASS%LIG .OR. jclass == CLASS%LIG) then
!                 cut = inligand%cutoff_exvol * inligand%cdist_rep12_llig
!                 if(iclass == CLASS%PRO .OR. jclass == CLASS%PRO) then
!                    cut = inligand%cutoff_exvol * inligand%cdist_rep12_lpro
!                 end if
!              else
!                 cut = inpro%cutoff_exvol * inpro%cdist_rep12
!              endif
!   
!              if (cut > max_cut) then
!                 max_cut = cut
!              endif
!           endif
!   
!           if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV6)) then
!              cut = inpro%cutoff_exvol * inpro%cdist_rep6
!              if (cut > max_cut) then
!                 max_cut = cut
!              endif
!           end if
!
!        else if (inmisc%i_residuenergy_radii == 1) then

           if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV12)) then
              do imp = lunit2mp(1,iunit), lunit2mp(2, iunit)
                 do jmp = lunit2mp(1,junit), lunit2mp(2, junit)
                    cut = inexv%exv12_cutoff * (exv_radius_mp(imp) + exv_radius_mp(jmp))
                    if (cut > max_cut) then
                       max_cut = cut
                    endif
                 enddo
              enddo
           endif

           if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV6)) then
              do imp = lunit2mp(1,iunit), lunit2mp(2, iunit)
                 do jmp = lunit2mp(1,junit), lunit2mp(2, junit)

                    if (iclass_unit(iunit) == CLASS%SOPSC .and. iclass_unit(junit) == CLASS%SOPSC) then
                       ! Backbone - Backbone
                       if (imp2type(imp) == MPTYPE%SOPBB .and. imp2type(jmp) == MPTYPE%SOPBB) then
                          cut = inexv%exv6_cutoff * (inexv%exv_rad_sopsc_BB_BB * 2)
                       ! Backbone - Sidechain
                       else if (imp2type(imp) == MPTYPE%SOPBB .and. imp2type(jmp) == MPTYPE%SOPSC) then
                          cut = inexv%exv6_cutoff * (inexv%exv_rad_sopsc_BB_SC + exv_radius_mp(jmp))
                       ! Sidechain - Backbone
                       else if (imp2type(imp) == MPTYPE%SOPSC .and. imp2type(jmp) == MPTYPE%SOPBB) then
                          cut = inexv%exv6_cutoff * (exv_radius_mp(imp) + inexv%exv_rad_sopsc_BB_SC)
                       ! Sidechain - Sidechain
                       else if (imp2type(imp) == MPTYPE%SOPSC .and. imp2type(jmp) == MPTYPE%SOPSC) then
                          cut = inexv%exv6_cutoff * (exv_radius_mp(imp) + exv_radius_mp(jmp))
                       else 
                          error_message = 'Error: logical defect in mloop_neigh_dist, imp2type of SOPSC'
                          call util_error(ERROR%STOP_ALL, error_message)
                       endif
                    else
                       cut = inexv%exv6_cutoff * (exv_radius_mp(imp) + exv_radius_mp(jmp))
                    endif

                    if (cut > max_cut) then
                       max_cut = cut
                    endif
                 enddo
              enddo
           endif

!        endif

        inmisc%rneighbordist2_unit(iunit, junit) = (max_cut + insimu%neigh_margin) ** 2
     enddo
  enddo


#ifdef MPI_PAR
  if (myrank == 0) then
#endif
  lunout = outfile%data
  ! ----------------------------------------------------------------------
  write(lunout, '(a)') '<neighbordistance> (=cutoff + margin)'
  write(lunout, '(a)') '--------------------------------------------'
  write(lunout, '(a6, a2, 200i6)') 'unit', '|', (i, i = 1, nunit_all)
  write(lunout, '(a)') '--------------------------------------------'

  do junit = 1, nunit_all
     write (lunout, '(i6, a2, 200f6.1)') &
          junit, '|', (sqrt(inmisc%rneighbordist2_unit(iunit, junit)), iunit = 1, junit)
  end do

  write (lunout, '(a)') '--------------------------------------------'
  write (lunout, '(a)')
  write (lunout, '(72(1H*))') 

#ifdef MPI_PAR
  end if
  !call MPI_Bcast(inmisc%rneighbordist2_unit, MXUNIT*MXUNIT, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
#endif

#ifdef _DEBUG
  write(*,*) '##### end mloop_neigh_dist'
#endif

end subroutine mloop_neigh_dist
