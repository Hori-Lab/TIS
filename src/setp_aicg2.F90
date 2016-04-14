! setp_aicg2
!> @brief Reset parameters of local and non-local interactions for &
!>        AICG2 simulation

subroutine setp_aicg2(iatomnum, xyz, cname_ha) 

  use const_maxsize
  use const_physical
  use const_index

  use var_struct, only : nbd, nba, ndih, ncon,  &
                         factor_bd, factor_aicg13, factor_aicg14, factor_go, &
                         coef_bd, coef_ba, coef_dih, &
                         coef_aicg13_gauss, coef_aicg14_gauss, coef_dih_gauss, coef_go, &
                         wid_aicg13_gauss, wid_aicg14_gauss,wid_dih_gauss, &
                         ibd2mp, iba2mp, idih2mp, &
                         icon2unit, imp2unit
  use var_inp,  only : outfile
  use var_setp, only : inaicg2, inmisc
#ifdef MPI_PAR
  use mpiconst
  use var_struct, only : nmp_all
#endif
  implicit none    
  real(PREC), intent(in) :: xyz(SPACE_DIM, MXATOM_MP, MXMP)
  integer,    intent(in) :: iatomnum(MXMP)
  character(4), intent(in) :: cname_ha(MXATOM_MP, MXMP)


  ! -----------------------------------------------------------------
  ! local variable
  integer :: ibd, iba, idih, icon, mba, mdih, mcon, iunit, junit
  integer :: ipara
  integer :: ncon_atm(MXCON,16), ncon_atm_13(MXMP,16), ncon_atm_14(MXMP,16)
  real(PREC) :: para_ms(17)
  real(PREC) :: e_ground, e_ground_local, e_ground_13, e_ground_14
  real(PREC) :: e_con(MXCON), e_loc_13(MXMP), e_loc_14(MXMP)
  integer :: lunout
!-------------------------------------------------------------------  
#ifdef MPI_PAR
  if (myrank == 0) then
#endif
  ncon_atm_13 = 0
  ncon_atm_14 = 0
  ncon_atm = 0
  e_con = 0.
  e_loc_13 = 0.
  e_loc_14 = 0.
  lunout = outfile%data
!------------------------------------------------------
! parameters from fitting the all atom contact energies
      para_ms(1) = -1.4247   !B-B hydrogen bonds
      para_ms(2) = -0.4921   !B-B donor-accetor contacts
      para_ms(3) = -0.2404   !B-B carbon-X contacts
      para_ms(4) = -0.1035   !B-B other
      para_ms(5) = -5.7267   !S-S hydrogen bonds
      para_ms(6) = -12.4878  !S-S salty bridge
      para_ms(7) = -0.0308   !S-S donor-accetor contacts
      para_ms(8) = -0.1113   !S-S carbon-X contacts
      para_ms(9) = -0.2168   !S-S charge-X contacts
      para_ms(10) = 0.2306   !S-S other
      para_ms(11) = -3.4819  !S-B hydrogen bonds
      para_ms(12) = -0.1809  !S-B donor-accetor contacts
      para_ms(13) = -0.1209  !S-B carbon-X contacts
      para_ms(14) = -0.2984  !S-B charge-X contacts
      para_ms(15) = -0.0487  !S-B other
      para_ms(16) = -0.0395  !long range contacts
      para_ms(17) = -0.1051  !offset

!------------------------------------------
!   counting the atomic contacts
         call contact_make_aicg2(iatomnum, xyz, cname_ha, ncon_atm_13, ncon_atm_14, ncon_atm)

!------------------------------------------
!  bond
        do ibd = 1, nbd
           iunit = imp2unit(ibd2mp(1, ibd))
           if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .OR. &
              inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then
              coef_bd(1,ibd) = factor_bd(ibd) * inaicg2%cbd_aicg2
           end if
        enddo
!-------------------
! local 13
         e_ground_local = 0.
         e_ground_13 = 0.
         mba = 0
         do iba = 1, nba
            iunit = imp2unit(iba2mp(1, iba))
            if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .OR. &
               inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then
              e_loc_13(iba) = para_ms(17)
              do ipara = 1, 16
                 e_loc_13(iba) = e_loc_13(iba) + para_ms(ipara) * ncon_atm_13(iba, ipara)
              end do
              if(e_loc_13(iba) >= -inaicg2%ecut_low_aicg2) e_loc_13(iba) = -inaicg2%ecut_low_aicg2
              if(e_loc_13(iba) <= -inaicg2%ecut_up_aicg2) e_loc_13(iba) = -inaicg2%ecut_up_aicg2
              e_ground_local = e_ground_local + e_loc_13(iba)
              e_ground_13 = e_ground_13 + e_loc_13(iba)
              mba = mba + 1
            end if
         end do
!-------------------
! local 14
         e_ground_14 = 0.
         mdih = 0
         do idih = 1, ndih
            iunit = imp2unit(idih2mp(1, idih))
            if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .OR. &
               inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then
              e_loc_14(idih) = para_ms(17)
              do ipara = 1, 16
                 e_loc_14(idih) = e_loc_14(idih) + para_ms(ipara) * ncon_atm_14(idih, ipara)
              end do
              if(e_loc_14(idih) >= -inaicg2%ecut_low_aicg2) e_loc_14(idih) = -inaicg2%ecut_low_aicg2
              if(e_loc_14(idih) <= -inaicg2%ecut_up_aicg2) e_loc_14(idih) = -inaicg2%ecut_up_aicg2
              e_ground_local = e_ground_local + e_loc_14(idih)
              e_ground_14 = e_ground_14 + e_loc_14(idih)
              mdih = mdih + 1
            end if
         end do
!-------------------
         e_ground_local = e_ground_local / (mba + mdih)
         e_ground_13 = e_ground_13 / mba
         e_ground_14 = e_ground_14 / mdih
!-------------------
!aicg13
         do iba = 1, nba
            iunit = imp2unit(iba2mp(1, iba))
            if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2))then
              if(inaicg2%iflag_scale_aicg2 == 0)then
                 coef_aicg13_gauss(iba) = inaicg2%ave_caicg2_loc * e_loc_13(iba) / e_ground_local
              else
                 coef_aicg13_gauss(iba) = -inaicg2%gen_caicg2_loc * e_loc_13(iba)
              end if
              coef_aicg13_gauss(iba) = factor_aicg13(iba) * coef_aicg13_gauss(iba)
              wid_aicg13_gauss(iba) = inaicg2%wid_aicg13
            else if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then
              if(inaicg2%iflag_scale_aicg2 == 0)then
                 coef_aicg13_gauss(iba) = inaicg2%ave_caicg2plus_13 * e_loc_13(iba) / e_ground_13
              else
                 coef_aicg13_gauss(iba) = -inaicg2%gen_caicg2plus_13 * e_loc_13(iba)
              end if
              coef_aicg13_gauss(iba) = factor_aicg13(iba) * coef_aicg13_gauss(iba)
              wid_aicg13_gauss(iba) = inaicg2%wid_aicg13
            end if
         end do
!aicg14
         do idih = 1, ndih
            iunit = imp2unit(idih2mp(1, idih))
            if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2))then 
              if(inaicg2%iflag_scale_aicg2 == 0)then
                 coef_aicg14_gauss(idih) = inaicg2%ave_caicg2_loc * e_loc_14(idih) / e_ground_local
              else
                 coef_aicg14_gauss(idih) = -inaicg2%gen_caicg2_loc * e_loc_14(idih)
              end if
              coef_aicg14_gauss(idih) = factor_aicg14(idih) * coef_aicg14_gauss(idih)
              wid_aicg14_gauss(idih) = inaicg2%wid_aicg14
            else if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then
              if(inaicg2%iflag_scale_aicg2 == 0)then
                 coef_dih_gauss(idih) = inaicg2%ave_caicg2plus_14 * e_loc_14(idih) / e_ground_14
              else
                 coef_dih_gauss(idih) = -inaicg2%gen_caicg2plus_14 * e_loc_14(idih)
              end if
              coef_dih_gauss(idih) = factor_aicg14(idih) * coef_dih_gauss(idih)
              wid_dih_gauss(idih) = inaicg2%wid_dih
            end if
         end do
!-------------------
! nlocal
         e_ground = 0.
         mcon = 0
         do icon = 1, ncon
            iunit = icon2unit(1, icon)
            junit = icon2unit(2, icon)
            if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG2))then
              e_con(icon) = para_ms(17)
              do ipara = 1, 16
                 e_con(icon) = e_con(icon) + para_ms(ipara) * ncon_atm(icon, ipara)
              end do
              if(e_con(icon) >= -inaicg2%ecut_low_aicg2) e_con(icon) = -inaicg2%ecut_low_aicg2
              if(e_con(icon) <= -inaicg2%ecut_up_aicg2) e_con(icon) = -inaicg2%ecut_up_aicg2
              e_ground = e_ground + e_con(icon)
              mcon = mcon + 1
            end if
         end do
         e_ground = e_ground / mcon

         do icon = 1, ncon
            iunit = icon2unit(1, icon)
            junit = icon2unit(2, icon)
            if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG2))then
            if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS) .or. &
               inmisc%flag_local_unit(junit, junit, LINTERACT%L_AICG2_PLUS)) then

              if(inaicg2%iflag_scale_aicg2 == 0)then
                 coef_go(icon) = inaicg2%ave_caicg2plus_nloc * e_con(icon) / e_ground
              else
                 coef_go(icon) = -inaicg2%gen_caicg2plus_nloc * e_con(icon)
              end if
              coef_go(icon) = factor_go(icon) * coef_go(icon)

            else

              if(inaicg2%iflag_scale_aicg2 == 0)then
                 coef_go(icon) = inaicg2%ave_caicg2_nloc * e_con(icon) / e_ground
              else
                 coef_go(icon) = -inaicg2%gen_caicg2_nloc * e_con(icon)
              end if
              coef_go(icon) = factor_go(icon) * coef_go(icon)
            end if
            end if
         end do

      write(lunout,*)
      write(lunout,'(A72)')'************************************************************************'
      write(lunout,'(A72)')'*************** aicg2 parameters by cafemol ****************************'
      write(lunout,'(A6,I8)')'bond: ', nbd
      write(lunout,'(10F10.3)')(coef_bd(1,ibd),ibd=1,nbd)
      write(lunout,'(A12,I8)')'bond angle: ', nba
      write(lunout,'(10F10.3)')(coef_ba(1,iba),iba=1,nba)
      write(lunout,'(A27,I8)')'1-3 contact of bond angle: ', nba
      write(lunout,'(10F10.3)')(coef_aicg13_gauss(iba),iba=1,nba)
      write(lunout,'(10F10.3)')(wid_aicg13_gauss(iba),iba=1,nba)
      write(lunout,'(A16,I8)')'dihedral angle: ', ndih
      write(lunout,'(10F10.3)')(coef_dih(1,idih),idih=1,ndih)
      write(lunout,'(10F10.3)')(coef_dih(2,idih),idih=1,ndih)
      write(lunout,'(A31,I8)')'1-4 contact of dihedral angle: ', ndih
      write(lunout,'(10F10.3)')(coef_aicg14_gauss(idih),idih=1,ndih)
      write(lunout,'(10F10.3)')(wid_aicg14_gauss(idih),idih=1,ndih)
      write(lunout,'(A51,I8)')'gaussian coefficient and width of dihedral angle: ', ndih
      write(lunout,'(10F10.3)')(coef_dih_gauss(idih),idih=1,ndih)
      write(lunout,'(10F10.3)')(wid_dih_gauss(idih),idih=1,ndih)
      write(lunout,'(A16,I8)')'native contact: ', ncon
      write(lunout,'(10F10.3)')(coef_go(icon),icon=1,ncon)
      write(lunout,*)
#ifdef MPI_PAR
        endif
  call MPI_Bcast(coef_bd,       2*MXMPBD *nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_aicg13_gauss,       MXMPBA *nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_aicg14_gauss,      MXMPDIH*nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_dih_gauss,      MXMPDIH*nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wid_aicg13_gauss,       MXMPBA *nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wid_aicg14_gauss,      MXMPDIH*nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wid_dih_gauss,      MXMPDIH*nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_go,         MXMPCON*nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)

#endif 

end subroutine setp_aicg2

!=============================================================================

! subroutine for making atomic contact

subroutine contact_make_aicg2(iatomnum, xyz, cname_ha, ncon_atm_13, ncon_atm_14, ncon_atm)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp, only : inpro, inmisc
  use var_struct, only : ncon, nba, ndih, nmp_all, cmp2seq, icon2mp, iba2mp, idih2mp, &
                         iclass_mp, icon2unit, imp2unit

  implicit none

  real(PREC), intent(in) :: xyz(SPACE_DIM, MXATOM_MP, MXMP)
  integer,    intent(in) :: iatomnum(MXMP)
  character(4), intent(in) :: cname_ha(MXATOM_MP, MXMP)
  integer,    intent(out) :: ncon_atm(MXCON,16), ncon_atm_13(MXMP,16), ncon_atm_14(MXMP,16)

  ! -----------------------------------------------------------------
  ! local variable
  integer :: icon, iba, idih
  integer :: imp, jmp, iunit, junit
  integer :: inum_end, inum, jnum_end, jnum
  integer :: ibackbone(MXMP, MXATOM_MP)
  integer :: idonor(MXMP, MXATOM_MP)
  integer :: iacceptor(MXMP, MXATOM_MP)
  integer :: ication(MXMP, MXATOM_MP)
  integer :: ianion(MXMP, MXATOM_MP)
  real(PREC) :: xyz_i(3), xyz_j(3)

  real(PREC) :: con_atm_cutoff, con_atm_cutoff2
  real(PREC) :: hb_cutoff, hb_cutoff2
  real(PREC) :: sb_cutoff, sb_cutoff2
  real(PREC) :: dist
  real(PREC) :: dfcontact2
  integer :: ncontact_atm

!--------------------------------------------------------------------
  dfcontact2 = inpro%dfcontact * inpro%dfcontact
  con_atm_cutoff = 5.0
  con_atm_cutoff2 = con_atm_cutoff * con_atm_cutoff
  hb_cutoff = 3.2
  hb_cutoff2 = hb_cutoff * hb_cutoff
  sb_cutoff = 3.5
  sb_cutoff2 = sb_cutoff * sb_cutoff

  do imp = 1, nmp_all
     if(iclass_mp(imp) == CLASS%PRO)then

     inum_end = iatomnum(imp)
     do inum = 1, inum_end
    ! ibackbone assignment

        if(cname_ha(inum, imp) == ' N  ' .or.    &
           cname_ha(inum, imp) == ' C  ' .or.    &
           cname_ha(inum, imp) == ' O  ' .or.    &
           cname_ha(inum, imp) == ' OXT' .or.    &
           cname_ha(inum, imp) == ' CA '      ) then
  
           ibackbone(imp, inum) =1
        else
           ibackbone(imp, inum) =0
        end if

    ! hydrogen bond iacceptor and idonor assignment
 
        if(cname_ha(inum, imp)(2:2) == 'N') then
           idonor(imp, inum) = 1
           iacceptor(imp, inum) = 0
        end if

        if(cname_ha(inum, imp)(2:2) == 'O') then
           idonor(imp, inum) = 0
           iacceptor(imp, inum) = 1

           if(cmp2seq(imp) == 'SER' .and. cname_ha(inum, imp) == ' OG ' .or. &
              cmp2seq(imp) == 'THR' .and. cname_ha(inum, imp) == ' OG1' .or. &
              cmp2seq(imp) == 'TYR' .and. cname_ha(inum, imp) == ' OH '     ) then

              idonor(imp, inum) = 1

           end if
        end if

        if(cname_ha(inum, imp)(2:2) == 'S') then
           if(cmp2seq(imp) == 'CYS') then
              idonor(imp, inum) = 1
           else
              idonor(imp, inum) = 0
           end if
           iacceptor(imp, inum) = 1
        end if

      ! ication and ianion assignment 

        if(cmp2seq(imp) == 'ARG' .and. cname_ha(inum, imp) == ' NH1' .or. &
           cmp2seq(imp) == 'ARG' .and. cname_ha(inum, imp) == ' NH2' .or. &
           cmp2seq(imp) == 'LYS' .and. cname_ha(inum, imp) == ' NZ '     ) then
   
           ication(imp, inum) = 1
           ianion(imp, inum) = 0
        end if

        if(cmp2seq(imp) == 'GLU' .and. cname_ha(inum, imp) == ' OE1' .or. &
           cmp2seq(imp) == 'GLU' .and. cname_ha(inum, imp) == ' OE2' .or. &
           cmp2seq(imp) == 'ASP' .and. cname_ha(inum, imp) == ' OD1' .or. &
           cmp2seq(imp) == 'ASP' .and. cname_ha(inum, imp) == ' OD2'     ) then
                  
           ication(imp, inum) = 0
           ianion(imp, inum) = 1
        end if

     end do  !inum

     end if
  end do   !imp

!------------------------------------------------------
! counting the atomic contacts
! 1-3 local residue pairs

  do iba = 1, nba
     iunit = imp2unit(iba2mp(1, iba))
     if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .OR. &
        inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then

     imp = iba2mp(1, iba)
     jmp = iba2mp(3, iba)
                 
     inum_end = iatomnum(imp)
     jnum_end = iatomnum(jmp)
                    
     ncontact_atm=0    
     do inum = 1, inum_end

     if (xyz(1,inum,imp) > INVALID_JUDGE) cycle 
     xyz_i(1:3) = xyz(1:3, inum, imp)

     do jnum = 1, jnum_end
                    
        if (xyz(1,jnum,jmp) > INVALID_JUDGE) cycle
        xyz_j(1:3) = xyz(1:3, jnum, jmp)
                     
        dist = (xyz_j(1) - xyz_i(1))**2  &
             + (xyz_j(2) - xyz_i(2))**2  &
             + (xyz_j(3) - xyz_i(3))**2

        if(dist < dfcontact2) ncon_atm_13(iba, 16) = ncon_atm_13(iba, 16) + 1
        if(dist < con_atm_cutoff2) then
           ncontact_atm = ncontact_atm + 1  !short range contacts
! ibackbone -ibackbone contacts
           if(ibackbone(imp, inum) == 1 .and. ibackbone(jmp, jnum) == 1)then
                 
              if(iacceptor(imp, inum) == 1 .and. idonor(jmp, jnum) == 1 .or. &
                 idonor(imp, inum) == 1 .and. iacceptor(jmp, jnum) == 1  ) then
                       
                 if(dist < hb_cutoff2) then
                    ncon_atm_13(iba, 1) = ncon_atm_13(iba, 1) + 1
                 else
                    ncon_atm_13(iba, 2) = ncon_atm_13(iba, 2) + 1
                 end if

              else if(cname_ha(inum, imp)(2:2) == 'C' .or. cname_ha(jnum, jmp)(2:2) == 'C') then
                      ncon_atm_13(iba, 3) = ncon_atm_13(iba, 3) + 1 

              else
                      ncon_atm_13(iba, 4) = ncon_atm_13(iba, 4) + 1
           
              end if

           end if    ! end of ibackbone -ibackbone contacts
! sidechain - sidechain contacts
  
           if(ibackbone(imp, inum) == 0 .and. ibackbone(jmp, jnum) == 0)then

              if(iacceptor(imp, inum) == 1 .and. idonor(jmp, jnum) == 1 .or. &
                 idonor(imp, inum) == 1 .and. iacceptor(jmp, jnum) == 1  ) then
     
                 if(ication(imp, inum) == 1 .and. ianion(jmp, jnum) == 1 .or.  &
                    ianion(imp, inum) == 1 .and. ication(jmp, jnum) == 1  )   then
     
                    if(dist < sb_cutoff2) then
                       ncon_atm_13(iba, 6) = ncon_atm_13(iba, 6) + 1
                    else 
                       ncon_atm_13(iba, 9) = ncon_atm_13(iba, 9) + 1
                    end if
        
                 else if(dist < hb_cutoff2) then
                     ncon_atm_13(iba, 5) = ncon_atm_13(iba, 5) + 1
        
                 else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                          ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                      ncon_atm_13(iba, 9) = ncon_atm_13(iba, 9) + 1

                 else
                      ncon_atm_13(iba, 7) = ncon_atm_13(iba, 7) + 1

                 end if

              else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                       ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                       ncon_atm_13(iba, 9)=ncon_atm_13(iba, 9) + 1

              else if(cname_ha(inum, imp)(2:2) == 'C' .or. cname_ha(jnum, jmp)(2:2) == 'C') then
                       ncon_atm_13(iba, 8) = ncon_atm_13(iba, 8) + 1

              else
                       ncon_atm_13(iba, 10) = ncon_atm_13(iba, 10) + 1

              end if

           end if  !end of sidechain -sidechain  contacts

!ibackbone -sidechain  contacts
           if(ibackbone(imp, inum) == 0 .and. ibackbone(jmp, jnum) == 1 .or. &
              ibackbone(imp, inum) == 1 .and. ibackbone(jmp, jnum) == 0   ) then

              if(iacceptor(imp, inum) == 1 .and. idonor(jmp, jnum) == 1 .or. &
                 idonor(imp, inum) == 1 .and. iacceptor(jmp, jnum) == 1  ) then

                 if(dist < hb_cutoff2) then
                    ncon_atm_13(iba, 11) = ncon_atm_13(iba, 11) + 1

                 else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                          ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                    ncon_atm_13(iba, 14) = ncon_atm_13(iba, 14) + 1

                 else
                    ncon_atm_13(iba, 12) = ncon_atm_13(iba, 12) + 1

                 end if

              else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                       ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                    ncon_atm_13(iba, 14) = ncon_atm_13(iba, 14) + 1

              else if(cname_ha(inum, imp)(2:2) == 'C' .or. cname_ha(jnum, jmp)(2:2) == 'C') then
                    ncon_atm_13(iba, 13) = ncon_atm_13(iba, 13) + 1

              else
                    ncon_atm_13(iba, 15) = ncon_atm_13(iba, 15) + 1

              end if

           end if  !ibackbone - sidechain  contacts
        end if !con_atom

     end do
     end do

!  defining the long range contacts
     if(dfcontact2 > con_atm_cutoff2) then
        ncon_atm_13(iba, 16) = ncon_atm_13(iba, 16) - ncontact_atm
     else
        ncon_atm_13(iba, 16) = 0
     end if

!  each pair of residues cannot form more than one salty bridges.
     if(ncon_atm_13(iba, 6) >= 2)then
        ncon_atm_13(iba, 9) = ncon_atm_13(iba, 9) + ncon_atm_13(iba, 6) -1
        ncon_atm_13(iba, 6) = 1
     end if

     end if
  end do

!--------------------------------
!1-4 residue pairs
  do idih = 1, ndih
     iunit = imp2unit(idih2mp(1, idih))

     if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .OR. &
        inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then

     imp = idih2mp(1, idih)
     jmp = idih2mp(4, idih)

     inum_end = iatomnum(imp)
     jnum_end = iatomnum(jmp)

     ncontact_atm=0
     do inum = 1, inum_end

     if (xyz(1,inum,imp) > INVALID_JUDGE) cycle 
     xyz_i(1:3) = xyz(1:3, inum, imp)

     do jnum = 1, jnum_end

        if (xyz(1,jnum,jmp) > INVALID_JUDGE) cycle 
        xyz_j(1:3) = xyz(1:3, jnum, jmp)

        dist = (xyz_j(1) - xyz_i(1))**2  &
             + (xyz_j(2) - xyz_i(2))**2  &
             + (xyz_j(3) - xyz_i(3))**2

        if(dist < dfcontact2) ncon_atm_14(idih, 16) = ncon_atm_14(idih, 16) + 1
        if(dist < con_atm_cutoff2) then
           ncontact_atm = ncontact_atm + 1  !short range contacts
! ibackbone -ibackbone contacts
           if(ibackbone(imp, inum) == 1 .and. ibackbone(jmp, jnum) == 1)then

              if(iacceptor(imp, inum) == 1 .and. idonor(jmp, jnum) == 1 .or. &
                 idonor(imp, inum) == 1 .and. iacceptor(jmp, jnum) == 1  ) then

                 if(dist < hb_cutoff2) then
                    ncon_atm_14(idih, 1) = ncon_atm_14(idih, 1) + 1
                 else
                    ncon_atm_14(idih, 2) = ncon_atm_14(idih, 2) + 1
                 end if

              else if(cname_ha(inum, imp)(2:2) == 'C' .or. cname_ha(jnum, jmp)(2:2) == 'C') then
                      ncon_atm_14(idih, 3) = ncon_atm_14(idih, 3) + 1

              else
                      ncon_atm_14(idih, 4) = ncon_atm_14(idih, 4) + 1

              end if

           end if    ! end of ibackbone -ibackbone contacts

! sidechain - sidechain contacts

           if(ibackbone(imp, inum) == 0 .and. ibackbone(jmp, jnum) == 0)then

              if(iacceptor(imp, inum) == 1 .and. idonor(jmp, jnum) == 1 .or. &
                 idonor(imp, inum) == 1 .and. iacceptor(jmp, jnum) == 1  ) then

                 if(ication(imp, inum) == 1 .and. ianion(jmp, jnum) == 1 .or.  &
                    ianion(imp, inum) == 1 .and. ication(jmp, jnum) == 1  )   then

                    if(dist < sb_cutoff2) then
                       ncon_atm_14(idih, 6) = ncon_atm_14(idih, 6) + 1
                    else
                       ncon_atm_14(idih, 9) = ncon_atm_14(idih, 9) + 1
                    end if

                 else if(dist < hb_cutoff2) then
                     ncon_atm_14(idih, 5) = ncon_atm_14(idih, 5) + 1

                 else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                          ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                      ncon_atm_14(idih, 9) = ncon_atm_14(idih, 9) + 1

                 else
                      ncon_atm_14(idih, 7) = ncon_atm_14(idih, 7) + 1

                 end if

              else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                       ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                       ncon_atm_14(idih, 9)=ncon_atm_14(idih, 9) + 1

              else if(cname_ha(inum, imp)(2:2) == 'C' .or. cname_ha(jnum, jmp)(2:2) == 'C') then
                       ncon_atm_14(idih, 8) = ncon_atm_14(idih, 8) + 1

              else
                       ncon_atm_14(idih, 10) = ncon_atm_14(idih, 10) + 1

              end if

           end if  !end of sidechain -sidechain  contacts

!ibackbone -sidechain  contacts
           if(ibackbone(imp, inum) == 0 .and. ibackbone(jmp, jnum) == 1 .or. &
              ibackbone(imp, inum) == 1 .and. ibackbone(jmp, jnum) == 0   ) then

              if(iacceptor(imp, inum) == 1 .and. idonor(jmp, jnum) == 1 .or. &
                 idonor(imp, inum) == 1 .and. iacceptor(jmp, jnum) == 1  ) then

                 if(dist < hb_cutoff2) then
                    ncon_atm_14(idih, 11) = ncon_atm_14(idih, 11) + 1

                 else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                          ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                    ncon_atm_14(idih, 14) = ncon_atm_14(idih, 14) + 1

                 else
                    ncon_atm_14(idih, 12) = ncon_atm_14(idih, 12) + 1

                 end if

              else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                       ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                    ncon_atm_14(idih, 14) = ncon_atm_14(idih, 14) + 1

              else if(cname_ha(inum, imp)(2:2) == 'C' .or. cname_ha(jnum, jmp)(2:2) == 'C') then
                    ncon_atm_14(idih, 13) = ncon_atm_14(idih, 13) + 1

              else
                    ncon_atm_14(idih, 15) = ncon_atm_14(idih, 15) + 1

              end if

           end if  !ibackbone - sidechain  contacts
        end if !con_atom

     end do
     end do

!  defining the long range contacts
     if(dfcontact2 > con_atm_cutoff2) then
        ncon_atm_14(idih, 16) = ncon_atm_14(idih, 16) - ncontact_atm
     else
        ncon_atm_14(idih, 16) = 0
     end if

!  each pair of residues cannot form more than one salty bridges.
     if(ncon_atm_14(idih, 6) >= 2)then
        ncon_atm_14(idih, 9) = ncon_atm_14(idih, 9) + ncon_atm_14(idih, 6) -1
        ncon_atm_14(idih, 6) = 1
     end if

     end if
  end do

!---------------------------------------------
! counting the atomic contacts
  do icon = 1, ncon

     iunit = icon2unit(1, icon)
     junit = icon2unit(2, icon)
     if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG2))then

     imp = icon2mp(1, icon)
     jmp = icon2mp(2, icon)

     inum_end = iatomnum(imp)
     jnum_end = iatomnum(jmp)

     ncontact_atm=0
     do inum = 1, inum_end

     if (xyz(1,inum,imp) > INVALID_JUDGE) cycle 
     xyz_i(1:3) = xyz(1:3, inum, imp)

     do jnum = 1, jnum_end

        if (xyz(1,jnum,jmp) > INVALID_JUDGE) cycle 
        xyz_j(1:3) = xyz(1:3, jnum, jmp)

        dist = (xyz_j(1) - xyz_i(1))**2  &
             + (xyz_j(2) - xyz_i(2))**2  &
             + (xyz_j(3) - xyz_i(3))**2

        if(dist < dfcontact2) ncon_atm(icon, 16) = ncon_atm(icon, 16) + 1
        if(dist < con_atm_cutoff2) then
           ncontact_atm = ncontact_atm + 1  !short range contacts
! ibackbone -ibackbone contacts
           if(ibackbone(imp, inum) == 1 .and. ibackbone(jmp, jnum) == 1)then

              if(iacceptor(imp, inum) == 1 .and. idonor(jmp, jnum) == 1 .or. &
                 idonor(imp, inum) == 1 .and. iacceptor(jmp, jnum) == 1  ) then

                 if(dist < hb_cutoff2) then
                    ncon_atm(icon, 1) = ncon_atm(icon, 1) + 1
                 else   
                    ncon_atm(icon, 2) = ncon_atm(icon, 2) + 1
                 end if
               
              else if(cname_ha(inum, imp)(2:2) == 'C' .or. cname_ha(jnum, jmp)(2:2) == 'C') then
                      ncon_atm(icon, 3) = ncon_atm(icon, 3) + 1

              else
                      ncon_atm(icon, 4) = ncon_atm(icon, 4) + 1
     
              end if
                              
           end if    ! end of ibackbone -ibackbone contacts

! sidechain - sidechain contacts

           if(ibackbone(imp, inum) == 0 .and. ibackbone(jmp, jnum) == 0)then 

              if(iacceptor(imp, inum) == 1 .and. idonor(jmp, jnum) == 1 .or. &
                 idonor(imp, inum) == 1 .and. iacceptor(jmp, jnum) == 1  ) then

                 if(ication(imp, inum) == 1 .and. ianion(jmp, jnum) == 1 .or.  &
                    ianion(imp, inum) == 1 .and. ication(jmp, jnum) == 1  )   then

                    if(dist < sb_cutoff2) then
                       ncon_atm(icon, 6) = ncon_atm(icon, 6) + 1
                    else 
                       ncon_atm(icon, 9) = ncon_atm(icon, 9) + 1
                    end if

                 else if(dist < hb_cutoff2) then
                     ncon_atm(icon, 5) = ncon_atm(icon, 5) + 1

                 else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                          ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then
  
                      ncon_atm(icon, 9) = ncon_atm(icon, 9) + 1

                 else
                      ncon_atm(icon, 7) = ncon_atm(icon, 7) + 1

                 end if

              else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                       ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                       ncon_atm(icon, 9)=ncon_atm(icon, 9) + 1
                       
              else if(cname_ha(inum, imp)(2:2) == 'C' .or. cname_ha(jnum, jmp)(2:2) == 'C') then
                       ncon_atm(icon, 8) = ncon_atm(icon, 8) + 1

              else 
                        ncon_atm(icon, 10) = ncon_atm(icon, 10) + 1

              end if

           end if  !end of sidechain -sidechain  contacts

!ibackbone -sidechain  contacts
           if(ibackbone(imp, inum) == 0 .and. ibackbone(jmp, jnum) == 1 .or. &
              ibackbone(imp, inum) == 1 .and. ibackbone(jmp, jnum) == 0   ) then

              if(iacceptor(imp, inum) == 1 .and. idonor(jmp, jnum) == 1 .or. &
                 idonor(imp, inum) == 1 .and. iacceptor(jmp, jnum) == 1  ) then

                 if(dist < hb_cutoff2) then
                    ncon_atm(icon, 11) = ncon_atm(icon, 11) + 1

                 else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                          ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                    ncon_atm(icon, 14) = ncon_atm(icon, 14) + 1

                 else
                    ncon_atm(icon, 12) = ncon_atm(icon, 12) + 1

                 end if

              else if (ication(imp, inum) == 1 .or. ication(jmp, jnum) == 1 .or. &
                       ianion(imp, inum) == 1 .or. ianion(jmp, jnum) == 1  ) then

                    ncon_atm(icon, 14) = ncon_atm(icon, 14) + 1

              else if(cname_ha(inum, imp)(2:2) == 'C' .or. cname_ha(jnum, jmp)(2:2) == 'C') then
                    ncon_atm(icon, 13) = ncon_atm(icon, 13) + 1

              else
                    ncon_atm(icon, 15) = ncon_atm(icon, 15) + 1

              end if

           end if  !ibackbone - sidechain  contacts
        end if !con_atom

     end do
     end do

!  defining the long range contacts
     if(dfcontact2 > con_atm_cutoff2) then
        ncon_atm(icon, 16) = ncon_atm(icon, 16) - ncontact_atm
     else
        ncon_atm(icon, 16) = 0
     end if

!  each pair of residues cannot form more than one salty bridges.

     if(ncon_atm(icon, 6) >= 2)then
        ncon_atm(icon, 9) = ncon_atm(icon, 9) + ncon_atm(icon, 6) -1
        ncon_atm(icon, 6) = 1
     end if

     end if
  end do

end subroutine contact_make_aicg2
