!simu_ele_set_allrep
!> @brief This subroutine is wrapper to call "simu_ele_set" for all of   &
!>        replicas.                                                      &
!>        Whenever temperature or ionic-strength are changed in REMD     &
!>        simulation, electrostatic parameters need to be reset          &
!>        by calling this subroutine.

subroutine simu_para2(tempk_in, ionic_strength_in)
  
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : indna, inmisc, inwind, inele
  use var_struct, only : nunit_real, lunit2mp, iclass_unit, &
                         lele, coef_ele, iele2mp, coef_charge, lmp2charge
  use var_replica,only : n_replica_all, n_replica_mpi, rep2val, flg_rep, inrep, irep2grep
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none
  ! ----------------------------------------------------------------------
  real(PREC), intent(in) :: tempk_in
  real(PREC), intent(in) :: ionic_strength_in

  ! ----------------------------------------------------------------------
  ! local variables
  integer    :: irep  ! local replica ID
  integer    :: grep  ! global replica ID
  integer    :: idna, n_nt, ln_nt, iunit, i_nt
  integer    :: ksta, kend, klen
  integer    :: iele, icharge, jcharge
  real(PREC) :: tempk
  real(PREC) :: ionic_strength
  real(PREC) :: pullforce
  character(CARRAY_MSG_ERROR) :: error_message

  ! ----------------------------------------------------------------------
  tempk = tempk_in
  ionic_strength = ionic_strength_in

  if (inmisc%class_flag(CLASS%DNA)) then
     n_nt = 0
     idna = 0
     i_nt = 0
     do iunit = 1, nunit_real
        if(iclass_unit(iunit) == CLASS%DNA) then
           ln_nt = (lunit2mp(2, iunit) - lunit2mp(1, iunit) + 3)/3
           if(idna == 0) then
              n_nt = ln_nt
              idna = 1
           else
              if(ln_nt /= n_nt) then 
                 i_nt = 1
              end if
           end if
        end if
     end do

     if(inmisc%iflag_solv_nt_dna == 0) then
        indna%e_n = indna%csolvmax_dna * (1.0 - 1.0/(1.40418 - 0.268231 * n_nt)) 
        if(i_nt == 1) then
           error_message = 'ERROR: should specify "iflag_solv_nt_dna = 1" when  using different lengths of DNA'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
     else
        indna%e_n = indna%csolvmax_dna
     end if
  end if

  ! ----------------------------------------------------------------------
  do grep = 1, n_replica_all
   
     if (flg_rep(REPTYPE%TEMP)) then
        tempk = rep2val(grep, REPTYPE%TEMP)
     endif

     if (flg_rep(REPTYPE%ION)) then
        ionic_strength = rep2val(grep, REPTYPE%ION)
     endif

     if (flg_rep(REPTYPE%WIND)) then
        inwind%iwind(grep) = int(rep2val(grep, REPTYPE%WIND))
     endif

     if (flg_rep(REPTYPE%WINZ)) then
        inwind%iwinz(grep) = int(rep2val(grep, REPTYPE%WINZ))
     endif

     
     if (flg_rep(REPTYPE%PULL)) then
        pullforce = rep2val(grep, REPTYPE%PULL)
        inmisc%pull_unravel_xyz(:,1:inrep%n_pull,grep) = pullforce * inrep%pull_direction(:, 1:inrep%n_pull)
     endif

     ! ----------------------------------------------------------------------
     if (inmisc%force_flag(INTERACT%ELE)) then
        call simu_ele_set(grep, tempk, ionic_strength)
     end if
  
     ! ----------------------------------------------------------------------
     if (inmisc%class_flag(CLASS%DNA)) then
        call simu_solv_set(grep, ionic_strength)
     end if
  
     ! ----------------------------------------------------------------------
     if (inmisc%force_flag(INTERACT%DTRNA)) then
        call simu_set_dtrna(grep, tempk)
     endif

  enddo


  if (inmisc%force_flag(INTERACT%ELE) .and. &
      inele%i_diele /= 0 .and. inele%i_calc_method == 0) then

     do irep = 1, n_replica_mpi
        grep = irep2grep(irep)
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH
        klen=(lele(irep)-1+npar_mpi)/npar_mpi
        ksta=1+klen*local_rank_mpi
        kend=min(ksta+klen-1, lele(irep))
#else
        ksta=1
        kend=lele(irep)
#endif
#else
        ksta=1
        kend=lele(irep)
#endif
        do iele=ksta,kend
           icharge = lmp2charge( iele2mp(1,iele,irep) )
           jcharge = lmp2charge( iele2mp(2,iele,irep) )
           coef_ele(iele,irep) = coef_charge(icharge,grep) * coef_charge(jcharge,grep) &
                                * inele%coef(grep)
        enddo
     enddo
  endif
end subroutine simu_para2
