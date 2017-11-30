! simu_ele_set
!> @brief Calculate parameters for electrostatic interaction

! ***********************************************************************
subroutine simu_ele_set(grep, tempk, ionic_strength)
  
  use const_maxsize
  use const_physical
  use const_index
  use var_setp, only : inele, inmisc, insimu, inperi, inpmf
  use var_struct, only : ncharge, coef_charge, icharge2mp, imp2type, cmp2seq, num_ion, nmp_all, lmp2charge
  use var_simu, only : ewld_d_coef, pmfdh_energy, pmfdh_force
  use var_io, only : outfile
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none
  ! ----------------------------------------------------------------------
  integer,    intent(in) :: grep
  real(PREC), intent(in) :: tempk
  real(PREC), intent(in) :: ionic_strength

  integer :: icharge, icharge_change, imp, i
  integer :: itype, ibin
  real(PREC) :: ek, Tc, lb, xi, Zp, theta, b, b3
  real(PREC) :: r, DH, e_add, Rbininv
  real(PREC) :: v1, v2, c1v1, c2v2
  real(PREC) :: conc_Mg
  character(CARRAY_MSG_ERROR) :: error_message
  real(PREC), parameter ::  MM_A=87.740e0_PREC, MM_B=-0.4008e0_PREC  ! i_diele=1
  real(PREC), parameter ::  MM_C=9.398e-4_PREC, MM_D=-1.410e-6_PREC  ! i_diele=1

  integer :: ifunc_seq2id
  logical, save :: flg_first = .True.

  Zp = 1.0 ! To suppress compiler warning
  b = inele%length_per_unit( ifunc_seq2id('P  '))

  ! -----------------------------------------------------------------------
  ! Dielectric constant
  ek = 0.0
  if(inele%i_diele == 0) then
     ek = inele%diele_water

  !Temperature dependent (Malmberg and Maryott, 1956)
  else if (inele%i_diele == 1) then

     if (inmisc%i_temp_independent <= 1) then
        Tc = tempk - 273.15e0_PREC
        ek =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc

     elseif (inmisc%i_temp_independent == 2) then
        Tc = insimu%tempk_ref - 273.15e0_PREC
        ek =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc
        inele%diele_dTcoef = 1.0e0_PREC + insimu%tempk_ref / ek   & 
                             * (MM_B + 2.0e0_PREC*MM_C*Tc + 3.0e0_PREC*MM_D*Tc*Tc)
     else
        error_message = 'Error: invalid i_temp_independent in simu_ele_set'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

  !Temperature dependent
  else if (inele%i_diele == 2) then

     if (inmisc%i_temp_independent <= 1) then
        ek = 1.0 / (0.03273600947 * tempk * BOLTZ_KCAL_MOL_ND - 0.00669625750)
     else if (inmisc%i_temp_independent == 2) then
        ek = 1.0 / (- 0.00669625750)
     else
        error_message = 'Error: invalid i_temp_independent in simu_ele_set'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

  else
     error_message = 'Error: logical defect in simu_ele_set'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  inele%diele = ek

  ! -----------------------------------------------------------------------
  ! For the dipoler term in Ewald
  if (inele%ewld_dipole == 1) then
     ewld_d_coef = 2.0 * F_PI / ((1.0+2.0) * inperi%psize(1)**3)
  else if (inele%ewld_dipole == 2) then
     ewld_d_coef = 2.0 * F_PI / ((1.0+2.0*ek) * inperi%psize(1)**3)
  else
     ewld_d_coef = 0.0e0_PREC
  endif

  ! -----------------------------------------------------------------------
  ! Charge value
  if(inele%i_charge == 1) then

     ! Bjerrum length
     if (inmisc%i_temp_independent <= 1) then
        lb = ELE * ELE / (4.0e0_PREC * F_PI * EPSI_0 * ek * BOLTZ_J * tempk) * 1.0e10_PREC
     elseif (inmisc%i_temp_independent == 2) then
        lb = ELE * ELE / (4.0e0_PREC * F_PI * EPSI_0 * ek * BOLTZ_J * insimu%tempk_ref) * 1.0e10_PREC
     else
        lb = 0.0
        error_message = 'Error: invalid i_temp_independent in simu_ele_set'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     ! Calculate phosphate charge taking into account monovalent salt condensation
     ! Semiexplicit ion model
     if (inmisc%i_dtrna_model == 2019) then
        xi = lb / b 
        b3 = b ** 3
      
        if (num_ion(IONTYPE%MG) > 0) then
           conc_Mg = num_ion(IONTYPE%MG) / (N_AVO * 1.0e-27_PREC * inperi%psize(1) * inperi%psize(2) * inperi%psize(3))
   
           ! Currently, only with divalent ion
   
           ! Formal derivation
           !v1 = 4 * F_PI * F_E * b3 * (xi - 1.0e0_PREC) * 2
           !v2 = 4 * F_PI * F_E * b3 * (xi - 0.5e0_PREC) * (2 + 1)
           !c1v1 = (N_AVO * 1.0e-27_PREC) * ionic_strength * v1
           !c2v2 = (N_AVO * 1.0e-27_PREC) * inele%conc_Mg  * v2
           !theta = (-c1v1**2 + sqrt(c1v1**4 + 8.0e0_PREC*F_E * c1v1**2 * c2v2 * (1.0e0_PREC - 1.0e0_PREC/xi))) &
           !       / (4.0e0_PREC * F_E * c2v2)
   
           ! Eliminate F_E (v1 is actually v1/e, v2 is v2/e [e is the base of the natural logarithm])
           v1 = 4 * F_PI * b3 * (xi - 1.0e0_PREC) * 2
           v2 = 4 * F_PI * b3 * (xi - 0.5e0_PREC) * (2 + 1)
           c1v1 = (N_AVO * 1.0e-27_PREC) * ionic_strength * v1
           c2v2 = (N_AVO * 1.0e-27_PREC) * conc_Mg  * v2
         
           theta = (-c1v1**2 + sqrt(c1v1**4 + 8.0e0_PREC * c1v1**2 * c2v2 * (1.0e0_PREC - 1.0e0_PREC/xi))) &
                  / (4.0e0_PREC * c2v2)
           Zp = 1.0 - theta
        else

           theta = 1.0 - 1.0/xi
           Zp = b / lb

        endif

     ! Implicit ion model
     else
        xi = lb / b   ! ~ 0.4
        theta = 1.0 - 1.0/xi
        !Zp = 1.0 - theta
        !Zp = 1.0 / xi
        Zp = b / lb

     endif

     ! Reflect Zp
     loop_charge:&
     do icharge = 1, ncharge
        imp = icharge2mp(icharge)
        do icharge_change = 1, inele%n_charge_change
           if (inele%charge_change_imp(icharge_change) == imp) then
              coef_charge(icharge,grep) = inele%charge_change_value(icharge_change)
              cycle loop_charge
           endif
        enddo
        if (imp2type(imp) == MPTYPE%RNA_PHOS) then
           coef_charge(icharge,grep) = - Zp
        else if (imp2type(imp) == MPTYPE%SOPSC) then
           coef_charge(icharge,grep) = inele%coef_charge_type( ifunc_seq2id( cmp2seq(imp) ) )
        endif
     enddo loop_charge

  endif
   
  ! ----------------------------------------------------------------------
  ! coef: j_kcal * eq**2 / (4.0e0_PREC * F_PI * e0 * ek * rij)
  !   =  332.063713019 / ek
  inele%coef(grep) = JOUL2KCAL_MOL * 1.0e10_PREC * ELE**2 &
                    / (4.0e0_PREC * F_PI * EPSI_0 * ek)
  !! Only for consistency to Denesyuk
  !inele%coef(grep) = 332.0637090 / ek

  ! Kd: sqrt(e0 * ek * RT / 2 * NA**2 * eq**2 * I)
  !   = sqrt(1.0e-3 * e0 * ek * kb / 2 * NA * eq**2) * sqrt(T(K) / I(M))
  if (inmisc%i_temp_independent == 0) then
     inele%cdist(grep) = 1.0e10_PREC                               &
                     * sqrt( (1.0e-3_PREC * EPSI_0 * ek * BOLTZ_J) &
                            / (2.0_PREC * N_AVO * ELE**2)  )       &
                     * sqrt(tempk / ionic_strength)
  else if (inmisc%i_temp_independent >= 1) then
     inele%cdist(grep) = 1.0e10_PREC                               &
                     * sqrt( (1.0e-3_PREC * EPSI_0 * ek * BOLTZ_J) &
                            / (2.0_PREC * N_AVO * ELE**2)  )       &
                     * sqrt(insimu%tempk_ref / ionic_strength)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Semiexplicit model
  if (inmisc%i_dtrna_model == 2019 .and. num_ion(IONTYPE%MG) > 0) then

     if(inele%i_charge /= 1) then
        error_message = 'i_charge must be 1 for semiexplicit ion model'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     itype = PMFTYPE%MG_P
   
     do ibin = 1, inpmf%Nbin(itype)
   
        r = inpmf%Rmin(itype) + inpmf%Rbin(itype) * real(ibin-1, kind=PREC)
                    
        DH = JOUL2KCAL_MOL * 1.0e10_PREC * ELE**2 / (4.0e0_PREC * F_PI * EPSI_0 * ek) &
            * 2.0     &  ! Mg charge
            * (- Zp)  &  ! Phosphate charge
            * exp(-r / inele%cdist(grep)) / r
   
        if (itype == PMFTYPE%MG_P) then
           e_add = (DH - inpmf%pmf(ibin, itype)) * exp(- (inpmf%Rmerge(itype) ** 2) / (r*r))
        else
           error_message = 'q must be 2 in simu_set_semiexplicit'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        
        pmfdh_energy(ibin, grep, itype) = inpmf%pmf(ibin, itype) + e_add
        !write(*,*) ibin, r, inpmf%pmf(ibin,itype), DH, e_add, pmfdh_energy(ibin,grep,itype)
     enddo
   
     Rbininv = 1.0 / inpmf%Rbin(itype)
     do ibin = 1, inpmf%Nbin(itype) - 1
        pmfdh_force(ibin,grep,itype) = - Rbininv * &
                   (pmfdh_energy(ibin+1, grep, itype) - pmfdh_energy(ibin, grep, itype))
     enddo

  endif ! Semiexplicit model
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (flg_first) then
#ifdef MPI_PAR
  if (myrank == 0) then
#endif
     
     write(outfile%data,'(a)') '<<<< simu_ele_set'
     write(outfile%data,'(a,i4)') 'i_dtrna_model: ', inmisc%i_dtrna_model

     if (inmisc%i_dtrna_model == 2019) then
        write(outfile%data,'(a,1x,i10)') '#Mg:', num_ion(IONTYPE%MG)
        write(outfile%data,'(a,1x,3(g12.5))') 'Box:', (inperi%psize(i), i=1,3)
        write(outfile%data,'(a,1x,f9.6)') '[M+]:', ionic_strength
        write(outfile%data,'(a,1x,f9.6)') '[Mg++]:', conc_Mg
        write(outfile%data,'(a,1x,g12.5)') 'ek:', ek
        write(outfile%data,'(a,1x,g12.5)') 'b:', b 
        write(outfile%data,'(a,1x,g12.5)') 'lb:', lb
        write(outfile%data,'(a,1x,g12.5)') 'theta:', theta
        write(outfile%data,'(a,1x,g12.5)') 'Zp:', Zp
   
        if (num_ion(IONTYPE%MG) > 0) then
           write(outfile%data,'(a)') ''
           write(outfile%data,'(a)') '<<<<<< Effective potential for Mg-P'
           write(outfile%data,'(a)') '<<    r     W(r)     PMF       DH      W-PMF'
   
           itype = PMFTYPE%MG_P
         
           do ibin = 1, inpmf%Nbin(itype)
              r = inpmf%Rmin(itype) + inpmf%Rbin(itype) * real(ibin-1, kind=PREC)
              DH = JOUL2KCAL_MOL * 1.0e10_PREC * ELE**2 / (4.0e0_PREC * F_PI * EPSI_0 * ek) &
                  * 2.0     &  ! Mg charge
                  * (- Zp)  &  ! Phosphate charge
                  * exp(-r / inele%cdist(grep)) / r
         
              e_add = (DH - inpmf%pmf(ibin, itype)) * exp(- (inpmf%Rmerge(itype) ** 2) / (r*r))
              
              pmfdh_energy(ibin, grep, itype) = inpmf%pmf(ibin, itype) + e_add
   
              write(outfile%data, '(6(1x,f8.3))') r, pmfdh_energy(ibin,grep,itype), inpmf%pmf(ibin,itype), DH, e_add , pmfdh_force(ibin,grep,itype)
           enddo
   
           write(outfile%data,*) '>>>>>>'
        endif
     endif

     write(outfile%data,'(a)') ''
     write(outfile%data,'(a)') '<<<<<< Charges'
     write(outfile%data,'(a)') '#  imp   coef_charge'
     do imp = 1, nmp_all
        if (lmp2charge(imp) /= 0) then
           write(outfile%data,'(i8,1x,f8.4)') imp, coef_charge(lmp2charge(imp),grep)
        endif
     enddo
     write(outfile%data,'(a)') '>>>>>>'

     write(outfile%data,'(a)') '>>>>'

#ifdef MPI_PAR
  end if
#endif

     flg_first = .False.
  endif

end subroutine simu_ele_set
