!simu_energy_nlocal_mgo
!> @brief Calculates the energy related to Go-type interctions in  &
!>        the multiple Go model. The Values are added into         &
!>        "e_exv(E_TYPE%GO)" and "e_exv_unit(,,E_TYPE%GO)".

! ************************************************************************
! formula of go1210
! ego = coef_go * {5*(go_nat/x)**12 -6*(go_nat/x)**10}
! coef_go = icon_dummy_mgo * factor_go * cgo1210
!
! parameter list
! cgo1210: constant of go energy
! factor_go: value of amino acid specifity (ex. 1.0, MJ)
! go_nat: distance of native contact 
! ************************************************************************
subroutine simu_energy_nlocal_mgo(irep, now_con, e_exv_unit, e_exv)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inpro
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, &
                          ncon, icon2mp, coef_go, go_nat
  use var_mgo,     only : ncontype_mgo, irefcon_mgo
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer,    intent(in)  :: irep
  integer,    intent(out) :: now_con(:,:)
  real(PREC), intent(inout) :: e_exv(:)
  real(PREC), intent(inout) :: e_exv_unit(:,:,:)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: imp1, imp2
  integer :: icon, iswitch, jcon
  integer :: iunit, junit, imirror
  real(PREC) :: efull
  real(PREC) :: cut_off, cut_off2, dist2, pre, pre5over6
  real(PREC) :: gorad2, gorad2j = 0.0
  real(PREC) :: sigma1, sigma2
  real(PREC) :: rovdist2, rovdist10, rovdist12, rjudge, rjudge_contact
  real(PREC) :: v21(SPACE_DIM)
  character(CARRAY_MSG_ERROR) :: error_message
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
!!$omp master

  ! ---------------------------------------------------------------------
  ! set parameter 
  cut_off = inpro%cutoff_go**2
  rjudge_contact = 1.2e0_PREC**2
  pre5over6 = 5.0e0_PREC / 6.0e0_PREC
  jcon = -999

  ! ---------------------------------------------------------------------
  ! zero clear
!  now_con(:,:) = 0
  pre = 0.0

  ! ---------------------------------------------------------------------
#ifdef MPI_PAR3
  klen=(ncon-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,ncon)
#else
  ksta = 1
  kend = ncon
#endif
!$omp do private(imp1,imp2,v21,dist2,gorad2,iswitch,jcon,gorad2j, &
!$omp&           sigma2,sigma1,pre,rovdist2, &
!$omp&           rovdist12,rovdist10,error_message,efull, &
!$omp&           rjudge,iunit,junit,imirror)
  do icon=ksta,kend

     imp1 = icon2mp(1, icon)
     imp2 = icon2mp(2, icon)
     
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
            
     gorad2 = go_nat(icon) * go_nat(icon)
     cut_off2 = gorad2 * cut_off
!this cut_off2 line will be deleted
     if(ncontype_mgo(icon) /= 2) then
        rjudge = gorad2 * rjudge_contact 
        if(dist2 < rjudge) then 
           now_con(1, icon) = 1
        else
           now_con(1, icon) = 0
        end if
        if(coef_go(icon) > ZERO_JUDGE) then
           now_con(2, icon) = 1
        else
           now_con(2, icon) = 0
        end if
     end if
     
     if(dist2 > cut_off2) cycle
!This cut off sentece will be deleted and switched to ones below commented out with double !s
     iswitch = 0
     if(ncontype_mgo(icon) == 0) then
        iswitch = 1

     else if(ncontype_mgo(icon) == 1 .or. &
             ncontype_mgo(icon) == 2) then
        jcon = irefcon_mgo(icon)
        gorad2j = go_nat(jcon)**2
        sigma2 = pre5over6 * gorad2j
        
        if(ncontype_mgo(icon) == 1) then
           sigma1 = pre5over6 * gorad2
           if(dist2 < sigma2) then
              iswitch = 2
           else if(dist2 >= sigma2 .and. dist2 <= sigma1) then
              iswitch = 3
           else
              iswitch = 1
           end if
              
        else if(ncontype_mgo(icon) == 2) then
           if(dist2 < sigma2) then
              iswitch = 2
           else
              iswitch = 3
           end if
        end if
     end if
        
        
     if(iswitch == 1) then
!!        cut_off2 = gorad2 * cut_off
        !add IF sentence for rna around here, if cutoff_go for rna is different from that for protein (Now such IF sentence is not added because by default both are same)
!!        if(dist2 > cut_off2) cycle
   
        pre = coef_go(icon)
        rovdist2 = gorad2 / dist2
        call dist12_10(rovdist2, rovdist12, rovdist10) 
     else if(iswitch == 2) then
!!        cut_off2 = gorad2j * cut_off
!!        if(dist2 > cut_off2) cycle
 
        pre = coef_go(jcon)
        rovdist2 = gorad2j / dist2
        call dist12_10(rovdist2, rovdist12, rovdist10) 
     else if(iswitch == 3) then
        cycle
     else
        error_message = 'Error: invalid value for iswitch in simu_energy_nlocal_mgo'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     
     efull = pre * (5.0e0_PREC * rovdist12 - 6.0e0_PREC * rovdist10)

     ! --------------------------------------------------------------------
     ! sum of the energy
     e_exv(E_TYPE%GO) = e_exv(E_TYPE%GO) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     e_exv_unit(iunit, junit, E_TYPE%GO) = e_exv_unit(iunit, junit, E_TYPE%GO) + efull
  end do
!$omp end do nowait
!!$omp end master

! ************************************************************************
contains 

! ************************************************************************
  subroutine dist12_10(dist2, dist12, dist10)
    
    implicit none

    ! -------------------------------------------------------------------
    real(PREC), intent(in) :: dist2
    real(PREC), intent(out) :: dist12, dist10

    ! -------------------------------------------------------------------
    ! local variables
    real(PREC) :: dist4, dist8

    ! -------------------------------------------------------------------
    dist4 = dist2 * dist2
    dist8 = dist4 * dist4
    dist12 = dist4 * dist8
    dist10 = dist2 * dist8

  end subroutine dist12_10
  
end subroutine simu_energy_nlocal_mgo
