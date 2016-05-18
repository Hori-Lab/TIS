!neighbor_list
!> @brief Constructs the neighboring list for non-local interactions. &
!>        In the runup to caluculations, this subroutine calls        &
!>        "neighbor_pre" to specify which unit pair(s) contain   &
!>        non-local contacts.

subroutine neighbor_list(irep, ineigh2mp, lmp2neigh)
  
  use if_neighbor
  use const_maxsize
  use const_index
  use var_setp, only : inmisc, inperi
  use var_struct, only : nmp_real, lunit2mp, xyz_mp_rep, pxyz_mp_rep, &
                         imp2unit, nmp_all
  use mpiconst

  implicit none

  ! -------------------------------------------------------------------
  integer, intent(in)  :: irep
  integer, intent(out) :: lmp2neigh((nmp_l+nthreads-1)/nthreads  ,0:nthreads-1)
  integer, intent(out) :: ineigh2mp(MXMPNEIGHBOR*nmp_all/nthreads,0:nthreads-1)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: n
  integer :: klen, ksta, kend
  integer :: imp, jmp, iunit, junit
  integer :: imirror
  integer :: ineighbor(0:nthreads-1)
  integer :: ineigh_unit(MXUNIT, MXUNIT)
  real(PREC) :: dist2, v21(3)
  character(CARRAY_MSG_ERROR) :: error_message
  integer :: imp_l

  ! -------------------------------------------------------------------
  ! calc neigh_unit
  call neighbor_pre(xyz_mp_rep(:,:,irep), ineigh_unit)

  ! -------------------------------------------------------------------
  ! calc ineigh2mp
  ineighbor(0:nthreads-1) = 0

!$omp parallel 
!$omp do private(klen,ksta,kend,imp_l,imp,iunit,jmp,junit,dist2,v21,imirror)
  do n = 0, nthreads-1
  klen=(nmp_l-1+nthreads)/nthreads
  ksta=1+klen*n
  kend=min(ksta+klen-1,nmp_l)
  
  do imp_l = ksta, kend
     imp = imp_l2g(imp_l)

     iunit = imp2unit(imp)
     jmp = imp + 1
     do while (jmp <= nmp_real)
        junit = imp2unit(jmp)

        if(ineigh_unit(iunit, junit) == 1) then

           if(inperi%i_periodic == 0) then
              v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
           else
              v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
              call util_pbneighbor(v21, imirror)
           end if

           dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

           if(dist2 < inmisc%rneighbordist2_unit(iunit, junit)) then
              ineighbor(n) = ineighbor(n) + 1
              ineigh2mp(ineighbor(n),n) = jmp
           end if
        else
           ! jump to last point of 'junit'
           jmp = lunit2mp(2, junit)
        endif

        jmp = jmp + 1
     end do
     lmp2neigh(imp_l-ksta+1,n) = ineighbor(n)
  end do
  end do
!$omp end do nowait
!$omp end parallel

  if(any(ineighbor > MXMPNEIGHBOR*nmp_all/nthreads)) then
     error_message = 'Error: too big ineighbor in neighbor_list'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

end subroutine neighbor_list
