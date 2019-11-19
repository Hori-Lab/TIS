subroutine neighbor_list_exv(irep)
  
  use if_neighbor
  use if_util
  use const_maxsize
  use const_index
  use var_setp, only : inmisc, inperi
  use var_struct, only : lunit2mp, xyz_mp_rep, pxyz_mp_rep, &
                         imp2unit, nmp_all, nexv_pairs, iexv_pairs, lexv, iexv2mp
  use mpiconst

  implicit none

  ! -------------------------------------------------------------------
  integer, intent(in)  :: irep

  ! -------------------------------------------------------------------
  ! local variables
  integer :: ipair, imp, jmp, iunit, junit
  integer :: n
  integer :: iexv(0:nthreads-1)
  integer :: iexv2mp_l(3,nexv_pairs,0:nthreads-1)
  integer :: imirror
  integer :: klen, ksta, kend
  real(PREC) :: dist2, v21(3)
  character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
  write(6,*) '####### start neighbor_list_exv'
#endif
  iexv(:) = 0
  lexv(1,:,irep) = 1
  lexv(2,:,irep) = 0

  klen = (nexv_pairs-1+nthreads) / nthreads

!$omp parallel
!$omp do private(n,ksta,kend,ipair,imp,jmp,v21,imirror,dist2,iunit,junit)
  do n = 0, nthreads-1
  ksta = 1 + klen*n
  kend = min(ksta+klen-1, nexv_pairs)

  !do ipair = 1, nexv_pairs
  do ipair = ksta, kend
     imp = iexv_pairs(1, ipair)
     jmp = iexv_pairs(2, ipair)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
        imirror = 1
     else
        v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
        call util_pbneighbor(v21, imirror)
     end if

     dist2 = dot_product(v21,v21)
     iunit = imp2unit(imp)
     junit = imp2unit(jmp)

     if(dist2 < inmisc%rneighbordist2_unit(iunit, junit)) then

        iexv(n) = iexv(n) + 1

        iexv2mp_l(1,iexv(n),n) = imp
        iexv2mp_l(2,iexv(n),n) = jmp
        iexv2mp_l(3,iexv(n),n) = imirror
     end if
  enddo
  enddo
!$omp end do nowait
!$omp end parallel

  !! Copy iexv2mp_l to iexv2mp
  iexv2mp(:,1:iexv(0),irep) = iexv2mp_l(:,1:iexv(0),0)
  do n = 1, nthreads-1
     ksta = iexv(0) + 1
     iexv(0) = iexv(0) + iexv(n)
     iexv2mp(:,ksta:iexv(0),irep) = iexv2mp_l(:,1:iexv(n),n)
  enddo

  lexv(2,E_TYPE%EXV_WCA,irep) = iexv(0)

#ifdef _DEBUG
  write(*,*) inmisc%rneighbordist2_unit(1,1)
  write(*,*) lexv(2,E_TYPE%EXV_WCA,irep)
  !write(*,*) '####'
  !do ipair = lexv(1,E_TYPE%EXV_WCA,irep), lexv(2,E_TYPE%EXV_WCA,irep)
  !   write(*,*) iexv2mp(1,ipair,irep), iexv2mp(2,ipair,irep), iexv2mp(3,ipair,irep)
  !enddo
  !write(*,*) '####'
  write(6,*) '####### end neighbor_list_exv'
#endif
end subroutine neighbor_list_exv
