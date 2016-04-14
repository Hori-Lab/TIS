subroutine dump_var_replica(lunout)

  use const_maxsize
  use var_replica

  implicit none

  integer, intent(in) :: lunout
  integer :: i, j, k, ni, nj, nk, ivar

  write(lunout, *) ''
  write(lunout, *) '### var_replica'

  write(lunout, *) 'n_replica_all,', n_replica_all
  write(lunout, *) 'inrep % n_step_replica,', inrep%n_step_replica

  do ivar = 1, REPTYPE%MAX
     write(lunout, *) 'inrep % lowest(',ivar,'),', inrep%lowest(ivar)
     write(lunout, *) 'inrep % highest(',ivar,'),', inrep%highest(ivar)
     !write(lunout, *) 'inrep % interval(',ivar,'),', inrep%interval(ivar)
     do i = 1, MXREPLICA
        write(lunout, *) 'inrep % var(',i,',',ivar,'),', inrep%var(i,ivar)
     enddo
  enddo

!  integer, allocatable :: hist_combination(:,:)
!  integer, allocatable :: hist_exchange(:)
!  integer              :: hist_attempt
  if (allocated(hist_combination)) then
     ni = ubound(hist_combination,1)
     nj = ubound(hist_combination,2)
     do j = 1, nj
        do i = 1, ni
           write(lunout, *) 'hist_combination(',i,',',j,'),', hist_combination(i,j)
        enddo
     enddo
  else
     write(lunout, *) 'hist_combination(:,:) is not allocated.'
  endif
  if (allocated(hist_exchange)) then
     ni = ubound(hist_exchange,1)
     nj = ubound(hist_exchange,2)
     nk = ubound(hist_exchange,3)
     do i = 1, ni
        do j = 1, nj
           do k = 1, nk
              write(lunout, *) 'hist_exchange(',i,',',j,',',k,'),', hist_exchange(i,j,k)
           enddo
        enddo
     enddo
  else
     write(lunout, *) 'hist_exchange(:,:,:) is not allocated.'
  endif
  write(lunout, *) 'hist_attempt,', hist_attempt

  ! changing table
!  integer,    save :: rep2lab(MXREPLICA)  ! replica => label
!  integer,    save :: lab2rep(MXREPLICA)  ! label   => replica
!  real(PREC), save :: lab2temp(MXREPLICA) ! label   => temperature
!  real(PREC), save :: lab2ion(MXREPLICA)  ! label   => ionic strength
  write(lunout, *) '# changing table'
  do i = 1, MXREPLICA
     write(lunout, *) 'rep2lab(',i,'),', rep2lab(i)
  enddo
  do i = 1, MXREPLICA
     write(lunout, *) 'lab2rep(',i,'),', lab2rep(i)
  enddo
  do ivar = 1, REPTYPE%MAX
     do i = 1, MXREPLICA
        write(lunout, *) 'lab2val(',i,',',ivar,'),', lab2val(i,ivar)
     enddo
  enddo

endsubroutine dump_var_replica
