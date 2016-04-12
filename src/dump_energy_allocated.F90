#ifdef _DUMP
#define _DUMP_ENERGY_ALLOCATED
#endif

#ifdef _DUMP_ENERGY_ALLOCATED

  ! about energy
  write(lundump, *) '# energy  step=',istep
!  real(PREC),allocatable :: pnlet(:,:)          ! (E_TYPE%MAX, replica)
  if (allocated(pnlet)) then
     ni_dump = ubound(pnlet,1)
     nj_dump = ubound(pnlet,2)
     do j_dump = 1, nj_dump
        do i_dump = 1, ni_dump
           write(lundump, *) 'pnlet(,',i_dump,',',j_dump,'),',pnlet(i_dump,j_dump)
        enddo
     enddo
  else 
     write(lundump, *) 'pnlet is not allocated.'
  endif

!  real(PREC),allocatable :: pnle_unit(:,:,:,:)  ! (unit, unit, E_TYPE%MAX, replica)
  if (allocated(pnle_unit)) then
     ni_dump = ubound(pnle_unit,1)
     nj_dump = ubound(pnle_unit,2)
     nk_dump = ubound(pnle_unit,3)
     nl_dump = ubound(pnle_unit,4)
     do l_dump = 1, nl_dump
           do j_dump = 1, nj_dump
              do i_dump = 1, ni_dump
        do k_dump = 1, nk_dump
                 write(lundump, *) 'pnle_unit(',i_dump,',',j_dump,',',k_dump,',',l_dump,'),', &
                 pnle_unit(i_dump,j_dump,k_dump,l_dump)
              enddo
           enddo
        enddo
     enddo
  else 
     write(lundump, *) 'pnle_unit is not allocated.'
  endif

#endif
