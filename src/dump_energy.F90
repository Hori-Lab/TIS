  ! about energy
  write(lundump, *) '# energy  step=',istep
!  real(PREC),allocatable :: pnlet(:,:)          ! (E_TYPE%MAX, replica)
     ni_dump = ubound(pnlet,1)
     nj_dump = ubound(pnlet,2)
     do j_dump = 1, nj_dump
        do i_dump = 1, ni_dump
           write(lundump, *) 'pnlet(,',i_dump,',',j_dump,'),',pnlet(i_dump,j_dump)
        enddo
     enddo

!  real(PREC),allocatable :: pnle_unit(:,:,:,:)  ! (unit, unit, E_TYPE%MAX, replica)
     ni_dump = ubound(pnle_unit,1)
     nj_dump = ubound(pnle_unit,2)
     nk_dump = ubound(pnle_unit,3)
     nl_dump = ubound(pnle_unit,4)
     do l_dump = 1, nl_dump
        do k_dump = 1, nk_dump
           do j_dump = 1, nj_dump
              do i_dump = 1, ni_dump
                 write(lundump, *) 'pnle_unit(',i_dump,',',j_dump,',',k_dump,',',l_dump,'),', &
                 pnle_unit(i_dump,j_dump,k_dump,l_dump)
              enddo
           enddo
        enddo
     enddo
