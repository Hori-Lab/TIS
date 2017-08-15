program test_tables

   implicit none

   integer, parameter :: nmp_real = 101
   integer, parameter :: npar_mpi = 4

   integer, allocatable :: imp_l2g(:)
   integer :: nmp_l
   integer :: n, lmp
   integer :: local_rank_mpi

   do local_rank_mpi = 0, npar_mpi-1

      n = nmp_real/npar_mpi + 1
      write(*,*) 'n=',n
   
      allocate( imp_l2g(n))
   
      call tables2( imp_l2g, nmp_l, nmp_real ) 
   
      write(*,*) 'nmp_l=', nmp_l 
      do lmp  = 1, nmp_l
         write(*,*) 'imp_l2g(',lmp,')=',imp_l2g(lmp)
      enddo

      deallocate(imp_l2g)

   enddo

contains
   subroutine tables2( il2g, n_l, n )
     integer,intent(in)  :: n
     integer,intent(out) :: n_l
     integer,intent(out) :: il2g(n/npar_mpi+1)
   
     logical,parameter :: cyclic = .true.
     integer :: incr, ip, len
     integer :: i, is, ie
   
   
     if( cyclic ) then
       n_l  =  0
       incr =  1
       ip   = -1
   
       do i = 1, n
         ip = ip + incr
         if( ip < 0 .or. ip > npar_mpi-1 ) then
           incr = -incr
           ip = ip + incr
         endif
   
         if( ip == local_rank_mpi ) then
           n_l = n_l + 1
           il2g(n_l) = i
         end if
       end do
   
       else
         len = (n-1+npar_mpi)/npar_mpi
         is  = 1+len*local_rank_mpi
         ie  = min(is+len-1,n)
   
         n_l = 0
         do i = is, ie
           n_l = n_l + 1
           il2g(n_l) = i
         end do 
   
       end if
   
       print *,"[mpi] n_l  = ",n_l
   !   print *,"[mpi] il2g = ",il2g(1:n_l)
   
   end subroutine tables2
endprogram test_tables

