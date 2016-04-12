program mt_cafe

  use mt_stream

  implicit none

  integer, parameter :: NSTREAM=10
  type(mt_state) :: mts(0:NSTREAM-1)
  integer :: iseed = 123456789
!  integer :: iseeda(4) = (/ Z'123', Z'234', Z'345', Z'456' /)
  integer :: i, k, is
  integer(8) :: jj


  call set_mt19937
  call new(mts(0))
  call init(mts(0),iseed)  ! init by scalar
!  call init(mts(0),iseeda)  ! init by array
  call print(mts(0))

  do is=1,NSTREAM-1
    call create_stream(mts(0),mts(is),is)
  enddo
  do is=0,NSTREAM-1
    call print(mts(is))
  enddo

  do is=0,NSTREAM-1
    write(*,'("@ first 10 random integers from stream(",I4,")")')is
    write(*,'("@ int32")')
    do i=1,10
       k = genrand_int32(mts(is))
      jj = k
      if (jj < 0) jj = jj + 2_8**32
      write(*,'("@ ",Z8,I11)')k,jj
    enddo
  enddo

  do is=0,NSTREAM-1
    call delete(mts(is))
  enddo

  stop

end program mt_cafe
