!force_dtrna_hbond13
!> @brief Calculate forces related to stacking interaction of RNA particles.
!
! See energy_dt13_hbond.F90 for detail

subroutine force_dtrna_hbond13(irep, force_mp) 

  use const_maxsize
  use const_physical
  use const_index
!  use var_setp,  only : inperi
  use var_struct,only : xyz_mp_rep, ndtrna_hb, idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb, nmp_all
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3,nmp_all)

  integer :: ihb, ksta, kend
  real(PREC) :: d, cos_theta, dih
  real(PREC) :: v12(3), v13(3), v53(3), v42(3), v46(3)
  real(PREC) :: a42, a13, a12
  real(PREC) :: d1212, d1313, d4242
  real(PREC) :: d1213, d1242, d4246, d1353
  real(PREC) :: d1213over1212, d1213over1313
  real(PREC) :: d1242over1212, d1242over4242
  real(PREC) :: d4246over4242
  real(PREC) :: d1353over1313
  real(PREC) :: m(3), n(3)
  real(PREC) :: c4212(3), c1213(3)
  real(PREC) :: c4212_abs2, c1213_abs2
  real(PREC) :: dnn, dmm
  real(PREC) :: pre
  real(PREC) :: for(3,6), f_i(3), f_k(3), f_l(3), ediv
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------

#ifdef MPI_PAR3
   klen=(ndtrna_hb-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,ndtrna_hb)
#else
   ksta = 1
   kend = ndtrna_hb
#endif

!$omp do private(f_i,f_k,f_l,for,pre,ediv,m,n,dmm,dnn,&
!$omp&           d,cos_theta,dih,&
!$omp&           v12,v13,v53,v42,v46,a12,a13,a42,d1212,d1313,d4242,d1213,d1242,d4246,d1353,&
!$omp&           d1213over1212,d1213over1313,d1242over1212,d1242over4242,&
!$omp&           d4246over4242,d1353over1313,c4212,c1213,c4212_abs2,c1213_abs2)
   do ihb=ksta,kend

      ediv = 1.0e0_PREC
      for(:,:) = 0.0e0_PREC
     
!      !===== calc vectors =====
!      if(inperi%i_periodic == 0) then
         v12 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
         v13 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
         v53 = xyz_mp_rep(1:3, idtrna_hb2mp(5,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
         v42 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
         v46 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(6,ihb), irep)
!      else
!         ! Add periodic boundary
!      end if

      d1212 = dot_product(v12,v12)
      d1313 = dot_product(v13,v13)
      d4242 = dot_product(v42,v42)
      a12 = sqrt(d1212)
      a13 = sqrt(d1313)
      a42 = sqrt(d4242)
      d1213 = dot_product(v13,v12)
      d1242 = dot_product(v12,v42)
      d4246 = dot_product(v42,v46)
      d1353 = dot_product(v13,v53)
      d1213over1212 = d1213 / d1212
      d1213over1313 = d1213 / d1313
      d1242over1212 = d1242 / d1212
      d4246over4242 = d4246 / d4242
      d1242over4242 = d1242 / d4242
      d1353over1313 = d1353 / d1313

      !===== Distance =====
      d = a12 - dtrna_hb_nat(1,ihb)
      ediv = ediv + coef_dtrna_hb(1, ihb) * d**2
      f_i(:) = 2.0e0_PREC * coef_dtrna_hb(1,ihb) * d * v12(:) / a12
      for(:,1) = + f_i(:)
      for(:,2) = - f_i(:)

      !===== Angle of 3-1=2  =====
      cos_theta = d1213 / (a13 * a12)
      d = acos(cos_theta) - dtrna_hb_nat(2,ihb)
      pre = 2.0e0_PREC * coef_dtrna_hb(2,ihb) * d / sqrt(d1313*d1212 - d1213**2)
      f_i(:) = pre * (v12(:) - (d1213over1313 * v13(:)))
      f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))
      for(:,3) = for(:,3) + f_i(:)
      for(:,2) = for(:,2) + f_k(:)
      for(:,1) = for(:,1) - f_i(:) - f_k(:)
      ediv = ediv + coef_dtrna_hb(2, ihb) * d**2

      !===== Angle of 1=2-4  =====
      cos_theta = d1242 / (a12 * a42)
      d = acos(cos_theta) - dtrna_hb_nat(3,ihb)
      pre = 2.0e0_PREC * coef_dtrna_hb(3,ihb) * d / sqrt(d1212*d4242 - d1242**2)
      f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1242over4242 * v42(:)))
      for(:,1) = for(:,1) + f_i(:)
      for(:,4) = for(:,4) + f_k(:)
      for(:,2) = for(:,2) - f_i(:) - f_k(:)
      ediv = ediv + coef_dtrna_hb(3, ihb) * d**2

     
      !===== Dihedral angle among 4-2=1=3 =====
      c4212(1) = v42(2)*v12(3) - v42(3)*v12(2)
      c4212(2) = v42(3)*v12(1) - v42(1)*v12(3)
      c4212(3) = v42(1)*v12(2) - v42(2)*v12(1)
      c1213(1) = v12(2)*v13(3) - v12(3)*v13(2)
      c1213(2) = v12(3)*v13(1) - v12(1)*v13(3)
      c1213(3) = v12(1)*v13(2) - v12(2)*v13(1)
      c4212_abs2 = dot_product(c4212,c4212)
      c1213_abs2 = dot_product(c1213,c1213)

      dih = atan2(dot_product(v42,c1213)*sqrt(d1212) , dot_product(c4212,c1213))
      d = dih - dtrna_hb_nat(4,ihb)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif
      ediv = ediv + coef_dtrna_hb(4,ihb) * d**2

      pre = 2.0e0_PREC * coef_dtrna_hb(4,ihb) * d * a12
      f_i(:) = + pre / c4212_abs2 * c4212(:)
      f_l(:) = - pre / c1213_abs2 * c1213(:)

      for(:,4) = for(:,4) + f_i(:)
      for(:,2) = for(:,2) + (-1.0e0_PREC + d1242over1212) * f_i(:) &
                          - (              d1213over1212) * f_l(:)
      for(:,1) = for(:,1) + (-1.0e0_PREC + d1213over1212) * f_l(:) &
                          - (              d1242over1212) * f_i(:)
      for(:,3) = for(:,3) + f_l(:)


      !===== Dihedral angle among 5-3-1=2 =====
      m(1) = v53(2) * v13(3) - v53(3) * v13(2)
      m(2) = v53(3) * v13(1) - v53(1) * v13(3)
      m(3) = v53(1) * v13(2) - v53(2) * v13(1)
      n(:) = -c1213(:)
      dmm = dot_product(m,m)
      dnn = c1213_abs2

      dih = atan2(dot_product(v53,n)*a13 , dot_product(m,n))
      d = dih - dtrna_hb_nat(5, ihb)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif
      ediv = ediv + coef_dtrna_hb(5,ihb) * d**2

      pre = 2.0e0_PREC * coef_dtrna_hb(5,ihb) * d * a13
      f_i(:) = + pre / dmm * m(:)
      f_l(:) = - pre / dnn * n(:)

      for(:,5) = for(:,5) + f_i(:)
      for(:,3) = for(:,3) + (-1.0e0_PREC + d1353over1313) * f_i(:) &
                          - (              d1213over1313) * f_l(:)
      for(:,1) = for(:,1) + (-1.0e0_PREC + d1213over1313) * f_l(:) &
                          - (              d1353over1313) * f_i(:)
      for(:,2) = for(:,2) + f_l(:)


      !===== Dihedral angle among 1=2-4-6 =====
      m(:) = -c4212(:)
      n(1) = v42(2) * v46(3) - v42(3) * v46(2)
      n(2) = v42(3) * v46(1) - v42(1) * v46(3)
      n(3) = v42(1) * v46(2) - v42(2) * v46(1)
      dmm = c4212_abs2
      dnn = dot_product(n,n)

      dih = atan2(dot_product(v12,n)*a42 , dot_product(m,n))
      d = dih - dtrna_hb_nat(6, ihb)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif
      ediv = ediv + coef_dtrna_hb(6,ihb) * d**2

      pre = 2.0e0_PREC * coef_dtrna_hb(6,ihb) * d * a42
      f_i(:) = + pre / dmm * m(:)
      f_l(:) = - pre / dnn * n(:)

      for(:,1) = for(:,1) + f_i(:)
      for(:,2) = for(:,2) + (-1.0e0_PREC + d1242over4242) * f_i(:) &
                          - (              d4246over4242) * f_l(:)
      for(:,4) = for(:,4) + (-1.0e0_PREC + d4246over4242) * f_l(:) &
                          - (              d1242over4242) * f_i(:)
      for(:,6) = for(:,6) + f_l(:)

      !===== Total =====
      for(:,:) = for(:,:) * coef_dtrna_hb(0,ihb) / ediv**2

      force_mp(1:3,idtrna_hb2mp(1,ihb)) = force_mp(1:3,idtrna_hb2mp(1,ihb)) + for(:,1)
      force_mp(1:3,idtrna_hb2mp(2,ihb)) = force_mp(1:3,idtrna_hb2mp(2,ihb)) + for(:,2)
      force_mp(1:3,idtrna_hb2mp(3,ihb)) = force_mp(1:3,idtrna_hb2mp(3,ihb)) + for(:,3)
      force_mp(1:3,idtrna_hb2mp(4,ihb)) = force_mp(1:3,idtrna_hb2mp(4,ihb)) + for(:,4)
      force_mp(1:3,idtrna_hb2mp(5,ihb)) = force_mp(1:3,idtrna_hb2mp(5,ihb)) + for(:,5)
      force_mp(1:3,idtrna_hb2mp(6,ihb)) = force_mp(1:3,idtrna_hb2mp(6,ihb)) + for(:,6)
   end do
!$omp end do nowait

end subroutine force_dtrna_hbond13
