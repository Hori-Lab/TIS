!force_dtrna_stack
!> @brief Calculate forces related to stacking interaction of RNA particles.
!
!   See energy_dtrna_stack.F90 for detail
!

 subroutine force_dtrna_stack(irep, force_mp) 
 
   use const_maxsize
   use const_physical
   use const_index
   use var_struct,  only : xyz_mp_rep,  ndtrna_st, idtrna_st2mp, dtrna_st_nat, coef_dtrna_st, nmp_all
   use var_simu,    only : st_status
   use var_replica, only : irep2grep
   use mpiconst
 
   implicit none
 
   integer,    intent(in)    :: irep
   real(PREC), intent(inout) :: force_mp(3,nmp_all)
 
   integer :: ist, grep
   integer :: ksta, kend
   real(PREC) :: dist, ddist, dih, d
   real(PREC) :: v21(3), v34(3), v54(3), v56(3), v76(3)
   real(PREC) :: abs54, d5454, abs56, d5656
   real(PREC) :: d5654
   real(PREC) :: d7656over5656, d5456over5656
   real(PREC) :: d3454over5454, d5654over5454
   real(PREC) :: m(3), n(3)
   real(PREC) :: dnn, dmm
   real(PREC) :: for(3,7), f_i(3), f_l(3), ediv
#ifdef MPI_PAR
   integer :: klen
#endif
 
   ! --------------------------------------------------------------------

   grep = irep2grep(irep)

#ifdef MPI_PAR
   klen=(ndtrna_st-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,ndtrna_st)
#else
   ksta = 1
   kend = ndtrna_st
#endif

!$omp do private(ediv, for, d, dist, ddist, dih, f_i, f_l,&
!$omp&           v21,v34,v54,v56,v76,&
!$omp&           m,n,dmm,dnn,d5454,d5656,d5654,abs54,abs56,&
!$omp&           d3454over5454,d5654over5454,d7656over5656,d5456over5656)
   do ist=ksta,kend

      if (.not. st_status(ist, irep)) then
         cycle
      endif

      ediv = 1.0e0_PREC
      for(:,:) = 0.0e0_PREC
     
      !===== calc vectors =====
      v21(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(2,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(1,ist), irep)
      v34(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(3,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(4,ist), irep)
      v54(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(5,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(4,ist), irep)
      v56(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(5,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(6,ist), irep)
      v76(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(7,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(6,ist), irep)

      !===== 1. Distance between 1 and 2 =====
      dist = norm2(v21)
      ddist = dist - dtrna_st_nat(1,ist)
      ediv = ediv + coef_dtrna_st(1, ist, grep) * ddist**2

      f_i(:) = 2.0e0_PREC * coef_dtrna_st(1,ist,grep) * ddist / dist * v21(:)
      for(:,1) = - f_i(:)
      for(:,2) = + f_i(:)

      !===== 2. Dihedral angle among 3-4-5-6 =====
      m(1) = v34(2)*v54(3) - v34(3)*v54(2)
      m(2) = v34(3)*v54(1) - v34(1)*v54(3)
      m(3) = v34(1)*v54(2) - v34(2)*v54(1)
      n(1) = v54(2)*v56(3) - v54(3)*v56(2)
      n(2) = v54(3)*v56(1) - v54(1)*v56(3)
      n(3) = v54(1)*v56(2) - v54(2)*v56(1)

      dmm = dot_product(m,m)
      dnn = dot_product(n,n)
      d5454 = dot_product(v54,v54)
      abs54 = sqrt(d5454)
      d5654 = dot_product(v56,v54)
      d3454over5454 = dot_product(v34,v54) / d5454
      d5654over5454 = d5654 / d5454

      dih = atan2(dot_product(v34,n)*sqrt(d5454) , dot_product(m,n))
      d = dih - dtrna_st_nat(2,ist)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif
      ediv = ediv + coef_dtrna_st(2, ist, grep) * d**2

      f_i(:) = + 2.0e0_PREC * coef_dtrna_st(2,ist,grep) * d * abs54 / dmm * m(:)
      f_l(:) = - 2.0e0_PREC * coef_dtrna_st(2,ist,grep) * d * abs54 / dnn * n(:)

      for(:,3) = f_i(:)
      for(:,4) = (-1.0e0_PREC + d3454over5454) * f_i(:) &
                -(             d5654over5454) * f_l(:)
      for(:,5) = (-1.0e0_PREC + d5654over5454) * f_l(:) &
                -(             d3454over5454) * f_i(:)
      for(:,6) = f_l(:)

      !===== 3. Dihedral angle among 7-6-5-4 =====
      m(1) = v76(2)*v56(3) - v76(3)*v56(2)
      m(2) = v76(3)*v56(1) - v76(1)*v56(3)
      m(3) = v76(1)*v56(2) - v76(2)*v56(1)
      !n(1) = v56(2)*v54(3) - v56(3)*v54(2)
      !n(2) = v56(3)*v54(1) - v56(1)*v54(3)
      !n(3) = v56(1)*v54(2) - v56(2)*v54(1)
      n(:) = -n(:)

      dmm = dot_product(m,m)
      !dnn = dot_product(n,n)  !! dnn does not change.
      d5656 = dot_product(v56, v56)
      abs56 = sqrt(d5656)
      d7656over5656 = dot_product(v76,v56) / d5656
      d5456over5656 = d5654 / d5656

      dih = atan2(dot_product(v76,n)*sqrt(d5656) , dot_product(m,n))
      d = dih - dtrna_st_nat(3,ist)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif
      ediv = ediv + coef_dtrna_st(3, ist, grep) * d**2

      f_i(:) = + 2.0e0_PREC * coef_dtrna_st(3,ist,grep) * d * abs56 / dmm * m(:)
      f_l(:) = - 2.0e0_PREC * coef_dtrna_st(3,ist,grep) * d * abs56 / dnn * n(:)

      for(:,7) = for(:,7) + f_i(:)
      for(:,6) = for(:,6) + (-1.0e0_PREC + d7656over5656) * f_i(:) &
                          - (              d5456over5656) * f_l(:)
      for(:,5) = for(:,5) - (              d7656over5656) * f_i(:) &
                          + (-1.0e0_PREC + d5456over5656) * f_l(:)
      for(:,4) = for(:,4) + f_l(:)


      !===== Total =====
      for(:,:) = coef_dtrna_st(0,ist,grep) / ediv**2 * for(:,:)

      force_mp(1:3,idtrna_st2mp(1,ist)) = force_mp(1:3,idtrna_st2mp(1,ist)) + for(1:3,1)
      force_mp(1:3,idtrna_st2mp(2,ist)) = force_mp(1:3,idtrna_st2mp(2,ist)) + for(1:3,2)
      force_mp(1:3,idtrna_st2mp(3,ist)) = force_mp(1:3,idtrna_st2mp(3,ist)) + for(1:3,3)
      force_mp(1:3,idtrna_st2mp(4,ist)) = force_mp(1:3,idtrna_st2mp(4,ist)) + for(1:3,4)
      force_mp(1:3,idtrna_st2mp(5,ist)) = force_mp(1:3,idtrna_st2mp(5,ist)) + for(1:3,5)
      force_mp(1:3,idtrna_st2mp(6,ist)) = force_mp(1:3,idtrna_st2mp(6,ist)) + for(1:3,6)
      force_mp(1:3,idtrna_st2mp(7,ist)) = force_mp(1:3,idtrna_st2mp(7,ist)) + for(1:3,7)
   end do
!$omp end do nowait

end subroutine force_dtrna_stack
