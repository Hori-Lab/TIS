program exact_coulomb
   implicit none

   integer, parameter :: PREC = 8

   ! Physical constants
   real(PREC), parameter :: N_AVO   = 6.022140857e23_PREC   !< Avogadro constant [/mol]
   real(PREC), parameter :: KCAL2JOUL = 4184.0              !< (kcal -> J)  [J/kcal]
   real(PREC), parameter :: JOUL2KCAL = 1.0/KCAL2JOUL   !< (J -> kcal)  [kcal/J]
   real(PREC), parameter :: JOUL2KCAL_MOL  = JOUL2KCAL * N_AVO  !< (J -> kcal/mol)
   !real(PREC), parameter :: JOUL2KCAL_MOL = 1.43862e20  ! (J -> kcal/mol)
   real(PREC), parameter :: FPI     = 3.14159265358979323846264338e0
   !real(PREC), parameter :: ELE     = 1.60217648740e-19 ! Elementary charge [C]
   real(PREC), parameter :: ELE     = 1.6021766208e-19_PREC !< Elementary charge [C]
   real(PREC), parameter :: EPSI_0  = 8.854187817e-12   ! Vacuum permittivity [F/m]
   real(PREC), parameter :: ORIGIN(3) = [0.0,0.0,0.0]

   ! Parameters of the lattice
   integer, parameter :: NLATTICE = 50
   !integer, parameter :: NLATTICE = 100
   integer, parameter :: MAXCELLS = (2*NLATTICE+1)**3 - 1
   integer, parameter :: FILEPDB = 10

   real(PREC), parameter :: CELLSIZE = 150.0

   ! Parameters of the input PDB
   !integer, parameter :: NMP = 202
   !character(100), parameter :: filename = 'BOX150_NaCl_0050mM_eq.pdb'
   integer, parameter :: NMP = 3904
   character(100), parameter :: filename = 'BOX150_NaCl_1000mM_eq.pdb'

   ! Temperature and dielectric constant
   real(PREC), parameter :: TK = 310.0
   real(PREC), parameter :: Tc = TK - 273.15
   real(PREC), parameter :: ek =  87.740 -0.4008*Tc + 9.398e-4*Tc*Tc -1.410e-6*Tc*Tc*Tc
   real(PREC), parameter :: coef = JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0e0_PREC * FPI * EPSI_0 * ek)

   integer :: i, il
   real(PREC) :: xyz(3,NMP), charge(NMP)
   real(PREC) :: ene, ene_total
   real(PREC), allocatable :: hs(:, :)
   real(PREC), allocatable :: norms(:)
   real(PREC) :: norm_pre
   integer :: ncells

   !###################
   !! Read PDB
   !###################
   open(FILEPDB, file=filename, status = 'old', action = 'read')
   do i = 1, NMP
      read(FILEPDB, '(6x,5x,1x,4x,1x,3x,2x,4x,4x,3f8.3,12x)') xyz(1,i), xyz(2,i), xyz(3,i) 
      if (i <= NMP/2) then
         charge(i) = 1.0
      else
         charge(i) = -1.0
      endif
   enddo

   write(*,*) '#coef=', coef

   allocate( hs(3, MAXCELLS) )
   allocate( norms(MAXCELLS) )

   !###################
   !! Generate reciprocal lattices (hs) and norms
   !###################
   !write(*, *) '# Calling generate_reciprocal_lattice'
   call generate_reciprocal_lattice(NLATTICE, ncells)

   !do i =1, NCELLS
   !   write(*,*) hs(1,i), hs(2,i), hs(3,i), norms(i)
   !enddo

   ene_total = 0.0e0_PREC

   !###################
   !! All surrounding cells
   !###################
   !norm_pre = 0.0
   norm_pre = norms(ncells)

   !do il = 1, NCELLS
   do il = ncells, 1, -1

      if (norms(il) /= norm_pre) then
         write(*,*) norm_pre, ene_total
         norm_pre = norms(il)
      endif

      ene = coef * coulomb( hs(:,il), .False. )
      ene_total = ene_total + ene

      !write(*,*) '#', il, norms(il), ene, ene_total
   enddo

   write(*,*) norm_pre, ene_total


   !###################
   !! Origin
   !###################
   ene = coef * coulomb( ORIGIN, .True.)
   ene_total = ene_total + ene
   write(*,*) '#', 0, 0.0, ene, ene_total
   write(*,*) 0.0, ene_total


contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Coulomb potential contributed from a cell "h"
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(PREC) function coulomb( h, flg_origin )

      real(PREC), intent(in) :: h(3)
      logical, intent(in) :: flg_origin
      real(PREC) :: xyzj(3)
      integer :: i, j

      coulomb = 0.0

      if (flg_origin) then
         !! For the origin cell, summation is all pairwise exept i=j.
         do i = 1, NMP
            do j = i+1, NMP
               coulomb = coulomb + pair_energy_coulomb(xyz(:,i), xyz(:,j), charge(i), charge(j))
            enddo
         enddo

         !! Note that there is not the factor 1/2 since we do not double-count.

      else
         do j = 1, NMP
            xyzj(:) = xyz(:,j) + h(1:3) * CELLSIZE
            do i = 1, NMP
               coulomb = coulomb + pair_energy_coulomb(xyz(:,i), xyzj, charge(i), charge(j))
            enddo
         enddo

         !! There is a factor 1/2 needed since a half of the energy is for the other cell.
         coulomb = 0.5e0_PREC * coulomb
      endif
   
   endfunction coulomb


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Pairwise Coulomb potential energy
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(PREC) function pair_energy_coulomb( r1, r2, c1, c2)
      real(PREC), intent(in) :: r1(3)
      real(PREC), intent(in) :: r2(3)
      real(PREC), intent(in) :: c1
      real(PREC), intent(in) :: c2
      real(PREC) :: v(3), d
      v(:) = r2(:) - r1(:)
      d = sqrt( dot_product(v, v))
      pair_energy_coulomb = c1 * c2 / d
      !write(*,*) c1, c2, d, pair_energy_coulomb
   endfunction pair_energy_coulomb


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine generate_reciprocal_lattice(m, ncells)
      integer, intent(in) :: m
      integer, intent(out) :: ncells
      !real(PREC), intent(out) :: hs(3, (m+1)**3 - 1)
      !real(PREC), intent(out) :: norms((m+1)**3 - 1)
      real(PREC) :: h(3), h_swp(3)
      real(PREC) :: r_swp, m_real, x
      integer :: i, ix, iy, iz
      integer :: icomb
      logical :: flg_swp

      m_real = real(m, kind=PREC)

      i = 0
      do ix = -m, m
         h(1) = real(ix)

         do iy = -m, m
            h(2) = real(iy)

            do iz = -m, m

               if (ix == 0 .and. iy == 0 .and. iz == 0) then
                  cycle
               endif
               
               h(3) = real(iz)

               x = sqrt(dot_product(h,h))
               if (x > m_real) then
                  cycle
               endif

               i = i+1
               hs(1:3,i) = h(1:3)
               norms(i) = x

            enddo
         enddo
      enddo
      
      ncells = i
      !if (i /= NCELLS) then
      !   write(*, *) 'Error: i /= NCELLS, i=', i, ' NCELLS=', NCELLS
      !   stop
      !endif

      !write(*,*) '# Calling comb11 sort'
      icomb = NCELLS
      flg_swp = .false.
      do while (icomb > 1 .or. flg_swp)
         icomb = (icomb * 10) / 13
         if (icomb < 1) then
            icomb = 1
         else if (icomb == 9 .or. icomb == 10) then 
            icomb = 11
         endif
         flg_swp = .false.
         do i = 1, NCELLS-icomb
            if (norms(i) > norms(i+icomb)) then
               r_swp          = norms(i) 
               norms(i)      = norms(i+icomb)
               norms(i+icomb)= r_swp
               h_swp(:)     = hs(:,i)
               hs(:,i)      = hs(:,i+icomb)
               hs(:,i+icomb)= h_swp(:)
               flg_swp = .true.
            endif
         enddo
      enddo

   endsubroutine generate_reciprocal_lattice

endprogram exact_coulomb
