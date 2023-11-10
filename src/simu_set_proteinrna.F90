module simu_set_proteinrna

    use var_struct, only: coef_LJ1210, LJ1210_nat, LJ1210_nat_2
    use const_maxsize

    implicit none

    LJ1210_nat = 3.6e0_PREC     ! This is what is used as the equilibrium bond length for protein-RNA aromatic interactions currently
    LJ1210_nat_2 = LJ1210_nat**2

    ! how do I reference to each aromatic interaction?

    do i = 1

