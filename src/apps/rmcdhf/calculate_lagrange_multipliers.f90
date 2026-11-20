!***********************************************************************
!                                                                      *
subroutine calculate_lagrange_multipliers(extended_optimal_level)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real32
    !                                                                      *
    !   Sets up the data structure  pertaining to the Lagrange multipli-   *
    !   ers  on the first entry;  on subsequent calls it  determines new   *
    !   estimates for the multipliers.                                     *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, QUAD, RINTI,                           *
    !               [RSCF92]: DACON, SETCOF, XPOT,  YPOT.                  *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 18 Dec 1992   *
    !   MPI version by Xinghong He              Last update: 03 Aug 1998   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: orbital_index_radix
    use array_memory_management, only: alloc
    use scf_orbital_order, only: iorder, max_orbitals
    use closed_core_subshells, only: closed_core_subshell_count
    use atomic_calculation_state, only: accy, max_radial_points, speed_of_light_au
    use fixed_orbital_flags, only: fixed_orbital_count, is_orbital_fixed, max_orbitals
    use radial_grid_state, only: radial_point_count, radial_derivative_over_radius
    use orbital_lagrange_multipliers, only: ecv, iecc, nec
    use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa, orbital_principal_n, &
        orbital_symmetry_labels, orbital_two_j
    use serial_mpi_context, only: mpi_rank
    use scf_radial_potentials, only: yp, xp, xq
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    use radial_quadrature_workspace, only: mtp, ta
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii

    ! Procedure interfaces
    use assemble_orbital_potential_coefficients_interface, only: assemble_orbital_potential_coefficients
    use build_direct_radial_potential_interface, only: build_direct_radial_potential
    use build_exchange_radial_potential_interface, only: build_exchange_radial_potential
    use apply_off_diagonal_integrals_interface, only: apply_off_diagonal_integrals

    ! Module dependencies
    use radial_quadrature, only: quad
    use one_electron_radial_integral, only: rinti
    implicit none (type, external)

    ! Arguments
    logical  :: extended_optimal_level

    ! Local constants
    real(real64), parameter :: p001 = 1.0e-1_real64
    integer, parameter :: orbital_label_radix = orbital_index_radix

    ! Local declarations
    integer :: itwice, liraw, li, lip1, naklj, orbital_index, lj, ieccli, l1, l2, &
        jlast, mlast, m, j, i
    real(real64) :: ypj(max_radial_points), ypm(max_radial_points), xpj(max_radial_points), xpm(max_radial_points), &
        xqj(max_radial_points), xqm(max_radial_points)
    real(real64) :: eps, ucfj, ucfm, result, rijm, qdif, obqdif, obqsum
    logical :: first, fixli, fixlj, fulli, fullj, vli, vlj

    ! Initial values
    !-----------------------------------------------
    !
    data first/ .true./
    !
    !-----------------------------------------------------------------------

    if (first) then

        !=======================================================================
        !   Determine the total number of Lagrange multipliers and store
        !   their indeces in IECC(1:NEC). Memories are allocated for IECC
        !   and ECV. The only outputs are NEC and IECC(). This is done only
        !   for the first time this routine is called.
        !
        !   Implementation changes:
        !   The "DO LIraw ..." is executed twice (if NEC obtained in the
        !   first time > 0) so the memory allocation is done only once.
        !
        !   This part is not distributed.
        !=======================================================================

        eps = accy*0.01_real64      ! criterion to see if an orb is occupied
        do itwice = 1, 2
            nec = 0
            do orbital_index = closed_core_subshell_count+1, orbital_count
                lj = iorder(orbital_index)
                naklj = orbital_kappa(lj)
                vlj = .not. is_orbital_fixed(lj)
                fullj = abs ( ucf(lj)-dble (orbital_two_j(lj)+1) ) < eps
                do liraw = 1, orbital_index-1
                    li = iorder(liraw)
                    vli = .not. is_orbital_fixed(li)
                    fulli = abs ( ucf(li)-dble (orbital_two_j(li)+1) ) < eps
                    if  (orbital_kappa(li) == naklj) then
                        if  (vli  .or. vlj ) then       ! at least one orbital varied
                            !                                                   but not both (varied and full)
                            if (.not. ((vli .and. vlj) .and. (fulli .and. fullj))) then
                                nec = nec + 1
                                !*** Encode index at 2nd round ***
                                if (itwice == 2) iecc(nec) = li + orbital_label_radix * lj
                            end if
                        end if
                    end if
                end do
            end do

            if (itwice==1 .and. nec>0) then
                call alloc (ecv, nec, "ECV", "SETLAG")
                call alloc (iecc, nec, "IECC", "SETLAG")
            else
                exit
            end if
        end do                                  !itwice

        !=======================================================================
        !   Print information about Lagrange multipliers
        !=======================================================================

        if (mpi_rank == 0) then
            if (nec == 0) then
                write (*, '(/,"Lagrange multipliers are not required")')
            else
                write (*, '(/,"Include Lagrange multipliers between:"/)')
                do li = 1, nec
                    !*** Decode index ***
                    ieccli = iecc(li)
                    l1 = ieccli/orbital_label_radix
                    l2 = ieccli - orbital_label_radix*l1
                    write (*, '(13x,2(2x,1i2,1a2))') &
                        orbital_principal_n(l2), orbital_symmetry_labels(l2), orbital_principal_n(l1), &
                        & orbital_symmetry_labels(l1)
                end do
            end if
        end if
        first = .false.
    end if

    !     Lagrange multipliers need to be computed also on the first call
    !     RETURN

    if (nec == 0) return
    if (mpi_rank == 0) write (*, '(/,"Lagrange multipliers:"/)')
    jlast = 0
    mlast = 0

    do li = 1, nec
        !*** Decode index ***
        ieccli = iecc(li)
        m = ieccli/orbital_label_radix
        j = ieccli - orbital_label_radix*m
        !
        if (j /= jlast) then
            ucfj = ucf(j)
            call assemble_orbital_potential_coefficients (extended_optimal_level, j)
            call build_direct_radial_potential (j)
            call build_exchange_radial_potential (j)
            call apply_off_diagonal_integrals
            ypj(:radial_point_count) = yp(:radial_point_count)
            xpj(:radial_point_count) = xp(:radial_point_count)
            xqj(:radial_point_count) = xq(:radial_point_count)
            jlast = j
        end if
        !
        if (m /= mlast) then
            ucfm = ucf(m)
            call assemble_orbital_potential_coefficients (extended_optimal_level, m)
            call build_direct_radial_potential (m)
            call build_exchange_radial_potential (m)
            call apply_off_diagonal_integrals
            ypm(:radial_point_count) = yp(:radial_point_count)
            xpm(:radial_point_count) = xp(:radial_point_count)
            xqm(:radial_point_count) = xq(:radial_point_count)
            mlast = m
        end if
        !
        mtp = max(orbital_last_grid_index(j),orbital_last_grid_index(m))
        !
        if (is_orbital_fixed(m)) then
            ta(1) = 0.0_real64
            do i = 2, mtp
                ta(i) = radial_derivative_over_radius(i)*((radial_large_components(i,m)*xqj(i) &
                    & -radial_small_components(i,m)*xpj(i))*speed_of_light_au+       &
                    (radial_large_components(i,m)*radial_large_components(i,j)+radial_small_components(i,m) &
                    & *radial_small_components(i,j))*ypj(i))
            end do

            call quad (result)
            ! note      ..rinti(m,j,1) is symmetric in (j,m)
            rijm = rinti(m,j,1)
            ecv(li) = (result - rijm)*ucfj

            ! start dbg
            !           WRITE (81,*)'1, RESULT, RIJM, UCFJ, ECV, TA' ! dbg
            !           WRITE (81,*)RESULT, RIJM, UCFJ, ECV ! dbg
            !           DO i = 1, MTP ! dbg
            !              WRITE (81,*) i, TA(i), r(i), rp(i) ! dbg
            !           ENDDO ! dbg
            ! end dbg

        else if (is_orbital_fixed(j)) then
            ta(1) = 0.0_real64
            do i = 2, mtp
                ta(i) = radial_derivative_over_radius(i)*((radial_large_components(i,j)*xqm(i) &
                    & -radial_small_components(i,j)*xpm(i))*speed_of_light_au+      &
                    (radial_large_components(i,j)*radial_large_components(i,m)+radial_small_components(i,j) &
                    & *radial_small_components(i,m))*ypm(i))
            end do

            !start dbg
            !           DO i = 1, MTP
            !              WRITE (81,*) i, TA(i)
            !               WRITE (83,*) i, r(i), rp(i), rpor(i)
            !               WRITE (84,*) i, pf(i,j), qf(i,j)
            !               write(85,*) i,ypm(i)
            !              write(86,*)i,xpm(i),xqm(i)
            !           ENDDO
            ! end dbg
            call quad (result)

            rijm = rinti(j,m,1)                  !/ nprocs
            ecv(li) = (result - rijm)*ucfm
            !start dbg
            !           WRITE (81,*)'2, RESULT, RIJM, UCFM, ECV, TA'
            !           WRITE (81,*)RESULT, RIJM, UCFJ, ECV, r(i), rp(i)
            !end dbg

        else

            qdif = abs((ucfj - ucfm)/max(ucfj,ucfm))

            if (qdif > p001) then
                obqdif = 1.0_real64/ucfj - 1.0_real64/ucfm
                ta(1) = 0.0_real64
                do i = 2, mtp
                    ta(i) = radial_derivative_over_radius(i)*((radial_large_components(i,m)*xqj(i) &
                        & -radial_small_components(i,m)*xpj(i)          &
                        -radial_large_components(i,j)*xqm(i)+radial_small_components(i,j)*xpm(i))* &
                        & speed_of_light_au       &
                        +(ypj(i)-ypm(i))*(radial_large_components(i,m)*radial_large_components(i,j)+ &
                        & radial_small_components(i,m)*radial_small_components(i,j)))
                end do

                call quad (result)
                ecv(li) = result/obqdif

                !start dbg
                !           WRITE (81,*)'3, RESULT, OBQDIF, ECV, TA'
                !           WRITE (81,*)RESULT, OBQDIF, ECV
                !           DO i = 1, MTP
                !              WRITE (81,*) i, TA(i), r(i), rp(i)
                !           ENDDO
                !end dbg

            else

                obqsum = 1.0_real64/ucfj + 1.0_real64/ucfm

                ta(1) = 0.0_real64
                do i = 2, mtp
                    ta(i) = radial_derivative_over_radius(i)*((radial_large_components(i,m)*xqj(i) &
                        & -radial_small_components(i,m)*xpj(i)         &
                        +radial_large_components(i,j)*xqm(i)-radial_small_components(i,j)*xpm(i))* &
                        & speed_of_light_au      &
                        +(ypj(i)+ypm(i))*(radial_large_components(i,m)*radial_large_components(i,j)+ &
                        & radial_small_components(i,m)*radial_small_components(i,j)))
                end do

                call quad (result)
                rijm = rinti(m,j,1)
                ecv(li) = (result - 2.0_real32*rijm)/obqsum


                !start dbg
                !           WRITE (81,*)'4, RESULT, RIUJM, OBQSUM, ECV, TA'
                !           WRITE (81,*)RESULT, RIUJM, OBQSUM, ECV
                !           DO i = 1, MTP
                !              WRITE (81,*) i, TA(i), r(i), rp(i)
                !           ENDDO
                !end dbg

            end if
        end if

        write (*, '(13x,2(2x,1i2,1a2),2x,1p,d16.9)') &
            orbital_principal_n(j), orbital_symmetry_labels(j), orbital_principal_n(m), &
            & orbital_symmetry_labels(m), ecv(li)

    end do


    return
end subroutine calculate_lagrange_multipliers
