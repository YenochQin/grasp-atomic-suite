module calculate_lagrange_multipliers_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_lagrange_multipliers_mpi
contains

    !***********************************************************************
    !                                                                      *
    subroutine calculate_lagrange_multipliers_mpi(extended_optimal_level)
        ! Module dependencies
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
        use grasp_parameters,    only: orbital_index_radix, max_radial_points
        use array_memory_management, only: alloc
        use scf_orbital_order, only: orbital_update_order
        use closed_core_subshells, only: closed_core_subshell_count
        use atomic_calculation_state, only: working_tolerance, max_radial_points, speed_of_light_au
        use fixed_orbital_flags, only: is_orbital_fixed
        use radial_grid_state, only: radial_point_count, radial_derivative_over_radius
        use orbital_lagrange_multipliers, only: lagrange_multiplier_values, lagrange_orbital_pair_labels, lagrange_multiplier_count
        use orbital_csf_basis, only: orbital_count, orbital_kappa, orbital_principal_n, orbital_symmetry_labels, orbital_two_j
        use mpi_state, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_process_count, &
            mpi_rank, mpi_sum
        use scf_radial_potentials, only: orbital_y_potential, large_exchange_source, small_exchange_source
        use scf_iteration_state, only: orbital_occupation_weights
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components

        ! Procedure interfaces
        use assemble_orbital_potential_coefficients_interface, only: assemble_orbital_potential_coefficients
        use build_direct_radial_potential_interface, only: build_direct_radial_potential
        use build_exchange_radial_potential_interface, only: build_exchange_radial_potential
        use apply_off_diagonal_integrals_interface, only: apply_off_diagonal_integrals

        ! Module dependencies
        use radial_quadrature, only: integrate_radial_grid
        use one_electron_radial_integral, only: integrate_one_electron_radial
        implicit none (type, external)

        ! Arguments
        logical  :: extended_optimal_level

        ! Local constants
        real(real64), parameter :: relative_occupation_difference_threshold = 1.0e-1_real64

        ! Local declarations
        integer :: itwice, orbital_index, li, lip1, nakli, ljraw, lj, ieccli, l1, l2, &
            jlast, mlast, m, j, i
        real(real64) :: ypj(max_radial_points), ypm(max_radial_points), xpj(max_radial_points), xpm(max_radial_points), &
            xqj(max_radial_points), xqm(max_radial_points)
        real(real64) :: &
            & eps,orbital_j_occupation,orbital_m_occupation,result,one_electron_radial_integral, &
            relative_occupation_difference, inverse_occupation_difference, inverse_occupation_sum, tmp
        logical :: first, fixli, fixlj, fulli, fullj

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

            eps = working_tolerance*0.01_real64      ! criterion to see if an orb is occupied
            do itwice = 1, 2
                lagrange_multiplier_count = 0
                !            IF (ITWICE /= 2) THEN
                !               DO LIRAW = 1, NW - 1
                !                  LI = IORDER(LIRAW)
                !                  LIP1 = MAX(NCORE,LIRAW) + 1
                !                  NAKLI = NAK(LI)
                !                  FIXLI = LFIX(LI)
                !                  FULLI = ABS(UCF(LI)-DBLE(NKJ(LI)+1)) < EPS
                !                  DO LJRAW = LIP1, NW
                !                     LJ = IORDER(LJRAW)
                !                     FIXLJ = LFIX(LJ)
                !                     FULLJ = ABS(UCF(LJ)-DBLE(NKJ(LJ)+1)) < EPS
                !                     IF (.NOT.(NAK(LJ)==NAKLI .AND. .NOT.(FIXLI .AND. FIXLJ)&
                !                         .AND. .NOT.(FULLI .AND. FULLJ))) CYCLE
                !                     NEC = NEC + 1
                !                     CYCLE
                !                  !*** Encode index at 2nd round ***
                !                  END DO
                !               END DO
                !            ELSE
                !               DO LIRAW = 1, NW - 1
                !                  LI = IORDER(LIRAW)
                !                  LIP1 = MAX(NCORE,LIRAW) + 1
                !                  NAKLI = NAK(LI)
                !                  FIXLI = LFIX(LI)
                !                  FULLI = ABS(UCF(LI)-DBLE(NKJ(LI)+1)) < EPS
                !                  !*** Encode index at 2nd round ***
                !                  DO LJRAW = LIP1, NW
                !                     LJ = IORDER(LJRAW)
                !                     FIXLJ = LFIX(LJ)
                !                     FULLJ = ABS(UCF(LJ)-DBLE(NKJ(LJ)+1)) < EPS
                !                     IF (.NOT.(NAK(LJ)==NAKLI .AND. .NOT.(FIXLI .AND. FIXLJ)&
                !                         .AND. .NOT.(FULLI .AND. FULLJ))) CYCLE
                !                     NEC = NEC + 1
                !                  !*** Encode index at 2nd round ***
                !                     IECC(NEC) = LI + KEY*LJ
                !                  END DO
                !               END DO
                !            ENDIF
                do orbital_index = 1, orbital_count - 1
                    li = orbital_update_order(orbital_index)
                    lip1 = max (closed_core_subshell_count, orbital_index) + 1
                    nakli = orbital_kappa(li)
                    fixli = is_orbital_fixed(li)
                    fulli = abs ( orbital_occupation_weights(li)-dble (orbital_two_j(li)+1) ) < eps
                    do ljraw = lip1, orbital_count
                        lj = orbital_update_order(ljraw)
                        fixlj = is_orbital_fixed(lj)
                        fullj = abs ( orbital_occupation_weights(lj)-dble (orbital_two_j(lj)+1) ) < eps
                        if ( (orbital_kappa(lj) == nakli) .and.            &
                            (.not. (fixli .and. fixlj)) .and.    &
                            (.not. (fulli .and. fullj)) ) then
                            lagrange_multiplier_count = lagrange_multiplier_count + 1
                            !*** Encode index at 2nd round ***
                            if (itwice == 2) lagrange_orbital_pair_labels(lagrange_multiplier_count) = li + orbital_index_radix * lj
                        end if
                    end do
                end do


                if (itwice==1 .and. lagrange_multiplier_count>0) then
                    call alloc (lagrange_multiplier_values, lagrange_multiplier_count, "lagrange_multiplier_values", &
                        "calculate_lagrange_multipliers_mpi")
                    call alloc (lagrange_orbital_pair_labels, lagrange_multiplier_count, "lagrange_orbital_pair_labels", &
                        "calculate_lagrange_multipliers_mpi")
                else
                    exit
                end if
            end do                                  !itwice

            !=======================================================================
            !   Print information about Lagrange multipliers
            !=======================================================================

            if (mpi_rank == 0) then
                if (lagrange_multiplier_count == 0) then
                    write (*, '(/,"Lagrange multipliers are not required")')
                else
                    write (*, '(/,"Include Lagrange multipliers between:"/)')
                    do li = 1, lagrange_multiplier_count
                        !*** Decode index ***
                        ieccli = lagrange_orbital_pair_labels(li)
                        l1 = ieccli/orbital_index_radix
                        l2 = ieccli - orbital_index_radix*l1
                        write (*, '(13x,2(2x,1i2,1a2))') &
                            orbital_principal_n(l2), orbital_symmetry_labels(l2), orbital_principal_n(l1), &
                            & orbital_symmetry_labels(l1)
                    end do
                end if
            end if
            first = .false.
        end if

        !FF+GG  12/07/05
        !     Lagrange multipliers need to be computed also on the first call
        !     RETURN

        if (lagrange_multiplier_count == 0) return
        if (mpi_rank == 0) write (*, '(/,"Lagrange multipliers:"/)')
        jlast = 0
        mlast = 0

        do li = 1, lagrange_multiplier_count
            !*** Decode index ***
            ieccli = lagrange_orbital_pair_labels(li)
            m = ieccli/orbital_index_radix
            j = ieccli - orbital_index_radix*m
            !
            if (j /= jlast) then
                orbital_j_occupation = orbital_occupation_weights(j)
                call assemble_orbital_potential_coefficients (extended_optimal_level, j)
                call build_direct_radial_potential (j)
                call build_exchange_radial_potential (j)
                call apply_off_diagonal_integrals
                ypj(:radial_point_count) = orbital_y_potential(:radial_point_count)
                xpj(:radial_point_count) = large_exchange_source(:radial_point_count)
                xqj(:radial_point_count) = small_exchange_source(:radial_point_count)
                jlast = j
            end if
            !
            if (m /= mlast) then
                orbital_m_occupation = orbital_occupation_weights(m)
                call assemble_orbital_potential_coefficients (extended_optimal_level, m)
                call build_direct_radial_potential (m)
                call build_exchange_radial_potential (m)
                call apply_off_diagonal_integrals
                ypm(:radial_point_count) = orbital_y_potential(:radial_point_count)
                xpm(:radial_point_count) = large_exchange_source(:radial_point_count)
                xqm(:radial_point_count) = small_exchange_source(:radial_point_count)
                mlast = m
            end if
            !
            integration_last_index = max(orbital_last_grid_index(j),orbital_last_grid_index(m))
            !
            if (is_orbital_fixed(m)) then
                radial_workspace_a(1) = 0.0_real64
                do i = 2, integration_last_index
                    radial_workspace_a(i) = radial_derivative_over_radius(i)*((radial_large_components(i,m)*xqj(i) &
                        & -radial_small_components(i,m)*xpj(i))*speed_of_light_au+(radial_large_components(i,m) &
                        & *radial_large_components(i&
                        ,j)+radial_small_components(i,m)*radial_small_components(i,j))*ypj(i))
                end do

                call integrate_radial_grid (result)
                one_electron_radial_integral = integrate_one_electron_radial(m,j,1)
                lagrange_multiplier_values(li) = (result - one_electron_radial_integral / mpi_process_count)*orbital_j_occupation
                ! start dbg
                !           WRITE (81,*)'1, RESULT, RIJM, UCFJ, ECV, TA' ! dbg
                !           WRITE (81,*)RESULT, RIJM, UCFJ, ECV ! dbg
                !           DO i = 1, MTP ! dbg
                !              WRITE (81,*) i, TA(i), r(i), rp(i) ! dbg
                !           ENDDO ! dbg
                ! end dbg

            else if (is_orbital_fixed(j)) then
                radial_workspace_a(1) = 0.0_real64
                do i = 2, integration_last_index
                    radial_workspace_a(i) = radial_derivative_over_radius(i)*((radial_large_components(i,j)*xqm(i) &
                        & -radial_small_components(i,j)*xpm(i))*speed_of_light_au+(radial_large_components(i,j) &
                        & *radial_large_components(i&
                        ,m)+radial_small_components(i,j)*radial_small_components(i,m))*ypm(i))
                    !GG             write(222,*)"XQM(I)",XQM(I)
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
                call integrate_radial_grid (result)

                one_electron_radial_integral = integrate_one_electron_radial(j,m,1)                  !/ nprocs
                lagrange_multiplier_values(li) = (result - one_electron_radial_integral / mpi_process_count)*orbital_m_occupation
                !start dbg
                !           WRITE (81,*)'2, RESULT, RIJM, UCFM, ECV, TA'
                !           WRITE (81,*)RESULT, RIJM, UCFJ, ECV, r(i), rp(i)
                !end dbg


            else
                relative_occupation_difference = abs((orbital_j_occupation - &
                    & orbital_m_occupation)/max(orbital_j_occupation,orbital_m_occupation))
                if (relative_occupation_difference > relative_occupation_difference_threshold) then
                    inverse_occupation_difference = 1.0_real64/orbital_j_occupation - 1.0_real64/orbital_m_occupation
                    radial_workspace_a(1) = 0.0_real64
                    do i = 2, integration_last_index
                        radial_workspace_a(i) = radial_derivative_over_radius(i)*((radial_large_components(i,m)*xqj(i) &
                            & -radial_small_components(i,m)*xpj(i)-radial_large_components(i,j)*xqm(i&
                            )+radial_small_components(i,j)*xpm(i))*speed_of_light_au+(ypj(i)-ypm(i))* &
                            & (radial_large_components(i,m)*radial_large_components(i,j)+radial_small_components(i,&
                            m)*radial_small_components(i,j)))
                    end do

                    call integrate_radial_grid (result)
                    lagrange_multiplier_values(li) = result/inverse_occupation_difference
                    !start dbg
                    !           WRITE (81,*)'3, RESULT, OBQDIF, ECV, TA'
                    !           WRITE (81,*)RESULT, OBQDIF, ECV
                    !           DO i = 1, MTP
                    !              WRITE (81,*) i, TA(i), r(i), rp(i)
                    !           ENDDO
                    !end dbg


                else
                    inverse_occupation_sum = 1.0_real64/orbital_j_occupation + 1.0_real64/orbital_m_occupation
                    radial_workspace_a(1) = 0.0_real64
                    do i = 2, integration_last_index
                        radial_workspace_a(i) = radial_derivative_over_radius(i)*((radial_large_components(i,m)*xqj(i) &
                            & -radial_small_components(i,m)*xpj(i)+radial_large_components(i,j)*xqm(i&
                            )-radial_small_components(i,j)*xpm(i))*speed_of_light_au+(ypj(i)+ypm(i))* &
                            & (radial_large_components(i,m)*radial_large_components(i,j)+radial_small_components(i,&
                            m)*radial_small_components(i,j)))
                    end do

                    call integrate_radial_grid (result)
                    one_electron_radial_integral = integrate_one_electron_radial(m,j,1)               !/ nprocs
                    lagrange_multiplier_values(li) = (result - 2.0_real64*one_electron_radial_integral / &
                        mpi_process_count)/inverse_occupation_sum
                    !start dbg
                    !           WRITE (81,*)'4, RESULT, RIUJM, OBQSUM, ECV, TA'
                    !           WRITE (81,*)RESULT, RIUJM, OBQSUM, ECV
                    !           DO i = 1, MTP
                    !              WRITE (81,*) i, TA(i), r(i), rp(i)
                    !           ENDDO
                    !end dbg

                end if
            end if
            !=======================================================================
            !  Collect contributions from all nodes.
            !  Another alternative is to modify mcpmpi to let every node
            !  have the same set of "LAB" and "NCONTR". So an easier aggregation
            !  can be done for arrays DA, XA, YA in setcof. And so the computation
            !  of YPOT, XPOT, DACON can be distributed.
            !=======================================================================
            call mpi_allreduce (lagrange_multiplier_values(li),tmp,1,mpi_double_precision, &
                mpi_sum, mpi_comm_world, mpi_error_code)
            lagrange_multiplier_values(li) = tmp

            if(mpi_rank == 0)write (*, '(13x,2(2x,1i2,1a2),2x,1p,d16.9)') &
                orbital_principal_n(j), orbital_symmetry_labels(j), orbital_principal_n(m), &
                & orbital_symmetry_labels(m), lagrange_multiplier_values(li)

        end do

        !db      close(81)
        !db      close(82)



        return
    end subroutine calculate_lagrange_multipliers_mpi
end module calculate_lagrange_multipliers_mpi_interface
