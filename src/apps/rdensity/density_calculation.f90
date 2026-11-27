! Accumulate one-body densities by level and form statistically averaged natural orbitals.
! MPI/serial migration of GRASP rdensity (upstream 488e875d).
! Accumulate one-body density matrices before contracting radial functions.
module density_calculation
    ! Module dependencies
    use calculation_state_files, only: natural_orbital_path
    use grasp_parameters, only: max_orbitals
    use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count, &
        orbital_kappa, orbital_principal_n
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, &
        final_orbital_last_grid_index, final_orbital_origin_parameters, initial_orbital_last_grid_index, &
        initial_orbital_origin_parameters, radial_large_components, radial_small_components, final_radial_large_components, &
        final_radial_small_components, initial_radial_large_components, initial_radial_small_components
    use radial_grid_state, only: radial_point_count, radial_grid
    use eigenstate_expansions, only: csf_mixing_coefficients
    use selected_level_indices, only: active_level_count, active_level_indices
    use level_symmetry_metadata, only: level_two_j_plus_one, level_parity_signs
    use debug_output_settings, only: matrix_element_screening_cutoff
    use angular_symmetry_labels, only: format_j_label, parity_labels
    use density_parallel, only: density_fail, mpi_rank, sum_density, mpi_process_count
    use natural_orbitals, only: diagonalize_density, transform_radial_block
    use one_electron_scalar_coefficients, only: calculate_one_electron_scalar_coefficients
    use csf_quantum_queries, only: csf_total_j_code, csf_parity_code
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_density
contains
    recursive subroutine calculate_density(name, mixing_file_kind, ordering)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        character(*), intent(in) :: name
        integer, intent(in) :: mixing_file_kind, ordering

        ! Local declarations
        real(real64), allocatable :: rho(:,:,:)
        real(real64) :: shell_coefficients(max_orbitals), coefficient
        integer, allocatable :: csf_two_j_plus_one(:), parity(:)
        integer :: csf_index, row_csf_index, bra_orbital_index, ket_orbital_index, orbital_index, state, offset, &
            & differing_subshell_count, occupation_difference

        allocate(rho(orbital_count,orbital_count,active_level_count), csf_two_j_plus_one(csf_count), parity(csf_count))
        rho = 0.0_real64
        do csf_index = 1,csf_count
            csf_two_j_plus_one(csf_index) = csf_total_j_code(csf_index)
            parity(csf_index) = csf_parity_code(csf_index)
        end do
        ! Same cyclic column ownership as RHFS/RIS. Each process has its
        ! own angular-library workspace, including ranks with no columns.
        do csf_index = mpi_rank+1,csf_count,mpi_process_count
            if (mpi_rank == 0 .and. mod((csf_index-1)/mpi_process_count,100) == 0) &
                write(*,"(a,i0)") "Processing root CSF column: ", csf_index
            do row_csf_index = csf_index,csf_count
                if (csf_two_j_plus_one(csf_index) /= csf_two_j_plus_one(row_csf_index) .or. parity(csf_index) /= &
                    & parity(row_csf_index)) cycle
                ! Exact one-body selection rules, also used by RIS_MPI.
                differing_subshell_count = 0
                do orbital_index = 1,orbital_count
                    occupation_difference = int(csf_subshell_occupations(orbital_index,csf_index))- &
                        & int(csf_subshell_occupations(orbital_index,row_csf_index))
                    if (occupation_difference == 0) cycle
                    differing_subshell_count = differing_subshell_count+1
                    if (abs(occupation_difference) > 1 .or. differing_subshell_count > 2) then
                        differing_subshell_count = 3
                        exit
                    end if
                end do
                if (differing_subshell_count /= 0 .and. differing_subshell_count /= 2) cycle
                if (differing_subshell_count == 0 .and. csf_index /= row_csf_index) cycle
                call calculate_one_electron_scalar_coefficients(csf_index, row_csf_index, bra_orbital_index, ket_orbital_index, &
                    & shell_coefficients)
                if (bra_orbital_index == 0) cycle
                if (bra_orbital_index == ket_orbital_index) then
                    do orbital_index = 1,orbital_count
                        if (abs(shell_coefficients(orbital_index)) <= matrix_element_screening_cutoff) cycle
                        do state = 1,active_level_count
                            offset = (state-1)*csf_count
                            coefficient = csf_mixing_coefficients(csf_index+offset)*csf_mixing_coefficients(row_csf_index+offset)* &
                                & shell_coefficients(orbital_index)
                            rho(orbital_index,orbital_index,state) = rho(orbital_index,orbital_index,state)+coefficient
                        end do
                    end do
                else
                    if (orbital_kappa(bra_orbital_index) /= orbital_kappa(ket_orbital_index) .or. abs(shell_coefficients(1)) <= &
                        & matrix_element_screening_cutoff) cycle
                    do state = 1,active_level_count
                        offset = (state-1)*csf_count
                        coefficient = csf_mixing_coefficients(csf_index+offset)*csf_mixing_coefficients(row_csf_index+offset)* &
                            & shell_coefficients(1)
                        rho(bra_orbital_index,ket_orbital_index,state) = rho(bra_orbital_index,ket_orbital_index,state)+coefficient
                        rho(ket_orbital_index,bra_orbital_index,state) = rho(ket_orbital_index,bra_orbital_index,state)+coefficient
                    end do
                end if
            end do
        end do
        do state = 1,active_level_count
            call sum_density(rho(:,:,state))
        end do
        if (mpi_rank == 0) then
            ! Use the original radial basis before rotating PF and QF.
            call write_density(name, mixing_file_kind, rho)
            call make_natural_orbitals(natural_orbital_path(name, mixing_file_kind == 0), ordering, rho)
        end if
    end subroutine calculate_density

    recursive subroutine write_density(name, mixing_file_kind, rho)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        character(*), intent(in) :: name
        integer, intent(in) :: mixing_file_kind
        real(real64), intent(in) :: rho(:,:,:)

        ! Local declarations
        real(real64), allocatable :: radial(:)
        integer :: unit, io_status, state, orbital_index, b, radial_index
        character(256) :: filename

        filename = trim(name)//".d"
        if (mixing_file_kind == 0) filename = trim(name)//".cd"
        open(newunit=unit, file=filename, status="new", action="write", iostat=io_status)
        if (io_status /= 0) call density_fail("Cannot create "//trim(filename)//" (remove an existing output first)")
        allocate(radial(radial_point_count))
        write(unit,*) "     r [au]          D(r)=r^2*rho(r)       rho(r)"
        write(unit,*)
        do state = 1,active_level_count
            radial = 0.0_real64
            do orbital_index = 1,orbital_count
                radial = radial+rho(orbital_index,orbital_index,state)*(radial_large_components(:radial_point_count, &
                    & orbital_index)**2+ &
                    & radial_small_components(:radial_point_count,orbital_index)**2)
                do b = orbital_index+1,orbital_count
                    if (orbital_kappa(orbital_index) /= orbital_kappa(b)) cycle
                    radial = radial+2.0_real64*rho(orbital_index,b,state)*(radial_large_components(:radial_point_count, &
                        & orbital_index)* &
                        & radial_large_components(:radial_point_count,b)+radial_small_components(:radial_point_count, &
                        & orbital_index)* &
                        & radial_small_components(:radial_point_count,b))
                end do
            end do
            write(unit,"(1x,i3,5x,2a)") active_level_indices(state), format_j_label(level_two_j_plus_one(state)), &
                parity_labels((level_parity_signs(state)+3)/2)
            do radial_index = 2,radial_point_count
                if (radial(radial_index) > 0.0_real64) write(unit,"(3es20.10e3)") radial_grid(radial_index), &
                    & radial(radial_index), radial(radial_index)/radial_grid(radial_index)**2
            end do
        end do
        close(unit, iostat=io_status)
        if (io_status /= 0) call density_fail("Error writing "//filename)
    end subroutine write_density

    recursive subroutine make_natural_orbitals(path, ordering, rho)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_wavefunction_format, only: write_radial_orbitals

        ! Arguments
        character(*), intent(in) :: path
        integer, intent(in) :: ordering
        real(real64), intent(in) :: rho(:,:,:)

        ! Local declarations
        real(real64), allocatable :: matrix(:,:), rotation(:,:), occupation(:)
        real(real64), allocatable :: p(:,:), q(:,:), origin(:)
        integer, allocatable :: orbits(:)
        integer :: kappa, ndim, a, b, orbital_index, state, points, unit, io_status, temp
        real(real64) :: weight
        character(160) :: error

        weight = sum(real(level_two_j_plus_one(:active_level_count), real64))
        do kappa = minval(orbital_kappa(:orbital_count)),maxval(orbital_kappa(:orbital_count))
            ndim = count(orbital_kappa(:orbital_count) == kappa)
            if (ndim == 0) cycle
            allocate(orbits(ndim), matrix(ndim,ndim), rotation(ndim,ndim), occupation(ndim))
            orbits = pack([(orbital_index,orbital_index=1,orbital_count)], orbital_kappa(:orbital_count) == kappa)
            ! Enumerate actual orbitals, so gaps in n and n>20 are supported.
            do a = 2,ndim
                b = a
                do while (b > 1)
                    if (orbital_principal_n(orbits(b-1)) <= orbital_principal_n(orbits(b))) exit
                    temp = orbits(b-1)
                    orbits(b-1) = orbits(b)
                    orbits(b) = temp
                    b = b-1
                end do
            end do
            matrix = 0.0_real64
            do state = 1,active_level_count
                matrix = matrix+real(level_two_j_plus_one(state),real64)*rho(orbits,orbits,state)/weight
            end do
            write(*,"(a,i5)") "KAPPA: ", kappa
            write(*,"(a,i5)") "NDIM: ", ndim
            write(*,*) "density matrix:"
            do a = 1,ndim
                write(*,"(i5,*(1x,es24.16))") orbital_principal_n(orbits(a)), matrix(:,a)
            end do
            call diagonalize_density(matrix, ordering, occupation, rotation, error)
            if (len_trim(error) /= 0) call density_fail(trim(error))
            points = maxval(orbital_last_grid_index(orbits))
            allocate(p(points,ndim), q(points,ndim), origin(ndim))
            p = radial_large_components(:points,orbits)
            q = radial_small_components(:points,orbits)
            origin = orbital_origin_parameters(orbits)
            call transform_radial_block(rotation, p, q, origin)
            radial_large_components(:points,orbits) = p
            radial_small_components(:points,orbits) = q
            orbital_origin_parameters(orbits) = origin
            orbital_last_grid_index(orbits) = points
            write(*,*) "occupations in output orbital order:"
            do a = 1,ndim
                write(*,"(i5,1x,es24.16)") orbital_principal_n(orbits(a)), occupation(a)
            end do
            write(*,*) "eigenvectors in output orbital order (columns printed as rows):"
            do a = 1,ndim
                write(*,"(i5,*(1x,es24.16))") orbital_principal_n(orbits(a)), rotation(:,a)
            end do
            deallocate(orbits, matrix, rotation, occupation, p, q, origin)
        end do
        open(newunit=unit, file=path, access="stream", form="unformatted", status="replace", &
            action="write", iostat=io_status)
        if (io_status /= 0) call density_fail("Cannot create natural-orbital file: "//path)
        call write_radial_orbitals(unit)
        close(unit, iostat=io_status)
        if (io_status /= 0) call density_fail("Error writing natural-orbital file: "//path)
        write(*,*) "Natural orbitals written to ", path
    end subroutine make_natural_orbitals
end module density_calculation
