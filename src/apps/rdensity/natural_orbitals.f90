! Density eigensystem and radial transformation, shared by serial/MPI targets.
! Algorithm derived from GRASP rdensity/natorbnew.f90 (488e875d).
module natural_orbitals
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none (type, external)

    ! Module interface
    private
    public :: diagonalize_density, transform_radial_block

    ! Procedure interfaces
    interface
        recursive subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(real64), intent(inout) :: a(lda,*), work(*)
            real(real64), intent(out) :: w(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface
contains
    recursive subroutine diagonalize_density(matrix, ordering, natural_occupations, natural_orbital_rotation, error)
        ! Arguments
        real(real64), intent(in) :: matrix(:,:)
        integer, intent(in) :: ordering
        real(real64), intent(out) :: natural_occupations(:), natural_orbital_rotation(:,:)
        character(*), intent(out) :: error

        ! Local declarations
        real(real64), allocatable :: ascending_occupations(:), density_eigenvectors(:,:), diagonalization_workspace(:)
        real(real64) :: workspace_size_query(1)
        integer, allocatable :: ordered_eigenvector_indices(:)
        logical, allocatable :: eigenvector_assigned(:)
        integer :: orbital_count, lapack_status, eigenvector_index, next_unused_eigenvector, orbital_slot
        integer :: diagonalization_workspace_size

        error = ""
        orbital_count = size(matrix,1)
        if (orbital_count < 1 .or. size(matrix,2) /= orbital_count .or. &
            size(natural_occupations) /= orbital_count .or. &
            any(shape(natural_orbital_rotation) /= [orbital_count, orbital_count])) then
            error = "Invalid density block dimensions"
            return
        end if
        if (.not.all(ieee_is_finite(matrix))) then
            error = "Non-finite density matrix"
            return
        end if
        allocate(ascending_occupations(orbital_count), density_eigenvectors(orbital_count, orbital_count), &
            ordered_eigenvector_indices(orbital_count), eigenvector_assigned(orbital_count))
        density_eigenvectors = matrix
        call dsyev("V", "U", orbital_count, density_eigenvectors, orbital_count, &
            ascending_occupations, workspace_size_query, -1, lapack_status)
        if (lapack_status /= 0) then
            error = "DSYEV workspace query failed"
            return
        end if
        diagonalization_workspace_size = max(1, ceiling(workspace_size_query(1)))
        allocate(diagonalization_workspace(diagonalization_workspace_size))
        call dsyev("V", "U", orbital_count, density_eigenvectors, orbital_count, &
            ascending_occupations, diagonalization_workspace, diagonalization_workspace_size, lapack_status)
        if (lapack_status /= 0) then
            error = "DSYEV failed to diagonalize density matrix"
            return
        end if
        ordered_eigenvector_indices = 0
        eigenvector_assigned = .false.
        select case(ordering)
        case(1)
            do eigenvector_index = 1,orbital_count
                orbital_slot = maxloc(abs(density_eigenvectors(:,eigenvector_index)), dim=1)
                if (ordered_eigenvector_indices(orbital_slot) /= 0) then
                    error = "Two eigenvectors have the same dominant component; use --order=2"
                    return
                end if
                ordered_eigenvector_indices(orbital_slot) = eigenvector_index
            end do
        case(2)
            ! Preserve pure orbitals first, then fill every other slot in
            ! descending occupation order. Unlike upstream, no slot can be
            ! overwritten or left unassigned by mixed/pure eigenvectors.
            do eigenvector_index = 1,orbital_count
                orbital_slot = maxloc(abs(density_eigenvectors(:,eigenvector_index)), dim=1)
                if (abs(density_eigenvectors(orbital_slot,eigenvector_index)) < 1.0_real64-1.0e-12_real64) cycle
                ordered_eigenvector_indices(orbital_slot) = eigenvector_index
                eigenvector_assigned(eigenvector_index) = .true.
            end do
            next_unused_eigenvector = orbital_count
            do orbital_slot = 1,orbital_count
                if (ordered_eigenvector_indices(orbital_slot) /= 0) cycle
                do while (eigenvector_assigned(next_unused_eigenvector))
                    next_unused_eigenvector = next_unused_eigenvector-1
                end do
                ordered_eigenvector_indices(orbital_slot) = next_unused_eigenvector
                eigenvector_assigned(next_unused_eigenvector) = .true.
            end do
        case default
            error = "Orbital ordering must be 1 or 2"
            return
        end select
        do eigenvector_index = 1,orbital_count
            natural_occupations(eigenvector_index) = ascending_occupations(ordered_eigenvector_indices(eigenvector_index))
            natural_orbital_rotation(:,eigenvector_index) = density_eigenvectors(:,ordered_eigenvector_indices(eigenvector_index))
        end do
    end subroutine diagonalize_density

    recursive subroutine transform_radial_block(natural_orbital_rotation, large_components, small_components, &
            large_component_origin_coefficients)
        ! Arguments
        real(real64), intent(inout) :: natural_orbital_rotation(:,:), large_components(:,:), small_components(:,:)
        real(real64), intent(inout) :: large_component_origin_coefficients(:)

        ! Local declarations
        integer :: orbital_index

        large_components = matmul(large_components, natural_orbital_rotation)
        small_components = matmul(small_components, natural_orbital_rotation)
        large_component_origin_coefficients = matmul(large_component_origin_coefficients, natural_orbital_rotation)
        ! Fix the phase after ALL origin coefficients have been transformed.
        do orbital_index = 1,size(large_component_origin_coefficients)
            if (large_component_origin_coefficients(orbital_index) >= 0.0_real64) cycle
            large_components(:,orbital_index) = -large_components(:,orbital_index)
            small_components(:,orbital_index) = -small_components(:,orbital_index)
            large_component_origin_coefficients(orbital_index) = -large_component_origin_coefficients(orbital_index)
            natural_orbital_rotation(:,orbital_index) = -natural_orbital_rotation(:,orbital_index)
        end do
    end subroutine transform_radial_block
end module natural_orbitals
