! Suite radial wavefunction streams.
module radial_wavefunction_format
    use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real64, file_storage_size, iostat_end
    use, intrinsic :: ieee_arithmetic, only: ieee_support_datatype
    implicit none (type, external)
    private
    public :: wavefunction_layout, write_radial_orbitals, open_wavefunction_input, read_wavefunction_header, &
        read_orbital_metadata, read_orbital_payload, write_radial_orbital_set

    character(8), parameter :: stream_signature = "RAD--WFN"
    integer, parameter :: header_bytes = 16, orbital_bytes = 32

    type :: wavefunction_layout
        private
        integer(int8), allocatable :: orbitals(:,:)
        real(real64), allocatable :: radial(:,:)
        integer :: orbital = 0
    end type
contains
    logical function supports_stream_encoding()
        integer(int8) :: byte_order(4)

        byte_order = transfer(1_int32, byte_order)
        supports_stream_encoding = file_storage_size == 8 .and. ieee_support_datatype(0.0_real64) .and. &
            all(byte_order == [1_int8, 0_int8, 0_int8, 0_int8])
    end function

    subroutine write_radial_orbitals(unit)
        use radial_grid_state, only: radial_grid, radial_point_count
        use orbital_csf_basis, only: orbital_count, orbital_principal_n, orbital_kappa, orbital_energy_parameters
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        integer, intent(in) :: unit
        call write_radial_orbital_set(unit, orbital_count, orbital_principal_n, orbital_kappa, &
            orbital_energy_parameters, orbital_last_grid_index, orbital_origin_parameters, &
            radial_large_components, radial_small_components, radial_grid(:radial_point_count))
    end subroutine

    subroutine write_radial_orbital_set(unit, orbital_count, principal_n, kappa, energies, point_counts, &
            origins, large_components, small_components, grid)
        integer, intent(in) :: unit, orbital_count, principal_n(:), kappa(:), point_counts(:)
        real(real64), intent(in) :: energies(:), origins(:), large_components(:,:), small_components(:,:), grid(:)
        integer :: orbital, points, longest
        integer(int8), allocatable :: orbitals(:,:)
        real(real64), allocatable :: radial(:,:)

        if (.not. supports_stream_encoding()) error stop "Wavefunction stream requires little-endian IEEE binary64"
        if (orbital_count < 1) error stop "Cannot write an empty wavefunction"
        if (any(point_counts(:orbital_count) < 1)) error stop "Invalid orbital length"
        longest = maxval(point_counts(:orbital_count))
        if (longest > size(grid)) error stop "Orbital exceeds the radial grid"

        allocate(orbitals(orbital_bytes,orbital_count), radial(longest,2*orbital_count+1))
        radial = 0.0_real64
        radial(:,1) = grid(:longest)
        do orbital = 1, orbital_count
            points = point_counts(orbital)
            orbitals(:16,orbital) = transfer([int(principal_n(orbital),int32), &
                int(kappa(orbital),int32), int(points,int32), 0_int32], [0_int8])
            orbitals(17:,orbital) = transfer([energies(orbital), origins(orbital)], [0_int8])
            radial(:points,2*orbital) = large_components(:points,orbital)
            radial(:points,2*orbital+1) = small_components(:points,orbital)
        end do
        ! Byte arrays keep the wire encoding independent of unformatted I/O conversion settings.
        write(unit) stream_signature, transfer(int([orbital_count,longest],int32), [0_int8]), &
            orbitals, transfer(radial, [0_int8])
    end subroutine

    subroutine open_wavefunction_input(unit, path, status)
        integer, intent(out) :: unit, status
        character(*), intent(in) :: path
        unit = -1
        open(newunit=unit, file=path, status="old", action="read", form="unformatted", access="stream", iostat=status)
    end subroutine

    subroutine read_wavefunction_header(unit, layout, status)
        integer, intent(in) :: unit
        type(wavefunction_layout), intent(out) :: layout
        integer, intent(out) :: status
        character(8) :: signature
        integer(int8) :: dimensions_bytes(8)
        integer(int32) :: dimensions(2)
        integer :: points, count, channels
        integer(int64) :: file_bytes, radial_bytes, channel_row_bytes
        integer(int8), allocatable :: packed_radial(:)

        status = 1
        if (.not. supports_stream_encoding()) return
        read(unit, iostat=status) signature, dimensions_bytes
        if (status /= 0) return
        dimensions = transfer(dimensions_bytes, dimensions)
        count = dimensions(1)
        points = dimensions(2)
        status = 1
        if (signature /= stream_signature) return
        if (min(count,points) < 1 .or. count > (huge(0)-1)/2) return

        inquire(unit=unit, size=file_bytes, iostat=status)
        if (status /= 0) return
        channels = 2*count+1
        radial_bytes = file_bytes - header_bytes - int(orbital_bytes,int64)*count
        channel_row_bytes = 8_int64*channels
        status = 1
        if (radial_bytes < 0 .or. radial_bytes/channel_row_bytes /= points .or. &
            mod(radial_bytes,channel_row_bytes) /= 0) return

        allocate(layout%orbitals(orbital_bytes,count), layout%radial(points,channels), packed_radial(radial_bytes), stat=status)
        if (status /= 0) return
        read(unit, iostat=status) layout%orbitals, packed_radial
        if (status /= 0) return
        layout%radial = reshape(transfer(packed_radial, [0.0_real64]), [points,channels])
    end subroutine

    subroutine read_orbital_metadata(layout, n, kappa, energy, points, status)
        type(wavefunction_layout), intent(inout) :: layout
        integer, intent(out) :: n, kappa, points, status
        real(real64), intent(out) :: energy
        integer(int32) :: integers(4)

        status = iostat_end
        if (layout%orbital == size(layout%orbitals,2)) return
        layout%orbital = layout%orbital+1
        integers = transfer(layout%orbitals(:16,layout%orbital), integers)
        n = integers(1)
        kappa = integers(2)
        points = integers(3)
        energy = transfer(layout%orbitals(17:24,layout%orbital), energy)
        status = 1
        if (n < 1 .or. kappa == 0 .or. points < 1 .or. points > size(layout%radial,1) .or. integers(4) /= 0) return
        status = 0
    end subroutine

    subroutine read_orbital_payload(layout, origin, large_component, small_component, radii)
        type(wavefunction_layout), intent(in) :: layout
        real(real64), intent(out) :: origin, large_component(:), small_component(:), radii(:)
        integer :: points, orbital

        points = size(radii)
        orbital = layout%orbital
        origin = transfer(layout%orbitals(25:32,orbital), origin)
        radii = layout%radial(:points,1)
        large_component = layout%radial(:points,2*orbital)
        small_component = layout%radial(:points,2*orbital+1)
    end subroutine

end module radial_wavefunction_format
