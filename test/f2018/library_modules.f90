! The public procedures support procedure arguments, pointers and keyword calls.
program library_modules
    use, intrinsic :: iso_fortran_env, only: int8, real64
    use grasp_parameters, only: max_radial_points
    use atomic_calculation_state, only: working_tolerance
    use radial_node_count_settings, only: relative_extremum_cutoff
    use radial_node_count, only: count_radial_nodes
    use csf_byte_packing, only: pack_csf_subshell_byte
    implicit none (type, external)

    procedure(count_radial_nodes), pointer :: count_nodes
    procedure(pack_csf_subshell_byte), pointer :: pack_shell

    count_nodes => count_radial_nodes
    pack_shell => pack_csf_subshell_byte
    if (.not. associated(count_nodes, count_radial_nodes)) error stop "Radial node count procedure association"
    if (.not. associated(pack_shell, pack_csf_subshell_byte)) error stop "CSF byte packing procedure association"
    call check_count(count_nodes)
    call check_pack(pack_shell)
contains
    subroutine check_count(entry)
        procedure(count_radial_nodes) :: entry
        real(real64) :: radial_values(max_radial_points), first_sign
        integer :: node_count

        radial_values = 0.0_real64
        radial_values(:3) = [1.0_real64, -1.0_real64, 1.0_real64]
        working_tolerance = 1.0e-10_real64
        relative_extremum_cutoff = 0.1_real64
        call entry(radial_values=radial_values, last_radial_index=3, node_count=node_count, &
            first_extremum_sign=first_sign)
        if (node_count /= 2 .or. first_sign /= 1.0_real64) error stop "Radial node count keyword call"
    end subroutine check_count

    subroutine check_pack(entry)
        procedure(pack_csf_subshell_byte) :: entry
        integer(int8) :: shells(3)

        shells = 0_int8
        call entry(unpacked_value=3, byte_index=2, packed_values=shells)
        if (any(shells /= [0_int8, 3_int8, 0_int8])) error stop "CSF byte packing keyword call"
    end subroutine check_pack
end program library_modules
