! Public names that coincide with intrinsics must remain usable as procedures.
program library_modules
    use, intrinsic :: iso_fortran_env, only: int8
    use radial_node_count, only: count
    use csf_byte_packing, only: pack
    implicit none (type, external)
    procedure(count), pointer :: count_nodes
    procedure(pack), pointer :: pack_shell
    integer(int8) :: shells(3)

    count_nodes => count
    pack_shell => pack
    if (.not. associated(count_nodes, count)) error stop "COUNT procedure association"
    shells = 0_int8
    call pack_shell(iunpkd=3, isubsh=2, ipackd=shells)
    if (any(shells /= [0_int8, 3_int8, 0_int8])) error stop "PACK keyword/procedure contract"
end program library_modules
