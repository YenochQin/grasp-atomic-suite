module transform_orbital_sets_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: transform_orbital_sets
contains

    !***********************************************************************
    !                                                                      *
    subroutine transform_orbital_sets(initial_large_components, initial_small_components, initial_orbitals_by_symmetry, &
                & initial_inactive_by_symmetry, final_large_components, final_small_components, final_orbitals_by_symmetry, &
                & final_inactive_by_symmetry, &
        radial_point_count, symmetry_count, workspace, workspace_size, global_trace_level, initial_orbital_rotations, &
                & initial_ci_rotations, final_orbital_rotations, final_ci_rotations)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        ! Generate Matrices for rotating radial functions and                  *
        ! for counter rotating CI coefficients                                 *
        !                                                                      *
        ! This subroutine is a modification of the corresponding MCHF          *
        ! routine                                                              *
        !                                                                      *
        ! Per Jonsson,     Lund June 1996                                      *
        !                                                                      *
        ! This subroutine has been modified to support transformations         *
        ! in cases like 3s2p1d 3s2p                                            *
        !                                                                      *
        ! Per Jonsson,     Lund Feb  1997                                      *
        !                                                                      *
        ! =====                                                                *
        ! Input                                                                *
        ! =====                                                                *
        !                                                                      *
        ! PI,QI   : Input radial shells for initial state                      *
        !         : (large and small component)                                *
        ! NLI     : Number of shells per kappa for initial state               *
        ! NINSHLI : Number of inactive shells  per kappa for initial state     *
        !                                                                      *
        ! PF,QF   : Input radial shells for final   state                      *
        !         : (large and small component)                                *
        ! NLF     : Number of shells per kappa for final   state               *
        ! NINSHLF : Number of inactive shells  per kappa for final state       *
        !                                                                      *
        ! NGRID   : Number of gridpoints                                       *
        ! MXL     : Number of kappa quantum numbers                            *
        ! SCR     :                                                            *
        ! LSCR    : Total length of scratch space.                             *
        ! NTESTG  : Global PRINT flag : = 0 => complete silence                *
        !                               = 1 => test and PRINT overlap matrices,*
        !                                      PRINT header and  t matrix      *
        !                            .gt. 1 => Hope you now what you are doing *
        !                                                                      *
        ! ======                                                               *
        ! Output                                                               *
        ! ======                                                               *
        !                                                                      *
        ! CISHL : Rotation matrix for initial state shells                     *
        ! CICI  : Rotation matrix for initial state CI coefficients            *
        ! CFSHL : Rotation matrix for final state shells                       *
        ! CFCI  : Rotation matrix for final state CI coefficients              *
        ! PI,QI : Initial shells in biorthogonal basis                         *
        !         biorthogonal basis                                           *
        ! PF,QF : Final   shells in biorthogonal basis                         *
        !                                                                      *
        !                                                                      *
        ! Inactive Shells : The functions of the inactive shells must be       *
        !                   be supplied to the program, and NLI,NLF            *
        !                   must refer to the total number of                  *
        !                   occupied shells ( inactive+active)                 *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64

        ! Procedure interfaces
        use find_integer_extremum_interface, only: find_integer_extremum
        use sum_integer_vector_interface, only: sum_integer_vector
        !      USE gets_I
        !      USE wrtmat_I
        !      USE copvec_I
        !      USE invmat_I
        !      USE ulla_I
        !      USE trpmat_I
        !      USE matml4_I
        !      USE scalve_I
        !      USE setvec_I
        !      USE pamtmt_I
        use copy_vector_interface, only: copy_vector
        use calculate_subshell_overlap_interface, only: calculate_subshell_overlap
        use invert_dense_submatrix_interface, only: invert_dense_submatrix
        use multiply_submatrices_interface, only: multiply_submatrices
        use build_subshell_transformation_matrix_interface, only: build_subshell_transformation_matrix
        use scale_vector_interface, only: scale_vector
        use fill_vector_interface, only: fill_vector
        use transpose_matrix_interface, only: transpose_matrix
        use combine_upper_lower_factors_interface, only: combine_upper_lower_factors
        use write_matrix_interface, only: write_matrix
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: radial_point_count,  symmetry_count
        integer, intent(in) :: workspace_size,  global_trace_level
        integer, intent(in) :: initial_orbitals_by_symmetry(symmetry_count), final_orbitals_by_symmetry(symmetry_count)
        integer :: initial_inactive_by_symmetry(symmetry_count), final_inactive_by_symmetry(symmetry_count)
        !     The following have dimensions (NGRID:*)
        !GG      REAL(DOUBLE), DIMENSION(:,:), pointer :: PI, QI, PF, QF
        real(real64) :: initial_large_components(radial_point_count, *), initial_small_components(radial_point_count, *), &
                & final_large_components(radial_point_count, *), final_small_components(radial_point_count, *)
        !     The following have fixed dimensions (NGRID, NLMAX)
        real(real64)  :: initial_orbital_rotations(*), final_orbital_rotations(*)
        !     The work array
        !     REAL(DOUBLE), DIMENSION(LWORK1) :: work
        !     Th following have fixed dimensions (20*NLMAX*NLMAX)
        real(real64) :: initial_ci_rotations(*), final_ci_rotations(*)
        real(real64), intent(inout) :: workspace(*)

        ! Local declarations
        integer :: local_trace_level, trace_level, first_initial_orbital, first_final_orbital, max_initial_symmetry_size, &
                & max_final_symmetry_size, max_symmetry_size, total_initial_orbitals, &
            total_final_orbitals, next_workspace_index, initial_components_start, final_components_start, &
                & rotated_components_start, full_overlap_start, overlap_block_start, inverse_overlap_start, &
            initial_rotation_start, final_rotation_start, auxiliary_workspace_start, symmetry_index, initial_rotation_offset, &
                & final_rotation_offset, initial_orbital_count, final_orbital_count, paired_orbital_count, initial_orbital,&
            final_orbital, larger_components_start, smaller_components_start, larger_orbital_count, smaller_orbital_count, &
                & larger_rotation_start, smaller_rotation_start, unpaired_orbital_count, rotation_column, &
            initial_ci_rotation_start, orbital_index, final_ci_rotation_start
        real(real64) :: diagonal_element, inverse_diagonal_element
        !-----------------------------------------------
        !
        local_trace_level = 0
        trace_level = max(local_trace_level,global_trace_level)
        if (trace_level >= 1) then
            write (output_unit, *)
            write (output_unit, *) "                   *************************"
            write (output_unit, *) "                   *   Entering BIOTR1     *"
            write (output_unit, *) "                   *************************"
            write (output_unit, *)
        end if
        !
        !. Scratch should at least be of length
        !
        first_initial_orbital = 1
        first_final_orbital = 1
        !
        !. Largest number of shells of a given symmetry
        !
        max_initial_symmetry_size = find_integer_extremum(initial_orbitals_by_symmetry,symmetry_count,1)
        max_final_symmetry_size = find_integer_extremum(final_orbitals_by_symmetry,symmetry_count,1)
        max_symmetry_size = max(max_initial_symmetry_size,max_final_symmetry_size)
        if (trace_level >= 10) write (output_unit, *) " NLIMX,NLFMX NLIFMX ", max_initial_symmetry_size, max_final_symmetry_size, &
            max_symmetry_size
        !
        !. Total numner of shells
        !
        total_initial_orbitals = sum_integer_vector(initial_orbitals_by_symmetry,symmetry_count)
        total_final_orbitals = sum_integer_vector(final_orbitals_by_symmetry,symmetry_count)
        if (trace_level >= 10) write (output_unit, *) " NLTI NLTF", total_initial_orbitals, total_final_orbitals
        !
        ! Scratch space for orbital rotations
        next_workspace_index = 1
        !
        initial_components_start = next_workspace_index
        next_workspace_index = next_workspace_index + max_initial_symmetry_size*radial_point_count
        !      WRITE (6, *) ' In biotrn: KLPI    = ', KLPI
        !      WRITE (6, *) '            KFREE   = ', KFREE
        !
        final_components_start = next_workspace_index
        next_workspace_index = next_workspace_index + max_final_symmetry_size*radial_point_count
        !      WRITE (6, *) ' In biotrn: KLPF    = ', KLPF
        !      WRITE (6, *) '            KFREE   = ', KFREE
        !
        rotated_components_start = next_workspace_index
        next_workspace_index = next_workspace_index + max_symmetry_size*radial_point_count
        !      WRITE (6, *) ' In biotrn: KLPIF   = ', KLPIF
        !      WRITE (6, *) '            KFREE   = ', KFREE
        !
        !. Total overlap matrix
        !
        full_overlap_start = next_workspace_index
        next_workspace_index = next_workspace_index + total_initial_orbitals*total_final_orbitals
        !      WRITE (6, *) ' In biotrn: KLSTOT  = ', KLSTOT
        !      WRITE (6, *) '            KFREE   = ', KFREE
        !
        overlap_block_start = next_workspace_index
        next_workspace_index = next_workspace_index + max_symmetry_size**2
        !      WRITE (6, *) ' In biotrn: KLSIF   = ', KLSIF
        !      WRITE (6, *) '            KFREE   = ', KFREE
        !
        inverse_overlap_start = next_workspace_index
        next_workspace_index = next_workspace_index + max_symmetry_size**2
        !      WRITE (6, *) ' In biotrn: KLSIFI  = ', KLSIFI
        !      WRITE (6, *) '            KFREE   = ', KFREE
        !
        initial_rotation_start = next_workspace_index
        next_workspace_index = next_workspace_index + max_symmetry_size**2
        !      WRITE (6, *) ' In biotrn: KLCI    = ', KLCI
        !      WRITE (6, *) '            KFREE   = ', KFREE
        !
        final_rotation_start = next_workspace_index
        next_workspace_index = next_workspace_index + max_symmetry_size**2
        !      WRITE (6, *) ' In biotrn: KLCF    = ', KLCF
        !      WRITE (6, *) '            KFREE   = ', KFREE
        !
        auxiliary_workspace_start = next_workspace_index
        next_workspace_index = next_workspace_index + max_symmetry_size**2 + max_symmetry_size*(max_symmetry_size + 1)
        !      WRITE (6, *) ' In biotrn: KLSCR   = ', KLSCR
        !      WRITE (6, *) '            KFREE   = ', KFREE
        !      WRITE (6, *) '         =>  FREE =     ', KFREE
        !
        !. Check length of scratch
        !
        if (workspace_size <= next_workspace_index - 1) then
            write (output_unit,*)" transform_orbital_sets in trouble ! "
            write (output_unit,*)" Increase dimension of scratch before call to BIOTR1"
            write (output_unit,*)" Current and required length (LSCR,KFREE-1)", workspace_size, &
                next_workspace_index - 1
            stop "Increase LWORK before call to BIOTR1"
        end if
        !
        !. Obtain overlap matrix
        !
        call calculate_subshell_overlap (workspace(full_overlap_start), total_initial_orbitals,total_final_orbitals)

        do symmetry_index = 1, symmetry_count
            if (trace_level >= 5) then
                write (output_unit, *) "   L = ", symmetry_index
                write (output_unit, *) "       Orbital rotation..."
            end if
            !
            !. Offsets for given L in shell matrices
            !
            if (symmetry_index == 1) then
                initial_rotation_offset = 1
                final_rotation_offset = 1
            else
                initial_rotation_offset = initial_rotation_offset + initial_orbitals_by_symmetry(symmetry_index-1)**2
                final_rotation_offset = final_rotation_offset + final_orbitals_by_symmetry(symmetry_index-1)**2
            end if
            if (trace_level >= 1) then
                write (output_unit, *)
                write (output_unit, *) &
                    " BIOTRN : Information on transformations of shells with L =",symmetry_index
                write (output_unit, *)
            end if
            !
            ! =========================================================
            ! 1 : Obtain Biorthogonal forms of initial and final shells
            ! =========================================================
            !
            !. The overlap matrix can be written
            !
            !     * * * * * * *
            !     *       *   *
            !     *   S   * X *
            !     *       *   *
            !     * * * * * * *
            !
            ! where S is the quadratic subblock.
            ! The basis functions corresponding to the quadratic
            ! subblock are made biorthogonal with an UL decomposition.
            ! The remaining basis
            ! functions are made biorthogonal by choosing  the
            ! transformation
            !            * * *
            !            *   *
            !            * Y *
            !            *   *
            !            * * *
            !            * 1 *
            !            * * *
            !
            !
            ! With Y = -S-1*X
            !
            initial_orbital_count = initial_orbitals_by_symmetry(symmetry_index)
            final_orbital_count = final_orbitals_by_symmetry(symmetry_index)
            !
            !.1.1 : Obtain shells of given L in proper order for
            !       biorthogonal treatment
            !
            if (symmetry_index == 1) then
                first_initial_orbital = 1
                first_final_orbital = 1
            else
                first_initial_orbital = first_initial_orbital + initial_orbitals_by_symmetry(symmetry_index-1)
                first_final_orbital = first_final_orbital + final_orbitals_by_symmetry(symmetry_index-1)
            end if

            !
            ! 1.2 obtain biorthogonal of the first min(ni,nf) shells
            !
            paired_orbital_count = min(initial_orbital_count,final_orbital_count)
            !
            !. Overlap matrix SIF = Integral (PI(I)*PF(J))
            !
            !ww Per change to support cases like 3s2p1d 3s2p. The belonging endif
            !ww Is just at the end

            if (paired_orbital_count <= 0) cycle

            do initial_orbital = 1, paired_orbital_count
                do final_orbital = 1, paired_orbital_count
                    workspace(overlap_block_start+(final_orbital-1)*paired_orbital_count+initial_orbital- &
                & 1) = workspace(full_overlap_start-1+(final_orbital+first_final_orbital-1-1)*total_initial_orbitals&
                        +initial_orbital+first_initial_orbital-1)
                end do
            end do
            if (trace_level >= 15) then
                write (output_unit, *) " Overlap matrix "
                call write_matrix (workspace(overlap_block_start), paired_orbital_count, paired_orbital_count, &
                & paired_orbital_count, paired_orbital_count)
            end if
            !
            ! Obtain upper triangular CI and CF so CI(T) S CF = 1
            ! or CF CI(T) = S-1, which corresponds to an UL decomposition
            !
            !. Invert S
            !
            call copy_vector (workspace(overlap_block_start), workspace(inverse_overlap_start), paired_orbital_count**2)
            call invert_dense_submatrix (workspace(inverse_overlap_start), workspace(initial_rotation_start), &
                & paired_orbital_count, paired_orbital_count)

            !. UL decompose
            call copy_vector (workspace(inverse_overlap_start), workspace(overlap_block_start), paired_orbital_count**2)
            call combine_upper_lower_factors (workspace(overlap_block_start), workspace(final_rotation_start), &
                & workspace(initial_rotation_start), paired_orbital_count, workspace(auxiliary_workspace_start))
            call transpose_matrix (workspace(initial_rotation_start), paired_orbital_count, paired_orbital_count, &
                & workspace(auxiliary_workspace_start))
            call copy_vector (workspace(auxiliary_workspace_start), workspace(initial_rotation_start), paired_orbital_count**2)
            !
            !. The transformation matrix between the first NIFMX
            !. shells is now known, biorthogonalize remaining orbitals
            !
            if (initial_orbital_count/=final_orbital_count .and. initial_orbital_count/=0 .and. final_orbital_count/=0) then
                if (initial_orbital_count > final_orbital_count) then
                    larger_components_start = initial_components_start
                    smaller_components_start = final_components_start
                    larger_orbital_count = initial_orbital_count
                    smaller_orbital_count = final_orbital_count
                    larger_rotation_start = initial_rotation_start
                    smaller_rotation_start = final_rotation_start
                else
                    larger_components_start = final_components_start
                    smaller_components_start = initial_components_start
                    larger_orbital_count = final_orbital_count
                    smaller_orbital_count = initial_orbital_count
                    larger_rotation_start = final_rotation_start
                    smaller_rotation_start = initial_rotation_start
                end if
                unpaired_orbital_count = larger_orbital_count - smaller_orbital_count
                !
                ! Y = -S-1 * X
                !. overlap X between remaining orbitals and the other set
                !
                if (initial_orbital_count > final_orbital_count) then
                    !
                    ! I columns F rows
                    !
                    do initial_orbital = smaller_orbital_count + 1, larger_orbital_count
                        do final_orbital = 1, final_orbital_count
                            workspace(overlap_block_start+(initial_orbital-smaller_orbital_count-1)*final_orbital_count+ &
                & final_orbital-1) = workspace(full_overlap_start-1+(final_orbital+first_final_orbital-1-&
                                1)*total_initial_orbitals+initial_orbital+first_initial_orbital-1)
                        end do
                    end do
                else if (final_orbital_count > initial_orbital_count) then
                    ! F columns I rows
                    do final_orbital = smaller_orbital_count + 1, larger_orbital_count
                        do initial_orbital = 1, initial_orbital_count
                            workspace(overlap_block_start+(final_orbital-smaller_orbital_count-1)*initial_orbital_count+ &
                & initial_orbital-1) = workspace(full_overlap_start-1+(final_orbital+first_final_orbital-1-&
                                1)*total_initial_orbitals+initial_orbital+first_initial_orbital-1)
                        end do
                    end do
                end if
                !
                if (initial_orbital_count > final_orbital_count) then
                    call transpose_matrix (workspace(inverse_overlap_start), smaller_orbital_count, smaller_orbital_count, &
                & workspace(auxiliary_workspace_start))
                    call copy_vector (workspace(auxiliary_workspace_start), workspace(inverse_overlap_start), &
                & smaller_orbital_count**2)
                end if
                call multiply_submatrices (workspace(auxiliary_workspace_start), workspace(inverse_overlap_start), &
                & workspace(overlap_block_start), smaller_orbital_count, unpaired_orbital_count, smaller_orbital_count, &
                    smaller_orbital_count, smaller_orbital_count, unpaired_orbital_count, 0)
                call scale_vector (workspace(auxiliary_workspace_start), -1.0_real64, smaller_orbital_count* &
                & unpaired_orbital_count)
                call copy_vector (workspace(auxiliary_workspace_start), workspace(overlap_block_start), smaller_orbital_count* &
                & unpaired_orbital_count)
                !
                ! Construct complete CMX
                !
                call fill_vector (workspace(auxiliary_workspace_start), 0.0_real64, larger_orbital_count**2)
                do rotation_column = 1, larger_orbital_count
                    if (rotation_column <= paired_orbital_count) then
                        call copy_vector (workspace(larger_rotation_start+(rotation_column-1)*paired_orbital_count), &
                & workspace(auxiliary_workspace_start+(rotation_column-1)*larger_orbital_count), &
                            smaller_orbital_count)
                    else
                        call copy_vector (workspace(overlap_block_start+(rotation_column-smaller_orbital_count-1)* &
                & smaller_orbital_count), workspace(auxiliary_workspace_start+(rotation_column-1)*larger_orbital_count), &
                            smaller_orbital_count)
                        workspace(auxiliary_workspace_start-1+(rotation_column-1)*larger_orbital_count+ &
                & rotation_column) = 1.0_real64
                    end if
                end do
                !
                call copy_vector (workspace(auxiliary_workspace_start), workspace(larger_rotation_start), larger_orbital_count**2)
            end if
            !ww Pertest
            !         ENDIF
            !
            !. The two upper triangular matrices CI and CF are now known
            !. Transfer to permanent arrays
            !
            call copy_vector (workspace(initial_rotation_start), initial_orbital_rotations(initial_rotation_offset), &
                & initial_orbital_count**2)
            call copy_vector (workspace(final_rotation_start), final_orbital_rotations(final_rotation_offset), &
                & final_orbital_count**2)
            !
            !. Rotate the large component of the shells
            !
            call copy_vector (initial_large_components(1,first_initial_orbital), workspace(initial_components_start), &
                & initial_orbital_count*radial_point_count)
            call copy_vector (final_large_components(1,first_final_orbital), workspace(final_components_start), &
                & final_orbital_count*radial_point_count)
            !
            !         WRITE (*, *) 'Transformation matrices initial'
            !         CALL WRTMAT (SCR(KLCI), NI, NI, NI, NI)
            call multiply_submatrices (workspace(rotated_components_start), workspace(initial_components_start), &
                & workspace(initial_rotation_start), radial_point_count, initial_orbital_count, radial_point_count, &
                & initial_orbital_count, &
                initial_orbital_count, initial_orbital_count, 0)
            call copy_vector (workspace(rotated_components_start), initial_large_components(1,first_initial_orbital), &
                & initial_orbital_count*radial_point_count)
            !         WRITE (*, *) 'Transformation matrices final'
            !         CALL WRTMAT (SCR(KLCF), NF, NF, NF, NF)
            call multiply_submatrices (workspace(rotated_components_start), workspace(final_components_start), &
                & workspace(final_rotation_start), radial_point_count, final_orbital_count, radial_point_count, &
                & final_orbital_count, &
                final_orbital_count, final_orbital_count, 0)
            call copy_vector (workspace(rotated_components_start), final_large_components(1,first_final_orbital), &
                & final_orbital_count*radial_point_count)
            !
            !. Rotate the small component of the shells
            !
            call copy_vector (initial_small_components(1,first_initial_orbital), workspace(initial_components_start), &
                & initial_orbital_count*radial_point_count)
            call copy_vector (final_small_components(1,first_final_orbital), workspace(final_components_start), &
                & final_orbital_count*radial_point_count)
            !
            call multiply_submatrices (workspace(rotated_components_start), workspace(initial_components_start), &
                & workspace(initial_rotation_start), radial_point_count, initial_orbital_count, radial_point_count, &
                & initial_orbital_count, &
                initial_orbital_count, initial_orbital_count, 0)
            call copy_vector (workspace(rotated_components_start), initial_small_components(1,first_initial_orbital), &
                & initial_orbital_count*radial_point_count)
            call multiply_submatrices (workspace(rotated_components_start), workspace(final_components_start), &
                & workspace(final_rotation_start), radial_point_count, final_orbital_count, radial_point_count, &
                & final_orbital_count, &
                final_orbital_count, final_orbital_count, 0)
            call copy_vector (workspace(rotated_components_start), final_small_components(1,first_final_orbital), &
                & final_orbital_count*radial_point_count)
            !
            if (trace_level >= 1) then
                write (output_unit, *) " Test of overlap of biorthonormal functions"
                ! F columns I rows
                do final_orbital = 1, final_orbital_count
                    do initial_orbital = 1, initial_orbital_count
                        workspace(overlap_block_start+(final_orbital-1)*initial_orbital_count+initial_orbital- &
                & 1) = workspace(full_overlap_start-1+(final_orbital+first_final_orbital-1-1)*total_initial_orbitals&
                            +initial_orbital+first_initial_orbital-1)
                    end do
                end do
                call multiply_submatrices (workspace(auxiliary_workspace_start), workspace(initial_rotation_start), &
                & workspace(overlap_block_start), initial_orbital_count, final_orbital_count, initial_orbital_count, &
                & initial_orbital_count, initial_orbital_count&
                    , final_orbital_count, 1)
                call multiply_submatrices (workspace(overlap_block_start), workspace(auxiliary_workspace_start), &
                & workspace(final_rotation_start), initial_orbital_count, final_orbital_count, initial_orbital_count, &
                & final_orbital_count, final_orbital_count&
                    , final_orbital_count, 0)
                write (output_unit, *) &
                    " new overlap matrix ( should be 1 on diag, 0 elsewhere )"
                call write_matrix (workspace(overlap_block_start), initial_orbital_count, final_orbital_count, &
                & initial_orbital_count, final_orbital_count)
            end if

            if (trace_level >= 1) then
                write (output_unit, *)
                write (output_unit, *) " Orbital Rotation matrix for I state"
                call write_matrix (initial_orbital_rotations(initial_rotation_offset), initial_orbital_count, &
                & initial_orbital_count, initial_orbital_count, initial_orbital_count)
                write (output_unit, *) " Orbital Rotation matrix for F state"
                call write_matrix (final_orbital_rotations(final_rotation_offset), final_orbital_count, final_orbital_count, &
                & final_orbital_count, final_orbital_count)
                write (output_unit, *)
            end if
            !
            !. Matrix for counterrotation of CI coefficients, initial state
            !
            initial_ci_rotation_start = overlap_block_start
            call build_subshell_transformation_matrix (workspace(initial_rotation_start), workspace(initial_ci_rotation_start), &
                & workspace(auxiliary_workspace_start), initial_orbital_count)
            do orbital_index = 1, initial_orbital_count
                diagonal_element = workspace(initial_ci_rotation_start-1+(orbital_index-1)*initial_orbital_count+orbital_index)
                inverse_diagonal_element = 1.0_real64/diagonal_element
                call scale_vector (workspace(initial_ci_rotation_start+(orbital_index-1)*initial_orbital_count), &
                & inverse_diagonal_element, orbital_index - 1)
            end do
            call copy_vector (workspace(initial_ci_rotation_start), initial_ci_rotations(initial_rotation_offset), &
                & initial_orbital_count*initial_orbital_count)
            !
            !. Matrix for counterrotation of CI coefficients, Final state
            !
            final_ci_rotation_start = overlap_block_start
            call build_subshell_transformation_matrix (workspace(final_rotation_start), workspace(final_ci_rotation_start), &
                & workspace(auxiliary_workspace_start), final_orbital_count)
            do orbital_index = 1, final_orbital_count
                diagonal_element = workspace(final_ci_rotation_start-1+(orbital_index-1)*final_orbital_count+orbital_index)
                inverse_diagonal_element = 1.0_real64/diagonal_element
                call scale_vector (workspace(final_ci_rotation_start+(orbital_index-1)*final_orbital_count), &
                & inverse_diagonal_element, orbital_index - 1)
            end do
            call copy_vector (workspace(final_ci_rotation_start), final_ci_rotations(final_rotation_offset), final_orbital_count* &
                & final_orbital_count)
            if (trace_level < 1) cycle
            write (output_unit, *)
            write (output_unit, *) " CI-Rotation matrix for I state"
            call write_matrix (initial_ci_rotations(initial_rotation_offset), initial_orbital_count, initial_orbital_count, &
                & initial_orbital_count, initial_orbital_count)
            write (output_unit, *) " CI-Rotation matrix for F state"
            call write_matrix (final_ci_rotations(final_rotation_offset), final_orbital_count, final_orbital_count, &
                & final_orbital_count, final_orbital_count)
            write (output_unit, *)

        end do
        !
        return
    end subroutine transform_orbital_sets
end module transform_orbital_sets_interface
