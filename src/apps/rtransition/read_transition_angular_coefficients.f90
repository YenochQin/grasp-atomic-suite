module read_transition_angular_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_transition_angular_coefficients
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_transition_angular_coefficients(operator_parity, transition_type_index, name, transition_angular_unit)
        ! Module dependencies
        !                                                                      *
        !   This routine loads  coefficients with parity and  rank specified   *
        !   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
        !   (+/- 1) and is determined from the sign of  KP(JKP).               *
        !                                                                      *
        !                                         Last revision: 28 Dec 1992   *
        !   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters, only: max_orbitals
        use orbital_csf_basis, only: max_orbitals, orbital_count
        use transition_matrix_workspace, only: transition_matrix_capacity, hb1, hb2, hc1, hc2, hm1, hm2, jja, jjb, &
            transition_coefficient_capacity, transition_angular_coefficients, transition_bra_csf_indices, &
            transition_ket_csf_indices, stored_transition_integral_count, transition_integral_count, &
            transition_integral_coefficient_offsets, packed_transition_integral_labels, requested_multipole_count

        ! Procedure interfaces
        use resize_transition_matrix_arrays_interface, only: resize_transition_matrix_arrays
        use resize_transition_integral_arrays_interface, only: resize_transition_integral_arrays
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: transition_angular_unit
        integer  :: operator_parity
        integer, intent(in) :: transition_type_index
        character  :: name(2)*24

        ! Local declarations
        integer :: m, k, initial_block_index, final_block_index, i, labl, ncsf, j
        !-----------------------------------------------
        !
        !   Read the data back as required by OSCL conventions
        !

        m = 0
        k = 0
        !
        read (transition_angular_unit) initial_block_index, final_block_index, orbital_count, requested_multipole_count
        read (transition_angular_unit) stored_transition_integral_count
        !
        do i = 1, stored_transition_integral_count
            !
            read (transition_angular_unit) labl, ncsf
            !
            m = m + 1
            !bieron
            if (m >= transition_matrix_capacity) call resize_transition_matrix_arrays (jja, jjb, hb1, hb2, hc1, &
                hc2, hm1, hm2, packed_transition_integral_labels, transition_integral_coefficient_offsets, &
                    transition_matrix_capacity, 2)
            packed_transition_integral_labels(m) = labl
            !
            !   Read configuration pairs and coefficients for this integral
            !
            do while (ncsf + k > transition_coefficient_capacity)
                call resize_transition_integral_arrays (transition_bra_csf_indices, transition_ket_csf_indices, &
                    transition_angular_coefficients, transition_coefficient_capacity, 2)
            end do
            transition_integral_coefficient_offsets(m) = k
            read (transition_angular_unit) (transition_bra_csf_indices(j + k),transition_ket_csf_indices(j + k), &
                transition_angular_coefficients(j + k),j=1,ncsf)
            !          write(*,*) (ISLDR(J+K),XSLDR(J+K),J = 1,NCSF)
            k = k + ncsf
            !
        end do
        !
        !   Close (and hence release) the scratch file
        !
        !     IF (AVAIL) THEN
        !
        !        CLOSE (unit=coefficient_scratch_unit,status="DELETE")
        !      ENDIF
        !
        transition_integral_coefficient_offsets(m+1) = k
        transition_integral_count = m
        !
        return
        !
        return
        !
    end subroutine read_transition_angular_coefficients
end module read_transition_angular_coefficients_interface
