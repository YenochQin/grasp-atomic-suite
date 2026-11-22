module read_ci_level_selection_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_ci_level_selection
contains

    !    Print  the  weights of the largest five CSF contributors to each
    !    ASF.
    !
    !                                           Last updated: 02 Nov 1992
    !
    ! Translated by Charlotte Froese Fischer
    !                        Gediminas Gaigalas  10/05/17
    subroutine read_ci_level_selection(iatjpo, iaspar)
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use eigenstate_expansions, only: csf_mixing_coefficients
        use angular_symmetry_labels, only: format_j_label, parity_labels
        use orbital_csf_basis, only: csf_count
        use selected_level_indices, only: active_level_indices, active_level_count
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: iatjpo, iaspar

        ! Local declarations
        integer :: iconf(5)
        real(real64) :: wt(5)

        real(real64) :: awt
        integer :: i, j, jj, ibegin, iend, nelt, iv, icf, ip
        !    Allocate storage for local arrays
        !
        ! These arrays are not needed
        write (summary_output_unit,'(/"Weights of major contributors to ASF:"//"Level J Parity      CSF contributions"/)')
        ! Initialized NELT for NCF >5
        nelt = 5
        if (csf_count < 5) nelt = csf_count
        do iv = 1,active_level_count
            icf = active_level_indices(iv)
            ! .., initialize for a new eigenvector
            wt = 0.0_real64;  iconf = 0
            ibegin = (iv-1)*csf_count
            iend   = iv*csf_count
            do i = ibegin+1, iend
                awt = abs(csf_mixing_coefficients(i))
                if (awt > wt(5)) then
                    ! ... the weight needs to be inserted
                    do j = 1,5
                        if ( awt > wt(j) ) then
                            if (j < 5) then
                                ! ... shift right
                                do jj = 5, j+1,-1
                                    wt(jj) = wt(jj-1)
                                    iconf(jj) = iconf(jj-1)
                                end do
                            end if
                            ! ...  insert the element
                            wt(j) = awt
                            iconf(j) = i
                            exit
                        end if
                    end do
                end if
            end do
            ! Print first five elements of list.
            ip = (iaspar + 3 ) / 2
            ! NELT = MIN (I,5)
            write (summary_output_unit,'(i3,2x,2a,5(3x,f8.5))') icf,format_j_label(iatjpo),parity_labels(ip),      &
                (csf_mixing_coefficients(iconf(j)),j = 1,nelt)
            write (summary_output_unit,'(13x,5(i11))') (iconf(j)-ibegin,j = 1,nelt)
        end do
        ! Deallocate storage for local arrays
        return

    end subroutine read_ci_level_selection
end module read_ci_level_selection_interface
