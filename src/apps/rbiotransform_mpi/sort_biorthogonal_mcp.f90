module sort_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: sort_biorthogonal_mcp
contains

    !***********************************************************************
    !                                                                      *
    subroutine sort_biorthogonal_mcp(nfile, number, kstart, name, kamax)
        ! Module dependencies
        use biorthogonal_subshell_layout, only: kappa_coefficient_units
        use biorthogonal_subshell_layout, only: biorthogonal_angular_unit
        !                                                                      *
        !     The list of unique integrals (j,i) is formed in the order of     *
        !     increasing symmetry, i.e. with j .le. i.                         *
        !     With each integral (INTGRL) there is a pointer (INTPTR)          *
        !     that points to the last element in the array of coefficients     *
        !     for that integral.                                               *
        !                                                                      *
        !     On Exit                                                          *
        !                                                                      *
        !     INTGRL -- array of integrals (packed form of j,i)                *
        !     NINTG  -- number of integrals                                    *
        !     INTPTR -- array of pointers to last element in list of coeff.    *
        !     CNN    -- array of coefficients                                  *
        !     NCOEFF -- number of coefficients                                 *
        !     JANN   -- row of matrix                                          *
        !     JBNN   -- column of matrix                                       *
        !                                                                      *
        !     Written by Per Jonsson                                           *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use array_memory_management, only: alloc, dalloc
        use orbital_csf_basis, only: csf_count, orbital_count
        use calculation_input_defaults, only: angular_output_dump_flag
        use mcp_coefficient_storage, only: mcp_coefficient_values, mcp_packed_orbital_pairs, mcp_coefficient_end_offsets, &
            mcp_bra_csf_indices, mcp_ket_csf_indices, mcp_coefficient_count, mcp_integral_count
        use mpi_state, only: mpi_process_count, mpi_rank
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nfile
        integer, intent(in) :: number
        integer  :: kstart
        integer, intent(in) :: kamax
        character  :: name*24

        ! Local constants

        ! Local declarations
        integer :: i, heap_root, heap_end, held_bra_csf_index, held_ket_csf_index, held_packed_orbital_pair, heap_child, &
            & current_packed_orbital_pair
        real(real64) :: held_coefficient
        !-----------------------------------------------
        !
        rewind (kappa_coefficient_units(nfile))

        mcp_coefficient_count = number
        !
        !   Sort the list
        !
        !   Allocate storage for all required arrays these arrays are then deallocated
        !   in mcp
        !
        call alloc (mcp_bra_csf_indices,   mcp_coefficient_count, "mcp_bra_csf_indices",   "sort_biorthogonal_mcp")
        call alloc (mcp_ket_csf_indices,   mcp_coefficient_count, "mcp_ket_csf_indices",   "sort_biorthogonal_mcp")
        call alloc (mcp_packed_orbital_pairs, mcp_coefficient_count, "mcp_packed_orbital_pairs", "sort_biorthogonal_mcp")
        call alloc (mcp_coefficient_values,    mcp_coefficient_count, "mcp_coefficient_values",    "sort_biorthogonal_mcp")
        call alloc (mcp_coefficient_end_offsets, mcp_coefficient_count, "mcp_coefficient_end_offsets", "sort_biorthogonal_mcp")
        !
        !   Read arrays into memory from NFILE
        !
        do i = 1, mcp_coefficient_count
            read (kappa_coefficient_units(nfile)) mcp_bra_csf_indices(i), mcp_ket_csf_indices(i), &
                mcp_packed_orbital_pairs(i), mcp_coefficient_values(i)
        end do
        !
        !   Sort INTGRL into ascending order using the heapsort algorithm;
        !   (Numerical recepies page 231.) move the associated members of
        !   of JANN and JBNN in the same
        !   manner;
        !
        if (mcp_coefficient_count > 1) then
            !
            heap_root = mcp_coefficient_count/2 + 1
            heap_end = mcp_coefficient_count
            do
                if (heap_root > 1) then
                    heap_root = heap_root - 1
                    held_bra_csf_index = mcp_bra_csf_indices(heap_root)
                    held_ket_csf_index = mcp_ket_csf_indices(heap_root)
                    held_packed_orbital_pair = mcp_packed_orbital_pairs(heap_root)
                    held_coefficient = mcp_coefficient_values(heap_root)
                else
                    held_bra_csf_index = mcp_bra_csf_indices(heap_end)
                    held_ket_csf_index = mcp_ket_csf_indices(heap_end)
                    held_packed_orbital_pair = mcp_packed_orbital_pairs(heap_end)
                    held_coefficient = mcp_coefficient_values(heap_end)
                    mcp_bra_csf_indices(heap_end) = mcp_bra_csf_indices(1)
                    mcp_ket_csf_indices(heap_end) = mcp_ket_csf_indices(1)
                    mcp_packed_orbital_pairs(heap_end) = mcp_packed_orbital_pairs(1)
                    mcp_coefficient_values(heap_end) = mcp_coefficient_values(1)
                    heap_end = heap_end - 1
                    if (heap_end == 1) then
                        mcp_bra_csf_indices(1) = held_bra_csf_index
                        mcp_ket_csf_indices(1) = held_ket_csf_index
                        mcp_packed_orbital_pairs(1) = held_packed_orbital_pair
                        mcp_coefficient_values(1) = held_coefficient
                        exit
                    end if
                end if
                i = heap_root
                heap_child = heap_root + heap_root
                do while (heap_child <= heap_end)
                    if (heap_child < heap_end) then
                        if (mcp_packed_orbital_pairs(heap_child) < mcp_packed_orbital_pairs(heap_child+1)) heap_child = &
                            & heap_child + 1
                    end if
                    if (held_packed_orbital_pair < mcp_packed_orbital_pairs(heap_child)) then
                        mcp_bra_csf_indices(i) = mcp_bra_csf_indices(heap_child)
                        mcp_ket_csf_indices(i) = mcp_ket_csf_indices(heap_child)
                        mcp_packed_orbital_pairs(i) = mcp_packed_orbital_pairs(heap_child)
                        mcp_coefficient_values(i) = mcp_coefficient_values(heap_child)
                        i = heap_child
                        heap_child = heap_child + heap_child
                    else
                        heap_child = heap_end + 1
                    end if
                end do
                mcp_bra_csf_indices(i) = held_bra_csf_index
                mcp_ket_csf_indices(i) = held_ket_csf_index
                mcp_packed_orbital_pairs(i) = held_packed_orbital_pair
                mcp_coefficient_values(i) = held_coefficient
            end do
        end if
        !
        !   Sorting complete
        !
        mcp_integral_count = 1
        current_packed_orbital_pair = mcp_packed_orbital_pairs(1)
        !
        do i = 1, mcp_coefficient_count
            if (mcp_packed_orbital_pairs(i) == current_packed_orbital_pair) cycle
            mcp_coefficient_end_offsets(mcp_integral_count) = i - 1
            mcp_integral_count = mcp_integral_count + 1
            current_packed_orbital_pair = mcp_packed_orbital_pairs(i)
        end do

        mcp_coefficient_end_offsets(mcp_integral_count) = mcp_coefficient_count
        !
        do i = 1, mcp_integral_count
            mcp_packed_orbital_pairs(i) = mcp_packed_orbital_pairs(mcp_coefficient_end_offsets(i))
        end do
        !
        !  If output option is set dump the data on file
        !
        if (angular_output_dump_flag == 1) then
            !
            !  If first set of data open the file and print
            !  some data to later be able to identify the file
            !
            !        IF (NFILE.EQ.1.AND.IBLK.EQ.1) THEN
            !         print *, ' open sorted ang. file .TB(biorthogonal_angular_unit)', biorthogonal_angular_unit
            !          J = INDEX(NAME,' ')
            !          OPEN (UNIT = biorthogonal_angular_unit,FILE=NAME(1:J-1)//'.TB',
            !     :          STATUS='UNKNOWN',FORM='UNFORMATTED')
            !          REWIND (biorthogonal_angular_unit)
            !        ENDIF
            if (nfile == 1) write (biorthogonal_angular_unit) csf_count,orbital_count,kamax,mpi_process_count,mpi_rank
            !
            !  Print out angular data for this kappa
            !
            write (biorthogonal_angular_unit) mcp_integral_count, mcp_coefficient_count
            do i = 1, mcp_integral_count
                write (biorthogonal_angular_unit) mcp_packed_orbital_pairs(i), mcp_coefficient_end_offsets(i)
            end do
            do i = 1, mcp_coefficient_count
                write (biorthogonal_angular_unit) mcp_coefficient_values(i), mcp_bra_csf_indices(i), mcp_ket_csf_indices(i)
            end do
        end if
        call dalloc (mcp_bra_csf_indices,   "mcp_bra_csf_indices",   "sort_biorthogonal_mcp")
        call dalloc (mcp_ket_csf_indices,   "mcp_ket_csf_indices",   "sort_biorthogonal_mcp")
        call dalloc (mcp_packed_orbital_pairs, "mcp_packed_orbital_pairs", "sort_biorthogonal_mcp")
        call dalloc (mcp_coefficient_values,    "mcp_coefficient_values",    "sort_biorthogonal_mcp")
        call dalloc (mcp_coefficient_end_offsets, "mcp_coefficient_end_offsets", "sort_biorthogonal_mcp")
        !
        !  Has all data been processed? If so close
        !  the file
        !
        !     IF (NFILE.EQ.KAMAX) CLOSE (biorthogonal_angular_unit)

        !  Debug output
        !
        !      IF (KSTART.EQ.1) THEN
        !        NNNN = 120
        !      ELSE
        !        NNNN = 140
        !      ENDIF
        !      WRITE(NNNN+NFILE,*) ' Data for '
        !      WRITE(NNNN+NFILE,*) NINTG, ' Integrals'
        !      DO I = 1,NINTG
        !        IA = INTGRL(I)/KEY
        !        IB = MOD(INTGRL(I),KEY)
        !        WRITE(NNNN+NFILE,302) NP(IA),NH(IA),NP(IB),NH(IB),INTPTR(I)
        !      ENDDO
        !      WRITE(NNNN+NFILE,*) NCOEFF,' Coefficients'
        !      WRITE(NNNN+NFILE,'(F12.8,2I6)')
        !     :     (CNN(I),JANN(I),JBNN(I),I=1,NCOEFF)
        !
        return
        !
        return
        !
    end subroutine sort_biorthogonal_mcp
end module sort_biorthogonal_mcp_interface
