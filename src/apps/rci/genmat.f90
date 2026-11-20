!***********************************************************************
!                                                                      *
      subroutine genmat(atwinv, symmetry_block_index, mpi_rank, mpi_process_count, elsto, irestart, slf_en)
      use calculation_file_units, only: summary_output_unit
!
!        Generate Hamiltonian matrix for all blocks
!        This routine calls setham to do the computation. It makes
!        sure that the hamiltonian matrix is complete.
!  Only node-0 has the correct, non-zero elsto. And in any case, elsto
!  is not added to the matrix elements in files *.res
!  This routine has been re-written to work in both single- and
!  multi- processors.
!
! Xinghong He 1998-06-23
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use csf_block_cutoffs, only: iccut, iccutblk
      use orbital_csf_basis, only: csf_count, orbital_count
      use where_c, only: restart_matrix_unit
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use file_positioning, only: posfile
      use setham_i, only: setham
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: symmetry_block_index
      integer, intent(in) :: mpi_rank
      integer, intent(in) :: mpi_process_count
      integer, intent(out) :: irestart
      real(real64)  :: atwinv
      real(real64)  :: elsto
      real(real64) :: slf_en(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iread, io_status, file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count, i, &
         ios2, nelc, ir, irowdum, ipos, row_count, j, icstrt
      real(real64) :: stoel, dum, eav0
!-----------------------------------------------
!
      stored_matrix_element_count = 0                      ! Counting continues in setham
      hamiltonian_energy_shift    = 0.0_real64
      elsto  = 0.0_real64

! See how much had been done (Hamiltonian matrix)
! irestart is set;
! iread accumulated;
! nelmnt, eav, elsto obtained (to be further modified in setham)

      iread = 0                       ! # of rows read, initialization necessary

      read (restart_matrix_unit, iostat=io_status) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
      irestart = io_status

      if (io_status == 0) then

         if (csf_count/=file_csf_count .or. iccut(1)/=iccutdum .or. mpi_rank/=file_mpi_rank .or. &
            mpi_process_count/= file_mpi_process_count) then
            write (stderr_unit, *) csf_count, file_csf_count, iccut(1), iccutdum, mpi_rank, file_mpi_rank, &
               mpi_process_count, file_mpi_process_count, "check"
            stop "genmat:1"
         end if

         do i = mpi_rank + 1, csf_count, mpi_process_count
            read (restart_matrix_unit, iostat=ios2) nelc, stoel, (dum,ir=2,nelc), eav0, (&
               irowdum,ir=1,nelc)
                            ! Lower triangle row-mode, diagonal last
            if (ios2 == 0) then
               iread = iread + 1
               stored_matrix_element_count = stored_matrix_element_count + nelc
               hamiltonian_energy_shift = hamiltonian_energy_shift + eav0
               elsto = stoel
            else
               exit
            end if
         end do
         ipos = 7 + orbital_count + orbital_count + iread + 1
      else
         ipos = 7 + orbital_count + orbital_count
      end if

! Find the maximum number of rows

      row_count = (csf_count - mpi_rank - 1 + mpi_process_count)/mpi_process_count
      if (csf_count < mpi_process_count) row_count = csf_count/(mpi_rank + 1)

! Report the number of rows read.
! A more suitable report on all nodes can be done here, but this will
! set a synchronization point.

      write (stderr_unit, *) iread, " (total ", row_count, ") rows read from .res"
      if (mpi_rank == 0) write (summary_output_unit, *) iread, " (total ", row_count, &
         ") rows read from .res"

! Position the file for the next record from setham

      do i = 1, symmetry_block_index - 1
         j = (csf_counts_by_block(i) - mpi_rank - 1 + mpi_process_count)/mpi_process_count
         if (csf_counts_by_block(i) < mpi_process_count) j = csf_counts_by_block(i)/(mpi_rank + 1)
         ipos = ipos + j + 1
      end do
      call posfile (0, restart_matrix_unit, ipos)

      if (io_status /= 0) write (restart_matrix_unit) csf_count, iccut(1), mpi_rank, mpi_process_count

      if (iread < row_count) then
         icstrt = iread*mpi_process_count + mpi_rank + 1
!     ...Generate the rest of the Hamiltonian matrix
         call setham (mpi_rank, mpi_process_count, symmetry_block_index, elsto, icstrt, stored_matrix_element_count, atwinv, &
            slf_en)
      else
         nelmnttmp = stored_matrix_element_count
         ncftmp = csf_count
      end if

      return
      end subroutine genmat
