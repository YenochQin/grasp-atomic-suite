!***********************************************************************
!                                                                      *
      subroutine strsum

!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use atomic_calculation_state, only: electron_count
      use csf_block_cutoffs, only: iccut, iccutblk
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use orbital_csf_basis, only: csf_count, orbital_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calendar_clock, only: calen
      use integer_format, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: lenth, i
      character :: ctime*10, cdate*8
      character :: record*256, cdata*26
!-----------------------------------------------
!
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
      call calen (ctime, cdate)
      write (24, *) "GENMCP run at ", ctime, " on ", cdate, "."
!
!   Write out the basic dimensions of the electron cloud
!
      write (24, *)
      call convrt (electron_count, record, lenth)
      write (24, *) "There are "//record(1:lenth)//" electrons in the cloud"
      call convrt (csf_count, record, lenth)
      write (24, *) " in "//record(1:lenth)//" relativistic CSFs"
      call convrt (orbital_count, record, lenth)
      write (24, *) " based on "//record(1:lenth)//" relativistic subshells."
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
      if (diagonal_only) then
         write (24, *)
         write (24, *) "Only diagonal matrix elements are computed."
      else
         if (use_first_order_restriction) then
            do i = 1, symmetry_block_count
               write (24, *)
               call convrt (iccut(i), record, lenth)
               write (24, *) "CSFs 1--"//record(1:lenth)//" constitute"//&
                  " the zero-order space."
            end do
         end if
      end if
!
      return
      end subroutine strsum
