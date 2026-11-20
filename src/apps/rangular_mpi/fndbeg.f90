!***********************************************************************
!                                                                      *
      subroutine fndbeg(first_bra_csf_index, first_ket_csf_index, matrix_element_position, one_electron_coefficient_count, &
         & two_electron_counts_by_rank)
      use mcp_generation_options, only: mcp_temporary_units
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use orbital_csf_basis, only: csf_count
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: first_bra_csf_index
      integer, intent(inout) :: first_ket_csf_index
      integer, intent(inout) :: matrix_element_position
      integer, intent(out) :: one_electron_coefficient_count
      integer, pointer :: two_electron_counts_by_rank(:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, io_status, icread, irread, indx, nrec, i, label
      integer, save :: sms_restart_unit = -1
      real(real64) :: coeff
      character :: srtlab*8, mcp_channel_label*3
!-----------------------------------------------
!
!   Begin by examining file 30; this is the last to be updated by
!   SUBROUTINE MCP
!
!   Read and check the character part of the file header
!
      rewind (mcp_temporary_units(30))
      read (mcp_temporary_units(30)) mcp_channel_label, srtlab
      if (srtlab == "  SORTED") then
         first_bra_csf_index = csf_count + 1
         first_ket_csf_index = csf_count + 1
         one_electron_coefficient_count = 0
         two_electron_counts_by_rank(:max_coulomb_multipole_rank) = 0
         goto 8
      end if
!
      read (mcp_temporary_units(30))
      read (mcp_temporary_units(30))
!
!   Read as many records as possible
!
    2 continue
      read (mcp_temporary_units(30), iostat=io_status) icread, irread, indx
!
      if (io_status == 0) then
!
!   No errors or end-of-file; keep reading
!
         first_bra_csf_index = icread
         first_ket_csf_index = irread
         matrix_element_position = indx
         goto 2
!
      else
!
         if (first_bra_csf_index==csf_count .and. first_ket_csf_index==csf_count) then
!
!   All coefficients have been generated; sorting may still
!   be necessary; force this option
!
            first_bra_csf_index = csf_count + 1
            first_ket_csf_index = csf_count + 1
!
         else
!
!   Some coefficients remain to be generated; reposition all files
!   for augmentation of lists by SUBROUTINE MCP; update JBSTRT and,
!   if appropriate, JASTRT
!
            do k = 31, 32 + max_coulomb_multipole_rank
               rewind (mcp_temporary_units(k))
               nrec = 3
               do i = 1, nrec
                  read (mcp_temporary_units(k))
               end do
    4          continue
               read (mcp_temporary_units(k), iostat=io_status) indx, label, coeff
               if (io_status==0 .and. indx<=matrix_element_position) then
                  nrec = nrec + 1
                  goto 4
               else
                  rewind (mcp_temporary_units(k))
                  do i = 1, nrec
                     read (mcp_temporary_units(k))
                  end do
                  if (k > 31) then
                     two_electron_counts_by_rank(k-32) = nrec - 3
                  else
                     one_electron_coefficient_count = nrec - 3
                  end if
               end if
            end do
!
!  Now, reposition the sms file. This file should contain the
!  same number of data records as file 33.

            if (sms_restart_unit == -1 .and. two_electron_counts_by_rank(1) > 0) then
               open(newunit=sms_restart_unit, file="fort.20", status="unknown", form="unformatted")
            end if
            if (sms_restart_unit /= -1) rewind (sms_restart_unit)
            do i = 1, two_electron_counts_by_rank(1)
               read (sms_restart_unit)
            end do

            first_ket_csf_index = first_ket_csf_index + 1
            if (first_ket_csf_index > csf_count) then
               first_bra_csf_index = first_bra_csf_index + 1
               first_ket_csf_index = first_bra_csf_index
            end if
!
         end if
!
      end if
!
    8 continue
      return
      end subroutine fndbeg
