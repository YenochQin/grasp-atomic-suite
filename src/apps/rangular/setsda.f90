!***********************************************************************
!                                                                      *
      subroutine setsda(write_coefficient_callback,nonzero_element_count,print_details,symmetry_block_index,mpi_rank, &
         & mpi_process_count,output_file_prefix)
      use mcp_generation_options, only: mcp_temporary_units
      use, intrinsic :: iso_fortran_env, only: output_unit
!                                                                      *
!   This routine examines lists                                        *
!                               (IC,IR,npos)                           *
!   to set up the array  IENDC  required by the Davidson eigensolver   *
!   of Stathopoulos and Fischer.                                       *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 10 Dec 1992   *
!   Modified for block interactions by C. F. Fischer        May 1997   *
!   Modified for multi-processors   by Xinghong He       03 Jul 1998   *
!
!  Currently shared by mcpblk, mcpmpi
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use array_memory_management, only: alloc, dalloc
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use orbital_csf_basis, only: csf_count
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use outsda_i, only: outsda
      implicit none (type, external)
      integer :: mcp_append_unit
      procedure(outsda) :: write_coefficient_callback
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nonzero_element_count
      integer, intent(in) :: symmetry_block_index
      integer, intent(in) :: mpi_rank
      integer, intent(in) :: mpi_process_count
      logical  :: print_details
      character, intent(in) :: output_file_prefix*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mb, iend, iclast, i, ic, npos, io_status
      integer, pointer :: iendc(:)
      integer, pointer :: irow(:)
      character :: mcp_channel_label*3
!-----------------------------------------------

      if (mpi_rank == 0) write (output_unit, *) &
         "Analysing sparse matrix array definition file ...", 30

      read (mcp_temporary_units(30)) mcp_channel_label, mb
      if (symmetry_block_index /= mb) then
         write (stderr_unit, *) "setsda: nb = ", symmetry_block_index, ".NE. mb (=", mb, ")"
         stop
      end if
!
!   Allocate storage for IENDC(0:NCF)
!
      call alloc (iendc, 0,  csf_count,"IENDC","SETSDA"  )
      call alloc (irow, nonzero_element_count, "IROW", "SETSDA" )
!
!   Analyse data on file 30; set up IENDC and IROW
! In multiprocessor environment, iendc of each node will have the
! same length (ncf+1); but will have its own part filled. irow is
! local, and its length is determined by the local parameter nnonz.

      iend = 0
      iclast = 0
      do i = 1, nonzero_element_count
         read (mcp_temporary_units(30)) ic, irow(i), npos
         if (ic /= iclast) then
            iendc(iclast) = iend
            iclast = ic
         end if
         iend = npos
      end do
!xhh - changed to suits MPI environment as well
!      IENDC(NCF) = IEND
      iendc(ic) = iend
!
!   Sorting complete; rewrite to mcpXXX.30 file
!
      open(newunit=mcp_append_unit, file=output_file_prefix//".30", status="OLD", form="UNFORMATTED", iostat=io_status&
         , position="APPEND")
      if (io_status /= 0) then
         write (stderr_unit, *) " Error when opening the file mcp.30"
         stop
      end if

      write (mcp_append_unit) "MCP", symmetry_block_index, csf_count
      write (mcp_append_unit) nonzero_element_count
      write (mcp_append_unit) (iendc(i),i=mpi_rank + 1,csf_count,mpi_process_count), (irow(i),i=1,nonzero_element_count)
      close(mcp_append_unit)
!cjb
!
!   Deallocate storage
!
      call dalloc (iendc, "IENDC", "SETSDA")
      call dalloc (irow, "IROW", "SETSDA")
!
      call write_coefficient_callback (print_details, nonzero_element_count, csf_count)
      return
      end subroutine setsda
