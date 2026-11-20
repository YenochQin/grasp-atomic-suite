!***********************************************************************

      subroutine setmcp(core_subshell_count, input_symmetry_block_count, block_labels, filehead)
      use mcp_file_metadata, only: mcp_input_units
!
!   Open, read, check the header of all .mcp files. Info for each
!   block is not accessed here.
!   File 30 stores the structure of H(DC) ;
!   file 31 stores the  T  coefficients;  files 32, 33, ...,
!   store V(0), V(1), ...
!
!   This version works in both serial/parallel(mpi) environment. The
!   difference comes from parameters filehead (in the arg list), myid
!   and nprocs (in common/mpi/).
!
!   The following items are read from mcp files. myid and nprocs
!   read from the files are checked against the ones in COMMON/mpi/
!
!       ncore, nblock, kmaxf             mcp.30 only
!       ncfblk(i), i=1, nblock)          mcp.30 only
!       idblk(i), i=1, nblock)             mcp.30 only
!
!       MCPLAB, nblock, myid, nprocs     all mcp files
!       nelec, ncf, nw               all mcp files
!       DIAG, ICCUT, LFORDR                all mcp files
!
!   Call(s) to: [LIB92]: CONVRT, OPENFL.
!
!   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
!   Modified by Xinghong He               Last revision: 06 Aug 1998   *
!   Modified by G. Gaigalas               Last revision: 07 Sep 2021   *
!      Spin-angular coefficiants have been put on memory               *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use array_memory_management, only: alloc
      use mpi_io_errors, only: fail_input_io => fail_input_io_mpi
      use atomic_calculation_state, only: electron_count
      use first_order_csf_cutoffs, only: iccut
      use mcp_file_metadata, only: max_mcp_file_multipole_rank
      use rmcdhf_mcp_options, only: diag, lfordr
      use mpi_state, only: mpi_process_count, mpi_rank
      use orbital_csf_basis, only: csf_count, orbital_count
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format, only: convrt
      use file_open, only: openfl
      use setmcp_mem_i, only: setmcp_mem
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: core_subshell_count
      integer , intent(in) :: input_symmetry_block_count
      character , intent(in) :: filehead*(*)
      character  :: block_labels(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: sparse_input_unit
      integer :: lenfh, io_status, ios, i, k, text_length, myidd, nprocss
      logical :: found, found1
      character :: number_text*2, mcp_channel_label*3
      character(len=len(filehead)+3) :: file_path
!-----------------------------------------------
!
      lenfh = len_trim(filehead)

      file_path = filehead(1:lenfh)//".30"
      sparse_input_unit = -1
      open(newunit=sparse_input_unit, file=file_path, form="UNFORMATTED", status="OLD", iostat=io_status, &
         position="asis")
      if (io_status /= 0) call fail_input_io(io_status, sparse_input_unit, "MCP setmcp open "//trim(file_path))

!  Check each record before using its contents.

      read (sparse_input_unit, iostat=ios) core_subshell_count, symmetry_block_count, max_mcp_file_multipole_rank
      if (ios /= 0) call fail_input_io(ios, sparse_input_unit, "MCP setmcp global header")
      if (symmetry_block_count>input_symmetry_block_count .or. symmetry_block_count<1) then
         write (stderr_unit, *) "setmcp: nblock = ", symmetry_block_count
         call fail_input_io(0, sparse_input_unit, "MCP setmcp invalid block count")
      end if

!cjb allocate ncfblk(0:*)
!cjb  CALL ALLOC (NCFBLK, NBLOCK + 1, 'NCFBLK', 'SETMCP')
      call alloc (csf_counts_by_block, 0, symmetry_block_count , "NCFBLK", "SETMCP")
!cjb
      csf_counts_by_block(0) = 0

      read (sparse_input_unit, iostat=ios) (csf_counts_by_block(i),i=1,symmetry_block_count)
      if (ios /= 0) call fail_input_io(ios, sparse_input_unit, "MCP setmcp block sizes")
      read (sparse_input_unit, iostat=ios) (block_labels(i),i=1,symmetry_block_count)
      if (ios /= 0) call fail_input_io(ios, sparse_input_unit, "MCP setmcp block labels")

!  Look for other mcp files

      found = .true.
      do k = 31, 32 + max_mcp_file_multipole_rank
         call convrt (k, number_text, text_length)
         file_path = filehead(1:lenfh)//"."//number_text(1:2)
         inquire(file=file_path, exist=found1)
         found = found .and. found1
      end do

      if (.not.found) then
         write (stderr_unit, *) "The mcp files do not exist"
         call fail_input_io(0, sparse_input_unit, "MCP setmcp missing channel file")
      end if

!  Open the files; check file headers

      allocate(mcp_input_units(30:32 + max_mcp_file_multipole_rank))
      mcp_input_units = -1
      mcp_input_units(30) = sparse_input_unit

      do k = 30, 32 + max_mcp_file_multipole_rank

         if (k /= 30) then
            call convrt (k, number_text, text_length)
            file_path = filehead(1:lenfh)//"."//number_text(1:2)
            call openfl (mcp_input_units(k), file_path, "UNFORMATTED", "OLD", io_status)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(k), "MCP setmcp open "//trim(file_path))
         end if

         read (mcp_input_units(k), iostat=ios) mcp_channel_label, symmetry_block_count, myidd, nprocss
         if (ios /= 0) call fail_input_io(ios, mcp_input_units(k), "MCP setmcp block or channel header")
         if (mpi_rank/=myidd .or. mpi_process_count/=nprocss) then
            write (stderr_unit, *) "mcp files were generated under different", &
               " processor configuration."
            call fail_input_io(0, mcp_input_units(k), "MCP setmcp incompatible rank layout")
         end if

         if (mcp_channel_label /= "MCP") then
            write (stderr_unit, *) "Not a sorted GRASP92 MCP File;"
            call fail_input_io(0, mcp_input_units(k), "MCP setmcp invalid file tag")
         end if

         read (mcp_input_units(k), iostat=ios) electron_count, csf_count, orbital_count
         if (ios /= 0) call fail_input_io(ios, mcp_input_units(k), "MCP setmcp basis dimensions")
         read (mcp_input_units(k), iostat=ios) diag, iccut, lfordr
         if (ios /= 0) call fail_input_io(ios, mcp_input_units(k), "MCP setmcp channel flags")
         call setmcp_mem(symmetry_block_count,k)
      end do
      if(mpi_rank == 0) write(*,"(A1)") ""

      return
      end subroutine setmcp
