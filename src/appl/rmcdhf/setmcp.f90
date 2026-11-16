!***********************************************************************

      subroutine setmcp(ncore, nblkin, idblk, filehead)
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
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use input_io, only: fail_input_io
      use def_c
      use foparm_c, only: iccut
      use mcpa_c, only: double, max_mcp_file_multipole_rank
      use mcpb_c, only: diag, double, lfordr
      use mpi_s
      use orb_c
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i, only: convrt
      use openfl_i, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncore
      integer , intent(in) :: nblkin
      character , intent(in) :: filehead*(*)
      character  :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: lenfh, ierror, ios, i, k, lck, myidd, nprocss
      logical :: found, found1
      character :: ck*2, mcplab*3
      character(len=120) :: filnam
!-----------------------------------------------
!
      lenfh = len_trim(filehead)

      filnam = filehead(1:lenfh)//".30"
      open(30, file=filnam, form="UNFORMATTED", status="OLD", iostat=ierror, &
         position="asis")
      if (ierror /= 0) call fail_input_io(ierror, 30, "MCP setmcp open "//trim(filnam))

!  Check each record before using its contents.

      read (30, iostat=ios) ncore, symmetry_block_count, max_mcp_file_multipole_rank
      if (ios /= 0) call fail_input_io(ios, 30, "MCP setmcp global header")
      if (symmetry_block_count>nblkin .or. symmetry_block_count<1) then
         write (stderr_unit, *) "setmcp: nblock = ", symmetry_block_count
         call fail_input_io(0, 30, "MCP setmcp invalid block count")
      end if

!cjb allocate ncfblk(0:*)
!cjb  CALL ALLOC (NCFBLK, NBLOCK + 1, 'NCFBLK', 'SETMCP')
      call alloc (csf_counts_by_block, 0, symmetry_block_count , "NCFBLK", "SETMCP")
!cjb
      csf_counts_by_block(0) = 0

      read (30, iostat=ios) (csf_counts_by_block(i),i=1,symmetry_block_count)
      if (ios /= 0) call fail_input_io(ios, 30, "MCP setmcp block sizes")
      read (30, iostat=ios) (idblk(i),i=1,symmetry_block_count)
      if (ios /= 0) call fail_input_io(ios, 30, "MCP setmcp block labels")

!  Look for other mcp files

      found = .true.
      do k = 31, 32 + max_mcp_file_multipole_rank
         call convrt (k, ck, lck)
         filnam = filehead(1:lenfh)//"."//ck(1:2)
         inquire(file=filnam, exist=found1)
         found = found .and. found1
      end do

      if (.not.found) then
         write (stderr_unit, *) "The mcp files do not exist"
         call fail_input_io(0, 30, "MCP setmcp missing channel file")
      end if

!  Open the files; check file headers

      do k = 30, 32 + max_mcp_file_multipole_rank

         if (k /= 30) then
            call convrt (k, ck, lck)
            filnam = filehead(1:lenfh)//"."//ck(1:2)
            call openfl (k, filnam, "UNFORMATTED", "OLD", ierror)
            if (ierror /= 0) call fail_input_io(ierror, k, "MCP setmcp open "//trim(filnam))
         end if

         read (k, iostat=ios) mcplab, symmetry_block_count, myidd, nprocss
         if (ios /= 0) call fail_input_io(ios, k, "MCP setmcp block or channel header")
         if (mpi_rank/=myidd .or. mpi_process_count/=nprocss) then
            write (stderr_unit, *) "mcp files were generated under different", &
               " processor configuration."
            call fail_input_io(0, k, "MCP setmcp incompatible rank layout")
         end if

         if (mcplab /= "MCP") then
            write (stderr_unit, *) "Not a sorted GRASP92 MCP File;"
            call fail_input_io(0, k, "MCP setmcp invalid file tag")
         end if

         read (k, iostat=ios) electron_count, csf_count, orbital_count
         if (ios /= 0) call fail_input_io(ios, k, "MCP setmcp basis dimensions")
         read (k, iostat=ios) diag, iccut, lfordr
         if (ios /= 0) call fail_input_io(ios, k, "MCP setmcp channel flags")

      end do

      return
      end subroutine setmcp
