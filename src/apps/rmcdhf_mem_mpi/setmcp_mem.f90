!*******************************************************************
!                                                                  *
      subroutine setmcp_mem(nblock,nfile)
      use mcp_file_metadata, only: mcp_input_units
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use array_memory_management, only: alloc
      use mpi_io_errors, only: fail_input_io => fail_input_io_mpi
      use rmcdhf_mem_c, only: coeff_31, coeff_32, coeff_33, coeff_34, coeff_35, coeff_36, coeff_37, coeff_38, &
         & coeff_39, coeff_40, coeff_41, coeff_42, coeff_43, coeff_44, coeff_45, coeff_46, coeff_47, coeff_48, &
         & coeff_49, coeff_50, iclmn_31, iclmn_32, iclmn_33, iclmn_34, iclmn_35, iclmn_36, iclmn_37, iclmn_38, &
         & iclmn_39, iclmn_40, iclmn_41, iclmn_42, iclmn_43, iclmn_44, iclmn_45, iclmn_46, iclmn_47, iclmn_48, &
         & iclmn_49, iclmn_50, iendc_30, imax1_30, imax2_30, imin1_30, imin2_30, indx_31, indx_32, indx_33, indx_34, &
         & indx_35, indx_36, indx_37, indx_38, indx_39, indx_40, indx_41, indx_42, indx_43, indx_44, indx_45, &
         & indx_46, indx_47, indx_48, indx_49, indx_50, irow_30, lab_31, lab_32, lab_33, lab_34, lab_35, lab_36, &
         & lab_37, lab_38, lab_39, lab_40, lab_41, lab_42, lab_43, lab_44, lab_45, lab_46, lab_47, lab_48, lab_49, &
         & lab_50, ncontr_31, ncontr_32, ncontr_33, ncontr_34, ncontr_35, ncontr_36, ncontr_37, ncontr_38, &
         & ncontr_39, ncontr_40, ncontr_41, ncontr_42, ncontr_43, ncontr_44, ncontr_45, ncontr_46, ncontr_47, &
         & ncontr_48, ncontr_49, ncontr_50
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use mpi_state, only: mpi_process_count, mpi_rank
      use sparse_hamiltonian_storage,          only: stored_matrix_element_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use systemmem_i, only: systemmem
      use systemfreemem_i, only: systemfreemem
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nblock, nfile
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character :: mcp_channel_label*3
      real(real64) :: coeff_dumm
      real(real64) :: sysmemt, sysmemf
      integer :: i, symmetry_block_index, jblockt, ncft, ncoeff, io_status, block_csf_count
      integer(int64) :: total_contribution_count, total_secondary_contribution_count
      integer :: iclmn_dumm, lab_dumm, ncontr_dumm, indx_dumm
      integer :: nelmntgg, iendcdum, irow_dumm
!-----------------------------------------------

      total_contribution_count = 0
      total_secondary_contribution_count = 0
      do symmetry_block_index = 1, nblock
        if(nfile == 30) then
          !*** Read in IROW from file mcp.30 ***
          read (mcp_input_units(30), iostat=io_status) mcp_channel_label, jblockt, block_csf_count
          if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(30), "MCP setmcp_mem block or channel header")
          if (jblockt /= symmetry_block_index) then
            write (stderr_unit, *) "setmcp_mem", "1: jblockt .NE. jblock"
            call fail_input_io(0, mcp_input_units(nfile), "MCP setmcp_mem invalid block or channel")
          end if
          read (mcp_input_units(30), iostat=io_status) nelmntgg
          if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(30), "MCP setmcp_mem sparse element count")
          stored_matrix_element_count = int(nelmntgg,kind=int64)
          if(symmetry_block_index == 1) then
            call alloc (imin1_30,nblock,"IMIN1_30","setmcp_mem")
            call alloc (imin2_30,nblock,"IMIN2_30","setmcp_mem")
            call alloc (imax1_30,nblock,"IMAX1_30","setmcp_mem")
            call alloc (imax2_30,nblock,"IMAX2_30","setmcp_mem")
            imin1_30(symmetry_block_index) = 0
            imax1_30(symmetry_block_index) = block_csf_count
            imin2_30(symmetry_block_index) = 0
            imax2_30(symmetry_block_index) = stored_matrix_element_count
          else
            imin1_30(symmetry_block_index) = imax1_30(symmetry_block_index-1) + 0
            imax1_30(symmetry_block_index) = imax1_30(symmetry_block_index-1) + block_csf_count
            imin2_30(symmetry_block_index) = imax2_30(symmetry_block_index-1) + 0
            imax2_30(symmetry_block_index) = imax2_30(symmetry_block_index-1) + stored_matrix_element_count
          end if
          read (mcp_input_units(30), iostat=io_status) (iendcdum,i=mpi_rank + 1,block_csf_count,mpi_process_count), &
             & &
                                                 (irow_dumm,i=1,stored_matrix_element_count)
          if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(30), "MCP setmcp_mem sparse indices")
        else
          read (mcp_input_units(nfile), iostat=io_status) mcp_channel_label, jblockt, ncft, ncoeff
          if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem block or channel header")
          if (jblockt /= symmetry_block_index) then
            write (stderr_unit, *) "setmcp_mem", ":2 jblockt .NE. jblock"
            call fail_input_io(0, mcp_input_units(nfile), "MCP setmcp_mem invalid block or channel")
          end if

          total_secondary_contribution_count = total_secondary_contribution_count + 1
          read (mcp_input_units(nfile), iostat=io_status) lab_dumm, ncontr_dumm
          if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem count-pass label/count")
          do while (lab_dumm /= 0 .or. ncontr_dumm /= 0)
            !*** Read column index, sparse matrix index, and
            !*** coefficient for all contributions from this integral
            read (mcp_input_units(nfile), iostat=io_status) (iclmn_dumm,indx_dumm,coeff_dumm,             &
                                                       i=1,ncontr_dumm)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem count-pass coefficients")
            total_contribution_count = total_contribution_count + ncontr_dumm
            total_secondary_contribution_count = total_secondary_contribution_count + 1
            read (mcp_input_units(nfile), iostat=io_status) lab_dumm, ncontr_dumm
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem count-pass label/count")
          end do
        end if
      end do

      if(nfile == 30) then
         ! Total physical memory
         sysmemt = systemmem()  ! GB
         sysmemf = systemfreemem()  ! GB
         if(mpi_rank == 0) write(*,"(A53)")                                &
         "-----------------------------------------------------"
         if(mpi_rank == 0) write(*,"(A53)")                                &
         "Spin-angular coefficient are putting into the memory:"
         if(mpi_rank == 0) write(*,"(A53)")                                &
         "-----------------------------------------------------"
         if(mpi_rank == 0) write(*,"(A1)") ""
         if(mpi_rank == 0) write(*,"(A30,F10.2,A3)")                       &
                          "Total memory on computer:",sysmemt," Gb"
         if(mpi_rank == 0) write(*,"(A30,F10.2,A3)")                       &
                          "Free  memory on computer:",sysmemf," Gb"
         if(mpi_rank == 0) write(*,"(A1)") ""
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.30:"
         call alloc (iendc_30,imax1_30(nblock),"IENDC_30","setmcp_mem")
         call alloc (irow_30, imax2_30(nblock),"IROW_30", "setmcp_mem")
      else if(nfile == 31) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.31:"
         call alloc (lab_31,   total_secondary_contribution_count,"LAB_31",   "setmcp_mem")
         call alloc (ncontr_31,total_secondary_contribution_count,"NCONTR_31","setmcp_mem")
         call alloc (iclmn_31, total_contribution_count, "ICLMN_31", "setmcp_mem")
         call alloc (indx_31,  total_contribution_count, "INDX_31",  "setmcp_mem")
         call alloc (coeff_31, total_contribution_count, "COEFF_31", "setmcp_mem")
      else if(nfile == 32) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.32:"
         call alloc (lab_32,   total_secondary_contribution_count,"LAB_32",   "setmcp_mem")
         call alloc (ncontr_32,total_secondary_contribution_count,"NCONTR_32","setmcp_mem")
         call alloc (iclmn_32, total_contribution_count, "ICLMN_32", "setmcp_mem")
         call alloc (indx_32,  total_contribution_count, "INDX_32",  "setmcp_mem")
         call alloc (coeff_32, total_contribution_count, "COEFF_32", "setmcp_mem")
      else if(nfile == 33) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.33:"
         call alloc (lab_33,   total_secondary_contribution_count,"LAB_33",   "setmcp_mem")
         call alloc (ncontr_33,total_secondary_contribution_count,"NCONTR_33","setmcp_mem")
         call alloc (iclmn_33, total_contribution_count, "ICLMN_33", "setmcp_mem")
         call alloc (indx_33,  total_contribution_count, "INDX_33",  "setmcp_mem")
         call alloc (coeff_33, total_contribution_count, "COEFF_33", "setmcp_mem")
      else if(nfile == 34) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.34:"
         call alloc (lab_34,   total_secondary_contribution_count,"LAB_34",   "setmcp_mem")
         call alloc (ncontr_34,total_secondary_contribution_count,"NCONTR_34","setmcp_mem")
         call alloc (iclmn_34, total_contribution_count, "ICLMN_34", "setmcp_mem")
         call alloc (indx_34,  total_contribution_count, "INDX_34",  "setmcp_mem")
         call alloc (coeff_34, total_contribution_count, "COEFF_34", "setmcp_mem")
      else if(nfile == 35) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.35:"
         call alloc (lab_35,   total_secondary_contribution_count,"LAB_35",   "setmcp_mem")
         call alloc (ncontr_35,total_secondary_contribution_count,"NCONTR_35","setmcp_mem")
         call alloc (iclmn_35, total_contribution_count, "ICLMN_35", "setmcp_mem")
         call alloc (indx_35,  total_contribution_count, "INDX_35",  "setmcp_mem")
         call alloc (coeff_35, total_contribution_count, "COEFF_35", "setmcp_mem")
      else if(nfile == 36) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.36:"
         call alloc (lab_36,   total_secondary_contribution_count,"LAB_36",   "setmcp_mem")
         call alloc (ncontr_36,total_secondary_contribution_count,"NCONTR_36","setmcp_mem")
         call alloc (iclmn_36, total_contribution_count, "ICLMN_36", "setmcp_mem")
         call alloc (indx_36,  total_contribution_count, "INDX_36",  "setmcp_mem")
         call alloc (coeff_36, total_contribution_count, "COEFF_36", "setmcp_mem")
      else if(nfile == 37) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.37:"
         call alloc (lab_37,   total_secondary_contribution_count,"LAB_37",   "setmcp_mem")
         call alloc (ncontr_37,total_secondary_contribution_count,"NCONTR_37","setmcp_mem")
         call alloc (iclmn_37, total_contribution_count, "ICLMN_37", "setmcp_mem")
         call alloc (indx_37,  total_contribution_count, "INDX_37",  "setmcp_mem")
         call alloc (coeff_37, total_contribution_count, "COEFF_37", "setmcp_mem")
      else if(nfile == 38) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.38:"
         call alloc (lab_38,   total_secondary_contribution_count,"LAB_38",   "setmcp_mem")
         call alloc (ncontr_38,total_secondary_contribution_count,"NCONTR_38","setmcp_mem")
         call alloc (iclmn_38, total_contribution_count, "ICLMN_38", "setmcp_mem")
         call alloc (indx_38,  total_contribution_count, "INDX_38",  "setmcp_mem")
         call alloc (coeff_38, total_contribution_count, "COEFF_38", "setmcp_mem")
      else if(nfile == 39) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.39:"
         call alloc (lab_39,   total_secondary_contribution_count,"LAB_39",   "setmcp_mem")
         call alloc (ncontr_39,total_secondary_contribution_count,"NCONTR_39","setmcp_mem")
         call alloc (iclmn_39, total_contribution_count, "ICLMN_39", "setmcp_mem")
         call alloc (indx_39,  total_contribution_count, "INDX_39",  "setmcp_mem")
         call alloc (coeff_39, total_contribution_count, "COEFF_39", "setmcp_mem")
      else if(nfile == 40) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.40:"
         call alloc (lab_40,   total_secondary_contribution_count,"LAB_40",   "setmcp_mem")
         call alloc (ncontr_40,total_secondary_contribution_count,"NCONTR_40","setmcp_mem")
         call alloc (iclmn_40, total_contribution_count, "ICLMN_40", "setmcp_mem")
         call alloc (indx_40,  total_contribution_count, "INDX_40",  "setmcp_mem")
         call alloc (coeff_40, total_contribution_count, "COEFF_40", "setmcp_mem")
      else if(nfile == 41) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.41:"
         call alloc (lab_41,   total_secondary_contribution_count,"LAB_41",   "setmcp_mem")
         call alloc (ncontr_41,total_secondary_contribution_count,"NCONTR_41","setmcp_mem")
         call alloc (iclmn_41, total_contribution_count, "ICLMN_41", "setmcp_mem")
         call alloc (indx_41,  total_contribution_count, "INDX_41",  "setmcp_mem")
         call alloc (coeff_41, total_contribution_count, "COEFF_41", "setmcp_mem")
      else if(nfile == 42) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.42:"
         call alloc (lab_42,   total_secondary_contribution_count,"LAB_42",   "setmcp_mem")
         call alloc (ncontr_42,total_secondary_contribution_count,"NCONTR_42","setmcp_mem")
         call alloc (iclmn_42, total_contribution_count, "ICLMN_42", "setmcp_mem")
         call alloc (indx_42,  total_contribution_count, "INDX_42",  "setmcp_mem")
         call alloc (coeff_42, total_contribution_count, "COEFF_42", "setmcp_mem")
      else if(nfile == 43) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.43:"
         call alloc (lab_43,   total_secondary_contribution_count,"LAB_43",   "setmcp_mem")
         call alloc (ncontr_43,total_secondary_contribution_count,"NCONTR_43","setmcp_mem")
         call alloc (iclmn_43, total_contribution_count, "ICLMN_43", "setmcp_mem")
         call alloc (indx_43,  total_contribution_count, "INDX_43",  "setmcp_mem")
         call alloc (coeff_43, total_contribution_count, "COEFF_43", "setmcp_mem")
      else if(nfile == 44) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.44:"
         call alloc (lab_44,   total_secondary_contribution_count,"LAB_44",   "setmcp_mem")
         call alloc (ncontr_44,total_secondary_contribution_count,"NCONTR_44","setmcp_mem")
         call alloc (iclmn_44, total_contribution_count, "ICLMN_44", "setmcp_mem")
         call alloc (indx_44,  total_contribution_count, "INDX_44",  "setmcp_mem")
         call alloc (coeff_44, total_contribution_count, "COEFF_44", "setmcp_mem")
      else if(nfile == 45) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.45:"
         call alloc (lab_45,   total_secondary_contribution_count,"LAB_45",   "setmcp_mem")
         call alloc (ncontr_45,total_secondary_contribution_count,"NCONTR_45","setmcp_mem")
         call alloc (iclmn_45, total_contribution_count, "ICLMN_45", "setmcp_mem")
         call alloc (indx_45,  total_contribution_count, "INDX_45",  "setmcp_mem")
         call alloc (coeff_45, total_contribution_count, "COEFF_45", "setmcp_mem")
      else if(nfile == 46) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.46:"
         call alloc (lab_46,   total_secondary_contribution_count,"LAB_46",   "setmcp_mem")
         call alloc (ncontr_46,total_secondary_contribution_count,"NCONTR_46","setmcp_mem")
         call alloc (iclmn_46, total_contribution_count, "ICLMN_46", "setmcp_mem")
         call alloc (indx_46,  total_contribution_count, "INDX_46",  "setmcp_mem")
         call alloc (coeff_46, total_contribution_count, "COEFF_46", "setmcp_mem")
      else if(nfile == 47) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.47:"
         call alloc (lab_47,   total_secondary_contribution_count,"LAB_47",   "setmcp_mem")
         call alloc (ncontr_47,total_secondary_contribution_count,"NCONTR_47","setmcp_mem")
         call alloc (iclmn_47, total_contribution_count, "ICLMN_47", "setmcp_mem")
         call alloc (indx_47,  total_contribution_count, "INDX_47",  "setmcp_mem")
         call alloc (coeff_47, total_contribution_count, "COEFF_47", "setmcp_mem")
      else if(nfile == 48) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.48:"
         call alloc (lab_48,   total_secondary_contribution_count,"LAB_48",   "setmcp_mem")
         call alloc (ncontr_48,total_secondary_contribution_count,"NCONTR_48","setmcp_mem")
         call alloc (iclmn_48, total_contribution_count, "ICLMN_48", "setmcp_mem")
         call alloc (indx_48,  total_contribution_count, "INDX_48",  "setmcp_mem")
         call alloc (coeff_48, total_contribution_count, "COEFF_48", "setmcp_mem")
      else if(nfile == 49) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.49:"
         call alloc (lab_49,   total_secondary_contribution_count,"LAB_49",   "setmcp_mem")
         call alloc (ncontr_49,total_secondary_contribution_count,"NCONTR_49","setmcp_mem")
         call alloc (iclmn_49, total_contribution_count, "ICLMN_49", "setmcp_mem")
         call alloc (indx_49,  total_contribution_count, "INDX_49",  "setmcp_mem")
         call alloc (coeff_49, total_contribution_count, "COEFF_49", "setmcp_mem")
      else if(nfile == 50) then
         if(mpi_rank == 0) write(*,"(A27)") "Allocation for mcp.50:"
         call alloc (lab_50,   total_secondary_contribution_count,"LAB_50",   "setmcp_mem")
         call alloc (ncontr_50,total_secondary_contribution_count,"NCONTR_50","setmcp_mem")
         call alloc (iclmn_50, total_contribution_count, "ICLMN_50", "setmcp_mem")
         call alloc (indx_50,  total_contribution_count, "INDX_50",  "setmcp_mem")
         call alloc (coeff_50, total_contribution_count, "COEFF_50", "setmcp_mem")
      else
         print*, "Error in setmcp_mem"
         call fail_input_io(0, mcp_input_units(nfile), "MCP setmcp_mem invalid block or channel")
      end if
      rewind (mcp_input_units(nfile), iostat=io_status)
      if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem rewind")
      if(nfile == 30) then
         read (mcp_input_units(nfile), iostat=io_status)
         if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem skip header or payload")
         read (mcp_input_units(nfile), iostat=io_status)
         if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem skip header or payload")
         read (mcp_input_units(nfile), iostat=io_status)
         if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem skip header or payload")
      end if
      read (mcp_input_units(nfile), iostat=io_status)
      if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem skip header or payload")
      read (mcp_input_units(nfile), iostat=io_status)
      if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem skip header or payload")
      read (mcp_input_units(nfile), iostat=io_status)
      if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem skip header or payload")
!
      total_contribution_count = 0
      total_secondary_contribution_count = 0
      do symmetry_block_index = 1, nblock
        if(nfile == 30) then
          !*** Read in IROW from file mcp.30 ***
          read (mcp_input_units(30), iostat=io_status) mcp_channel_label, jblockt, block_csf_count
          if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(30), "MCP setmcp_mem block or channel header")
          if (jblockt /= symmetry_block_index) then
            write (stderr_unit, *) "setmcp_mem", "3: jblockt .NE. jblock"
            call fail_input_io(0, mcp_input_units(nfile), "MCP setmcp_mem invalid block or channel")
          end if
        else
          read (mcp_input_units(nfile), iostat=io_status) mcp_channel_label, jblockt, ncft, ncoeff
          if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem block or channel header")
          if (jblockt /= symmetry_block_index) then
            write (stderr_unit, *) "setmcp_mem", ":4 jblockt .NE. jblock"
            call fail_input_io(0, mcp_input_units(nfile), "MCP setmcp_mem invalid block or channel")
          end if
          total_secondary_contribution_count = total_secondary_contribution_count + 1
        end if

         if(nfile == 30) then
            read (mcp_input_units(30), iostat=io_status) nelmntgg
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(30), "MCP setmcp_mem sparse element count")
            stored_matrix_element_count = int(nelmntgg,kind=int64)
            read (mcp_input_units(30), iostat=io_status) (iendc_30(i),i=                                  &
                     imin1_30(symmetry_block_index)+mpi_rank+1,imax1_30(symmetry_block_index),mpi_process_count), &
                     (irow_30(i),i=imin2_30(symmetry_block_index)+1,imax2_30(symmetry_block_index))
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(30), "MCP setmcp_mem sparse indices")
         else if(nfile == 31) then
            read (mcp_input_units(nfile), iostat=io_status) lab_31(total_secondary_contribution_count),              &
                                  ncontr_31(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_31(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_31(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_31(total_contribution_count+i),                   &
                              indx_31(total_contribution_count+i),                   &
                             coeff_31(total_contribution_count+i),                   &
                        i=1,ncontr_31(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_31(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_31(total_secondary_contribution_count),           &
                                     ncontr_31(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 32) then
            read (mcp_input_units(nfile), iostat=io_status) lab_32(total_secondary_contribution_count),              &
                                  ncontr_32(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_32(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_32(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_32(total_contribution_count+i),                   &
                              indx_32(total_contribution_count+i),                   &
                             coeff_32(total_contribution_count+i),                   &
                        i=1,ncontr_32(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_32(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_32(total_secondary_contribution_count),           &
                                     ncontr_32(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 33) then
            read (mcp_input_units(nfile), iostat=io_status) lab_33(total_secondary_contribution_count),              &
                                  ncontr_33(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_33(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_33(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_33(total_contribution_count+i),                   &
                              indx_33(total_contribution_count+i),                   &
                             coeff_33(total_contribution_count+i),                   &
                        i=1,ncontr_33(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_33(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_33(total_secondary_contribution_count),           &
                                     ncontr_33(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 34) then
            read (mcp_input_units(nfile), iostat=io_status) lab_34(total_secondary_contribution_count),              &
                                  ncontr_34(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_34(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_34(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_34(total_contribution_count+i),                   &
                              indx_34(total_contribution_count+i),                   &
                             coeff_34(total_contribution_count+i),                   &
                        i=1,ncontr_34(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_34(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_34(total_secondary_contribution_count),           &
                                     ncontr_34(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 35) then
            read (mcp_input_units(nfile), iostat=io_status) lab_35(total_secondary_contribution_count),              &
                                  ncontr_35(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_35(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_35(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_35(total_contribution_count+i),                   &
                              indx_35(total_contribution_count+i),                   &
                             coeff_35(total_contribution_count+i),                   &
                        i=1,ncontr_35(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_35(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_35(total_secondary_contribution_count),           &
                                     ncontr_35(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 36) then
            read (mcp_input_units(nfile), iostat=io_status) lab_36(total_secondary_contribution_count),              &
                                  ncontr_36(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_36(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_36(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_36(total_contribution_count+i),                   &
                              indx_36(total_contribution_count+i),                   &
                             coeff_36(total_contribution_count+i),                   &
                        i=1,ncontr_36(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_36(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_36(total_secondary_contribution_count),           &
                                     ncontr_36(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 37) then
            read (mcp_input_units(nfile), iostat=io_status) lab_37(total_secondary_contribution_count),              &
                                  ncontr_37(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_37(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_37(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_37(total_contribution_count+i),                   &
                              indx_37(total_contribution_count+i),                   &
                             coeff_37(total_contribution_count+i),                   &
                        i=1,ncontr_37(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_37(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_37(total_secondary_contribution_count),           &
                                     ncontr_37(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 38) then
            read (mcp_input_units(nfile), iostat=io_status) lab_38(total_secondary_contribution_count),              &
                                  ncontr_38(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_38(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_38(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_38(total_contribution_count+i),                   &
                              indx_38(total_contribution_count+i),                   &
                             coeff_38(total_contribution_count+i),                   &
                        i=1,ncontr_38(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_38(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_38(total_secondary_contribution_count),           &
                                     ncontr_38(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 39) then
            read (mcp_input_units(nfile), iostat=io_status) lab_39(total_secondary_contribution_count),              &
                                  ncontr_39(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_39(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_39(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_39(total_contribution_count+i),                   &
                              indx_39(total_contribution_count+i),                   &
                             coeff_39(total_contribution_count+i),                   &
                        i=1,ncontr_39(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_39(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_39(total_secondary_contribution_count),           &
                                     ncontr_39(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 40) then
            read (mcp_input_units(nfile), iostat=io_status) lab_40(total_secondary_contribution_count),              &
                                  ncontr_40(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_40(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_40(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_40(total_contribution_count+i),                   &
                              indx_40(total_contribution_count+i),                   &
                             coeff_40(total_contribution_count+i),                   &
                        i=1,ncontr_40(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_40(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_40(total_secondary_contribution_count),           &
                                     ncontr_40(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 41) then
            read (mcp_input_units(nfile), iostat=io_status) lab_41(total_secondary_contribution_count),              &
                                  ncontr_41(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_41(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_41(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_41(total_contribution_count+i),                   &
                              indx_41(total_contribution_count+i),                   &
                             coeff_41(total_contribution_count+i),                   &
                        i=1,ncontr_41(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_41(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_41(total_secondary_contribution_count),           &
                                     ncontr_41(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 42) then
            read (mcp_input_units(nfile), iostat=io_status) lab_42(total_secondary_contribution_count),              &
                                  ncontr_42(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_42(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_42(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_42(total_contribution_count+i),                   &
                              indx_42(total_contribution_count+i),                   &
                             coeff_42(total_contribution_count+i),                   &
                        i=1,ncontr_42(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_42(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_42(total_secondary_contribution_count),           &
                                     ncontr_42(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 43) then
            read (mcp_input_units(nfile), iostat=io_status) lab_43(total_secondary_contribution_count),              &
                                  ncontr_43(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_43(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_43(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_43(total_contribution_count+i),                   &
                              indx_43(total_contribution_count+i),                   &
                             coeff_43(total_contribution_count+i),                   &
                        i=1,ncontr_43(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_43(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_43(total_secondary_contribution_count),           &
                                     ncontr_43(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 44) then
            read (mcp_input_units(nfile), iostat=io_status) lab_44(total_secondary_contribution_count),              &
                                  ncontr_44(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_44(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_44(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_44(total_contribution_count+i),                   &
                              indx_44(total_contribution_count+i),                   &
                             coeff_44(total_contribution_count+i),                   &
                        i=1,ncontr_44(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_44(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_44(total_secondary_contribution_count),           &
                                     ncontr_44(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 45) then
            read (mcp_input_units(nfile), iostat=io_status) lab_45(total_secondary_contribution_count),              &
                                  ncontr_45(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_45(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_45(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_45(total_contribution_count+i),                   &
                              indx_45(total_contribution_count+i),                   &
                             coeff_45(total_contribution_count+i),                   &
                        i=1,ncontr_45(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_45(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_45(total_secondary_contribution_count),           &
                                     ncontr_45(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 46) then
            read (mcp_input_units(nfile), iostat=io_status) lab_46(total_secondary_contribution_count),              &
                                  ncontr_46(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_46(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_46(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_46(total_contribution_count+i),                   &
                              indx_46(total_contribution_count+i),                   &
                             coeff_46(total_contribution_count+i),                   &
                        i=1,ncontr_46(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_46(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_46(total_secondary_contribution_count),           &
                                     ncontr_46(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 47) then
            read (mcp_input_units(nfile), iostat=io_status) lab_47(total_secondary_contribution_count),              &
                                  ncontr_47(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_47(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_47(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_47(total_contribution_count+i),                   &
                              indx_47(total_contribution_count+i),                   &
                             coeff_47(total_contribution_count+i),                   &
                        i=1,ncontr_47(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_47(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_47(total_secondary_contribution_count),           &
                                     ncontr_47(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 48) then
            read (mcp_input_units(nfile), iostat=io_status) lab_48(total_secondary_contribution_count),              &
                                  ncontr_48(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_48(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_48(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_48(total_contribution_count+i),                   &
                              indx_48(total_contribution_count+i),                   &
                             coeff_48(total_contribution_count+i),                   &
                        i=1,ncontr_48(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_48(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_48(total_secondary_contribution_count),           &
                                     ncontr_48(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 49) then
            read (mcp_input_units(nfile), iostat=io_status) lab_49(total_secondary_contribution_count),              &
                                  ncontr_49(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_49(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_49(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_49(total_contribution_count+i),                   &
                              indx_49(total_contribution_count+i),                   &
                             coeff_49(total_contribution_count+i),                   &
                        i=1,ncontr_49(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_49(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_49(total_secondary_contribution_count),           &
                                     ncontr_49(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else if(nfile == 50) then
            read (mcp_input_units(nfile), iostat=io_status) lab_50(total_secondary_contribution_count),              &
                                  ncontr_50(total_secondary_contribution_count)
            if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            do while (lab_50(total_secondary_contribution_count) /= 0 .or.                    &
                   ncontr_50(total_secondary_contribution_count) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (mcp_input_units(nfile), iostat=io_status) (iclmn_50(total_contribution_count+i),                   &
                              indx_50(total_contribution_count+i),                   &
                             coeff_50(total_contribution_count+i),                   &
                        i=1,ncontr_50(total_secondary_contribution_count))
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem coefficient payload")
               total_contribution_count  = total_contribution_count + ncontr_50(total_secondary_contribution_count)
               total_secondary_contribution_count = total_secondary_contribution_count + 1
               read (mcp_input_units(nfile), iostat=io_status) lab_50(total_secondary_contribution_count),           &
                                     ncontr_50(total_secondary_contribution_count)
               if (io_status /= 0) call fail_input_io(io_status, mcp_input_units(nfile), "MCP setmcp_mem load-pass label/count")
            end do
         else
            print*, "Error in setmcp_mem"
            call fail_input_io(0, mcp_input_units(nfile), "MCP setmcp_mem invalid block or channel")
         end if

      end do
      sysmemf = systemfreemem()  ! GB
      if(mpi_rank == 0)                                                    &
      write(*,"(A30,F10.2,A3)") "Free memory on computer",sysmemf," Gb"
      close(mcp_input_units(nfile))

      return
      end subroutine setmcp_mem
