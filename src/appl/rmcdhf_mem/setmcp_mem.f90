!*******************************************************************
!                                                                  *
      subroutine setmcp_mem(nblock,nfile)
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use memory_man
      use input_io, only: fail_input_io
      use rmcdhf_mem_c
      use iounit_c
      use mpi_s
      use hmat_c,          only: stored_matrix_element_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use systemmem_i
      use systemfreemem_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nblock, nfile
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character :: mcplab*3
      real(real64) :: coeff_dumm
      real(real64) :: sysmemt, sysmemf
      integer :: i, jblock, jblockt, ncft, ncoeff, ios, ncf
      integer(int64) :: ncontr_tot, ncontr_tot2
      integer :: iclmn_dumm, lab_dumm, ncontr_dumm, indx_dumm
      integer :: nelmntgg, iendcdum, irow_dumm
!-----------------------------------------------

      ncontr_tot = 0
      ncontr_tot2 = 0
      do jblock = 1, nblock
        if(nfile == 30) then
          !*** Read in IROW from file mcp.30 ***
          read (30, iostat=ios) mcplab, jblockt, ncf
          if (ios /= 0) call fail_input_io(ios, 30, 'MCP setmcp_mem block or channel header')
          if (jblockt /= jblock) then
            write (stderr_unit, *) 'setmcp_mem', '1: jblockt .NE. jblock'
            call fail_input_io(0, nfile, 'MCP setmcp_mem invalid block or channel')
          endif
          read (30, iostat=ios) nelmntgg
          if (ios /= 0) call fail_input_io(ios, 30, 'MCP setmcp_mem sparse element count')
          stored_matrix_element_count = int(nelmntgg,kind=int64)
          if(jblock == 1) then
            call alloc (imin1_30,nblock,'IMIN1_30','setmcp_mem')
            call alloc (imin2_30,nblock,'IMIN2_30','setmcp_mem')
            call alloc (imax1_30,nblock,'IMAX1_30','setmcp_mem')
            call alloc (imax2_30,nblock,'IMAX2_30','setmcp_mem')
            imin1_30(jblock) = 0
            imax1_30(jblock) = ncf
            imin2_30(jblock) = 0
            imax2_30(jblock) = stored_matrix_element_count
          else
            imin1_30(jblock) = imax1_30(jblock-1) + 0
            imax1_30(jblock) = imax1_30(jblock-1) + ncf
            imin2_30(jblock) = imax2_30(jblock-1) + 0
            imax2_30(jblock) = imax2_30(jblock-1) + stored_matrix_element_count
          end if
          read (30, iostat=ios) (iendcdum,i=mpi_rank + 1,ncf,mpi_process_count),                  &
                                                 (irow_dumm,i=1,stored_matrix_element_count)
          if (ios /= 0) call fail_input_io(ios, 30, 'MCP setmcp_mem sparse indices')
        else
          read (nfile, iostat=ios) mcplab, jblockt, ncft, ncoeff
          if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem block or channel header')
          if (jblockt /= jblock) then
            write (stderr_unit, *) 'setmcp_mem', ':2 jblockt .NE. jblock'
            call fail_input_io(0, nfile, 'MCP setmcp_mem invalid block or channel')
          endif

          ncontr_tot2 = ncontr_tot2 + 1
          read (nfile, iostat=ios) lab_dumm, ncontr_dumm
          if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem count-pass label/count')
          do while (lab_dumm /= 0 .or. ncontr_dumm /= 0)
            !*** Read column index, sparse matrix index, and
            !*** coefficient for all contributions from this integral
            read (nfile, iostat=ios) (iclmn_dumm,indx_dumm,coeff_dumm,             &
                                                       i=1,ncontr_dumm)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem count-pass coefficients')
            ncontr_tot = ncontr_tot + ncontr_dumm
            ncontr_tot2 = ncontr_tot2 + 1
            read (nfile, iostat=ios) lab_dumm, ncontr_dumm
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem count-pass label/count')
          end do
        end if
      end do

      if(nfile == 30) then
         ! Total physical memory
         sysmemt = systemmem()  ! GB
         sysmemf = systemfreemem()  ! GB
         write(*,'(A53)')                                              &
         '-----------------------------------------------------'
         write(*,'(A53)')                                              &
         'Spin-angular coefficient are putting into the memory:'
         write(*,'(A53)')                                              &
         '-----------------------------------------------------'
         write(*,'(A1)') ''
         write(*,'(A30,F10.2,A3)')                                     &
                          'Total memory on computer:',sysmemt,' Gb'
         write(*,'(A30,F10.2,A3)')                       &
                          'Free  memory on computer:',sysmemf,' Gb'
         write(*,'(A1)') ''
         write(*,'(A27)') 'Allocation for mcp.30:'
         call alloc (iendc_30,imax1_30(nblock),'IENDC_30','setmcp_mem')
         call alloc (irow_30, imax2_30(nblock),'IROW_30', 'setmcp_mem')
      else if(nfile == 31) then
         write(*,'(A27)') 'Allocation for mcp.31:'
         call alloc (lab_31,   ncontr_tot2,'LAB_31',   'setmcp_mem')
         call alloc (ncontr_31,ncontr_tot2,'NCONTR_31','setmcp_mem')
         call alloc (iclmn_31, ncontr_tot, 'ICLMN_31', 'setmcp_mem')
         call alloc (indx_31,  ncontr_tot, 'INDX_31',  'setmcp_mem')
         call alloc (coeff_31, ncontr_tot, 'COEFF_31', 'setmcp_mem')
      else if(nfile == 32) then
         write(*,'(A27)') 'Allocation for mcp.32:'
         call alloc (lab_32,   ncontr_tot2,'LAB_32',   'setmcp_mem')
         call alloc (ncontr_32,ncontr_tot2,'NCONTR_32','setmcp_mem')
         call alloc (iclmn_32, ncontr_tot, 'ICLMN_32', 'setmcp_mem')
         call alloc (indx_32,  ncontr_tot, 'INDX_32',  'setmcp_mem')
         call alloc (coeff_32, ncontr_tot, 'COEFF_32', 'setmcp_mem')
      else if(nfile == 33) then
         write(*,'(A27)') 'Allocation for mcp.33:'
         call alloc (lab_33,   ncontr_tot2,'LAB_33',   'setmcp_mem')
         call alloc (ncontr_33,ncontr_tot2,'NCONTR_33','setmcp_mem')
         call alloc (iclmn_33, ncontr_tot, 'ICLMN_33', 'setmcp_mem')
         call alloc (indx_33,  ncontr_tot, 'INDX_33',  'setmcp_mem')
         call alloc (coeff_33, ncontr_tot, 'COEFF_33', 'setmcp_mem')
      else if(nfile == 34) then
         write(*,'(A27)') 'Allocation for mcp.34:'
         call alloc (lab_34,   ncontr_tot2,'LAB_34',   'setmcp_mem')
         call alloc (ncontr_34,ncontr_tot2,'NCONTR_34','setmcp_mem')
         call alloc (iclmn_34, ncontr_tot, 'ICLMN_34', 'setmcp_mem')
         call alloc (indx_34,  ncontr_tot, 'INDX_34',  'setmcp_mem')
         call alloc (coeff_34, ncontr_tot, 'COEFF_34', 'setmcp_mem')
      else if(nfile == 35) then
         write(*,'(A27)') 'Allocation for mcp.35:'
         call alloc (lab_35,   ncontr_tot2,'LAB_35',   'setmcp_mem')
         call alloc (ncontr_35,ncontr_tot2,'NCONTR_35','setmcp_mem')
         call alloc (iclmn_35, ncontr_tot, 'ICLMN_35', 'setmcp_mem')
         call alloc (indx_35,  ncontr_tot, 'INDX_35',  'setmcp_mem')
         call alloc (coeff_35, ncontr_tot, 'COEFF_35', 'setmcp_mem')
      else if(nfile == 36) then
         write(*,'(A27)') 'Allocation for mcp.36:'
         call alloc (lab_36,   ncontr_tot2,'LAB_36',   'setmcp_mem')
         call alloc (ncontr_36,ncontr_tot2,'NCONTR_36','setmcp_mem')
         call alloc (iclmn_36, ncontr_tot, 'ICLMN_36', 'setmcp_mem')
         call alloc (indx_36,  ncontr_tot, 'INDX_36',  'setmcp_mem')
         call alloc (coeff_36, ncontr_tot, 'COEFF_36', 'setmcp_mem')
      else if(nfile == 37) then
         write(*,'(A27)') 'Allocation for mcp.37:'
         call alloc (lab_37,   ncontr_tot2,'LAB_37',   'setmcp_mem')
         call alloc (ncontr_37,ncontr_tot2,'NCONTR_37','setmcp_mem')
         call alloc (iclmn_37, ncontr_tot, 'ICLMN_37', 'setmcp_mem')
         call alloc (indx_37,  ncontr_tot, 'INDX_37',  'setmcp_mem')
         call alloc (coeff_37, ncontr_tot, 'COEFF_37', 'setmcp_mem')
      else if(nfile == 38) then
         write(*,'(A27)') 'Allocation for mcp.38:'
         call alloc (lab_38,   ncontr_tot2,'LAB_38',   'setmcp_mem')
         call alloc (ncontr_38,ncontr_tot2,'NCONTR_38','setmcp_mem')
         call alloc (iclmn_38, ncontr_tot, 'ICLMN_38', 'setmcp_mem')
         call alloc (indx_38,  ncontr_tot, 'INDX_38',  'setmcp_mem')
         call alloc (coeff_38, ncontr_tot, 'COEFF_38', 'setmcp_mem')
      else if(nfile == 39) then
         write(*,'(A27)') 'Allocation for mcp.39:'
         call alloc (lab_39,   ncontr_tot2,'LAB_39',   'setmcp_mem')
         call alloc (ncontr_39,ncontr_tot2,'NCONTR_39','setmcp_mem')
         call alloc (iclmn_39, ncontr_tot, 'ICLMN_39', 'setmcp_mem')
         call alloc (indx_39,  ncontr_tot, 'INDX_39',  'setmcp_mem')
         call alloc (coeff_39, ncontr_tot, 'COEFF_39', 'setmcp_mem')
      else if(nfile == 40) then
         write(*,'(A27)') 'Allocation for mcp.40:'
         call alloc (lab_40,   ncontr_tot2,'LAB_40',   'setmcp_mem')
         call alloc (ncontr_40,ncontr_tot2,'NCONTR_40','setmcp_mem')
         call alloc (iclmn_40, ncontr_tot, 'ICLMN_40', 'setmcp_mem')
         call alloc (indx_40,  ncontr_tot, 'INDX_40',  'setmcp_mem')
         call alloc (coeff_40, ncontr_tot, 'COEFF_40', 'setmcp_mem')
      else if(nfile == 41) then
         write(*,'(A27)') 'Allocation for mcp.41:'
         call alloc (lab_41,   ncontr_tot2,'LAB_41',   'setmcp_mem')
         call alloc (ncontr_41,ncontr_tot2,'NCONTR_41','setmcp_mem')
         call alloc (iclmn_41, ncontr_tot, 'ICLMN_41', 'setmcp_mem')
         call alloc (indx_41,  ncontr_tot, 'INDX_41',  'setmcp_mem')
         call alloc (coeff_41, ncontr_tot, 'COEFF_41', 'setmcp_mem')
      else if(nfile == 42) then
         write(*,'(A27)') 'Allocation for mcp.42:'
         call alloc (lab_42,   ncontr_tot2,'LAB_42',   'setmcp_mem')
         call alloc (ncontr_42,ncontr_tot2,'NCONTR_42','setmcp_mem')
         call alloc (iclmn_42, ncontr_tot, 'ICLMN_42', 'setmcp_mem')
         call alloc (indx_42,  ncontr_tot, 'INDX_42',  'setmcp_mem')
         call alloc (coeff_42, ncontr_tot, 'COEFF_42', 'setmcp_mem')
      else if(nfile == 43) then
         write(*,'(A27)') 'Allocation for mcp.43:'
         call alloc (lab_43,   ncontr_tot2,'LAB_43',   'setmcp_mem')
         call alloc (ncontr_43,ncontr_tot2,'NCONTR_43','setmcp_mem')
         call alloc (iclmn_43, ncontr_tot, 'ICLMN_43', 'setmcp_mem')
         call alloc (indx_43,  ncontr_tot, 'INDX_43',  'setmcp_mem')
         call alloc (coeff_43, ncontr_tot, 'COEFF_43', 'setmcp_mem')
      else if(nfile == 44) then
         write(*,'(A27)') 'Allocation for mcp.44:'
         call alloc (lab_44,   ncontr_tot2,'LAB_44',   'setmcp_mem')
         call alloc (ncontr_44,ncontr_tot2,'NCONTR_44','setmcp_mem')
         call alloc (iclmn_44, ncontr_tot, 'ICLMN_44', 'setmcp_mem')
         call alloc (indx_44,  ncontr_tot, 'INDX_44',  'setmcp_mem')
         call alloc (coeff_44, ncontr_tot, 'COEFF_44', 'setmcp_mem')
      else if(nfile == 45) then
         write(*,'(A27)') 'Allocation for mcp.45:'
         call alloc (lab_45,   ncontr_tot2,'LAB_45',   'setmcp_mem')
         call alloc (ncontr_45,ncontr_tot2,'NCONTR_45','setmcp_mem')
         call alloc (iclmn_45, ncontr_tot, 'ICLMN_45', 'setmcp_mem')
         call alloc (indx_45,  ncontr_tot, 'INDX_45',  'setmcp_mem')
         call alloc (coeff_45, ncontr_tot, 'COEFF_45', 'setmcp_mem')
      else if(nfile == 46) then
         write(*,'(A27)') 'Allocation for mcp.46:'
         call alloc (lab_46,   ncontr_tot2,'LAB_46',   'setmcp_mem')
         call alloc (ncontr_46,ncontr_tot2,'NCONTR_46','setmcp_mem')
         call alloc (iclmn_46, ncontr_tot, 'ICLMN_46', 'setmcp_mem')
         call alloc (indx_46,  ncontr_tot, 'INDX_46',  'setmcp_mem')
         call alloc (coeff_46, ncontr_tot, 'COEFF_46', 'setmcp_mem')
      else if(nfile == 47) then
         write(*,'(A27)') 'Allocation for mcp.47:'
         call alloc (lab_47,   ncontr_tot2,'LAB_47',   'setmcp_mem')
         call alloc (ncontr_47,ncontr_tot2,'NCONTR_47','setmcp_mem')
         call alloc (iclmn_47, ncontr_tot, 'ICLMN_47', 'setmcp_mem')
         call alloc (indx_47,  ncontr_tot, 'INDX_47',  'setmcp_mem')
         call alloc (coeff_47, ncontr_tot, 'COEFF_47', 'setmcp_mem')
      else if(nfile == 48) then
         write(*,'(A27)') 'Allocation for mcp.48:'
         call alloc (lab_48,   ncontr_tot2,'LAB_48',   'setmcp_mem')
         call alloc (ncontr_48,ncontr_tot2,'NCONTR_48','setmcp_mem')
         call alloc (iclmn_48, ncontr_tot, 'ICLMN_48', 'setmcp_mem')
         call alloc (indx_48,  ncontr_tot, 'INDX_48',  'setmcp_mem')
         call alloc (coeff_48, ncontr_tot, 'COEFF_48', 'setmcp_mem')
      else if(nfile == 49) then
         write(*,'(A27)') 'Allocation for mcp.49:'
         call alloc (lab_49,   ncontr_tot2,'LAB_49',   'setmcp_mem')
         call alloc (ncontr_49,ncontr_tot2,'NCONTR_49','setmcp_mem')
         call alloc (iclmn_49, ncontr_tot, 'ICLMN_49', 'setmcp_mem')
         call alloc (indx_49,  ncontr_tot, 'INDX_49',  'setmcp_mem')
         call alloc (coeff_49, ncontr_tot, 'COEFF_49', 'setmcp_mem')
      else if(nfile == 50) then
         write(*,'(A27)') 'Allocation for mcp.50:'
         call alloc (lab_50,   ncontr_tot2,'LAB_50',   'setmcp_mem')
         call alloc (ncontr_50,ncontr_tot2,'NCONTR_50','setmcp_mem')
         call alloc (iclmn_50, ncontr_tot, 'ICLMN_50', 'setmcp_mem')
         call alloc (indx_50,  ncontr_tot, 'INDX_50',  'setmcp_mem')
         call alloc (coeff_50, ncontr_tot, 'COEFF_50', 'setmcp_mem')
      else
         print*, "Error in setmcp_mem"
         call fail_input_io(0, nfile, 'MCP setmcp_mem invalid block or channel')
      end if
      rewind (nfile, iostat=ios)
      if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem rewind')
      if(nfile == 30) then
         read (nfile, iostat=ios)
         if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem skip header or payload')
         read (nfile, iostat=ios)
         if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem skip header or payload')
         read (nfile, iostat=ios)
         if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem skip header or payload')
      end if
      read (nfile, iostat=ios)
      if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem skip header or payload')
      read (nfile, iostat=ios)
      if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem skip header or payload')
      read (nfile, iostat=ios)
      if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem skip header or payload')
!
      ncontr_tot = 0
      ncontr_tot2 = 0
      do jblock = 1, nblock
        if(nfile == 30) then
          !*** Read in IROW from file mcp.30 ***
          read (30, iostat=ios) mcplab, jblockt, ncf
          if (ios /= 0) call fail_input_io(ios, 30, 'MCP setmcp_mem block or channel header')
          if (jblockt /= jblock) then
            write (stderr_unit, *) 'setmcp_mem', '3: jblockt .NE. jblock'
            call fail_input_io(0, nfile, 'MCP setmcp_mem invalid block or channel')
          endif
        else
          read (nfile, iostat=ios) mcplab, jblockt, ncft, ncoeff
          if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem block or channel header')
          if (jblockt /= jblock) then
            write (stderr_unit, *) 'setmcp_mem', ':4 jblockt .NE. jblock'
            call fail_input_io(0, nfile, 'MCP setmcp_mem invalid block or channel')
          endif
          ncontr_tot2 = ncontr_tot2 + 1
        endif

         if(nfile == 30) then
            read (30, iostat=ios) nelmntgg
            if (ios /= 0) call fail_input_io(ios, 30, 'MCP setmcp_mem sparse element count')
            stored_matrix_element_count = int(nelmntgg,kind=int64)
            read (30, iostat=ios) (iendc_30(i),i=                                  &
                     imin1_30(jblock)+mpi_rank+1,imax1_30(jblock),mpi_process_count), &
                     (irow_30(i),i=imin2_30(jblock)+1,imax2_30(jblock))
            if (ios /= 0) call fail_input_io(ios, 30, 'MCP setmcp_mem sparse indices')
         else if(nfile == 31) then
            read (nfile, iostat=ios) lab_31(ncontr_tot2),              &
                                  ncontr_31(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_31(ncontr_tot2) /= 0 .or.                    &
                   ncontr_31(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_31(ncontr_tot+i),                   &
                              indx_31(ncontr_tot+i),                   &
                             coeff_31(ncontr_tot+i),                   &
                        i=1,ncontr_31(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_31(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_31(ncontr_tot2),           &
                                     ncontr_31(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 32) then
            read (nfile, iostat=ios) lab_32(ncontr_tot2),              &
                                  ncontr_32(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_32(ncontr_tot2) /= 0 .or.                    &
                   ncontr_32(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_32(ncontr_tot+i),                   &
                              indx_32(ncontr_tot+i),                   &
                             coeff_32(ncontr_tot+i),                   &
                        i=1,ncontr_32(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_32(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_32(ncontr_tot2),           &
                                     ncontr_32(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 33) then
            read (nfile, iostat=ios) lab_33(ncontr_tot2),              &
                                  ncontr_33(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_33(ncontr_tot2) /= 0 .or.                    &
                   ncontr_33(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_33(ncontr_tot+i),                   &
                              indx_33(ncontr_tot+i),                   &
                             coeff_33(ncontr_tot+i),                   &
                        i=1,ncontr_33(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_33(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_33(ncontr_tot2),           &
                                     ncontr_33(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 34) then
            read (nfile, iostat=ios) lab_34(ncontr_tot2),              &
                                  ncontr_34(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_34(ncontr_tot2) /= 0 .or.                    &
                   ncontr_34(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_34(ncontr_tot+i),                   &
                              indx_34(ncontr_tot+i),                   &
                             coeff_34(ncontr_tot+i),                   &
                        i=1,ncontr_34(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_34(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_34(ncontr_tot2),           &
                                     ncontr_34(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 35) then
            read (nfile, iostat=ios) lab_35(ncontr_tot2),              &
                                  ncontr_35(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_35(ncontr_tot2) /= 0 .or.                    &
                   ncontr_35(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_35(ncontr_tot+i),                   &
                              indx_35(ncontr_tot+i),                   &
                             coeff_35(ncontr_tot+i),                   &
                        i=1,ncontr_35(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_35(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_35(ncontr_tot2),           &
                                     ncontr_35(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 36) then
            read (nfile, iostat=ios) lab_36(ncontr_tot2),              &
                                  ncontr_36(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_36(ncontr_tot2) /= 0 .or.                    &
                   ncontr_36(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_36(ncontr_tot+i),                   &
                              indx_36(ncontr_tot+i),                   &
                             coeff_36(ncontr_tot+i),                   &
                        i=1,ncontr_36(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_36(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_36(ncontr_tot2),           &
                                     ncontr_36(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 37) then
            read (nfile, iostat=ios) lab_37(ncontr_tot2),              &
                                  ncontr_37(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_37(ncontr_tot2) /= 0 .or.                    &
                   ncontr_37(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_37(ncontr_tot+i),                   &
                              indx_37(ncontr_tot+i),                   &
                             coeff_37(ncontr_tot+i),                   &
                        i=1,ncontr_37(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_37(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_37(ncontr_tot2),           &
                                     ncontr_37(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 38) then
            read (nfile, iostat=ios) lab_38(ncontr_tot2),              &
                                  ncontr_38(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_38(ncontr_tot2) /= 0 .or.                    &
                   ncontr_38(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_38(ncontr_tot+i),                   &
                              indx_38(ncontr_tot+i),                   &
                             coeff_38(ncontr_tot+i),                   &
                        i=1,ncontr_38(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_38(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_38(ncontr_tot2),           &
                                     ncontr_38(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 39) then
            read (nfile, iostat=ios) lab_39(ncontr_tot2),              &
                                  ncontr_39(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_39(ncontr_tot2) /= 0 .or.                    &
                   ncontr_39(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_39(ncontr_tot+i),                   &
                              indx_39(ncontr_tot+i),                   &
                             coeff_39(ncontr_tot+i),                   &
                        i=1,ncontr_39(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_39(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_39(ncontr_tot2),           &
                                     ncontr_39(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 40) then
            read (nfile, iostat=ios) lab_40(ncontr_tot2),              &
                                  ncontr_40(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_40(ncontr_tot2) /= 0 .or.                    &
                   ncontr_40(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_40(ncontr_tot+i),                   &
                              indx_40(ncontr_tot+i),                   &
                             coeff_40(ncontr_tot+i),                   &
                        i=1,ncontr_40(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_40(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_40(ncontr_tot2),           &
                                     ncontr_40(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 41) then
            read (nfile, iostat=ios) lab_41(ncontr_tot2),              &
                                  ncontr_41(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_41(ncontr_tot2) /= 0 .or.                    &
                   ncontr_41(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_41(ncontr_tot+i),                   &
                              indx_41(ncontr_tot+i),                   &
                             coeff_41(ncontr_tot+i),                   &
                        i=1,ncontr_41(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_41(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_41(ncontr_tot2),           &
                                     ncontr_41(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 42) then
            read (nfile, iostat=ios) lab_42(ncontr_tot2),              &
                                  ncontr_42(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_42(ncontr_tot2) /= 0 .or.                    &
                   ncontr_42(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_42(ncontr_tot+i),                   &
                              indx_42(ncontr_tot+i),                   &
                             coeff_42(ncontr_tot+i),                   &
                        i=1,ncontr_42(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_42(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_42(ncontr_tot2),           &
                                     ncontr_42(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 43) then
            read (nfile, iostat=ios) lab_43(ncontr_tot2),              &
                                  ncontr_43(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_43(ncontr_tot2) /= 0 .or.                    &
                   ncontr_43(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_43(ncontr_tot+i),                   &
                              indx_43(ncontr_tot+i),                   &
                             coeff_43(ncontr_tot+i),                   &
                        i=1,ncontr_43(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_43(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_43(ncontr_tot2),           &
                                     ncontr_43(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 44) then
            read (nfile, iostat=ios) lab_44(ncontr_tot2),              &
                                  ncontr_44(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_44(ncontr_tot2) /= 0 .or.                    &
                   ncontr_44(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_44(ncontr_tot+i),                   &
                              indx_44(ncontr_tot+i),                   &
                             coeff_44(ncontr_tot+i),                   &
                        i=1,ncontr_44(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_44(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_44(ncontr_tot2),           &
                                     ncontr_44(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 45) then
            read (nfile, iostat=ios) lab_45(ncontr_tot2),              &
                                  ncontr_45(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_45(ncontr_tot2) /= 0 .or.                    &
                   ncontr_45(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_45(ncontr_tot+i),                   &
                              indx_45(ncontr_tot+i),                   &
                             coeff_45(ncontr_tot+i),                   &
                        i=1,ncontr_45(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_45(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_45(ncontr_tot2),           &
                                     ncontr_45(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 46) then
            read (nfile, iostat=ios) lab_46(ncontr_tot2),              &
                                  ncontr_46(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_46(ncontr_tot2) /= 0 .or.                    &
                   ncontr_46(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_46(ncontr_tot+i),                   &
                              indx_46(ncontr_tot+i),                   &
                             coeff_46(ncontr_tot+i),                   &
                        i=1,ncontr_46(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_46(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_46(ncontr_tot2),           &
                                     ncontr_46(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 47) then
            read (nfile, iostat=ios) lab_47(ncontr_tot2),              &
                                  ncontr_47(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_47(ncontr_tot2) /= 0 .or.                    &
                   ncontr_47(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_47(ncontr_tot+i),                   &
                              indx_47(ncontr_tot+i),                   &
                             coeff_47(ncontr_tot+i),                   &
                        i=1,ncontr_47(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_47(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_47(ncontr_tot2),           &
                                     ncontr_47(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 48) then
            read (nfile, iostat=ios) lab_48(ncontr_tot2),              &
                                  ncontr_48(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_48(ncontr_tot2) /= 0 .or.                    &
                   ncontr_48(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_48(ncontr_tot+i),                   &
                              indx_48(ncontr_tot+i),                   &
                             coeff_48(ncontr_tot+i),                   &
                        i=1,ncontr_48(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_48(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_48(ncontr_tot2),           &
                                     ncontr_48(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 49) then
            read (nfile, iostat=ios) lab_49(ncontr_tot2),              &
                                  ncontr_49(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_49(ncontr_tot2) /= 0 .or.                    &
                   ncontr_49(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_49(ncontr_tot+i),                   &
                              indx_49(ncontr_tot+i),                   &
                             coeff_49(ncontr_tot+i),                   &
                        i=1,ncontr_49(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_49(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_49(ncontr_tot2),           &
                                     ncontr_49(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else if(nfile == 50) then
            read (nfile, iostat=ios) lab_50(ncontr_tot2),              &
                                  ncontr_50(ncontr_tot2)
            if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            do while (lab_50(ncontr_tot2) /= 0 .or.                    &
                   ncontr_50(ncontr_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               read (nfile, iostat=ios) (iclmn_50(ncontr_tot+i),                   &
                              indx_50(ncontr_tot+i),                   &
                             coeff_50(ncontr_tot+i),                   &
                        i=1,ncontr_50(ncontr_tot2))
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem coefficient payload')
               ncontr_tot  = ncontr_tot + ncontr_50(ncontr_tot2)
               ncontr_tot2 = ncontr_tot2 + 1
               read (nfile, iostat=ios) lab_50(ncontr_tot2),           &
                                     ncontr_50(ncontr_tot2)
               if (ios /= 0) call fail_input_io(ios, nfile, 'MCP setmcp_mem load-pass label/count')
            end do
         else
            print*, "Error in setmcp_mem"
            call fail_input_io(0, nfile, 'MCP setmcp_mem invalid block or channel')
         end if

      end do
      sysmemf = systemfreemem()  ! GB
      write(*,'(A30,F10.2,A3)') 'Free memory on computer',sysmemf,' Gb'
      close(nfile)

      return
      end subroutine setmcp_mem
