!***********************************************************************
!                                                                      *
      subroutine sortmem(mcp_channel_index,sort_capacity_per_channel,coefficient_count,distinct_integral_count,print_details, &
         & symmetry_block_index,output_file_prefix)
!                                                                      *
!   This routine sorts lists                                           *
!                                                                      *
!                     (ICLMN,INDEX,LABEL,COEFF)                        *
!   into lists                                                         *
!                     (LABEL,ICLMN,INDEX,COEFF)                        *
!                                                                      *
!   using Heapsort. File NFILE is closed by this routine.  NCOEFF is   *
!   the number of triads (INDEX, ...) and is an input. NTGRAL is the   *
!   number of different values of LABEL, and is an output. If LPRINT   *
!   is  .TRUE. , the contents of the sorted file are interpreted and   *
!   printed to unit 99.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: orbital_index_radix
      use array_memory_management, only: alloc, dalloc
      use orbital_csf_basis,           only: orbital_principal_n, csf_count, orbital_symmetry_labels
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use angular_coefficient_channels, only: channel_coefficients, channel_packed_orbital_labels, &
         & channel_csf_columns, channel_matrix_positions
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mcp_channel_index
      integer, intent(in) :: coefficient_count
      integer, intent(out)  :: distinct_integral_count
      integer, intent(in)   :: symmetry_block_index,sort_capacity_per_channel
      logical, intent(in)   :: print_details
      character, intent(in) :: output_file_prefix*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER, DIMENSION(:), pointer :: ICLMN, INDEX, LABEL, NSWAP
      integer, pointer :: nswap(:)
!      REAL(DOUBLE), DIMENSION(:), pointer :: COEFF
      integer :: text_length, lck, io_status, mb, i, l, ir, icl, ind, lab, nsw, j,&
         last, ibeg, iend, ncontr, ia, ib, k, id, ic, igg
      real(real64) :: cof
      character :: number_text*20, srtlab*8, mcp_channel_label*3, ck*2
      character (len = len (output_file_prefix) + 3):: file_path
!-----------------------------------------------
!
      if(mcp_channel_index == 31) then
         igg = 1
      else
         igg = mcp_channel_index - 30
      end if
!      IGG = IGGMAX_K*(IGG-1)
!   Message
!
      call convrt (coefficient_count, number_text, text_length)
      if (mcp_channel_index > 31) then
         call convrt (mcp_channel_index - 32, ck, lck)
         write (6, *) "Sorting "//number_text(1:text_length)//" V(k="//ck(1:lck)//&
            ") coefficients ...", mcp_channel_index
      else
         write (6, *) "Sorting "//number_text(1:text_length)//" T coefficients ...", mcp_channel_index
      end if
!
      call convrt (mcp_channel_index, ck, lck)
      if (lck > 2) then
         write (stderr_unit, *) "sort: nfile > 99; check fullname"
         stop
      end if

      file_path = output_file_prefix//"."//ck(1:2)

      open(29,file=file_path,status="OLD",form="UNFORMATTED", &
         iostat=io_status, position="APPEND")
      if (io_status /= 0) then
         write (stderr_unit, *) " Error when opening the file ", file_path
         stop
      end if
!
!   Sort the list
!
      if (coefficient_count > 0) then
!
!   Allocate storage for all required arrays
!
         if (mcp_channel_index == 33) call alloc (nswap, coefficient_count,"NSWAP", "SORTMEM")
      end if
!
!   Sort LABEL into ascending order using the heapsort algorithm;
!   move the associated members of COEFF and INDEX in the same
!   manner; the code below is adapted from Press et al.
!
      if (mcp_channel_index==33 .and. coefficient_count>1) then

         l = coefficient_count/2 + 1
         ir = coefficient_count
  234    continue
         if (l > 1) then
            l = l - 1
            cof = channel_coefficients(l,igg)
            icl = channel_csf_columns(l,igg)
            ind = channel_matrix_positions(l,igg)
            lab = channel_packed_orbital_labels(l,igg)
            nsw = nswap(l)
         else
            cof = channel_coefficients(ir,igg)
            icl = channel_csf_columns(ir,igg)
            ind = channel_matrix_positions(ir,igg)
            lab = channel_packed_orbital_labels(ir,igg)
            nsw = nswap(ir)
            channel_coefficients(ir,igg) = channel_coefficients(1,igg)
            channel_csf_columns(ir,igg) = channel_csf_columns(1,igg)
            channel_matrix_positions(ir,igg) = channel_matrix_positions(1,igg)
            channel_packed_orbital_labels(ir,igg) = channel_packed_orbital_labels(1,igg)
            nswap(ir) = nswap(1)
            ir = ir - 1
            if (ir == 1) then
               channel_coefficients(1,igg) = cof
               channel_csf_columns(1,igg) = icl
               channel_matrix_positions(1,igg) = ind
               channel_packed_orbital_labels(1,igg) = lab
               nswap(1) = nsw
               goto 456
            end if
         end if
         i = l
         j = l + l
  345    continue
         if (j <= ir) then
            if (j < ir) then
               if (channel_packed_orbital_labels(j,igg) < channel_packed_orbital_labels(j+1,igg)) j = j + 1
            end if
            if (lab < channel_packed_orbital_labels(j,igg)) then
               channel_coefficients(i,igg) = channel_coefficients(j,igg)
               channel_csf_columns(i,igg) = channel_csf_columns(j,igg)
               channel_matrix_positions(i,igg) = channel_matrix_positions(j,igg)
               channel_packed_orbital_labels(i,igg) = channel_packed_orbital_labels(j,igg)
               nswap(i) = nswap(j)
               i = j
               j = j + j
            else
               j = ir + 1
            end if
            goto 345
         end if
         channel_coefficients(i,igg) = cof
         channel_csf_columns(i,igg) = icl
         channel_matrix_positions(i,igg) = ind
         channel_packed_orbital_labels(i,igg) = lab
         nswap(i) = nsw
         goto 234

      else if (mcp_channel_index/=33 .and. coefficient_count>1) then
!
!   Sort LABEL into ascending order using the heapsort algorithm;
!   move the associated members of COEFF and INDEX in the same
!   manner; the code below is adapted from Press et al.
!
         l = coefficient_count/2 + 1
         ir = coefficient_count
   92    continue
         if (l > 1) then
            l = l - 1
            cof = channel_coefficients(l,igg)
            icl = channel_csf_columns(l,igg)
            ind = channel_matrix_positions(l,igg)
            lab = channel_packed_orbital_labels(l,igg)
         else
            cof = channel_coefficients(ir,igg)
            icl = channel_csf_columns(ir,igg)
            ind = channel_matrix_positions(ir,igg)
            lab = channel_packed_orbital_labels(ir,igg)
            channel_coefficients(ir,igg) = channel_coefficients(1,igg)
            channel_csf_columns(ir,igg) = channel_csf_columns(1,igg)
            channel_matrix_positions(ir,igg) = channel_matrix_positions(1,igg)
            channel_packed_orbital_labels(ir,igg) = channel_packed_orbital_labels(1,igg)
            ir = ir - 1
            if (ir == 1) then
               channel_coefficients(1,igg) = cof
               channel_csf_columns(1,igg) = icl
               channel_matrix_positions(1,igg) = ind
               channel_packed_orbital_labels(1,igg) = lab
               goto 456
            end if
         end if
         i = l
         j = l + l
   93    continue
         if (j <= ir) then
            if (j < ir) then
               if (channel_packed_orbital_labels(j,igg) < channel_packed_orbital_labels(j+1,igg)) j = j + 1
            end if
            if (lab < channel_packed_orbital_labels(j,igg)) then
               channel_coefficients(i,igg) = channel_coefficients(j,igg)
               channel_csf_columns(i,igg) = channel_csf_columns(j,igg)
               channel_matrix_positions(i,igg) = channel_matrix_positions(j,igg)
               channel_packed_orbital_labels(i,igg) = channel_packed_orbital_labels(j,igg)
               i = j
               j = j + j
            else
               j = ir + 1
            end if
            goto 93
         end if
         channel_coefficients(i,igg) = cof
         channel_csf_columns(i,igg) = icl
         channel_matrix_positions(i,igg) = ind
         channel_packed_orbital_labels(i,igg) = lab
         goto 92

      end if
!
!   Sorting complete; rewrite the file header
!
!
  456 continue
      write (29) "MCP", symmetry_block_index, csf_count, coefficient_count
!GG      WRITE (9999,*) 'MCP', NB, NCF, NCOEFF, NFILE
!
!   Write the sorted list to mcp.xx
!
      if (coefficient_count > 0) then
!
         last = channel_packed_orbital_labels(1,igg)
         ibeg = 1
         iend = 1
         distinct_integral_count = 1
!
         do i = 2, coefficient_count
            if (channel_packed_orbital_labels(i,igg) == last) then
               iend = iend + 1
            else
               write (29) last, iend - ibeg + 1
!GG               WRITE (9999,*) LAST, IEND - IBEG + 1
               write (29) (channel_csf_columns(j,igg),channel_matrix_positions(j,igg),channel_coefficients(j,igg),     &
                                                          j=ibeg,iend)
!cjb-GG              DO J = IBEG,IEND
!GG               WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
!cjb-GG              END DO

!              IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,IEND)
               distinct_integral_count = distinct_integral_count + 1
               last = channel_packed_orbital_labels(i,igg)
               ibeg = iend + 1
               iend = ibeg
            end if
         end do
!
         if (ibeg <= coefficient_count) then
            write (29) last, coefficient_count - ibeg + 1
!GG            WRITE (9999,*) LAST, NCOEFF - IBEG + 1
            write (29) (channel_csf_columns(j,igg),channel_matrix_positions(j,igg),channel_coefficients(j,igg),       &
                                                       j=ibeg,coefficient_count)
!GG          DO J=IBEG,NCOEFF
!GG            WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
!GG          END DO

!           IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,NCOEFF)
         end if
!
      else
!
         distinct_integral_count = 0
!
      end if
!
!     write the terminator record for this block
!
      write (29) 0, 0
!GG      WRITE (9999,*) 0, 0
      close(29)
!
!   Completion message
!
      !PRINT *, ' ... sort complete; ', ntgral, ' integrals;'
!
!   Debug printout
!
      if (print_details) then
         write (99, '(/,"From SORT:")')
         write (6, '(/,"From SORT:")')
         if (coefficient_count > 0) then
!
            last = channel_packed_orbital_labels(1,igg)
            ibeg = 1
            iend = 1
!
            do i = 2, coefficient_count
               if (channel_packed_orbital_labels(i,igg) == last) iend = iend + 1
  567          continue
               if (channel_packed_orbital_labels(i,igg)==last .and. i/=coefficient_count) cycle
               lab = last
               ncontr = iend - ibeg + 1
               if (mcp_channel_index == 31) then
                  ia = mod(lab,orbital_label_radix)
                  ib = lab/orbital_label_radix
                  write (99, '(" I(",1i2,1a2,",",1i2,1a2,"):")') &
                     orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                     & orbital_symmetry_labels(ib)
                  do j = ibeg, iend
                     write (99, '("  T_[",1i2,",",1i4,"] = ",1p,d19.12)') channel_csf_columns(j,igg), &
                        & channel_matrix_positions(j,igg),       &
                                                        channel_coefficients(j,igg)
                  end do
               else
                  k = mcp_channel_index - 32
                  id = mod(lab,orbital_label_radix)
                  lab = lab/orbital_label_radix
                  ib = mod(lab,orbital_label_radix)
                  lab = lab/orbital_label_radix
                  ic = mod(lab,orbital_label_radix)
                  ia = lab/orbital_label_radix
                  write (99, '(" R^[(",1i2,")] (",1i2,1a2,",",1i2,1a2,";",1i2,1a2,",",1i2,1a2,"):")') &
                     k, orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                     & orbital_symmetry_labels(ib), orbital_principal_n(ic), orbital_symmetry_labels&
                     (ic), orbital_principal_n(id), orbital_symmetry_labels(id)
                  do j = ibeg, iend
                     write (99, '("  V^[(",1i2,")]_[",1i8,",",1i8,"] = ",1p,d19.12)') &
                        k,channel_csf_columns(j,igg),channel_matrix_positions(j,igg),     &
                                                          channel_coefficients(j,igg)
                  end do
               end if
               last = channel_packed_orbital_labels(i,igg)
               ibeg = iend + 1
               iend = ibeg
               if (iend == coefficient_count) goto 567
            end do
         end if
         write (99, '("  Number of integrals is ",1i4)') distinct_integral_count
      end if
!
!   Deallocate storage
!
      if (coefficient_count > 0) then
         if (mcp_channel_index == 33) call dalloc (nswap, "NSWAP", "SORTMEM")
      end if


      return
      end subroutine sortmem
