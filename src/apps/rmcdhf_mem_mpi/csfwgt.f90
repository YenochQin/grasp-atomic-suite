!***********************************************************************
      subroutine csfwgt(write_standard_output)
!                                                                      *
!   Print  the  weights of the largest five CSF contributors to each   *
!   ASF.                                                               *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, ISPAR, ITJPO.          *
!                                                                      *
!                                          Last updated: 21 Dec 1992   *
!                                          Last updated: 24 Feb 1997   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: iccmin, ncmin
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use angular_symmetry_labels, only: right_justified_j_labels => jlbr, parity_labels => jlbp
      use level_symmetry_metadata,          only: jpgg
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
          csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
      use transition_block_indices, only: idxblk
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) :: write_standard_output
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iconf(5)
      integer :: iv, symmetry_block_index, block_csf_count, csf_offset, eigenvalue_offset, block_mixing_offset, &
         & level_mixing_offset, reported_coefficient_count, &
         icf, ivtjpo, ivspar, i, j, itemp, irem, ip
      real(real64) :: temp, w
      real(real64) :: coeff(5)
      character :: record*256, number_text*8
!-----------------------------------------------
!

      if (write_standard_output) then
         write (stdout_unit, '(/,"Weights of major contributors to ASF:"/,/,&
            &"Block Level J Parity      CSF contributions"/)')
      else
         write (24, '(/,"Weights of major contributors to ASF:"/,/,"Block Level J Parity      CSF contributions"/)')
      end if

      do iv = 1, ncmin
!         loop over eigenvectors

         symmetry_block_index = idxblk(iv)
         block_csf_count = csf_counts_by_block(symmetry_block_index)
         csf_offset = csf_offsets_by_block(symmetry_block_index)
         eigenvalue_offset = eigenvalue_offsets_by_block(symmetry_block_index)
         block_mixing_offset = mixing_offsets_by_block(symmetry_block_index)
         level_mixing_offset = block_mixing_offset + (iv - eigenvalue_offset - 1)*block_csf_count
         reported_coefficient_count = min(5,block_csf_count)                       !Find maximum 5 ... within block
         icf = iccmin(iv)
!GGGG
         ivtjpo = iabs(jpgg(symmetry_block_index))
         if(jpgg(symmetry_block_index) >= 0) then
            ivspar = 2
         else
            ivspar = 1
         end if
!GG         IVTJPO = IATJPO(IV)                     ! j-value related
!GG         IVSPAR = IASPAR(IV)                     ! parity related


         do i = 1, reported_coefficient_count
            coeff(i) = csf_mixing_coefficients(level_mixing_offset + i)
            iconf(i) = i
         end do

!         sort the first nelt in decreasing order
         do i = 1, reported_coefficient_count
            do j = i + 1, reported_coefficient_count
               if (abs(coeff(j)) <= abs(coeff(i))) cycle
               temp = coeff(i)
               coeff(i) = coeff(j)
               coeff(j) = temp
               itemp = iconf(i)
               iconf(i) = iconf(j)
               iconf(j) = itemp
            end do
         end do

         l20: do i = reported_coefficient_count + 1, block_csf_count
            w = csf_mixing_coefficients(level_mixing_offset + i)
            if (w==0.0_real64 .or. abs(w)<=abs(coeff(reported_coefficient_count))) cycle  l20
!               we have a non-zero value larger than the largest so far
            do j = 1, reported_coefficient_count
               if (abs(w) <= abs(coeff(j))) cycle
               coeff(reported_coefficient_count:1+j:(-1)) = coeff(reported_coefficient_count-1:j:(-1))
               iconf(reported_coefficient_count:1+j:(-1)) = iconf(reported_coefficient_count-1:j:(-1))
               coeff(j) = w
               iconf(j) = i
               cycle  l20
            end do
         end do l20

!GG         IP = (IASPAR(IV) + 3)/2
           ip =  ivspar

         if (write_standard_output) then
            write (stdout_unit, '(i3,1x,i5,2x,2a4,5(3x,f8.4))') symmetry_block_index, icf, right_justified_j_labels(ivtjpo), &
               & parity_labels(ip), (coeff(i),i&
               =1,reported_coefficient_count)
            write (stdout_unit, '(19x,5(3x,i8))') (iconf(i),i=1,reported_coefficient_count)
         else

            write (24, '(i3,1x,i5,2x,2a4,5(3x,f8.4))') symmetry_block_index, icf, right_justified_j_labels(ivtjpo), &
               & parity_labels(ip), (coeff(i),i=1,&
               reported_coefficient_count)
            write (24, '(19x,5(3x,i8))') (iconf(i),i=1,reported_coefficient_count)
         end if
      end do


      return
      end subroutine csfwgt
