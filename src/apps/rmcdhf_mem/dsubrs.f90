!***********************************************************************
!                                                                      *
      real(real64) function dsubrs (extended_optimal_level, bra_csf_index, ket_csf_index, symmetry_block_index)
!                                                                      *
!   The coefficients d   for  I = r, J = s  are calculated here.       *
!                     rs                                               *
!                                                                      *
!                         NCMIN                                        *
!                          Sum     w     c          c                  *
!                         t = 1     t     r Gamma    s Gamma           *
!                                                t          t          *
!    EOL:          d   = -------------------------------------         *
!                   rs             NCMIN                               *
!                                   Sum     w                          *
!                                  t = 1     t                         *
!                                                                      *
!                                                                      *
!                                           /  NCF                     *
!    EAL:             d   =  delta     w   /   Sum   w                 *
!                      rs         rs    r /   t = 1   t                *
!                                                                      *
!                                                                      *
!   Written by Farid A Parpia               Last update: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!   This routine assumes (and does not check) that I and J belong
!   the same block and the job is done within the block specified
!   by parameter jblock
!
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: wt
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
          csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: bra_csf_index
      integer , intent(in) :: ket_csf_index
      integer  :: symmetry_block_index
      logical , intent(in) :: extended_optimal_level
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: level_index, kcmin, ik, jk
!-----------------------------------------------

      if (extended_optimal_level) then
         dsubrs = 0.0_real64
         do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)
            kcmin = level_index + eigenvalue_offsets_by_block(symmetry_block_index)
            ik = bra_csf_index + mixing_offsets_by_block(symmetry_block_index) + (level_index - 1)*csf_counts_by_block( &
               & symmetry_block_index)
            jk = ket_csf_index + mixing_offsets_by_block(symmetry_block_index) + (level_index - 1)*csf_counts_by_block( &
               & symmetry_block_index)
            dsubrs = dsubrs + csf_mixing_coefficients(ik)*csf_mixing_coefficients(jk)*wt(kcmin)
         end do
      else if (bra_csf_index == ket_csf_index) then
         dsubrs = wt(bra_csf_index)
      else
         dsubrs = 0.0_real64
      end if

      return
      end function dsubrs
