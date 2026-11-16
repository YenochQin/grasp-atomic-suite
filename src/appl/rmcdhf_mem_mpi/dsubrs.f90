!***********************************************************************
!                                                                      *
      real(real64) function dsubrs (eol, i, j, jblock)
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
      use def_c
      use eigv_c
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use pos_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: i
      integer , intent(in) :: j
      integer  :: jblock
      logical , intent(in) :: eol
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, kcmin, ik, jk
!-----------------------------------------------

      if (eol) then
         dsubrs = 0.0_real64
         do k = 1, eigenvector_counts_by_block(jblock)
            kcmin = k + eigenvalue_offsets_by_block(jblock)
            ik = i + mixing_offsets_by_block(jblock) + (k - 1)*csf_counts_by_block(jblock)
            jk = j + mixing_offsets_by_block(jblock) + (k - 1)*csf_counts_by_block(jblock)
            dsubrs = dsubrs + csf_mixing_coefficients(ik)*csf_mixing_coefficients(jk)*wt(kcmin)
         end do
      else if (i == j) then
         dsubrs = wt(i)
      else
         dsubrs = 0.0_real64
      end if

      return
      end function dsubrs
