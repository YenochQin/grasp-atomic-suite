!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function dsubrs (eol, i, j, jblock)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
      use eigv_c
      use hblock_c
      use pos_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i
      integer , intent(in) :: j
      integer  :: jblock
      logical , intent(in) :: eol
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, kcmin, ik, jk
!-----------------------------------------------

      if (eol) then
         dsubrs = 0.d0
         do k = 1, nevblk(jblock)
            kcmin = k + ncminpast(jblock)
            ik = i + nevecpast(jblock) + (k - 1)*ncfblk(jblock)
            jk = j + nevecpast(jblock) + (k - 1)*ncfblk(jblock)
            dsubrs = dsubrs + evec(ik)*evec(jk)*wt(kcmin)
         end do
      else if (i == j) then
         dsubrs = wt(i)
      else
         dsubrs = 0.d0
      endif

      return
      end function dsubrs
