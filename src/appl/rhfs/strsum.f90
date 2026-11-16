!***********************************************************************
!                                                                      *
      subroutine strsum
!                                                                      *
!   Generates the first part of  hfs92.sum  (on stream 24 and 29).     *
!                                                                      *
!   Call(s) to: ENGOUTH                                                *
!                                                                      *
!   Written by P. Jonsson                 Last revision: 20 Oct 1999   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:03   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/01/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use eigv_c
      use nsmdat_c,        only: sqn, dmomnm, qmomb,                  &
                                 hfsi=>sqn, hfsd=>dmomnm, hfsq=>qmomb
      use prnt_c
      use syma_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use engouth_i
      implicit none (type, external)
!-----------------------------------------------
!
!    Write nuclear data
!
      write (24, 302) hfsi
      write (24, 303) hfsd
      write (24, 304) hfsq
      write (29, 302) hfsi
      write (29, 303) hfsd
      write (29, 304) hfsq
!
!    Write the list of eigenpair indices
!
!      CALL ENGOUTH (EAV, EVAL, IATJPO, IASPAR, IVEC, NVEC, 3)
!
      return
!
  302 format('Nuclear spin                        ',1p,d22.15,' au')
  303 format('Nuclear magnetic dipole moment      ',1p,d22.15,' n.m.')
  304 format('Nuclear electric quadrupole moment  ',1p,d22.15,' barns')
      return
!
      end subroutine strsum
