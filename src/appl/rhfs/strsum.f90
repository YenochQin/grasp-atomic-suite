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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use eigv_c
      use nsmdat_c,        only: sqn, dmomnm, qmomb,                  &
                                 hfsi=>sqn, hfsd=>dmomnm, hfsq=>qmomb
      use prnt_c, only: byte, double, ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use syma_c, only: double, iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use engouth_i, only: engouth
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
  302 format("Nuclear spin                        ",1p,d22.15," au")
  303 format("Nuclear magnetic dipole moment      ",1p,d22.15," n.m.")
  304 format("Nuclear electric quadrupole moment  ",1p,d22.15," barns")
      return
!
      end subroutine strsum
