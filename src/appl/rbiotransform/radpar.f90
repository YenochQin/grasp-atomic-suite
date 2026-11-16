!***********************************************************************
!                                                                      *
      subroutine radpar
!                                                                      *
!   This subroutine sets the parameters controlling the radial grid    *
!                                                                      *
!   Last revision: June 1996                                           *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnp
      use def_c,           only: c, cvac, z, accy
      use default_c
      use grid_c
      use npar_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical   :: yes
      character :: answer
!-----------------------------------------------
!
      if (nparm == 0) then
         rnt = exp((-65.0d00/16.0d00))/z
         h = 0.5d00**4
         n = min(220,nnnp)
      else
!cff     .. Should be Z-dependent
         rnt = 2.0d-06/z
         h = 5.0d-02
         n = nnnp
      endif
      hp = 0.0d00
      if (ndef /= 0) then
         write (6, *) 'The default radial grid parameters'
         write (6, *) ' for this case are:'
         write (6, *) ' RNT = ', rnt, ';'
         write (6, *) ' H = ', h, ';'
         write (6, *) ' HP = ', hp, ';'
         write (6, *) ' N = ', n, ';'
         write (6, *) ' revise these values?'
         yes = getyn()
         if (yes) then
            write (6, *) 'Enter RNT:'
            read (5, *) rnt
            write (6, *) 'Enter H:'
            read (5, *) h
            write (6, *) 'Enter HP:'
            read (5, *) hp
            write (6, *) 'Enter N:'
            read (5, *) n
         endif
      endif
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = h**6
!

      if (ndef /= 0) then
         write (6, *) 'The physical speed of light in'
         write (6, *) ' atomic units is', cvac, ';'
         write (6, *) ' revise this value?'
         yes = getyn()
         if (yes) then
            write (6, *) 'Enter the revised value:'
            read (5, *) c
         else
            c = cvac
         endif
      else
         c = cvac
      endif

      return
      end subroutine radpar
