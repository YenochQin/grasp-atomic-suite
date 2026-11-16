

!***********************************************************************
!                                                                      *
      subroutine getinf
!                                                                      *
!   Interactively determines data useful for generating estimates of   *
!   the subshell radial wavefunctions.                                 *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
!               [RCI92]: SETISO.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use iounit_c
      use def_c
      use default_c
      use grid_c, only: rnt, h, hp, n
      use npar_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setiso_i
      use setqic_i
      use radgrd_i
      use nucpot_i
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: yes
!-----------------------------------------------
!
!
!
!   Open, check, load data from, and close the  .iso  file
!
      call setiso ('isodata')
!
!   Set defaults
!
      c = cvac
!
      if (nparm == 0) then
         rnt = exp((-65.0d00/16.0d00))/z
         h = 0.5d00**4
         n = min(220,nnnp)
      else
!CFF     .. should be Z-dependent
         rnt = 2.0d-06/z
         h = 5.0d-02
         n = nnnp
      endif
      hp = 0.0d00
!
      if (ndef /= 0) then
         write (istde, *) 'Change the default speed of light ', &
            'or radial grid parameters?'
!
         yes = getyn()
         if (yes) then
!
!   Modify the speed of light
!
            write (istde, *) 'The physical speed of light in ', &
               'atomic units is', cvac, ';'
            write (istde, *) 'revise this value?'
            yes = getyn()
            if (yes) then
               write (istde, *) 'Enter the revised value:'
               read (5, *) c
            endif
!
!   Modify the parameters controlling the radial grid
!
            write (istde, *) 'The default radial grid parameters ', &
               'for this case are:'
            write (istde, *) ' RNT = ', rnt, ';'
            write (istde, *) ' H = ', h, ';'
            write (istde, *) ' HP = ', hp, ';'
            write (istde, *) ' N = ', n, ';'
            write (istde, *) 'revise these values?'
            yes = getyn()
            if (yes) then
               write (istde, *) 'Enter RNT:'
               read (5, *) rnt
               write (istde, *) 'Enter H:'
               read (5, *) h
               write (istde, *) 'Enter HP:'
               read (5, *) hp
               write (istde, *) 'Enter N:'
               read (5, *) n
            endif
!
         endif
      endif
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = h**6
!
!   Set up the coefficients for the numerical procedures
!
      call setqic
!
!   Generate the radial grid and all associated arrays
!
      call radgrd
!
!   Generate $- r \times V_ (r)$
!
      call nucpot
!
      return
      end subroutine getinf
