!***********************************************************************
!                                                                      *
      subroutine getsmd(name)
!                                                                      *
!   Interactively determines the data governing the SMS problem.       *
!                                                                      *
!   Call(s) to: [LIB92]: NUCPOT, RADGRD, SETQIC.                       *
!               [RCI92]: SETISO, SETRWF.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnp
      use decide_c
      use def_c
      use default_c
      use foparm_c
      use grid_c
      use npar_c
      use npot_c
      use orb_c
      use wave_c
      use wfac_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setiso_i
      use setqic_i
      use radgrd_i
      use nucpot_i
      use setrwfa_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character, intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: yes
!-----------------------------------------------
!
!
!   Open, check, load data from, and close the  .iso  file
!
      call setiso ('isodata')
!
!   Determine the physical effects specifications
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
!
      if (ndef /= 0) then
!
         write (6, *) 'Treat contributions of some CSFs'
         write (6, *) ' as first-order perturbations?'
         yes = getyn()
         if (yes) then
            lfordr = .true.
            write (6, *) 'The contribution of CSFs'
            write (6, *) ' 1 -- ICCUT will be treated'
            write (6, *) ' variationally; the remainder'
            write (6, *) ' perturbatively; enter ICCUT:'
            read (5, *) iccut
         else
            lfordr = .false.
            iccut = 0
         endif
      else
         lfordr = .false.
         iccut = 0
      endif
!
!   Determine the parameters controlling the radial grid
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
!   Set up the coefficients for the numerical procedures
!
      call setqic
!
!   Generate the radial grid and all associated arrays
!
      call radgrd
!
!   Generate $- r \times V_nuc (r)$
!
      call nucpot
!
!   Load the radial wavefunctions
!
!      CALL SETRWFA(NAME)
      call setrwfa (trim(name)//'.w')
!
      return
      end subroutine getsmd
