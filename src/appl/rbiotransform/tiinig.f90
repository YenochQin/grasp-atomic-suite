!***********************************************************************
!                                                                      *
      subroutine tiinig(ciin, ncsf, nciv, i, l, const, ciout, ntestg)
!                                                                      *
!   Calculates the action of the operator                              *
!   Const ** E(li,li) on a set of vectors                              *
!                                                                      *
!   Adapted for GRASP, daughter of TIINI, born February 1996           *
!                                                                      *
!   =====                                                              *
!   Input                                                              *
!   =====                                                              *
!   CIIN      : Input CI vectors                                       *
!   NCSCF     : Length of CI expansion                                 *
!   NCIV      : Number of CI vectors                                   *
!   I         : Shell number                                           *
!   L         : symmetry                                               *
!   NSHLP(L,K): Gives the shell number as defined in getcsl for        *
!               the K:th shell with symmetry L                         *
!   CONST     : The constant                                           *
!                                                                      *
!   ======                                                             *
!   Output                                                             *
!   ======                                                             *
!                                                                      *
!   CIOUT : List of output CI vectors                                  *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:41:42   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix
      use mcpdata_c
      use sbdat1_c

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setvec_i
      use wrtmat_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ncsf
      integer, intent(in) :: nciv
      integer, intent(in) :: i
      integer, intent(in) :: l
      integer, intent(in) :: ntestg
      real(real64), intent(in) :: const
      real(real64), intent(in) :: ciin(ncsf,nciv)
      real(real64)  :: ciout(ncsf,nciv)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntestl, ntest, nfound, k, ia, ib, ifirst, ielmnt, ival, ileft&
         , ivec
      real(real64) :: constn
!-----------------------------------------------
!
!
      ntestl = 0
      ntest = max(ntestl,ntestg)


      if (ntest >= 10) write (6, *) ' Entering TIINI'

      call setvec (ciout, 0.0_real64, ncsf*nciv)
!
!.  Obtain address of first coupling coefficient for h(il,il) :  IFIRST
!.  Obtain number of        coupling coefficient for h(il,il) :  NFOUND
!.
!     NFOUND : Number of coefficients obtained
!     IVAL   : actual RACAH coefficient <CSF(L)!E(il,il)!CSF(R)>
!     ILEFT  = CSF(L) ?
!
      nfound = 0
      do k = 1, nintg
         ia = intgrl(k)/key
         ib = mod(intgrl(k),key)
         if (nshlp(l,i)/=ia .or. ia/=ib) cycle
         if (k == 1) then
            ifirst = 1
         else
            ifirst = intptr(k - 1) + 1
         endif
         nfound = intptr(k) - ifirst + 1
         exit
      end do
      do ielmnt = 1, nfound
! Bug 2011-08-18 Per Jonsson        IVAL = CNN(IFIRST-1+IELMNT)
         ival = idnint(cnn(ifirst-1+ielmnt))
         constn = const**ival
         ileft = jann(ifirst - 1 + ielmnt)
         ciout(ileft,:nciv) = constn*ciin(ileft,:nciv)
      end do
!
!.  The previous provided us with all
!   terms with nonvanishing occupation.
!   For terms with vanishing occupation of il,
!   just copy coefficients, since (x) ** 0 = 1
!
      where (ciout(:ncsf,:nciv) == 0.0_real64)
         ciout(:ncsf,:nciv) = ciin(:ncsf,:nciv)
      end where
!
      if (ntest >= 100) then
         write (6, *)
         write (6, *) ' Input and output vectors from TIINI I,L', i, l
         call wrtmat (ciin, ncsf, nciv, ncsf, nciv)
         write (6, *)
         call wrtmat (ciout, ncsf, nciv, ncsf, nciv)
         write (6, *)
      endif
!
      if (ntest >= 10) write (6, *) ' Leaving  TIINI'

      return
      end subroutine tiinig
