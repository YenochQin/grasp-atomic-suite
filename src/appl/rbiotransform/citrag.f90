!***********************************************************************
!                                                                      *
!                                                                      *
      subroutine citrag(ciin, ncsf, nciv, l, nshl, t, nin, ntestg, ciout, scr)
!                                                                      *
!   Calculate the action of the  operator                              *
!                                                                      *
!    T(L) = T(nl) T(n-1 L ) .... T(1 l)                                *
!                                                                      *
!   where                                                              *
!                                                                      *
!    T(i l ) = ( sum(k=0,2l+1) t(i l ) ** k ) (t_ilil ** E(il,il))     *
!    t(i l ) = sum(j.lt.i) t(j,i) E(lj,li)                             *
!    E(lj,li) total symmetric shell excitation operator                *
!                                                                      *
!   Modification of the CITRA routine for MCHF                         *
!                                                                      *
!   =====                                                              *
!   Input                                                              *
!   =====                                                              *
!   CIIN      : Input CI vectors ( Destroyed in the process )          *
!   NCSCF     : Length of CI expansion                                 *
!   NCIV      : Number of CI vectors                                   *
!   L         : L value of shells in excitations                       *
!   NSHL      : Number of shells with this L                           *
!   T         : List of coefficients                                   *
!   NIN       : Number of inactive shells                              *
!                                                                      *
!   ======                                                             *
!   Output                                                             *
!   ======                                                             *
!                                                                      *
!   CIOUT : List of output CI vectors                                  *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use sbdat1_c
      use orb_c
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use wrtmat_i
      use scalve_i
      use tiinig_i
      use ti1tv_i
      use vecsum_i
      use copvec_i, only: copvec
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ncsf
      integer, intent(in) :: nciv
      integer, intent(in) :: l
      integer, intent(in) :: nshl
      integer, intent(in) :: nin
      integer, intent(in) :: ntestg
      real(real64)  :: ciin(ncsf,nciv)
      real(real64), intent(in) :: t(nshl,nshl)
      real(real64)  :: ciout(ncsf,nciv)
      real(real64)  :: scr(ncsf,nciv)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntestl, ntest, iin, ipot, i, n
      real(real64) :: factor, tii, xnfaci
!-----------------------------------------------
!

      ntestl = 0
      ntest = max(ntestg,ntestl)
      ntest = 0

!
      if (ntest >= 10) then
         write (6, *)
         write (6, *) ' ***************'
         write (6, *) ' Entering CITRAG'
         write (6, *) ' ***************'
         write (6, *)
      endif
      if (ntest >= 100) then
         write (6, *) ' Input CI vectors'
         call wrtmat (ciin, ncsf, nciv, ncsf, nciv)
         write (6, *) ' Transformation matrix T'
         call wrtmat (t, nshl, nshl, nshl, nshl)
      endif
!
!. Factor from inactive shells
!
      if (nin /= 0) then
         factor = 1.0_real64
         do iin = 1, nin
            factor = factor*t(iin,iin)
         end do
!
!       IPOT = 2*(2*L+1)  (number of m_lm_s. This should be replaced
!                          by (2j+1) corresponding to L)
!
         ipot = 2*iabs(orbital_kappa(nshlp(l,iin)))
         factor = factor**ipot
         call scalve (ciin, factor, nciv*ncsf)
      endif
      if (nin == nshl) call copvec (ciin, ciout, nciv*ncsf)
!
      do i = nin + 1, nshl
         if (ntest >= 100) write (6, *) ' Loop I,L = ', i, l
!
!. The diagonal contribution
!
         tii = t(i,i)
         call tiinig (ciin, ncsf, nciv, i, l, tii, ciout, ntestg)
!        IF (LWORK2.LT.NCIV*NCSF) THEN
!          WRITE(*,*) 'In CITRAG: Dimension of LWORK2 must be',
!     &               'increased to at least',NCIV*NCSF
!        ENDIF
         call copvec (ciout, scr, nciv*ncsf)
!
!.  Off diagonal contributions
!
         xnfaci = 1.0_real64
         do n = 1, 2*iabs(orbital_kappa(nshlp(l,i)))
            if (ntest >= 100) write (6, *) ' Loop N = ', n
!
!   T ** (N-1) is supposed to be in SCR, copy to CIIN
!   and apply S
!
            call copvec (scr, ciin, nciv*ncsf)
            call ti1tv (ciin, ncsf, nciv, i, l, t(1,i), nshl, scr, ntestg)
            xnfaci = xnfaci/float(n)
            call vecsum (ciout, ciout, scr, 1.0_real64, xnfaci, nciv*ncsf)
         end do
         call copvec (ciout, ciin, nciv*ncsf)
!
      end do
!
      if (ntest >= 100) then
         write (6, *) ' Output CI vectors L = ', l
         call wrtmat (ciout, ncsf, nciv, ncsf, nciv)
      endif
!
      return
      end subroutine citrag
