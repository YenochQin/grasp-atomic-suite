!***********************************************************************
!                                                                      *
      real(real64) function fco (k, ir, ia, ib)
!                                                                      *
!   This routine evaluates a coefficient                               *
!                                                                      *
!                                K                                     *
!                               f   (IA,IB)                            *
!                                IR                                    *
!                                                                      *
!   Here  K  is the multipolarity, IR  is the sequence number of the   *
!   configuration, and  IA  and  IB  are orbital  sequence  numbers.   *
!   ( I P Grant,  B J McKenzie,  P H Norrington, D F Mayers, and N C   *
!   Pyper,  Computer Phys Commun 21 (1980) 207-231, Eqs (6). )         *
!                                                                      *
!   Call(s) to: [LIB92]: CLRX, IQ.                                     *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 21 Dec 1992   *
!                                                                      *
!***********************************************************************
! . Two calls to IQ are physically inlined here to reduce the overhead.
! . In the inlined text, the checks for the integer range (ref.
!    lib92/iq.f for details) are removed since FCO is called only by
!    SETCOF and SETHAM where the input arguments always fall in the
!    suitable range.
!XHH 1997.03.05
!
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_c
      use orb_c, iiqa=>csf_subshell_occupations
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use clrx_i, only: clrx
      use iq_i, only: iq
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: k
      integer  :: ir
      integer , intent(in) :: ia
      integer , intent(in) :: ib
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqa, iqf, kappa
      real(real64) :: fac
!-----------------------------------------------
!
      if (ia == ib) then
!
         iqa = iq (ia,ir)
!         IQA = IBITS(IIQA((IA - 1)/4 + 1,IR),8*MOD(IA - 1,4),8)
!
         if (k == 0) then
            fco = dble((iqa*(iqa - 1))/2)
         else
            iqf = orbital_two_j(ia) + 1
            if (iqa == iqf) then
               kappa = orbital_kappa(ia)
               fac = clrx(kappa,k,kappa)*dble(iqa)
               fco = -0.5_real64*fac*fac
            else
               fco = 0.0_real64
            end if
         end if
!
      else
!
         if (k == 0) then
            fco = dble (iq (ia,ir)*iq (ib,ir))
!            FCO = DBLE(IBITS(IIQA((IA - 1)/4 + 1,IR),8*MOD(IA - 1,4),8)*IBITS(&
!               IIQA((IB - 1)/4 + 1,IR),8*MOD(IB - 1,4),8))
         else
            fco = 0.0_real64
         end if
!
      end if
      return
!*
!  300 FORMAT (/'  ',1I2
!     :        /' f    (',1I2,1A2,',',1I2,1A2,') = ',1PD21.14,
!     :        /'  ',1I3/)
!*
      end function fco
