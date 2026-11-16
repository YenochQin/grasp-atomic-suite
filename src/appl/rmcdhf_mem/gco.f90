!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function gco (k, ir, ia, ib)
!                                                                      *
!   This routine evaluates a coefficient                               *
!                                                                      *
!                                K                                     *
!                               g   (IA,IB)                            *
!                                IR                                    *
!                                                                      *
!                                                                      *
!   Here  K  is the multipolarity, IR  is the sequence number of the   *
!   configuration, and  IA and IB are orbital sequence  numbers. See   *
!   I P Grant,  B J McKenzie,  P H Norrington,  D F  Mayers, and N C   *
!   Pyper, Computer Phys Commun 21 (1980) 207-231, Eq (7).             *
!                                                                      *
!   Call(s) to: [LIB92]: CLRX.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use debug_c
      use orb_c, iiqa=>csf_subshell_occupations
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use clrx_i
      use iq_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: k
      integer  :: ir
      integer , intent(in) :: ia
      integer , intent(in) :: ib
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqa, iqb
      real(double) :: fac
!      LOGICAL :: FULLA, FULLB
!-----------------------------------------------
!
         iqa = iq (ia,ir)
         iqb = iq (ib,ir)
!GG      IQA = IBITS(IIQA((IA - 1)/4 + 1,IR),8*MOD(IA - 1,4),8)
!GG      IQB = IBITS(IIQA((IB - 1)/4 + 1,IR),8*MOD(IB - 1,4),8)
      if (iqa==orbital_two_j(ia) + 1 .or. iqb==orbital_two_j(ib)+1) then
         fac = clrx(orbital_kappa(ia),k,orbital_kappa(ib))
         gco = -dble(iqa*iqb)*fac*fac
      else
         gco = 0.0d00
      endif

      return
      end function gco
