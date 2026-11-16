!***********************************************************************
!                                                                      *
      subroutine hmout(imcdf)
!                                                                      *
!   Routine for printing the Hamiltonian matrix.  File IMCDF must be   *
!   positioned correctly before a call is made to this module.         *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT.                                       *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use hmat_c
      use orb_c,           only: ncf, nw, iqa
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: imcdf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ic, nelc, ir, lenthr, lenthc
      real(double) :: elsto
      character    :: cir*8, cic*8
!-----------------------------------------------
!
      do ic = 1, ncf
         read (imcdf) nelc, elsto, (emt(ir),ir=1,nelc), (irow(ir),ir=1,nelc)
         do ir = 1, nelc
            call convrt (irow(ir), cir, lenthr)
            call convrt (ic, cic, lenthc)
            write (99, 300) cir(1:lenthr), cic(1:lenthc), emt(ir)
         end do
      end do
!
  300 format(' H (',a,',',a,') = ',1p,1d19.12)
      return
!
      end subroutine hmout
