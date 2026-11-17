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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use hmat_c, only: double, long, stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use orb_c,           only: csf_count, orbital_count, csf_subshell_occupations
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: imcdf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ic, nelc, ir, lenthr, lenthc
      real(real64) :: elsto
      character    :: cir*8, cic*8
!-----------------------------------------------
!
      do ic = 1, csf_count
         read (imcdf) nelc, elsto, (hamiltonian_values(ir),ir=1,nelc), (matrix_partner_indices(ir),ir=1,nelc)
         do ir = 1, nelc
            call convrt (matrix_partner_indices(ir), cir, lenthr)
            call convrt (ic, cic, lenthc)
            write (99, 300) cir(1:lenthr), cic(1:lenthc), hamiltonian_values(ir)
         end do
      end do
!
  300 format(" H (",a,",",a,") = ",1p,1d19.12)
      return
!
      end subroutine hmout
