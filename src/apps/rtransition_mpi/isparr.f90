!***********************************************************************
!                                                                      *
      integer function isparr (csf_index)
      use, intrinsic :: iso_fortran_env, only: output_unit
!                                                                      *
!   ISPARR is the value of P for CSF number ICSF.                      *
!                                                                      *
!   Call(s) to: [LIB92]: IUNPCK.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grasp_parameters, only: max_orbitals
      use csf_coupling_data,        only: jcupar
      use standard_io_units,      only: stderr_unit
      use orbital_csf_basis,         only: ncfr
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: csf_index
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!
      if (csf_index>=1 .and. csf_index<=ncfr) then
         isparr = jcupar(max_orbitals,csf_index)

         if (isparr > 127) isparr = isparr - 256
         isparr = sign(1,isparr)
      else
         write (output_unit, *) "ISPARR: Argument ICSF is out of range."
         stop
      end if
!
      return
      end function isparr
