!***********************************************************************
!                                                                      *
      subroutine smsread(vint,vint2)
!
!   Call(s) to: [LIB92]: ALCBUF                                        *
!                                                                      *
!   WRITTEN  by C. Naz\'e  Oct. 2011                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,    only: orbital_index_radix, max_orbitals
      use prnt_c
      use ris_c
      use orb_c
      use eigv_c
      use buffer_c
      use debug_c,          only: cutoff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcbuf_i
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: contri, contrik1, coeffsms
      integer :: ic, ir, i, j, iia, iib, iic, iid, loc, lab, ios
!-----------------------------------------------
!
      call alcbuf (1)

      rewind(51)
  16  read (51,iostat = ios) ic,ir
!C      WRITE(*,*) IC,IR
      if (ios .eq. 0) then
!cc        DO J = 1,NVEC
          read(51) coeffsms,lab
          iid = mod (lab, key)
          lab = lab/key
          iib  = mod (lab, key)
          lab = lab/key
          iic = mod (lab, key)
          iia = lab/key
        do j = 1,nvec
          loc = (j-1)*csf_count
            contrik1 = - csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)                     &
               * coeffsms                                              &
               * vint (iia,iic)*vint(iib,iid)
            contri = - csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)                       &
               * coeffsms                                              &
               * ( vint2(iia,iic)*vint(iib,iid)                        &
               + vint2(iib,iid)*vint(iia,iic))/2.0_real64
          if (ir.ne.ic) then
            contri = 2.0_real64 * contri
            contrik1 = 2.0_real64 * contrik1
          endif
          smsc1(j) = smsc1(j) + contrik1
          smsc2(j) = smsc2(j) + contri
        enddo

      goto 16
      endif
      call alcbuf (3)
      return
      end subroutine smsread
