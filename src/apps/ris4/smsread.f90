!***********************************************************************
!                                                                      *
      subroutine smsread(vint,vint2)
      use calculation_file_units, only: two_body_cache_unit
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
      use grasp_parameters,    only: orbital_index_radix, max_orbitals
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use ris_c, only: dens1, dens2, dens3, dens4, dens5, dens6, dens7, smsc1, smsc2
      use orbital_csf_basis, only: csf_count, max_orbitals
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use angular_coefficient_buffer_state, only: angular_buffer_capacity, angular_coefficient_count, &
          angular_coefficient_labels, angular_coefficients
      use debug_output_settings,          only: cutoff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use angular_coefficient_buffer, only: alcbuf
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: contri, contrik1, coeffsms
      integer :: ic, ir, i, j, iia, iib, iic, iid, loc, lab, io_status
!-----------------------------------------------
!
      call alcbuf (1)

      rewind(two_body_cache_unit)
  16  read (two_body_cache_unit,iostat = io_status) ic,ir
!C      WRITE(*,*) IC,IR
      if (io_status == 0) then
!cc        DO J = 1,NVEC
          read(two_body_cache_unit) coeffsms,lab
          iid = mod (lab, orbital_label_radix)
          lab = lab/orbital_label_radix
          iib  = mod (lab, orbital_label_radix)
          lab = lab/orbital_label_radix
          iic = mod (lab, orbital_label_radix)
          iia = lab/orbital_label_radix
        do j = 1,nvec
          loc = (j-1)*csf_count
            contrik1 = - csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)                     &
               * coeffsms                                              &
               * vint (iia,iic)*vint(iib,iid)
            contri = - csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)                       &
               * coeffsms                                              &
               * ( vint2(iia,iic)*vint(iib,iid)                        &
               + vint2(iib,iid)*vint(iia,iic))/2.0_real64
          if (ir/=ic) then
            contri = 2.0_real64 * contri
            contrik1 = 2.0_real64 * contrik1
          end if
          smsc1(j) = smsc1(j) + contrik1
          smsc2(j) = smsc2(j) + contri
        end do

      goto 16
      end if
      call alcbuf (3)
      return
      end subroutine smsread
