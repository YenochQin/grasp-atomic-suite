!***********************************************************************
!                                                                      *
!JE   SUBROUTINE DENSREAD(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
      subroutine densread_seltz(dint1,dint2,dint3,                     &
                          dint4,dint5,dint6,                           &
                          dint7,dint1vec,dens1vec,nrnuc)
      use calculation_file_units, only: one_body_cache_unit
!                                                                      *
!   IF angular coefficients already exist                              *
!   This routine controls combines the radial and angular parts for the*
!   calculation of the NMS parameter, the electron density at the      *
!   origin and radial expectation values.                              *
!                                                                      *
!   Written by C\'edric Naz\'e                                         *
!                                                                      *
!                                         Last revision: Feb. 2011     *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: orbital_index_radix, max_orbitals, max_radial_points
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use ris_c, only: dens1, dens2, dens3, dens4, dens5, dens6, dens7, smsc1, smsc2
      use orbital_csf_basis, only: csf_count, max_orbitals
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
      integer :: ncount
      integer, intent(in) :: nrnuc
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals),  &
                                          dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                                             dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals),  &
                                          dint7(max_orbitals, max_orbitals)
      real(real64), intent(out)     :: dens1vec(nvec, nrnuc) !  JE ADD
      real(real64), intent(in) :: dint1vec(max_orbitals, max_orbitals, nrnuc) !  JE ADD
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: tshell_r(max_orbitals)
      real(real64) :: contri1vec(nrnuc), elemnt1vec(nrnuc)
      real(real64) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5
      real(real64) :: elemnt6, elemnt7
      real(real64) :: contri1, contri2, contri3, contri4, contri5
      real(real64) :: contri6, contri7
      integer :: icold, irold, io_status, ia, ib, ic, ir, i, j, l, loc, lab
!-----------------------------------------------
!
! DINT1 contain the density
! DINT2 contain the uncorrected NMS parameter: K^1_NMS
! DINT3 contain the expect. value <r>
! DINT4 contain the expect. value <r2>
! DINT5 contain the expect. value <r-1>
! DINT6 contain the expect. value <r-2>
! DINT7 contain the sum of NMS parameters: K^1_NMS+K^2_NMS+K^3_NMS
!
      dens1vec(:,:) = 0.0_real64                                     ! JE ADD

      icold = 0
      irold = 0

      rewind (one_body_cache_unit)
   16 read (one_body_cache_unit,iostat = io_status) ic,ir,ncount

      if (io_status == 0) then
!*      print*, 'ic', IC,IR,NCOUNT
!
!   Read successful; decode the labels of I(ab)
!
!      Initialise
        if ((ic/=icold).or.(ir/=irold)) then
!*            ISPARC = ISPAR (IC)
!*            ITJPOC = ITJPO (IC)
!*            ITJPOR = ITJPO (IR)
!*            IDIFF = ITJPOC - ITJPOR
            icold = ic
            irold = ir
!*        DO I = 1,NELMNT
          elemnt1 = 0.0_real64
          elemnt2 = 0.0_real64
          elemnt3 = 0.0_real64
          elemnt4 = 0.0_real64
          elemnt5 = 0.0_real64
          elemnt6 = 0.0_real64
          elemnt7 = 0.0_real64
          elemnt1vec(:) = 0.0_real64                               ! JE ADD
!*        ENDDO
        end if
        do i= 1,ncount
           read(one_body_cache_unit) tshell_r(i),lab
           ia = lab/orbital_label_radix
           ib = mod(lab,orbital_label_radix)
!*           LOC = INDEX(I)
           elemnt1 = elemnt1 + dint1(ia,ib)*tshell_r(i)
           do 21 l = 2,nrnuc                                                  ! JE ADD
              elemnt1vec(l) = elemnt1vec(l) + dint1vec(ia,ib,l)*       &     ! JE ADD
                   tshell_r(i)                                               ! JE ADD
   21      continue                                                          ! JE ADD
           elemnt2 = elemnt2 + dint2(ia,ib)*tshell_r(i)
           elemnt3 = elemnt3 + dint3(ia,ib)*tshell_r(i)
           elemnt4 = elemnt4 + dint4(ia,ib)*tshell_r(i)
           elemnt5 = elemnt5 + dint5(ia,ib)*tshell_r(i)
           elemnt6 = elemnt6 + dint6(ia,ib)*tshell_r(i)
           elemnt7 = elemnt7 + dint7(ia,ib)*tshell_r(i)
        end do
!
!   Return to the start of the loop
!

!*      ICI = 0
!*      DO 21 I = 1,NELMNT
!*        IRI = IROW(I)
!*         IF (I .GT. IENDC(ICI)) ICI = ICI+1
         do 22 j = 1,nvec
            loc = (j-1)*csf_count
            contri1 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt1
            contri1vec(:) = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt1vec(:)           ! JE ADD
            contri2 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt2
            contri3 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt3
            contri4 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt4
            contri5 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt5
            contri6 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt6
            contri7 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt7
            if (ir/=ic) then
               contri1 = 2.0_real64 * contri1
               contri1vec(:) = 2.0_real64 * contri1vec(:)                        ! JE ADD
               contri2 = 2.0_real64 * contri2
               contri3 = 2.0_real64 * contri3
               contri4 = 2.0_real64 * contri4
               contri5 = 2.0_real64 * contri5
               contri6 = 2.0_real64 * contri6
               contri7 = 2.0_real64 * contri7
            end if
            dens1(j) = dens1(j) + contri1
            do 23 l = 2,nrnuc                                                  ! JE ADD
               dens1vec(j,l) = dens1vec(j,l) + contri1vec(l)                  ! JE ADD
   23       continue                                                          ! JE ADD
            dens2(j) = dens2(j) + contri2
            dens3(j) = dens3(j) + contri3
            dens4(j) = dens4(j) + contri4
            dens5(j) = dens5(j) + contri5
            dens6(j) = dens6(j) + contri6
            dens7(j) = dens7(j) + contri7
   22    continue
!*   21 CONTINUE

      goto 16
      end if
      return
      end subroutine densread_seltz
