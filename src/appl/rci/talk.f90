!***********************************************************************
!                                                                      *
      subroutine talk(ja, jb, nu, ia, ib, ic, id, itype, coef)
!                                                                      *
!   Print  coefficients  and  integral  parameters  if IBUG1 > 0 and   *
!   write to disk.                                                     *
!                                                                      *
!                                           Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix
      use buffer_c
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use orb_c
      use cons_c, only: double, eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcbuf_i, only: alcbuf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ja, jb, nu, ia, ib, ic, id, itype
      real(real64), intent(in) :: coef
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   Print coefficient if requested
!
      if (ibug1 /= 0) write (99, 300)ja,jb,orbital_principal_n(ia),orbital_symmetry_labels(ia),orbital_principal_n(ib), &
         & orbital_symmetry_labels(ib),&
         orbital_principal_n(ic),orbital_symmetry_labels(ic),orbital_principal_n(id),orbital_symmetry_labels(id),nu,itype,coef
!
!   Increment counter
!
      if(dabs(coef) > eps) then
         angular_coefficient_count = angular_coefficient_count + 1
!
!   Ensure that arrays are of adequate size; reallocate if necessary
!
         if (angular_coefficient_count > angular_buffer_capacity) call alcbuf (2)
!
!   Store integral indices and coefficient in COMMON/BUFFER/
!
         angular_coefficient_labels(1,angular_coefficient_count) = ia
         angular_coefficient_labels(2,angular_coefficient_count) = ib
         angular_coefficient_labels(3,angular_coefficient_count) = ic
         angular_coefficient_labels(4,angular_coefficient_count) = id
         angular_coefficient_labels(5,angular_coefficient_count) = nu
         angular_coefficient_labels(6,angular_coefficient_count) = itype
         angular_coefficients(angular_coefficient_count)   = coef
      end if
!
      return
!
  300 format(2(1x,1i2),4(1x,i2,a2),1x,1i2,1x,1i2,1x,1p,d19.12)
      return
!
      end subroutine talk
