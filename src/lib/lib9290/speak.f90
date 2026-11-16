module speak_i
    implicit none (type, external)
    private
    public :: speak
contains
!***********************************************************************
!                                                                      *
      subroutine speak(ja, jb, ia1, ib1, ia2, ib2, k, x)
!                                                                      *
!   Output MCP  coefficients and integral parameters to COMMON block   *
!   /BUFFER/. Also print these if  IBUG1 = 1 .                         *
!                                                                      *
!                                           Last Update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use buffer_c,        only: angular_coefficient_labels, angular_coefficients, angular_buffer_capacity, &
         & angular_coefficient_count
      use debug_c ,        only: ibug1
      use orb_c,           only: orbital_principal_n, orbital_symmetry_labels
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcbuf_i, only: alcbuf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ja
      integer, intent(in) :: jb
      integer, intent(in) :: ia1
      integer, intent(in) :: ib1
      integer, intent(in) :: ia2
      integer, intent(in) :: ib2
      integer, intent(in) :: k
      real(real64), intent(in) :: x
!-----------------------------------------------
!
!
      if (ibug1 /= 0) write (99, 300) ja, jb, orbital_principal_n(ia1), orbital_symmetry_labels(ia1), &
         & orbital_principal_n(ib1), orbital_symmetry_labels(ib1&
         ), orbital_principal_n(ia2), orbital_symmetry_labels(ia2), orbital_principal_n(ib2), orbital_symmetry_labels(ib2), k, x
!
!   Increment counter
!
      angular_coefficient_count = angular_coefficient_count + 1
!
!   Ensure that arrays are of adequate size; reallocate if necessary
!
      if (angular_coefficient_count > angular_buffer_capacity) call alcbuf (2)
!
!   Store integral indices and coefficient in COMMON/BUFFER/
!
      angular_coefficient_labels(1,angular_coefficient_count) = ia1
      angular_coefficient_labels(2,angular_coefficient_count) = ib1
      angular_coefficient_labels(3,angular_coefficient_count) = ia2
      angular_coefficient_labels(4,angular_coefficient_count) = ib2
      angular_coefficient_labels(5,angular_coefficient_count) = k
      angular_coefficients(angular_coefficient_count) = x
!
      return
!
  300 format(2(1x,1i2),4(1x,i2,a2),1x,i2,1x,1p,d19.12)
      return
!
      end subroutine speak
end module speak_i
