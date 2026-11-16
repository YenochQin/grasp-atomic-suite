      module alcnsa_i
      implicit none (type, external)
      interface
      subroutine alcnsa(jja, jjb, hb1, hb2, hc1, hc2, hm1, &
          hm2, lab, nptr, nsdim, imode)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer  :: nsdim
      integer, intent(in) :: imode
      integer, pointer :: jja(:), jjb(:), lab(:), nptr(:)
      real(real64), pointer :: hb1(:), hb2(:), hc1(:), hc2(:), hm1(:), hm2(:)
      end subroutine
      end interface
      end module
