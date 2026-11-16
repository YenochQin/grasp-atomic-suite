      module bilst_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer :: ndtpa(6), ntpi(6)
      integer, pointer :: indtp1(:), indtp2(:), indtp3(:), indtp4(:), &
                                        indtp5(:), indtp6(:)
      real(real64), pointer :: valtp1(:), valtp2(:), valtp3(:), valtp4(:), &
                                        valtp5(:), valtp6(:)
      logical :: first(6)
      end module bilst_c
