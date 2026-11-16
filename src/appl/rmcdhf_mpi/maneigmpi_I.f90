      module maneigmpi_i
      implicit none (type, external)
      interface
      subroutine maneigmpi (dvdfirst,lprint, jblock, ncfpat, ncminpat, nevecpat, ncftot)
      implicit none (type, external)
      logical, intent(in) ::  dvdfirst
      logical :: lprint
!VAST...Dummy argument LPRINT is not referenced in this routine.
      integer, intent(in) :: jblock
      integer, intent(in) :: ncfpat
      integer, intent(in) :: ncminpat
      integer, intent(in) :: nevecpat
      integer :: ncftot
      end subroutine
      end interface
      end module
