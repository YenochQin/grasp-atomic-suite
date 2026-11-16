      module rkco_gg_i
      implicit none (type, external)
      interface
!
      subroutine rkco_gg (ja,jb,coefficient_consumer,incor,icolbrei)
      use cord_i, only: cord_callback => cord
      implicit none (type, external)
      procedure(cord_callback) :: coefficient_consumer
      integer, intent(in) :: ja,jb,incor,icolbrei
      end subroutine
      end interface
      end module
