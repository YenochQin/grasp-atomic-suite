      module rkco_gg_i
      interface
!
      subroutine rkco_gg (ja,jb,coefficient_consumer,incor,icolbrei)
!      USE vast_kind_param, ONLY:  DOUBLE
      use cord_i, only: cord_callback => cord
      procedure(cord_callback) :: coefficient_consumer
      integer, intent(in) :: ja,jb,incor,icolbrei
      end subroutine
      end interface
      end module
