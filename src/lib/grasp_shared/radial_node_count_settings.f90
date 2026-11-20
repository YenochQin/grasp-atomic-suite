! Stores the radial node-counting threshold and diagnostic output settings.
      module radial_node_count_settings
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      real(real64) :: thresh
      integer :: count_context = 0
      ! The current rank is zero in serial. Only rank zero writes the node trace.
      integer :: count_trace_rank = 0
      end module radial_node_count_settings
