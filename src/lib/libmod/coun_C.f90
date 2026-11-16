      module coun_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      real(real64) :: thresh
      integer :: count_context = 0
!     Every rank replicates the radial solve and shares one rmcdhf.log, so
!     the node trace is written by one rank only.  Serial callers leave this
!     at 0 and keep tracing.
      integer :: count_trace_rank = 0
      end module coun_c
