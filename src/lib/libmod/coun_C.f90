      module coun_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double) :: thresh
      integer :: count_context = 0
!     Every rank replicates the radial solve and shares one rmcdhf.log, so
!     the node trace is written by one rank only.  Serial callers leave this
!     at 0 and keep tracing.
      integer :: count_trace_rank = 0
      end module coun_c
