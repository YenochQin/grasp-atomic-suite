! Stores list offsets by multipole rank for the two Breit integral groups.
      module breit_rank_offsets
      implicit none (type, external)
      integer, parameter :: kmax = 20
      integer :: kstartbreit1(0:kmax),kstartbreit2(0:kmax)
      end module breit_rank_offsets
