! Intersect the triangle bounds of two pairs of angular momenta.
module shared_tensor_rank_range
    implicit none (type, external)
    private
    public :: itrexg
contains
      integer function itrexg(two_j_1,two_j_2,two_j_3,two_j_4,two_rank_range_size) result(min_two_rank)
      implicit none (type, external)
      integer, intent(in)  :: two_j_1, two_j_2, two_j_3, two_j_4
      integer, intent(out) :: two_rank_range_size
      integer :: lower_two_rank
      lower_two_rank=max0(iabs(two_j_1-two_j_2),iabs(two_j_3-two_j_4))
      two_rank_range_size=min0(iabs(two_j_1+two_j_2),iabs(two_j_3+two_j_4))-lower_two_rank+1
      min_two_rank=lower_two_rank
      return
      end function itrexg
end module shared_tensor_rank_range
