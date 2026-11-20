! Test the four triangle conditions of a Wigner 6j symbol.
module wigner_6j_triangle_selection
    implicit none (type, external)
    private
    public :: ixjtik
contains
      integer function ixjtik (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, two_j_6) result(triangles_allowed)
      use angular_triangle_and_parity, only: ittk
      implicit none (type, external)
      integer  :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, two_j_6
      triangles_allowed = 0
      if (ittk(two_j_1,two_j_2,two_j_3) == 0) return
      if (ittk(two_j_1,two_j_5,two_j_6) == 0) return
      if (ittk(two_j_4,two_j_2,two_j_6) == 0) return
      if (ittk(two_j_4,two_j_5,two_j_3) == 0) return
      triangles_allowed = 1
      return
      end function ixjtik
end module wigner_6j_triangle_selection
