! Test the four triangle conditions of a Wigner 6j symbol.
module wigner_6j_triangle_selection
    implicit none (type, external)

    ! Module interface
    private
    public :: allows_wigner_6j_triangles
contains
    integer function allows_wigner_6j_triangles (two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, two_j_6) result(triangles_allowed)
        ! Module dependencies
        use angular_triangle_and_parity, only: satisfies_angular_triangle_and_parity
        implicit none (type, external)

        ! Arguments
        integer  :: two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, two_j_6
        triangles_allowed = 0
        if (satisfies_angular_triangle_and_parity(two_j_1,two_j_2,two_j_3) == 0) return
        if (satisfies_angular_triangle_and_parity(two_j_1,two_j_5,two_j_6) == 0) return
        if (satisfies_angular_triangle_and_parity(two_j_4,two_j_2,two_j_6) == 0) return
        if (satisfies_angular_triangle_and_parity(two_j_4,two_j_5,two_j_3) == 0) return
        triangles_allowed = 1
        return
    end function allows_wigner_6j_triangles
end module wigner_6j_triangle_selection
