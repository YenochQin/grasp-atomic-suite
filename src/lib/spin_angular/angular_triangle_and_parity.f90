! Test triangle bounds and even parity for three integer angular labels.
module angular_triangle_and_parity
    implicit none (type, external)

    ! Module interface
    private
    public :: satisfies_angular_triangle_and_parity
contains
    integer function satisfies_angular_triangle_and_parity (first_angular_label, second_angular_label, &
        & third_angular_label) result(triangle_allowed)
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_angular_label, second_angular_label, third_angular_label
        triangle_allowed = 0
        if (iabs(first_angular_label - second_angular_label) > third_angular_label) return
        if (first_angular_label + second_angular_label < third_angular_label) return
        if (mod(first_angular_label + second_angular_label + third_angular_label,2) /= 0) return
        triangle_allowed = 1
        return
    end function satisfies_angular_triangle_and_parity
end module angular_triangle_and_parity
