! Fill a real vector at the requested stride, including reverse traversal for a negative stride.
module real_vector_initialization
    implicit none (type, external)

    ! Module interface
    private
    public :: dinit
contains
    subroutine dinit(element_count, fill_value, values, stride)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: element_count
        integer, intent(in) :: stride
        real(real64), intent(in) :: fill_value
        real(real64), intent(out) :: values(*)

        ! Local declarations
        integer :: first_index, i
        if (stride == 1) then
            values(:element_count) = fill_value
        else
            first_index = 1
            if (stride < 0) first_index = ((-element_count) + 1)*stride + 1
            values(first_index:(element_count-1)*stride+first_index:stride) = fill_value
        end if
        return
    end subroutine dinit
end module real_vector_initialization
