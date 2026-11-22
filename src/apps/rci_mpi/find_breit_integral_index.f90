! Find the position of a packed Breit integral label in its sorted table.
module find_breit_integral_index_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: find_breit_integral_index
contains

    ! Written by Farid A Parpia. Last revision: 09 Oct 1992
    integer function find_breit_integral_index (itype, i)
        ! Module dependencies
        use breit_integral_tables, only: breit_integral_indices_1, breit_integral_indices_2, breit_integral_indices_3, &
            breit_integral_indices_4, breit_integral_indices_5, breit_integral_indices_6
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: itype
        integer  :: i

        select case (itype)
        case (1)
            find_breit_integral_index = breit_integral_indices_1(i)
        case (2)
            find_breit_integral_index = breit_integral_indices_2(i)
        case (3)
            find_breit_integral_index = breit_integral_indices_3(i)
        case (4)
            find_breit_integral_index = breit_integral_indices_4(i)
        case (5)
            find_breit_integral_index = breit_integral_indices_5(i)
        case (6)
            find_breit_integral_index = breit_integral_indices_6(i)
        end select
        return
    end function find_breit_integral_index
end module find_breit_integral_index_interface
