!***********************************************************************
!                                                                      *
integer function find_breit_integral_index (itype, i)
    ! Module dependencies
    !                                                                      *
    !   Written by Farid A Parpia             Last revision: 09 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
        valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: itype
    integer  :: i
    !-----------------------------------------------
    !
    select case (itype)
    case (1)
        find_breit_integral_index = indtp1(i)
    case (2)
        find_breit_integral_index = indtp2(i)
    case (3)
        find_breit_integral_index = indtp3(i)
    case (4)
        find_breit_integral_index = indtp4(i)
    case (5)
        find_breit_integral_index = indtp5(i)
    case (6)
        find_breit_integral_index = indtp6(i)
    end select
    !
    return
end function find_breit_integral_index
