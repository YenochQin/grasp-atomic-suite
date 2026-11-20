module read_ci_level_selection_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine read_ci_level_selection(iatjpo, iaspar)
            implicit none (type, external)

            ! Arguments
            integer, intent(in)  :: iatjpo, iaspar
        end subroutine read_ci_level_selection
    end interface
end module read_ci_level_selection_interface
