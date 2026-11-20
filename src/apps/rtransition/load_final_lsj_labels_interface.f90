module load_final_lsj_labels_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine load_final_lsj_labels (name)
            implicit none (type, external)

            ! Arguments
            character(len=24) :: name
        end subroutine load_final_lsj_labels
    end interface
end module load_final_lsj_labels_interface
