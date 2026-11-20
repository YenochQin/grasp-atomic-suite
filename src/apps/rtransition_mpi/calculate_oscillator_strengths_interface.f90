module calculate_oscillator_strengths_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine calculate_oscillator_strengths (name,fullname,tmpdir,startup_directory,block_label_text)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(inout) :: name(2)
            character(len=3)  ::  block_label_text
            character (len = 128), intent(inout) :: fullname(2)
            character(len=128) :: startup_directory,tmpdir
        end subroutine calculate_oscillator_strengths
    end interface
end module calculate_oscillator_strengths_interface
