module apply_subshell_transfer_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine apply_subshell_transfer(input_mixing_coefficients,csf_count,level_count,i,l,t,subshell_count, &
            & transformed_mixing_coefficients,ntestg)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: csf_count, level_count, subshell_count, i, l, ntestg
            real(real64), intent(in) :: input_mixing_coefficients(csf_count, level_count)
            real(real64), intent(in) :: t(subshell_count)
            real(real64), intent(out) :: transformed_mixing_coefficients(csf_count, level_count)
        end subroutine apply_subshell_transfer
    end interface
end module apply_subshell_transfer_interface
