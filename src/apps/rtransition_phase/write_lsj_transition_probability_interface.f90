module write_lsj_transition_probability_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_lsj_transition_probability (inum_ii,inum_ff,asfa,asfb,i,j,omega,factor)
            ! Module dependencies
            !...Translated by Charlotte Froese Fischer
            !                       Gediminas Gaigalas  10/05/17
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer  :: i, j, inum_ii, inum_ff
            real(real64), intent(in) :: asfa, asfb
            real(real64), intent(in) :: omega, factor
        end subroutine write_lsj_transition_probability
    end interface
end module write_lsj_transition_probability_interface
