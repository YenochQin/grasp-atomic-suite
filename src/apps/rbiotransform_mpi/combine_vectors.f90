module combine_vectors_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: combine_vectors
contains

    !
    !     ------------------------------------------------------------------
    !     V E C S U M
    !     ------------------------------------------------------------------
    !
    subroutine combine_vectors(combined_vector, first_vector, second_vector, first_vector_weight, second_vector_weight, &
        & element_count)
        ! Module dependencies
        !
        !     CACLULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: element_count
        real(real64), intent(in) :: first_vector_weight
        real(real64), intent(in) :: second_vector_weight
        real(real64), intent(in)  :: first_vector(element_count)
        real(real64), intent(in)  :: second_vector(element_count)
        real(real64), intent(out) :: combined_vector(element_count)

        ! Local declarations
        !-----------------------------------------------
        !
        if (first_vector_weight/=0.0_real64 .and. second_vector_weight/=0.0_real64) then
            combined_vector(:element_count) = first_vector_weight*first_vector(:element_count) + second_vector_weight* &
                & second_vector(:element_count)
        !
        else if (first_vector_weight==0.0_real64 .and. second_vector_weight/=0.0_real64) then
            combined_vector(:element_count) = second_vector_weight*second_vector(:element_count)
        !
        else if (first_vector_weight/=0.0_real64 .and. second_vector_weight==0.0_real64) then
            combined_vector(:element_count) = first_vector_weight*first_vector(:element_count)
        !
        else if (first_vector_weight==0.0_real64 .and. second_vector_weight==0.0_real64) then
            combined_vector(:element_count) = 0.0_real64
        end if
        !
        return
    end subroutine combine_vectors
end module combine_vectors_interface
