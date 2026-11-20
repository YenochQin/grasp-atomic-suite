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
    subroutine combine_vectors(c, a, b, faca, facb, ndim)
        ! Module dependencies
        !
        !     CACLULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: ndim
        real(real64), intent(in) :: faca
        real(real64), intent(in) :: facb
        real(real64), intent(in)  :: a(ndim)
        real(real64), intent(in)  :: b(ndim)
        real(real64), intent(out) :: c(ndim)

        ! Local declarations
        !-----------------------------------------------
        !
        if (faca/=0.0_real64 .and. facb/=0.0_real64) then
            c(:ndim) = faca*a(:ndim) + facb*b(:ndim)
        !
        else if (faca==0.0_real64 .and. facb/=0.0_real64) then
            c(:ndim) = facb*b(:ndim)
        !
        else if (faca/=0.0_real64 .and. facb==0.0_real64) then
            c(:ndim) = faca*a(:ndim)
        !
        else if (faca==0.0_real64 .and. facb==0.0_real64) then
            c(:ndim) = 0.0_real64
        end if
        !
        return
    end subroutine combine_vectors
end module combine_vectors_interface
