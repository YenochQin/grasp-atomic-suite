! Store the electronic isotope-shift options selected for a RIS4 run.
module isotope_shift_options
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    !  Module state
    ! Translated by Gediminas Gaigalas 11/18/19
    real(real64), pointer :: contact_density_by_level(:), total_nms_by_level(:), radius_by_level(:),  &
        radius_squared_by_level(:), inverse_radius_by_level(:), inverse_radius_squared_by_level(:),  &
        leading_nms_by_level(:), leading_sms_by_level(:), total_sms_by_level(:)
end module isotope_shift_options
