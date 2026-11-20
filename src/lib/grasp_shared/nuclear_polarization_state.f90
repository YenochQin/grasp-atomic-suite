! Stores the Coulomb baseline and the published nuclear-polarization configuration and potentials.
module nuclear_polarization_state
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters, only: max_radial_points
    implicit none (type, external)

    ! Module interface
    private

    ! Module constants
    integer, parameter, public :: npol_message_len = 512
    character(*), parameter, public :: npol_model = "flambaum2021_scalar_e1_e2"

    ! Derived types
    type, public :: npol_config
        ! Components
        integer :: schema_version = 1, isotope_z = 0, isotope_a = 0
        logical :: enabled = .false., use_e1 = .false., use_e2 = .false.
        real(real64) :: k1 = -1, b1 = -1, k2 = -1, b2 = -1
        character(16) :: units = ""
        character(64) :: model_id = ""
        character(1024) :: parameter_source = ""
    end type npol_config

    ! Module state
    type(npol_config), protected, public :: npol_settings
    logical, protected, public :: npol_ready = .false.
    integer, protected, public :: npol_n = 0
    real(real64), protected, public :: coulomb_potential_baseline(max_radial_points) = 0
    real(real64), protected, public :: vnp_e1(max_radial_points) = 0, vnp_e2(max_radial_points) = 0, &
        & vnp_total(max_radial_points) = 0

    ! Module interface
    public :: npol_reset, npol_capture, npol_publish
contains
    subroutine npol_reset()
        npol_settings = npol_config()
        npol_ready = .false.
        npol_n = 0
        coulomb_potential_baseline = 0
        vnp_e1 = 0
        vnp_e2 = 0
        vnp_total = 0
    end subroutine npol_reset

    ! Capture the unmodified -r*V_Coulomb before composing the nuclear potential.
    subroutine npol_capture(n, coul)
        ! Arguments
        integer, intent(in) :: n
        real(real64), intent(in) :: coul(n)
        if (n < 1 .or. n > max_radial_points) return
        npol_n = n
        coulomb_potential_baseline(:n) = coul
    end subroutine npol_capture

    ! Publish the first n grid values after the caller validates and constructs the potentials.
    subroutine npol_publish(cfg, n, e1, e2, total)
        ! Arguments
        type(npol_config), intent(in) :: cfg
        integer, intent(in) :: n
        real(real64), intent(in) :: e1(n), e2(n), total(n)
        if (n /= npol_n .or. n < 1 .or. n > max_radial_points) return
        npol_settings = cfg
        vnp_e1(:n) = e1
        vnp_e2(:n) = e2
        vnp_total(:n) = total
        npol_ready = .true.
    end subroutine npol_publish
end module nuclear_polarization_state
