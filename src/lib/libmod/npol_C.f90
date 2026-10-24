! NP state is published only after parsing, construction and output succeed.
module npol_C
    use vast_kind_param, only: double
    use parameter_def, only: NNNP
    implicit none
    private
    integer, parameter, public :: NPOL_MESSAGE_LEN = 512
    character(*), parameter, public :: NPOL_MODEL = 'flambaum2021_scalar_e1_e2'
    type, public :: npol_config
        integer :: schema_version = 1, isotope_z = 0, isotope_a = 0
        logical :: enabled = .false., use_e1 = .false., use_e2 = .false.
        real(double) :: k1 = -1, b1 = -1, k2 = -1, b2 = -1
        character(16) :: units = ''
        character(64) :: model_id = ''
        character(1024) :: parameter_source = ''
    end type
    type(npol_config), protected, public :: npol_settings
    logical, protected, public :: npol_ready = .false.
    integer, protected, public :: npol_n = 0
    real(double), protected, public :: zz_coul(NNNP) = 0
    real(double), protected, public :: vnp_e1(NNNP) = 0, vnp_e2(NNNP) = 0, vnp_total(NNNP) = 0
    public :: npol_reset, npol_capture, npol_publish
contains
    subroutine npol_reset()
        npol_settings = npol_config()
        npol_ready = .false.
        npol_n = 0
        zz_coul = 0
        vnp_e1 = 0
        vnp_e2 = 0
        vnp_total = 0
    end subroutine

    ! Initialization calls this once, immediately after NUCPOT, with checked N.
    ! Composition never captures its output as a new Coulomb baseline.
    subroutine npol_capture(n, coul)
        integer, intent(in) :: n
        real(double), intent(in) :: coul(n)
        if (n < 1 .or. n > NNNP) return
        npol_n = n
        zz_coul(:n) = coul
    end subroutine

    ! Internal publication seam; callers must finish all validation first.
    subroutine npol_publish(cfg, n, e1, e2, total)
        type(npol_config), intent(in) :: cfg
        integer, intent(in) :: n
        real(double), intent(in) :: e1(n), e2(n), total(n)
        if (n /= npol_n .or. n < 1 .or. n > NNNP) return
        npol_settings = cfg
        vnp_e1(:n) = e1
        vnp_e2(:n) = e2
        vnp_total(:n) = total
        npol_ready = .true.
    end subroutine
end module npol_C
