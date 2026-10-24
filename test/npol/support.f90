module npol_test_support
    use vast_kind_param, only: double
    use npol_C, only: npol_config, NPOL_MODEL
    implicit none
contains
    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (condition) return
        write(*,'(A)') 'FAIL: '//label
        error stop 2
    end subroutine

    subroutine near(actual, expected, label)
        real(double), intent(in) :: actual, expected
        character(*), intent(in) :: label
        ! Arithmetic tolerance, not a physical accuracy threshold.
        call require(abs(actual-expected) <= 64*epsilon(actual)*max(1.0_double,abs(expected)), label)
    end subroutine

    function fixture() result(cfg)
        type(npol_config) :: cfg
        ! Software-only constants. These do not describe an isotope model.
        cfg%enabled = .true.
        cfg%use_e1 = .true.
        cfg%use_e2 = .true.
        cfg%units = 'atomic'
        cfg%model_id = NPOL_MODEL
        cfg%parameter_source = 'artificial software fixture; E2 convention: test only'
        cfg%isotope_z = 2
        cfg%isotope_a = 4
        cfg%k1 = 16
        cfg%b1 = 2
        cfg%k2 = 128
        cfg%b2 = 2
    end function

    subroutine write_fixture(path)
        character(*), intent(in) :: path
        integer :: u
        open(newunit=u, file=path, status='replace')
        write(u,'(A)') '&npol schema_version=1, enabled=.true., use_e1=.true., use_e2=.true.,', &
            "units='atomic', model_id='flambaum2021_scalar_e1_e2',", &
            "parameter_source='artificial software fixture; E2 convention: test only',", &
            'isotope_z=2, isotope_a=4, k1=16, b1=2, k2=128, b2=2 /'
        close(u)
    end subroutine
end module npol_test_support
