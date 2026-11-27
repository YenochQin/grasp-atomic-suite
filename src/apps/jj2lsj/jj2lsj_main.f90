! Load the jj-coupled basis and run the LSJ conversion for the selected states.
! When both relativistic partners occupy a shell, list j = l - 1/2 before j = l + 1/2.
! The LS basis follows the same shell order.
! The LS–jj coefficient tables are in the ls_jj_coefficients modules.
! Written by G. Gaigalas,
!     NIST                                                  May 2011
!     VILNIUS                                               May 2017
!
! Modified by G. Gaigalas and C. Cychen. 2021
! Modified by G. Gaigalas. 2022
! Modified by G. Gaigalas. Jan 2024
! Modified by G. Gaigalas. Jan 2025
program jj2lsj_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int64
    use jj2lsj_profile, only: profile_finish, profile_start
    use jj2lsj_conversion, only: run_ls_coupling_conversion
    use log_factorial_table, only: initialize_log_factorial_table
    use isotope_data_input, only: read_isotope_data
    use run_timer, only: begin_run_timer, finish_run_timer
    use physical_constants, only: initialize_physical_constants
    implicit none (type, external)
    integer(int64) :: start_clock_count
    call begin_run_timer (start_clock_count, "jj2lsj_omp")
    call profile_start()
    print *, " "
    print *, "jj2lsj_omp: Transformation of ASFs from a jj-coupled CSF basis"
    print *, "        into an LS-coupled CSF basis  (Fortran 95 version)"
    print *, "        (C) Copyright by   G. Gaigalas and Ch. F. Fischer,"
    print *, "        (2025)."
    print *, "        Input files: name_desc.parquet, name_header.toml, name.mix/name.cimix"
    print *, "         (optional)  name.lsj.transform.cache"
    print *, "        Output files: name.mcdhf.lsj.labels or name.ci.lsj.labels,"
    print *, "         (optional)  With the same name.mcdhf.lsj or name.ci.lsj stem:"
    print *, "                     .csf-list, .expansion, .unique-labels, .unique-labels.summary"
    print *, "         (optional)  Transformation cache:"
    print *, "                     name.lsj.transform.cache"
    print *, "         (optional)  .profile.csv with the same name.mcdhf.lsj or name.ci.lsj stem"
    print *, " "
    ! Initialize the log-factorial table for angular algebra.
    call initialize_physical_constants
    call initialize_log_factorial_table
    call read_isotope_data("isodata")
    call run_ls_coupling_conversion
    call profile_finish()
    call finish_run_timer (start_clock_count, "jj2lsj_omp")
    stop
end program jj2lsj_main
