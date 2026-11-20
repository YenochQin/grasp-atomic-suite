! Keep memory reporting deterministic independently of the test host.
module query_total_system_memory_interface
    implicit none (type, external)
    private
    public :: query_total_system_memory
contains
    function query_total_system_memory() result(value)
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)
        real(real64) :: value
        value = 1.0_real64
    end function query_total_system_memory
end module query_total_system_memory_interface

module query_available_system_memory_interface
    implicit none (type, external)
    private
    public :: query_available_system_memory
contains
    function query_available_system_memory() result(value)
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)
        real(real64) :: value
        value = 1.0_real64
    end function query_available_system_memory
end module query_available_system_memory_interface
