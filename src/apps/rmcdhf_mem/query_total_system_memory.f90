! Query the installed memory capacity used by MCP cache diagnostics.
module query_total_system_memory_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: query_total_system_memory
contains

    ! Written by Chongyang Chen, Fudan University, Shanghai, China
    !                                                 Oct. 2020
    !
    ! Translated by Gediminas Gaigalas  May 2021
    function query_total_system_memory()
        !  Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: memory_info_unit
        ! MemTotal:       264566412 kB
        ! MemTotal:       1056570236 kB
        real(real64) :: query_total_system_memory
        integer      :: summem, leng
        character(len=128) :: line,cmem

        open(newunit=memory_info_unit,file="/proc/meminfo", status="OLD", form="FORMATTED")
        read(memory_info_unit,"(A)")line
        close(memory_info_unit)
        leng = len_trim(line)
        cmem = line(10:leng-3)
        read(cmem,*)summem  ! KB
        ! SYSTEMMem = SUMMem / 1.0D6  ! GB
        query_total_system_memory = summem / 1048576  ! GB

        return
    end function query_total_system_memory
end module query_total_system_memory_interface
