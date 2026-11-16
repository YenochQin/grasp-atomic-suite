! Compatibility entry point for downstream users of the installed modules.
module vast_kind_param
    use, intrinsic :: iso_fortran_env, only: byte => int8, short => int16, long => int64, double => real64, &
        byte_log => int8, short_log => int16, long_log => int64, dble_complex => real64
    use memory_man, only: extended, double_ext, ext_complex
    implicit none (type, external)
end module vast_kind_param
