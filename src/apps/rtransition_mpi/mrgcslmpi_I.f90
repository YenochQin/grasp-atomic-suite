      module mrgcsl_i
      implicit none (type, external)
      interface
      subroutine mrgcsl (state_name)
      implicit none (type, external)
      character (len = 128) :: state_name(2)
      end subroutine mrgcsl
      end interface
      end module mrgcsl_i
