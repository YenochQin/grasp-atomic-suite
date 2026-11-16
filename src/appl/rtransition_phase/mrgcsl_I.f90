      module mrgcsl_i
      implicit none (type, external)
      interface
      subroutine mrgcsl (name)
      implicit none (type, external)
      character (len = 24) :: name(2)
      end subroutine mrgcsl
      end interface
      end module mrgcsl_i
