module generator_I
    implicit none
    interface
        subroutine rcsfexcitation()
        end subroutine
        subroutine jjgen15()
        end subroutine
        subroutine rcsfsym()
        end subroutine
        subroutine rcsfblock()
        end subroutine
        subroutine rcsfsym_cyc()
        end subroutine
        subroutine removeblank(line1, m, line2, n)
            character(300) :: line1, line2
            integer :: m, n
        end subroutine
        subroutine indexa(n, a, ldown, indx)
            integer :: n, indx(n)
            character(*) :: a(n)
            logical :: ldown
        end subroutine
    end interface
end module
