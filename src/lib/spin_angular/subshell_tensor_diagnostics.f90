! Report unsupported occupation or subshell term encoding and stop.
module subshell_tensor_diagnostics
    implicit none (type, external)
    private
    public :: mes
contains
      subroutine mes(diagnostic_code)
      implicit none (type, external)
      integer, intent(in) :: diagnostic_code
      integer :: message_index
      character(len=10) :: routine_names(6)
      data routine_names/"   I T L S"," I T L S 2"," I T L S 3", &
                   "AW P 1 L S","WA P 1 L S","       W 1"/
      if(diagnostic_code > 50) then
        message_index=diagnostic_code-50
        write(6,"(A)") " error in func./sub. "
        write(6,"(20X,A10)") routine_names(message_index)
        write(6,"(A)") " susimaise f sluoksnio termu kodavimas  "
      else
        write(6,"(A)") " yra daugiau nei 2 ele. sluoks. f,g,h,i,k,l,m"
        write(6,"(3X,I5)") diagnostic_code
        if(diagnostic_code == 1) then
          write(6,"(A)") " error in Subroutine   W 1 G "
        else if(diagnostic_code == 2) then
          write(6,"(A)") " error in Subroutine   W 1 "
        else if(diagnostic_code == 11) then
          write(6,"(A)") " error in Function     I T L S  "
        else if(diagnostic_code == 12) then
          write(6,"(A)") " error in Function     I T L S 2  "
        else if(diagnostic_code == 13) then
          write(6,"(A)") " error in Function     I T L S 3  "
        else if(diagnostic_code == 30) then
          write(6,"(A)") " error in Subroutine   A 1 A 2 A 3 A 4 L S "
        else if(diagnostic_code == 31) then
          write(6,"(A)") " error in Subroutine   A 1 A 2 L S "
        else if(diagnostic_code == 32) then
          write(6,"(A)") " error in Subroutine   A 1 A 2 W 3 L S "
        else if(diagnostic_code == 33) then
          write(6,"(A)") " error in Subroutine   A 1 A W 2 L S "
        else if(diagnostic_code == 34) then
          write(6,"(A)") " error in Subroutine   W A 1 A 2 L S "
        else if(diagnostic_code == 35) then
          write(6,"(A)") " error in Subroutine   W 1 W 2 L S "
        else
          write(6,"(A)") " error in unknown Subroutine  "
        end if
      end if
      write(6,"(A)") " Contact to   G. Gediminas please ! ! ! "
      stop
      end subroutine mes
end module subshell_tensor_diagnostics
