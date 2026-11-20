module getinf_i
   use csf_descriptor_reader, only: descriptor_source
   implicit none (type, external)
   interface
      subroutine getinf(source, status)
         import descriptor_source
         implicit none (type, external)
         type(descriptor_source), intent(in) :: source
         integer, intent(out) :: status
      end subroutine getinf
   end interface
end module getinf_i
