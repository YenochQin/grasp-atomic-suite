integer function ispar(icsf)
   use rmcdhf_descriptor_input, only: rmcdhf_csf_symmetry
   implicit none (type, external)
   integer, intent(in) :: icsf
   ispar = sign(1,rmcdhf_csf_symmetry(icsf))
end function ispar
