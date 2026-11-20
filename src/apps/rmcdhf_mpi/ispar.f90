integer function ispar(csf_index)
   use rmcdhf_csf_input, only: rmcdhf_csf_symmetry
   implicit none (type, external)
   integer, intent(in) :: csf_index
   ispar = sign(1,rmcdhf_csf_symmetry(csf_index))
end function ispar
