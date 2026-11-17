integer function itjpo(icsf)
   use rmcdhf_descriptor_input, only: rmcdhf_csf_symmetry
   implicit none (type, external)
   integer, intent(in) :: icsf
   itjpo = abs(rmcdhf_csf_symmetry(icsf))
end function itjpo
