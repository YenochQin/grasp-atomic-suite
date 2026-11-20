! Stores intermediate-term index ranges by occupation parity for jj-coupled subshells with 2j < 9.
      module jj_reduced_coefficients
      implicit none (type, external)
      integer :: imptjj(63), imgtjj(63), impnjj(63), imgnjj(63)
      data impnjj/2,1,4,2*3,3*9,3*6,6*18,8*12,20*46,18*26/
      data imgnjj/2,1,5,2*3,3*11,3*8,6*25,8*17,20*63,18*45/
      data imptjj/1,2,3,2*4,3*6,3*9,6*12,8*18,20*26,18*46/
      data imgtjj/1,2,3,2*5,3*8,3*11,6*17,8*25,20*45,18*63/
      end module jj_reduced_coefficients
