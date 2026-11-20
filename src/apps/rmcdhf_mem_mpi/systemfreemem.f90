!*******************************************************************
!                                                                  *
      function systemfreemem()
      use, intrinsic :: iso_fortran_env, only: real32
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
      integer :: memory_info_unit
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!MemTotal:       264566412 kB
!MemTotal:       1056570236 kB
      real(real64) :: systemfreemem
      integer      :: summem, leng
      character(len=128) :: line,cmem
!-----------------------------------------------
!
      open(newunit=memory_info_unit,file="/proc/meminfo", status="OLD", form="FORMATTED")
      read(memory_info_unit,"(A)")line
      read(memory_info_unit,"(A)")line
      close(memory_info_unit)
      leng = len_trim(line)
      cmem = line(10:leng-3)
      read(cmem,*)summem  ! KB
      systemfreemem = summem * 1.0_real32 / 1048576.0_real32  ! GB

      return
      end function systemfreemem
