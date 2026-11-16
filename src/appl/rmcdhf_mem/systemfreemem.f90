!*******************************************************************
!                                                                  *
      function systemfreemem()
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!MemTotal:       264566412 kB
!MemTotal:       1056570236 kB
      real(double) :: systemfreemem
      integer      :: summem, leng
      character(len=128) :: line,cmem
!-----------------------------------------------
!
      open(5321,file='/proc/meminfo', status='OLD', form='FORMATTED')
      read(5321,'(A)')line
      read(5321,'(A)')line
      close(5321)
      leng = len_trim(line)
      cmem = line(10:leng-3)
      read(cmem,*)summem  ! KB
      systemfreemem = summem * 1.0 / 1048576.0  ! GB

      return
      end function systemfreemem
