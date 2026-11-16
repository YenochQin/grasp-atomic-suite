!***********************************************************************
!                                                                      *
      function systemmem()
!                                                                      *
!***********************************************************************
!...Written by Chongyang Chen, Fudan University, Shanghai, China       *
!                                                Oct. 2020             *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!MemTotal:       264566412 kB
!MemTotal:       1056570236 kB
      real(real64) :: systemmem
      integer      :: summem, leng
      character(len=128) :: line,cmem
!-----------------------------------------------
!
      open(5321,file="/proc/meminfo", status="OLD", form="FORMATTED")
      read(5321,"(A)")line
      close(5321)
      leng = len_trim(line)
      cmem = line(10:leng-3)
      read(cmem,*)summem  ! KB
!GG      SYSTEMMem = SUMMem / 1.0D6  ! GB
      systemmem = summem / 1048576  ! GB

      return
      end function systemmem
