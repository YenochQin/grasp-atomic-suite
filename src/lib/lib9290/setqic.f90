!***********************************************************************
!                                                                      *
      subroutine setqic
!                                                                      *
!   This  subroutine sets up the coefficients for SUBROUTINEs DPBDT,   *
!   QUAD, RINTI, START, YZK, ZKF.                                      *
!                                                                      *
!   No SUBROUTINES called.                                             *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 05 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:39   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cnc_c
      use grid_c
      use lic13_c, a13=>a
      use ncc_c
      use sbc_c
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      real(real64) :: b13(13, 13)
      real(real64) :: cg(6)
      real(real64) :: c5num(5, 2:5)
      real(real64) :: c6num(6, 2:6)
      real(real64) :: b13den, denom, c5den, c6den, factor
      logical :: first
!-----------------------------------------------
!
!
!
!----------------------------------------------------------------------*
!                                                                      *
!   Thirteen-point  Lagrange  interpolation  coefficients for first    *
!   derivative                                                         *
!
      data (b13(1,i),i=1,13)/ -1486442880.0_real64, 5748019200.0_real64, &
         -15807052800.0_real64, 35126784000.0_real64, -59276448000.0_real64, &
         75873853440.0_real64, -73766246400.0_real64, 54195609600.0_real64, &
         -29638224000.0_real64, 11708928000.0_real64, -3161410560.0_real64, 522547200.0_real64&
         ,  - 39916800.0_real64/
      data (b13(2,i),i=1,13)/ -39916800.0_real64, -967524480.0_real64, 2634508800.0_real64&
         , -4390848000.0_real64, 6586272000.0_real64, -7903526400.0_real64, 7376624640.0_real64&
         , -5269017600.0_real64, 2822688000.0_real64, -1097712000.0_real64, 292723200.0_real64&
         , -47900160.0_real64, 3628800.0_real64/
      data (b13(3,i),i=1,13)/ 3628800.0_real64, -87091200.0_real64, -684478080.0_real64, &
         1596672000.0_real64, -1796256000.0_real64, 1916006400.0_real64, -1676505600.0_real64, &
         1149603840.0_real64, -598752000.0_real64, 228096000.0_real64, -59875200.0_real64, &
         9676800.0_real64,  - 725760.0_real64/
      data (b13(4,i),i=1,13)/ -725760.0_real64, 13063680.0_real64, -143700480.0_real64, &
         -476910720.0_real64, 1077753600.0_real64, -862202880.0_real64, 670602240.0_real64, &
         -431101440.0_real64, 215550720.0_real64, -79833600.0_real64, 20528640.0_real64, &
         -3265920.0_real64, 241920.0_real64/
      data (b13(5,i),i=1,13)/ 241920.0_real64, -3870720.0_real64, 31933440.0_real64, &
         -212889600.0_real64, -303937920.0_real64, 766402560.0_real64, -447068160.0_real64, &
         255467520.0_real64, -119750400.0_real64, 42577920.0_real64, -10644480.0_real64, &
         1658880.0_real64,  - 120960.0_real64/
      data (b13(6,i),i=1,13)/ -120960.0_real64, 1814400.0_real64, -13305600.0_real64, &
         66528000.0_real64, -299376000.0_real64, -148262400.0_real64, 558835200.0_real64, &
         -239500800.0_real64, 99792000.0_real64, -33264000.0_real64, 7983360.0_real64, &
         -1209600.0_real64, 86400.0_real64/
      data (b13(7,i),i=1,13)/ 86400.0_real64, -1244160.0_real64, 8553600.0_real64, &
         -38016000.0_real64, 128304000.0_real64, -410572800.0_real64, 0.0_real64, &
         410572800.0_real64, -128304000.0_real64, 38016000.0_real64, -8553600.0_real64, &
         1244160.0_real64,  - 86400.0_real64/
      data (b13(8,i),i=1,13)/ -86400.0_real64, 1209600.0_real64, -7983360.0_real64, &
         33264000.0_real64, -99792000.0_real64, 239500800.0_real64, -558835200.0_real64, &
         148262400.0_real64, 299376000.0_real64, -66528000.0_real64, 13305600.0_real64, &
         -1814400.0_real64, 120960.0_real64/
      data (b13(9,i),i=1,13)/ 120960.0_real64, -1658880.0_real64, 10644480.0_real64, &
         -42577920.0_real64, 119750400.0_real64, -255467520.0_real64, 447068160.0_real64, &
         -766402560.0_real64, 303937920.0_real64, 212889600.0_real64, -31933440.0_real64, &
         3870720.0_real64,  - 241920.0_real64/
      data (b13(10,i),i=1,13)/ -241920.0_real64, 3265920.0_real64, -20528640.0_real64, &
         79833600.0_real64, -215550720.0_real64, 431101440.0_real64, -670602240.0_real64, &
         862202880.0_real64, -1077753600.0_real64, 476910720.0_real64, 143700480.0_real64, &
         -13063680.0_real64, 725760.0_real64/
      data (b13(11,i),i=1,13)/ 725760.0_real64, -9676800.0_real64, 59875200.0_real64, &
         -228096000.0_real64, 598752000.0_real64, -1149603840.0_real64, 1676505600.0_real64, &
         -1916006400.0_real64, 1796256000.0_real64, -1596672000.0_real64, 684478080.0_real64, &
         87091200.0_real64,  - 3628800.0_real64/
      data (b13(12,i),i=1,13)/ -3628800.0_real64, 47900160.0_real64, -292723200.0_real64, &
         1097712000.0_real64, -2822688000.0_real64, 5269017600.0_real64, -7376624640.0_real64, &
         7903526400.0_real64, -6586272000.0_real64, 4390848000.0_real64, -2634508800.0_real64, &
         967524480.0_real64, 39916800.0_real64/
      data (b13(13,i),i=1,13)/ 39916800.0_real64, -522547200.0_real64, 3161410560.0_real64&
         , -11708928000.0_real64, 29638224000.0_real64, -54195609600.0_real64, &
         73766246400.0_real64, -75873853440.0_real64, 59276448000.0_real64, &
         -35126784000.0_real64, 15807052800.0_real64, -5748019200.0_real64, 1486442880.0_real64&
         /
!
      data b13den/ 479001600.0_real64/
!
!----------------------------------------------------------------------*
!                                                                      *
!   Coefficients for Sienkiewicz-Baylis formula                        *
!
      data cg/ 1771.0_real64, 9235.0_real64, 5890.0_real64, 4610.0_real64, 35.0_real64, 59.0_real64/
!
      data denom/ 5760.0_real64/
!
!----------------------------------------------------------------------*
!                                                                      *
!   Five-point Newton-Cotes coefficients for closed integration. Ex-   *
!   pressed as rational numbers                                        *
!
      data (c5num(i,2),i=1,5)/ 251.0_real64, 646.0_real64, -264.0_real64, 106.0_real64,  - &
         19.0_real64/
      data (c5num(i,3),i=1,5)/ 232.0_real64, 992.0_real64, 192.0_real64, 32.0_real64,  - 8.0_real64&
         /
      data (c5num(i,4),i=1,5)/ 243.0_real64, 918.0_real64, 648.0_real64, 378.0_real64,  - &
         27.0_real64/
      data (c5num(i,5),i=1,5)/ 224.0_real64, 1024.0_real64, 384.0_real64, 1024.0_real64, &
         224.0_real64/
!
      data c5den/ 720.0_real64/
!
!----------------------------------------------------------------------*
!                                                                      *
!   Six-point Newton-Cotes coefficients for closed integration.  Ex-   *
!   pressed as rational numbers                                        *
!
      data (c6num(i,2),i=1,6)/ 475.0_real64, 1427.0_real64, -798.0_real64, 482.0_real64, &
         -173.0_real64, 27.0_real64/
      data (c6num(i,3),i=1,6)/ 448.0_real64, 2064.0_real64, 224.0_real64, 224.0_real64, &
         -96.0_real64, 16.0_real64/
      data (c6num(i,4),i=1,6)/ 459.0_real64, 1971.0_real64, 1026.0_real64, 1026.0_real64, &
         -189.0_real64, 27.0_real64/
      data (c6num(i,5),i=1,6)/ 448.0_real64, 2048.0_real64, 768.0_real64, 2048.0_real64, &
         448.0_real64, 0.0_real64/
      data (c6num(i,6),i=1,6)/ 475.0_real64, 1875.0_real64, 1250.0_real64, 1250.0_real64, &
         1875.0_real64, 475.0_real64/
!
      data c6den/ 1440.0_real64/
!
!----------------------------------------------------------------------*
!
      data first/ .true./
!
!   Lagrange interpolation coefficients
!
!   Do this initialization once per run only
!   Do this initialization once per run only
!
      if (first) then
!
!   Thirteen-point coefficients for DPBDT
!
         factor = 1.0_real64/b13den
         do 2 j = 1,13
            do 1 i = 1,13
               a13(i,j) = b13(i,j)*factor
    1       continue
    2    continue
!
         first = .false.
!
      endif
!
!   Sienkiewicz-Baylis coefficients for SBSTEP
!
      c(1) = cg(1)/denom
      factor = transformed_grid_step/denom
      do 3 i = 2,6
         c(i) = cg(i)*factor
    3 continue
!
!   Newton-Cotes coefficients for YZK and QUAD
!
      factor = transformed_grid_step/c5den
      do 5 j = 2,4
         do 4 i = 2,5
            cnc5c(i,j) = factor*c5num(i,j)
    4    continue
    5 continue
!
      c1 = factor*c5num(1,5)
      c2 = factor*c5num(2,5)
      c3 = factor*c5num(3,5)
      c4 = c1+c1
!
!   Newton-Cotes coefficients for START
!
      factor = transformed_grid_step/c6den
      do 7 j = 2,6
         do 6 i = 1,6
            cnc6c(i,j) = factor*c6num(i,j)
    6    continue
    7 continue
!
!
      return
!
      end subroutine setqic
