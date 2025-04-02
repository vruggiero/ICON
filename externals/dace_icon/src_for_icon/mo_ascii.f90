!
!+ Constant definitions for non-printable ASCII characters
!
! $Id$
!
MODULE mo_ascii
!
! Description:
!   Constant definitions for non-printable ASCII characters
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_8         2009/12/09 Harald Anlauf
!  Add "implicit none"
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  add some comments
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD  1999-2007
!------------------------------------------------------------------------------
   implicit none
   private
   public :: NUL, BEL, HT
!------------------------------------------------------------------------------
!                          mnemonic  decimal ! octal !
   character, parameter :: NUL = achar ( 00) ! 000   ! NUL
!  character, parameter :: SOH = achar ( 01) ! 001   !
!  character, parameter :: STX = achar ( 02) ! 002   !
!  character, parameter :: ETX = achar ( 03) ! 003   !
!  character, parameter :: EOT = achar ( 04) ! 004   !
!  character, parameter :: ENQ = achar ( 05) ! 005   !
!  character, parameter :: ACK = achar ( 06) ! 006   !
   character, parameter :: BEL = achar ( 07) ! 007   ! Bell
!  character, parameter :: BS  = achar ( 08) ! 010   !
   character, parameter :: HT  = achar ( 09) ! 011   ! Horizontal TAB
!  character, parameter :: LF  = achar ( 10) ! 012   !
!  character, parameter :: VT  = achar ( 11) ! 013   !
!  character, parameter :: FF  = achar ( 12) ! 014   !
!  character, parameter :: CR  = achar ( 13) ! 015   !
!  character, parameter :: SO  = achar ( 14) ! 016   !
!  character, parameter :: SI  = achar ( 15) ! 017   !
!  character, parameter :: DLE = achar ( 16) ! 020   !
!  character, parameter :: DC1 = achar ( 17) ! 021   !
!  character, parameter :: DC2 = achar ( 18) ! 022   !
!  character, parameter :: DC3 = achar ( 19) ! 023   !
!  character, parameter :: DC4 = achar ( 20) ! 024   !
!  character, parameter :: NAK = achar ( 21) ! 025   !
!  character, parameter :: SYN = achar ( 22) ! 026   !
!  character, parameter :: ETB = achar ( 23) ! 027   !
!  character, parameter :: CAN = achar ( 24) ! 030   !
!  character, parameter :: EM  = achar ( 25) ! 031   !
!  character, parameter :: SUB = achar ( 26) ! 032   !
!  character, parameter :: ESC = achar ( 27) ! 033   !
!  character, parameter :: FS  = achar ( 28) ! 034   !
!  character, parameter :: GS  = achar ( 29) ! 035   !
!  character, parameter :: RS  = achar ( 30) ! 036   !
!  character, parameter :: US  = achar ( 31) ! 037   !
!  character, parameter :: DEL = achar (127) ! 177   !
end module mo_ascii
