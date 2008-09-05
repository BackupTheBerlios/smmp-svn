! **************************************************************
!
! This file contains the subroutines: difang,addang
!
! Copyright 2003       Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      real*8 function difang(a1,a2)

! ......................................................
!  PURPOSE:  difang = a2 - a1  with:  -pi < difang <= pi
!            
!  INPUT:    a1,a2-two angles [rad.]
!
!  CALLS: none
!
! ......................................................

      implicit real*8 (a-h,o-z)

      parameter (pi=3.141592653589793d0,
     &           pi2=2.d0*pi)

      d=mod((a2-a1),pi2)
      if (abs(d).le.pi) then
        difang=d
      else
        difang=d-sign(pi2,d)
      endif

      return
      end
! *********************************
      real*8 function addang(a1,a2)

! ......................................................
!  PURPOSE:  addang = a1 + a2  with:  -pi < addang <= pi
!            
!  INPUT:    a1,a2-two angles [rad.]
!
!  CALLS: none
!
! ......................................................

      implicit real*8 (a-h,o-z)

      parameter (pi=3.141592653589793d0,
     &           pi2=2.d0*pi)

      d=mod((a1+a2),pi2)
      if (abs(d).le.pi) then
        addang=d
      else
        addang=d-sign(pi2,d)
      endif

      return
      end

