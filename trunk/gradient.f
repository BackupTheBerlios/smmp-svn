! **************************************************************
!
! This file contains the subroutines: gradient
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************


      subroutine gradient()

! -------------------------------------------
! PURPOSE: calculate energy & gradients
!
! CALLS:   opeflx,opereg,opeshe,opesol,setvar
! -------------------------------------------

      include 'INCL.H'


      double precision esm

      integer i, ivr1, ivr2, j

      esm = 0.d0

      do i = 1,ntlml  ! molecules

        call setvar(i,vlvr)  ! set variables & rebuild

        if (flex) then
          call opeflx(i)
        else
          call opeshe(i)
        endif

        esm = esm + eysm

        if (itysol.lt.0) then

          call opesol(i)
          esm = esm + eysl

          ivr1=ivrml1(i)
          ivr2=ivr1+nvrml(i)-1

          do j=ivr1,ivr2
            gdeyvr(j) = gdeyvr(j)+gdeysl(j)
          enddo

        else if (itysol.eq.0) then

          eysl = 0.d0
        else

          write(*,*)  'gradient> Set itysol < 0'
          stop
        endif

        if (ireg.eq.1)  call opereg(i)

      enddo

      eysm = esm

      return
      end
