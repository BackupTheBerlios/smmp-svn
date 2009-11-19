!**************************************************************
!
! This file contains the subroutines: nursvr, nursat
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************
      integer*4 function nursvr(ivr)

! ...........................................................
!  PURPOSE: defines index of residue for given variable 'ivr'
!
!  CALLS: none
!
! ...........................................................
      include 'INCL.H'

      integer i, ifirs, ivr, j

      do i=ntlml,1,-1
        ifirs=irsml1(i)
        if (ivr.ge.ivrrs1(ifirs).and.nvrml(i).gt.0) then
          do j=irsml2(i),ifirs,-1
            if (ivr.ge.ivrrs1(j).and.nvrrs(j).gt.0) then
              nursvr=j
              return
            endif
          enddo
        endif
      enddo

      write (logString, '(a,i5)') ' nursvr > Cannot find variable # '
     &   ,ivr
      stop

      end

! **********************************
      integer*4 function nursat(iat)

! .......................................................
!  PURPOSE: defines index of residue for given atom 'iat'
! .......................................................

      include 'INCL.H'

      integer i, ifirs, ilars, iat, j

      do i=1,ntlml

        ifirs=irsml1(i)
        ilars=irsml2(i)

        if (iat.ge.iatrs1(ifirs).and.iat.le.iatrs2(ilars)) then

          do j=ifirs,ilars

            if (iat.ge.iatrs1(j).and.iat.le.iatrs2(j)) then

              nursat=j

              return
            endif

          enddo

        endif
      enddo

      write (logString, '(a,i5)') ' nursat > Cannot find atom # ',iat
      stop

      end
