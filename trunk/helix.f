! **************************************************************
!
! This file contains the subroutines: helix
!
! Copyright 2003       Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************


      subroutine helix(nhel,mhel,nbet,mbet) 
!---------------------------------------------------------------
!
!   PURPOSE: simple identification of secondary structure content
!
!   CALLS: none
!
! ---------------------------------------------------------------
      include 'INCL.H'

!f2py intent(out) nhel
!f2py intent(out) mhel
!f2py intent(out) nbet
!f2py intent(out) mbet
           
      logical lhel,lbet
      double precision hlim, philim, psilim, hlim2, philim2, psilim2
      double precision xphi, xpsi

      integer nhel, mhel, nbet, mbet, i, iv

      parameter(hlim=30.0d0,philim=-70.0d0,psilim=-37.0d0)
      parameter(hlim2=30.0d0,philim2=-150.0d0,psilim2=150.0d0)

      nhel = 0 ! Number of helical residues
      mhel = 0 ! Number of helical segments
      nbet = 0 ! Number of sheet-like residues
      mbet = 0 ! Number of sheet-like segments
      lhel = .false.
      lbet = .false.
      do i=1,nvr
       iv=idvr(i)
       if(nmvr(iv).eq.'phi') then
         xphi = vlvr(iv)*crd
         xpsi = vlvr(idvr(i+1))*crd
! Helicity
         if(abs(xphi-philim).le.hlim) then
          lbet=.false.
          if(abs(xpsi-psilim).le.hlim) then
            nhel = nhel+1
            if(.not.lhel) mhel = mhel + 1
            lhel = .true.
          else
            lhel = .false.
          end if
! Sheetness
         else if(abs(xphi-philim2).le.hlim2) then
          lhel = .false.
          if(abs(xpsi-psilim2).le.hlim2) then
            nbet = nbet + 1
            if(.not.lbet) mbet = mbet + 1
             lbet = .true.
           else
             lbet = .false.
           end if
         else
          lhel=.false.
          lbet=.false.
         end if
       end if
      end do

      return
      end

