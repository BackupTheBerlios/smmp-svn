!**************************************************************
!
! This file contains the subroutines: setvar
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      
      subroutine setvar(nml,vlvrx)

! ..............................................................
! PURPOSE: Reset variables in molecule 'nml' to new values given
!          in 'vlvrx' and rebuild molecule
!          
! ! assure constant PHASE angles for branches from same atom
!
! CALLS: bldmol,difang
! ....................................................

      include 'INCL.H'
Cf2py intent(in) nml, vlvrx
      integer nml
      double precision vlvrx
      dimension vlvrx(mxvr)
      
!     functions
      double precision difang

      integer i, iat,ity, iow, j, jat, i1vr
      double precision ba, to, tsh
      

      i1vr=ivrml1(nml)
      do i=i1vr,i1vr+nvrml(nml)-1  ! __________________ Variables of 'nml'

        iat=iatvr(i)
        ity=ityvr(i)

        if (ity.eq.3) then   ! torsion (assure phase=const.)

          tsh=difang(toat(iat),vlvrx(i))
          iow=iowat(iat)   ! (cannot be 1st atom of 'nml')

          do j=1,nbdat(iow)
            jat=ibdat(j,iow)
            if (iowat(jat).eq.iow) then ! excl. ring
              to=toat(jat)+tsh
              if (abs(to).gt.pi) to=to-sign(pi2,to)
              toat(jat)=to
              sntoat(jat)=sin(to)
              cstoat(jat)=cos(to)
            endif
          enddo

        elseif (ity.eq.2) then     ! valence angle
          ba=vlvrx(i)
          baat(iat)=ba
          snbaat(iat)=sin(ba)
          csbaat(iat)=cos(ba)
        elseif (ity.eq.1) then     ! valence length
          blat(iat)=vlvrx(i)
        endif

      enddo  ! ... Variables

      call bldmol(nml)

      return
      end

