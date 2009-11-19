! **************************************************************
!
! This file contains the subroutines: minim,move
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************


      subroutine minim(imin, maxit, eps)

! ......................................................................
! PURPOSE: Use minimizers
!
!          imin = 1:  use Quasi-Newton
!          imin = 2:  use Conjugated Gradients
!
!          @param maxit maximum number of iterations 
!          @param eps acceptance criterium
!Institute
! CALLS: difang,energy,gradient, mincjg,minqsn, nursvr
! ......................................................................

      include 'INCL.H'
!f2py intent(in) imin
!f2py intent(in) maxit
!f2py intent(in) eps
      integer msvmx

      integer imin

      integer maxit
      
      integer nursvr
      double precision energy, difang
            
      integer ngbvr, msv, i, j, mxop, n6, n3, n, n1, n2, n5, n4
      integer ntlvr, nop, nv
      
      double precision eps, acc, esm, gdey2, gbpro, gdvr, gdrg2, scl
      double precision vlvrn, w, vlvro, vr      
      
      parameter (msvmx=mxvr*(mxvr+5)/(2*(2*mxvr+1)),   
     &            msv  = 50 )                             

      dimension w(mxvr*(mxvr+13)/2)

      dimension vlvrn(mxvr),vlvro(mxvr),gdvr(mxvr),scl(mxvr)

! --------------------------- new
      dimension gbpro(6,mxml)

      mxop=maxit
      acc=eps   

      if (msv.gt.msvmx) then
        write (logString, '(a,i5)') ' minreg> parameter MSV > ',msvmx
        stop
      endif

! ----------------------- energy & gradient

      call gradient()

      write (logString, '(/,a,/)') ' Energy BEFORE minimization:' 

      if (ireg.eq.0) then

        write (logString, '(a,e12.5,/,3(a,e11.4),/,2(a,e11.4),/)') 
     &    ' Total: ', eysm,
     &    '   Coulomb: ',eyel,' Lennard-Jones: ',eyvw,' HB: ',eyhb,
     &    '   Variables: ',eyvr,'  Solvatation: ',eysl

       else


        write (logString, '(a,e12.5,/,3(a,e11.4),/,3(a,e11.4),/)') 
     &    ' Total: ', wtey*eysm + wtrg*eyrg,
     &    '   Coulomb: ',eyel,' Lennard-Jones: ',eyvw,' HB: ',eyhb,
     &    '   Variables: ',eyvr,'  Solvatation: ',eysl,
     &    ' Regularization: ',eyrg

       endif

! --------------------------------------- variables

      ntlvr=ivrml1(ntlml)+nvrml(ntlml)-1
      nv=0

      gdey2=0.d0
      gdrg2=0.d0

      do i=1,ntlvr
        if (.not.fxvr(i)) then

          nv=nv+1

          vlvrn(nv) = vlvr(i)

          scl(nv) = 1.d0 ! scale grad_i

          if (ireg.eq.0) then
            gdvr(nv) = gdeyvr(i)
          else
            gdvr(nv) = wtey*gdeyvr(i)+wtrg*gdeyrg(i)
            gdrg2 = gdrg2+gdeyrg(i)**2
          endif

          gdey2=gdey2+gdeyvr(i)**2

        endif

        vlvro(i)=vlvr(i)
      enddo

      n=nv

      if (ireg.eq.0) then

        esm=eysm

      else

        ngbvr=0

        do i=1,ntlml
          do j=1,6

            n=n+1
            ngbvr=ngbvr+1

            vlvrn(n) = gbpr(j,i)
            scl(n) = 1.d0 ! Scale

            gbpro(j,i) = gbpr(j,i)   ! save

            gdvr(n) = wtrg*gdeygb(ngbvr)

            gdrg2=gdrg2+gdeygb(ngbvr)**2
          enddo
        enddo

!        if (abs(wtrg-1.d0).gt.1.d-4.and.abs(wtey-1.d0).gt.1.d-4) then
!          gdey2 = max(acc,gdey2)
!          gdrg2 = max(acc,gdrg2)
!          wtrg = wtrg * sqrt(gdey2/gdrg2)
!          write (logString, *)  ' -->    Wt_energy = ',wtey,'  Wt_regul. = ',wtrg
!          write (logString, *)  '  '
!        endif

        esm=wtey*eysm+wtrg*eyrg

      endif


      if (imin.eq.1) then ! Quasi-Newton

        n1=1+(n*(n+1))/2
        n2=n1+n
        n3=n2+n
        n4=n3+n
        n5=n4+n
        n6=n5+n

        call minqsn(n,mxvr,vlvrn,esm,gdvr,scl,acc,w,w(n1),w(n2),
     &              w(n3),w(n4),w(n5),w(n6),mxop,nop)

      elseif (imin.eq.2) then ! Conjugated Gradients

        n1=1+n
        n2=n1+n
        n3=n2+n
        n4=n3+n
        n5=n4+n

        call mincjg(n,mxvr,vlvrn,esm,gdvr,acc,w,w(n1),w(n2),  ! no 'scl'
     &              w(n3),w(n4),w(n5),mxop,nop)

      endif


      if (nop.lt.mxop) then
        write (logString, '(a)') ' ---- CONVERGENCE ----'
      else
        write (logString, '(a)')  '---- STEP LIMIT ----'
      endif


      write (logString, '(/,2a,/)') ' Final energies ',
     & '__________________________________________________'

      eysm = energy()

      if (ireg.eq.0) then

        write (logString, '(a,e12.5,/,3(a,e11.4),/,2(a,e11.4))') 
     &  ' Total: ',eysm,'   Coulomb: ',eyel,' Lennard-Jones: ',eyvw,
     &  ' HB: ',eyhb, '   Variables: ',eyvr,'  Solvatation: ',eysl

      else

        write (logString, '(a,e12.5,/,3(a,e11.4),/,3(a,e11.4))') 
     &  ' Total: ', wtey*eysm + wtrg*eyrg,
     &  '   Coulomb: ',eyel,' Lennard-Jones: ',eyvw,' HB: ',eyhb,
     &  '   Variables: ',eyvr,'  Solvatation: ',eysl,
     &  ' Regularization: ',eyrg

      endif

      write (logString, '(/,a,/)') ' Variables _________________'

      nv = 0
      do i=1,ntlvr  !! nvr
        if (.not.fxvr(i)) then

          nv=nv+1
          vr=vlvrn(nv)
          if (abs(vr).gt.pi) vr=vr-sign(pi2,vr)

          write (logString, '(1x,a,1x,i4,f8.1,a,f5.1,a)') nmvr(i),
     &        nursvr(i), vr*crd,'    (',abs(difang(vr,vlvro(i)))*crd,')'

          vlvr(i) = vr
        endif
      enddo


      if (ireg.ne.0) then

        write (logString,'(/,a,/)') ' Global Variables ___________'

        do i=1,ntlml
          write (logString, *) ' Molecule #',i,' old          new'
          do j=1,3
            write (logString, *) gbpro(j,i),' ',gbpr(j,i)
          enddo
          do j=4,6
            write (logString, *) gbpro(j,i)*crd,' ',gbpr(j,i)*crd
          enddo
        enddo

      endif

      write (logString,'(/,2a)') ' Gradient ',
     & '______________________________________________________________'

      write (logString,'(8(1x,f8.3))') (gdvr(i),i=1,nv)

      if (ireg.ne.0) then

        write (logString, *) 
     &             ' -------------- global variables ------------'
        write (logString, '(6(1x,f8.3))') (gdvr(i+nv),i=1,ngbvr)

      endif

      return
      end
! ********************************************
      subroutine move(nop,nvr1,esm,vlvrn,gdvr)
!
! CALLS: gradient
!
      include 'INCL.H'
      integer ngbvr

      integer nop

      integer nvr1

      double precision esm

      double precision vlvrn

      double precision gdvr, gdsmey, gdsmrg

      integer i, ii, j, ntlvr, n
      
      dimension vlvrn(mxvr),gdvr(mxvr)


! ------------------------ compile & new variables

      ntlvr=ivrml1(ntlml)+nvrml(ntlml)-1
      n=0

      do i=1,ntlvr
        if (.not.fxvr(i)) then
          n=n+1
          vlvr(i)=vlvrn(n)
        endif
      enddo

      if (ireg.ne.0) then

        ii=0
        do i=1,ntlml
          do j=1,6   ! global vars.
            ii=ii+1
            gbpr(j,i)=vlvrn(ii+n)
          enddo
        enddo

      endif

! -------------------------- new minimz. gradient

      call gradient()

      gdsmey=0.d0
      gdsmrg=0.d0

      n=0

      do i=1,ntlvr
        if (.not.fxvr(i)) then
          n=n+1

          if (ireg.eq.0) then
            gdvr(n) = gdeyvr(i)
          else
            gdvr(n) = wtey*gdeyvr(i) + wtrg*gdeyrg(i)
            gdsmrg = gdsmrg + gdeyrg(i)**2
          endif

          gdsmey = gdsmey + gdeyvr(i)**2

        endif
      enddo


      if (ireg.eq.0) then

        esm=eysm

        write (logString, '(a,i5,a,2(e13.6,a))') ' Step ',nop,
     &                    ': energy ', esm, '  (',gdsmey,' )'

      else

        ii=0
        do i=1,ntlml  ! global vars.
          do j=1,6
            n=n+1
            ii=ii+1

            gdvr(n) = wtrg*gdeygb(ii)
            gdsmrg = gdsmrg+gdeygb(ii)**2

          enddo
        enddo

        esm=wtey*eysm+wtrg*eyrg

        write (logString, '(a,i5,a,3(e13.6,a))') ' Step ',nop,
     &             ': energy ',esm,'  (',gdsmey,',',gdsmrg,' )'

      endif

      return
      end

