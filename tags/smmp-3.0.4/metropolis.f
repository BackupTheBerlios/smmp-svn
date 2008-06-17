c **************************************************************
c
c This file contains the subroutines:  metropolis
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************


      subroutine metropolis(eol1,acz,dummy)
C
C SUBROUTINE FOR METROPOLIS UPDATE OF CONFIGURATIONS
C
C CALLS: energy,addang,grnd,dummy (function provided as argument)
C
      include 'INCL.H'
      include 'INCP.H'
      include 'incl_lund.h'
      external dummy
!      common/bet/beta
      common/updstats/ncalls(5),nacalls(5)
      integer updtch1, updtch2,bgs
      double precision vrol(mxvr), gbol
cf2py intent(in, out) eol
cf2py intent(in, out) acz
      eol = energy()
c     external rand
      do nsw=1,nvr
c Loop over dihedrals
c
         iupstate=0
         iupt=1
         if (rndord) then 
            ivar = 1+int(nvr*grnd())
         else
            ivar=nsw
         endif
         jv = idvr(ivar)
         if (nmvr(jv)(1:1).eq.'x') then 
            iupt=1
         else 
            if (upchswitch.eq.0) then
               iupt=1
            else if (upchswitch.eq.1) then
               iupt=updtch1(ivar,beta)
            else 
               iupt=updtch2(ivar,beta)
            endif
         endif
         if (iupt.eq.2) then
            iupstate=bgs(eol,dummy)
         else
c           Simple twist of 
c           Get Proposal configuration
            vrol=vlvr!(jv)  
            dv=axvr(jv)*(grnd()-0.5)
            vlvr(jv)=addang(vrol(jv),dv)
c
c           Get dummy of proposal configuration
c
            enw = energy()
c
            delta =  dummy(enw) - dummy(eol)
c           ___________________________ check acceptance criteria
            if (delta.le.0.0d0) then
               eol=enw  
               iupstate=1
            else 
               rd=grnd()
               delta = min(delta,100.0d0)
               ex=exp( -delta )
               if ( ex .ge. rd ) then
                  eol=enw
                  iupstate=1
               else
                  vlvr=vrol
               endif
            endif
         endif
         enw1=energy()
!          if (enw1.ne.eol) then
!             write (*,*) 'metropolis: variable, enw1,eol,acptd',
!      #           ivar,enw1,eol,acptd,(enw1-eol),enw
!          endif
         acz=acz+iupstate*1.0
         call accanalyze(iupt,iupstate)
      end do

c Updates on relative position of different molecules when there are many      
      if (ntlml.gt.1) then
         do iml=1, ntlml
            do i=1,3
               iupstate=0
               gbol = gbpr(i, iml)
               gbpr(i, iml) = gbpr(i, iml) + (grnd()-0.5)
               if (boxsize.gt.0) then
                  if (gbpr(i, iml).gt.boxsize) then 
                     gbpr(i, iml)=boxsize
                  elseif (gbpr(i, iml).lt.-boxsize) then 
                     gbpr(i, iml)=-boxsize
                  endif
               endif
c     
c              Get dummy of proposal configuration
c     
               enw = energy()
c
               delta =  dummy(enw) - dummy(eol)
c     
c              ____________________________ check acceptance criteria
c     
               if (delta.le.0.0d0) then
                  eol=enw
                  iupstate=1
               else
                  rd=grnd()
                  delta = min(delta,100.0d0)
                  ex=exp( -delta )
                  if ( ex .ge. rd ) then
                     eol=enw
                     iupstate=1
                  else
                     gbpr(i, iml)=gbol
                  endif
               endif
               acz=acz+iupstate*1.0
               call accanalyze(3,iupstate)
            enddo
            do i=4,6
               iupstate=0
               gbol = gbpr(i, iml)
               if (i.eq.5) then ! global psi value must be in [-pi/2, pi/2]
                  gbpr(i, iml) = (grnd()-0.5) * pi
               else               
                  gbpr(i, iml) = (grnd()-0.5) * pi2
               endif
c     
c              Get dummy of proposal configuration
c     
               enw = energy()
c     
               delta =  dummy(enw) - dummy(eol)
c     
c              ____________________________ check acceptance criteria
c     
               if (delta.le.0.0d0) then
                  eol=enw
                  iupstate=1
               else
                  rd=grnd()
                  delta = min(delta,100.0d0)
                  ex=exp( -delta )
                  if ( ex .ge. rd ) then
                     eol=enw
                     iupstate=1
                  else
                     gbpr(i, iml)=gbol
                  endif
               endif
               acz=acz+iupstate*1.0
               call accanalyze(4,iupstate)
            enddo
         enddo
      endif
c     
c     Re-calculate energy
c     
      enw = energy()
      if(abs(eol-enw).gt.0.000001)  then
         write(*,*) 'metropolis: eol,enw,difference:', eol, enw, eol-enw
         if (eol.lt.100000) then 
            stop 'metropolis: eol and enw difference unacceptable'
         endif
      endif
c     
      eol1 = eol
      return
c     
      end
      
      subroutine accanalyze(iuptype,iupdstate)
      common/updstats/ncalls(5),nacalls(5)
      ncalls(5)=ncalls(5)+1
      nacalls(5)=nacalls(5)+iupdstate
      ncalls(iuptype)=ncalls(iuptype)+1
      nacalls(iuptype)=nacalls(iuptype)+iupdstate
      end      
      

      integer function updtch1(iiii,bbbb)
      logical rndord
      integer upchswitch, iiii
      double precision bgsprob, grnd, bbbb
      common /updchois/rndord,upchswitch,bgsprob
      if (grnd().lt.bgsprob) then
         updtch1=2
      else 
         updtch1=1
      endif
      return
      end

      integer function updtch2(iiii,bbbb)
      integer iiii
      double precision bbbb,curprob, grnd,up2bmax,up2bmin
      common/updtparam/up2bmax,up2bmin
      curprob=(bbbb-up2bmin)/(up2bmax-up2bmin)
      if (grnd().lt.curprob) then
         updtch2=2
      else
         updtch2=1
      endif
      end
      block data updtchs
      double precision up2bmax, up2bmin
      common/updtparam/up2bmax,up2bmin
      data up2bmax,up2bmin/4.0d0,0.5d0/
      common/updstats/ncalls(5),nacalls(5)
      data ncalls/5*0/
      data nacalls/5*0/
      end
