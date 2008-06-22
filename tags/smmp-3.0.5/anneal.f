c **************************************************************
c
c This file contains the subroutines:  anneal
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c $Id: anneal.f 334 2007-08-07 09:23:59Z meinke $
c **************************************************************

      subroutine  anneal(nequi, nswp, nmes, tmax, tmin, lrand)

C --------------------------------------------------------------
C PURPOSE: SIMULATED ANNEALING SEARCH OF LOWEST-POTENTIAL-ENERGY
C          CONFORMATIONS OF PROTEINS
C
C CALLS: addang,energy,metropolis,outvar,outpdb,rgyr,setvar,zimmer
C
C ---------------------------------------------------------------
 
      include 'INCL.H'

cf2py intent(in) nequi
cf2py intent(in) nswp
cf2py intent(in) nmes
cf2py intent(in) Tmax
cf2py intent(in) Tmin
cf2py logical optional, intent(in):: lrand = 1

c     external rand
      external can_weight
c      parameter(lrand=.true.)
c      parameter(nequi=100, nswp=100000,nmes=1000)
c      parameter(tmax=1000.0,tmin=100.0)
C     lrand=.true.: creates random start configuration 
      logical lrand
C     nequi: Number of sweeps for equilibrisation of system
      integer nequi
C     nswp:  Number of sweeps for simulation run
      integer nswp
c     nmes:  Number of sweeps between measurments
      integer nmes
C     tmax: Start temperature
      double precision tmax
C     tmin: Final temperature
      double precision tmin
     
      
!      common/bet/beta
C
      dimension vlvrm(mxvr)

     
 
c     Define files for output:
      open(14,file='time.d')
      write(14, *) '# $Id: anneal.f 334 2007-08-07 09:23:59Z meinke $'
      write(14, *) '# nsw, temp, eol, eysl, eyslh, eyslp, asa, rgy, ',
     &             '# rgyh, rgyp, eyhb, eyvw, eyel, eyvr, zimm'
      bmin=1.0/ ( tmax * 1.98773d-3 )
      bmax=1.0/ ( tmin * 1.98773d-3 )
      db = exp(log(bmax/bmin)/nswp)

c     nresi: Number of residues
c FIXME: Should loop over all proteins
      nresi=irsml2(ntlml)-irsml1(1)+1
c _________________________________ random start
      if(lrand) then
       do i=1,nvr
        iv=idvr(i)  
        dv=axvr(iv)*(grnd()-0.5)
        vr=addang(pi,dv)
        vlvr(iv)=vr
       enddo
      end if

      eol=energy()
      write (*,'(a,e12.5,/)')  'energy of start configuration: ',eol

C Write start configuration in pdb-format into file
        call outpdb(0, "start.pdb")

c =====================Equilibration by  Metropolis
      beta =  bmin
      do nsw=1,nequi
         call metropolis(eol,acz,can_weight)
      end do
      write(*,*) 'Energy after  equilibration:',eol

C======================Simulation by simulated annealing
      acz = 0.0d0
      ymin = eol
      do nsw=0,nswp
        beta = bmin*db**nsw
        call metropolis(eol,acz,can_weight)
c Store lowest-energy conformation
        if(eol.lt.ymin) then
         ymin = eol
         nemin = nsw
         call outvar(0,'global.var')
C     Output of lowest-energy conformation as pdb-file
         call outpdb(0,"global.pdb")
         do j=1,nvr
          iv=idvr(j)
          vlvrm(j) = vlvr(iv)
         end do
        end if
c
        if(mod(nsw,nmes).eq.0) then
C Measure radius of gyration and end-to-end distance
         call rgyr(1, rgy, ee)
C Determine Zimmerman code of actual conformation
         call zimmer(nresi)
C Write down information on actual conformation
         temp =  1.0d0/beta/0.00198773
         write(14,'(i6,13f12.3,1x,a)')  
     &   nsw, temp, eol, eysl, eyslh, eyslp, asa, 
     &   rgy, rgyh, rgyp,
     &   eyhb, eyvw, eyel, eyvr, zimm(1:nresi)
        end if
C
      end do

      acz = acz/dble(nsw*nvr)
      write(*,*) 'acceptance rate:',acz
      write(*,*)
c ------------ Output Dihedreals of final configuration
      write(*,*) 'last energy',eol
      call outvar(0,' ')
C     Output final conformation as pdb-file
      call outpdb(0,"final.pdb")
      write(*,*)

c ------------ Output Dihedreals of conformation with lowest energy
      write(*,*) 'lowest energy ever found:',nemin,ymin
      close(14)
c =====================


       end


