! **************************************************************
!
! This file contains the subroutines:  anneal
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! $Id: anneal.f 334 2007-08-07 09:23:59Z meinke $
! **************************************************************

      subroutine  anneal(nequi, nswp, nmes, tmax, tmin, lrand)

! --------------------------------------------------------------
! PURPOSE: SIMULATED ANNEALING SEARCH OF LOWEST-POTENTIAL-ENERGY
!          CONFORMATIONS OF PROTEINS
!
! CALLS: addang,energy,metropolis,outvar,outpdb,rgyr,setvar,zimmer
!
! ---------------------------------------------------------------
 
      include 'INCL.H'

!f2py intent(in) nequi
!f2py intent(in) nswp
!f2py intent(in) nmes
!f2py intent(in) Tmax
!f2py intent(in) Tmin
!f2py logical optional, intent(in):: lrand = 1

!     external rand
      external can_weight
      
      double precision bmin, bmax, db, dv, grnd, vr, addang, eol, energy
      double precision acz, ymin, vlvrm, rgy, ee, temp

      integer nresi, i, iv, nsw, nemin, j
      
!      parameter(lrand=.true.)
!      parameter(nequi=100, nswp=100000,nmes=1000)
!      parameter(tmax=1000.0,tmin=100.0)
!     lrand=.true.: creates random start configuration 
      logical lrand
!     nequi: Number of sweeps for equilibrisation of system
      integer nequi
!     nswp:  Number of sweeps for simulation run
      integer nswp
!     nmes:  Number of sweeps between measurments
      integer nmes
!     tmax: Start temperature
      double precision tmax
!     tmin: Final temperature
      double precision tmin
     
      
!      common/bet/beta
!
      dimension vlvrm(mxvr)

     
 
!     Define files for output:
      open(14,file='time.d')
      write(14, *) '# $Id: anneal.f 334 2007-08-07 09:23:59Z meinke $'
      write(14, *) '# nsw, temp, eol, eysl, eyslh, eyslp, asa, rgy, ',
     &             '# rgyh, rgyp, eyhb, eyvw, eyel, eyvr, zimm'
      bmin=1.0/ ( tmax * 1.98773d-3 )
      bmax=1.0/ ( tmin * 1.98773d-3 )
      db = exp(log(bmax/bmin)/nswp)

!     nresi: Number of residues
! FIXME: Should loop over all proteins
      nresi=irsml2(ntlml)-irsml1(1)+1
! _________________________________ random start
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

! Write start configuration in pdb-format into file
        call outpdb(0, "start.pdb")

! =====================Equilibration by  Metropolis
      beta =  bmin
      do nsw=1,nequi
         call metropolis(eol,acz,can_weight)
      end do
      write(*,*) 'Energy after  equilibration:',eol

!======================Simulation by simulated annealing
      acz = 0.0d0
      ymin = eol
      do nsw=0,nswp
        beta = bmin*db**nsw
        call metropolis(eol,acz,can_weight)
! Store lowest-energy conformation
        if(eol.lt.ymin) then
         ymin = eol
         nemin = nsw
         call outvar(0,'global.var')
!     Output of lowest-energy conformation as pdb-file
         call outpdb(0,"global.pdb")
         do j=1,nvr
          iv=idvr(j)
          vlvrm(j) = vlvr(iv)
         end do
        end if
!
        if(mod(nsw,nmes).eq.0) then
! Measure radius of gyration and end-to-end distance
         call rgyr(1, rgy, ee)
! Determine Zimmerman code of actual conformation
         call zimmer(nresi)
! Write down information on actual conformation
         temp =  1.0d0/beta/0.00198773
         write(14,'(i6,13f12.3,1x,a)')  
     &   nsw, temp, eol, eysl, eyslh, eyslp, asa, 
     &   rgy, rgyh, rgyp,
     &   eyhb, eyvw, eyel, eyvr, zimm(1:nresi)
        end if
!
      end do

      acz = acz/dble(nsw*nvr)
      write(*,*) 'acceptance rate:',acz
      write(*,*)
! ------------ Output Dihedreals of final configuration
      write(*,*) 'last energy',eol
      call outvar(0,' ')
!     Output final conformation as pdb-file
      call outpdb(0,"final.pdb")
      write(*,*)

! ------------ Output Dihedreals of conformation with lowest energy
      write(*,*) 'lowest energy ever found:',nemin,ymin
      close(14)
! =====================


       end


