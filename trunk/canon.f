! **************************************************************
!
! This file contains the subroutines: canon,can_weight 
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      subroutine  canon(nequi, nswp, nmes, temp, lrand)
! -----------------------------------------------------------------
! PURPOSE: CANONICAL SIMULATION OF PROTEINS USING METROPOLIS UPDATES
!
! CALLS:  addang,energy,metropolis,hbond,helix,outvar,outpdb,rgyr
!
!-------------------------------------------------------------------
      include 'INCL.H'

!f2py intent(in) nequi
!f2py intent(in) nswp
!f2py intent(in) nmes
!f2py intent(in) temp
!f2py logical optional, intent(in):: lrand = 1

!     external rand
      external can_weight
      
      double precision dv, grnd, vr, addang, eol, energy, acz, rgy, ee

      integer i, iv, nsw, nhel, mhel, nbet, mbet, mhb

     
      logical lrand
!      parameter(lrand=.false.)
!      parameter(nequi=10, nswp=1000,nmes=10)
!      parameter(temp=300.0)
!     lrand=.true.: creates random start configuration 
!     nequi: Number of sweeps for equilibrisation of system
      integer nequi
!     nswp:  Number of sweeps for simulation run
      integer nswp
!     nmes:  Number of sweeps between measurments
      integer nmes
!     temp:  Temperature of simulation
      double precision temp, e_min
!
!      common/bet/beta

      character*80 file

!     Define files for output:
      open(13,file='time.d')

 
      beta=1.0/ ( temp * 1.98773d-3 )

! _________________________________ random start
      if(lrand) then
       do i=1,nvr
        iv=idvr(i)  ! provides index of non-fixed variable
        dv=axvr(iv)*(grnd()-0.5)
        vr=addang(pi,dv)
        vlvr(iv)=vr
       enddo
      end if

      eol = energy()
      e_min = eol
      write (logString, '(a,e12.5,/)')  'energy of start configuration:',eol

! Write start configuration in pdb-format into file
      call outpdb(0,'start.pdb')

! =====================Equilibration by  Metropolis
      acz = 0.0d0
      do nsw=1,nequi
         call metropolis(eol,acz,can_weight)
      end do
      write (logString, *) 'Energy after equilibration:',eol

!======================Simulation in canonical ensemble
      acz = 0.0d0
      do nsw=0,nswp
        call metropolis(eol,acz,can_weight)
        if (eol.lt.e_min) then
           write (logString, *) "New minimum energy:", eol, "t=", nsw
           write (logString, *) eyel,eyhb,eyvr,eysl
           e_min = eol
           call outpdb(0, "minen.pdb")
        endif
!
        if(mod(nsw,nmes).eq.0) then
! Measure radius of gyration and end-to-end distance
! rgy: radius of gyration
! ee:  end-to-end distance
         call rgyr(1,rgy,ee)
! Measure helicity 
! nhel: number of helical residues
! mhel: number of helical segments
! nbet: number of sheet-like residues
! mbet: number of sheet-like segments
         call helix(nhel,mhel,nbet,mbet) 
! Measure number of hydrogen bonds (mhb)
        do i=1,ntlml
         call hbond(i,mhb,0)
        end do
! Write down information on actual conformation
         write(13,'(i7,2f15.3,5i7,4f15.3)')  nsw,  eol, rgy,
     &                              nhel,mhel,nbet,mbet,mhb,
     &                              eyel,eyhb,eyvr,eysl
        end if
!
      end do

      acz = acz/dble(nsw*nvr)
      write (logString, *) 'acceptance rate:',acz
      write (logString, *)
! ------------ Output Dihedreals of final configuration
      write (logString, *) 'last energy',eol
      call outvar(0,'lastconf.var')
!     Output final conformation as pdb-file
      call outpdb(0,'final.pdb')

      close(11)
      close(12)
      close(13)
! =====================


       end

! ********************************************************
      real*8 function can_weight(x)
!
! CALLS: none
!
      implicit none
      double precision beta, x


      common/bet/beta

      can_weight = beta*x

      return

      end
