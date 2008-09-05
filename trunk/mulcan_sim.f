! **************************************************************
!
! This file contains the subroutines: mulcan_sim,muca_weight2
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      subroutine  mulcan_sim
!
! PURPOSE: PERFORM A MULTICANONICAL SIMULATION
! REQUIRES AS INPUT THE MULTICANONICAL PARAMETER AS CALCULATED
! BY THE SUBROUTINE mulcan_par
!
! CALLS: addang, contacts,energy,metropolis
!
      include 'INCL.H'

!     external rand
      external muca_weight2

      logical restart

      parameter(restart=.false.)
      parameter(kmin=-12,kmax=20,ebin=1.0d0)
      parameter(nsweep=100000,nequi=100)
      Parameter(nsave=1000,nmes=10)
!
!     restart: .true. =  restart of simulation
!              .false. = start of simulation with random configuration
!     kmin,kmax: Range of multicanonical parameter
!     ebin:      bin size for multicanonical parameter
!     nequi: Number of sweeps for equilibrisation
!     nsweep:  Number of sweeps for simulation run
!     nsave:  Number of sweeps after which actual configuration is saved
!             for re-starts
!     nmes: Number of sweeps between measurments
!

      dimension xhist(kmin:kmax),ihist(kmin:kmax)
      common/muca2/b(kmin:kmax),alpha(kmin:kmax)

 
! FILE with last conformation (for re-starts)
      open(8,file='EXAMPLES/start.d')
! File with contact map of reference configuration
      open(9,file='EXAMPLES/enkefa.ref')
! File with multicanonical parameter
      open(10,file='EXAMPLES/muca.d')
! Result file: Time series of certain quantities
      open(11, file='EXAMPLES/time.d')


      do j=kmin,kmax
       ihist(j)= 0
      end do

      nresi=irsml2(1)-irsml1(1) + 1
!     nresi:  Number of residues

! READ  REFERENCE CONTACT MAP
      nci = 0
      do i=1,nresi
       read(9,*) (iref(i,j) , j=1,nresi)
      end do
      do i=1,nresi
       do j=nresi,i+3,-1
        if(iref(i,j).eq.1) nci = nci + 1
       end do
      end do
      write(*,*) 'Number of contacts in reference conformation:',nci

! READ IN FIELDS WITH MULTICANONICAL PARAMETER
      Do j=kmin,kmax
       read(10,*) i,b(i),alpha(i)
      end do
!

      if(restart) then
       read(8,*) nswm, eol_old
       read(8,*) (xhist(j), j=kmin,kmax)
       do i=1,nvr
        read(8,*) j, x
        iv = idvr(j)
        vlvr(iv) = x
       end do
       write(*,*) 'Last iteration, energy:',nswm,eol_old
      else
! _________________________________ random start
       do i=1,nvr
        iv=idvr(i)  ! provides index of non-fixed variable
        dv=axvr(iv)*(grnd()-0.5)
        vr=addang(pi,dv)
        vlvr(iv)=vr
       enddo
      end if
!
      eol = energy()
      write (*,'(e12.5,/)')  eol
      call contacts(nhy,nhx,dham)
      write(*,*) 'Number of contacts in start configuration:',nhy
      write(*,*) 'Number of native contacts in start configuration:',
     &            nhx
      do i=1,nresi
       write(*,'(62I1)') (ijcont(i,j), j=1,nresi)
      end do
      write(*,*)
!

      
      if(.not.restart) then
! =====================Equilibrization by  Metropolis
       do nsw=1,nequi
        call metropolis(eol,acz,muca_weight2)
       end do
       do i=kmin,kmax
        ihist(i) = 0
       end do 
       nswm = 1
      end if

!======================Simulation
      acz = 0.0d0
! LOOP OVER SWEEPS
      do nsw=nswm,nsweep
!
! METROPOLIS UPDATE
       call metropolis(eol,acz,muca_weight2)
       muold = min(kmax,max(kmin,int(eol/ebin+sign(0.5d0,eol))))
       ihist(muold) = ihist(muold) + 1
!
!  SAVE ACTUAL CONFORMATIONS FOR RE-STARTS:
       if(mod(nsw,nsave).eq.0) then
        rewind 8
        write(8,*) nswm, eol
        write(8,*) (xhist(j), j=kmin,kmax)
        do i=1,nvr
         iv = idvr(i)
         write(8,*) i,vlvr(iv)
        end do
       end if
! Measurements after NMES sweeps
       if(mod(nsw,nmes).eq.0) then
! Take a histogram of energy
        do i=kmin,kmax
         xhist(i) = xhist(i) + ihist(i)
         ihist(i) = 0
        end do
! Calculate contacts in actual configuartion and compare with reference
! configuration
!       call contacts(nhx,nhy,dham)
! nhx : Number of contcats in actual conformation
! nhy : Number of contacts which are identical in actual and reference 
!       configuration
! dham:  Hamming distance between actual and reference configuration      
!
        write(11,'(i7,f12.2,2i8,f12.4)')  nsw,eol,nhx,nhy,dham
       end if
      end do
! END OF SIMULATION

! FINAL OUTPUT:
      acz = acz/dble(nsw*nvr)
      write(*,*) 'last energy',eol
      write(*,*) 'aczeptance rate:',acz

! WRITE DOWN (UN-REWEIGHTED) HISTOGRAM OF MULTICANONICAL SIMULATION
      do i=kmin,kmax
       if(xhist(i).gt.0.0d0) then
        write(*,*) i,xhist(i)
       end if
      end do
! =====================
      close(8)
      close(9)
      close(10)
      close(11)

      return
      end

! ************************************************************
      real*8 function muca_weight2(x)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      Parameter(kmin=-12,kmax=20,ebin=1.0d0)
 
      common/muca2/b(kmin:kmax),alpha(kmin:kmax)

      muold = min(kmax,max(kmin,int(x/ebin+sign(0.5d0,x))))
      muca_weight2 = b(muold)*x + alpha(muold)

      return

      end


