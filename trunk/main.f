! **************************************************************
! This file contains the:  main (SINGLE PROCESSOR JOBS ONLY,
!                                FOR PARALLEL JOBS USE pmain)
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! CALLS: init_energy,init_molecule
! CALLS TASK SUBROUTINE: anneal,canon,elp,minim,mulcan_par,
!                        mulcan_sim,partem_s, or regul
! CAN ALSO CALL MEASUREMENT ROUTINES: cnteny,contacts,helix,hbond,
!                                    outpdb,outvar,rgyr,
!                                    rmsinit and rsmdfun,zimmer
! $Id: main.f 334 2007-08-07 09:23:59Z meinke $
! **************************************************************
      
      program main

      include 'INCL.H'
      include 'INCP.H'
      common/updstats/ncalls(5),nacalls(5)
      character*80 libdir, seqfile, varfile
      character grpn*4,grpc*4
      logical lrand,bgsposs

! =================================================== Energy setup

!            Directory for SMMP libraries
!     Change the following directory path to where you want to put SMMP
!     libraries of residues. 
      libdir='./SMMP/'

!      The switch in the following line is now not used.
      flex=.false.        ! .true. for Flex  / .false. for ECEPP

!     Choose energy type with the following switch instead ...
      ientyp = 0
!        0  => ECEPP2 or ECEPP3 depending on the value of sh2
!        1  => FLEX 
!        2  => Lund force field
!        3  => ECEPP with Abagyan corrections
!

      sh2=.false.         ! .true. for ECEPP/2; .false. for ECEPP3
      epsd=.false.        ! .true. for  distance-dependent  dielectric
                          !  permittivity

      itysol= 0    !  0: vacuum
                   ! >0: numerical solvent energy
                   ! <0: analytical solvent energy & gradients

      call init_energy(libdir)

! ================================================= Structure setup

      grpn = 'nh2' ! N-terminal group
      grpc = 'cooh'! C-terminal group

      iabin = 1  ! =0: read from PDB-file
                 ! =1: ab Initio from sequence (& variables)
      seqfile='EXAMPLES/enkefa.seq'
      varfile='EXAMPLES/enkefa.var'
!       varfile = ' '
      
      ntlml = 0
      write (*,*) 'Solvent: ', itysol
!     Initialize random number generator.
      call sgrnd(31433)
      
      if (itysol.eq.0.and.ientyp.eq.3) then
         print *,'Can not use Abagyan entropic corrections without '
         print *,'solvent term. '
         stop
      endif

      call init_molecule(iabin,grpn,grpc,seqfile,varfile)

! Decide if and when to use BGS, and initialize Lund data structures 
      bgsprob=0.75   ! Prob for BGS, given that it is possible
! upchswitch= 0 => No BGS 1 => BGS with probability bgsprob 
! 2 => temperature dependent choice 
      upchswitch=1
      rndord=.true.
      call init_lund
      if (ientyp.eq.2) call init_lundff
      if (ientyp.eq.3) call init_abgn
      

! ========================================  Add your task down here

      imin = 1 ! Quasi-Newton
      maxit = 15000 ! maximum number of iterations in minimization
      eps = 1.0d-7 ! requested precision
      call minim(imin, maxit, eps)
!     To do a canonical Monte Carlo simulation uncomment the lines below
!       nequi = 100
!       nsweep = 50000
!       nmes = 10
!       temp = 300.0
!       lrand = .true.
!      Canonical Monte Carlo
!       call canon(nequi, nsweep, nmes, temp, lrand)

!      For simulated annealing uncomment the lines below
!      tmin = 200.0
!      tmax = 500.0
!      call anneal(nequi, nsweep, nmes, tmax, tmin, lrand);
! ========================================  End of main      
       end
