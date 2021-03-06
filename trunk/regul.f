!**************************************************************
!
! This file contains the subroutines: regul
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************


      subroutine regul(nml, iter, nsteps, acc)

! ----------------------------------------------------------
! PURPOSE: regularization of PDB-structure into SMMP geometry
!
!          @param nml molecule to be regularized
!          @param iter number of iterations during regularization
!          @param nsteps maximum number of steps in minimization
!          @param acc acceptance criterium for minimization
!
! CALLS:   minim, cnteny, outvar,rmsdopt
! ----------------------------------------------------------

      include 'INCL.H'
      include 'INCP.H'

!f2py intent(in) nml
!f2py intent(in) iter
!f2py intent(in) nsteps
!f2py intent(in) acc

      double precision acc, rm, av1, av2, rmsd, dn

      integer nsteps, nml, nrs, i, n, iter, it

      dimension rm(3,3),av1(3),av2(3)
      logical ishy(mxvr),fxvro(mxvr)



      wtrg = 1.d0
      wtey = 0.d0

      write (logString, '(/,a,2(a,f4.2),/)')
     &  ' ====================== Regularization only',
     &  '   Wt(energy) = ',wtey,'  Wt(regul.) = ',wtrg

      call minim(2, nsteps, acc)

      write (logString, *) ' '
      write (logString, *) ' -------- contacts after 1st regularization'
      write (logString, *) ' '
      call cnteny(nml)
      write (logString, *) ' '

      nrs = irsml2(nml)-irsml1(nml)+1
      call rmsdopt(nml,1,nrs,ixatp,xatp,yatp,zatp,0,rm,av1,av2,rmsd)

      write (logString, *) ' RMSD = ',rmsd

! --------------------------------------- fix vars. defined in PDB


      do i = ivrml1(nml),nvrml(nml)
        fxvro(i) = fxvr(i)  ! save
        if (isrfvr(i)) fxvr(i) = .true.  ! fix vars. defined in ref.str.
      enddo  ! vars.
      ireg = 0

      write (logString, '(/,a,2(a,f4.2),/)')
     &  ' ====================== Internal Energy for Hydrogens only',
     &  '   Wt(energy) = ',wtey,'  Wt(regul.) = ',wtrg

      call minim(1, nsteps, acc)

      write (logString, *) ' '
      write (logString, *) ' -------- contacs after Emin. for hydrogens'
      write (logString, *) ' '
      call cnteny(nml)

      do i = ivrml1(nml),nvrml(nml)
        fxvr(i) = fxvro(i)  ! restore
      enddo  ! vars.
      ireg = 1


      wtrg = 1.d0
      wtey = 0.d0

      n=iter
      dn=1.d0/dble(n)

      do it = 1,n

        wtrg = 1.d0 - dn*dble(it)
        wtey = 1.d0 - wtrg

        write (logString, '(/,a,i2,2(a,e11.3),/)')
     &    ' ================ Minimization #',it,
     &        '   Wt(energy) = ',wtey,'  Wt(regul.) = ',wtrg

        call minim(1, nsteps, acc)

        nrs = irsml2(nml)-irsml1(nml)+1
        call rmsdopt(nml,1,nrs,ixatp,xatp,yatp,zatp,0,rm,av1,av2,rmsd)

        write (logString, *) ' '
        write (logString, *) ' RMSD = ',rmsd

      enddo

      write (logString, *) ' '
      write (logString, *) ' ------- contacts after full regularization'
      write (logString, *) ' '
      call cnteny(nml)

!      call outpdb(nml,12)

! Output of dihedral angles of the regularized structure
      write (logString, *) 
     &   'Dihedral angles of the regularized structure;'
      call outvar(nml, 'regd.var')

      ireg = 0

      return
      end

