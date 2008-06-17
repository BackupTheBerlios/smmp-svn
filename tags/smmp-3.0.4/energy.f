c **************************************************************
c
c This file contains the subroutines: energy, enyinternal
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      real*8 function energy()
c ------------------------------------------
c PURPOSE: calculate the *total* energy of the system
c
c ientyp = 0 for ECEPP3, 1 for FLEX, 2 for Lund
c          3 for ECEPP3 with Abagyan corrections
c
c CALLS:  enyflx,enyreg,enyshe,enysol,esolan,setvar, eninteract
c
c -------------------------------------------------

      include 'INCL.H'
      double precision teysm, teyel, teyvw, teyhb, teyvr
c      print *,'energy function with ientyp  = ',ientyp   
      esm = 0.d0
      teysm = 0.d0
      teyel = 0.d0
      teyvw = 0.d0
      teyhb = 0.d0
      teyvr = 0.d0
      teysl = 0.d0

      do i = 1,ntlml  
         eysm=0
         eyel=0
         eyvr=0
         eyhb=0
         eyvw=0
         eysl=0

        call setvar(i,vlvr)  ! set variables & rebuild

        if (ientyp.eq.0.or.ientyp.eq.3) then
          esm=esm+enyshe(i)
        else if (ientyp.eq.1) then 
          esm=esm+enyshe(i)
        else if (ientyp.eq.2) then
           esm=enylun(i)
        endif

        teysm = teysm + eysm
        teyel = teyel + eyel
        teyhb = teyhb + eyhb
        teyvr = teyvr + eyvr
        if (ientyp.eq.2) then
c     The Lund term stores the hydrophobicity energy in eysl
           teysl = teysl + eysl
        else 
c     .. and the excluded volume term in eyvw, which is calculated once.
           teyvw = teyvw + eyvw
        endif

        if (ireg.eq.1)  eyrg=enyreg(i)        

      enddo

      if (ientyp.ne.2) then 
c     Don't touch eysl if using Lund potential, as enylun stores 
c     its hydrophobicity term there.
         if (itysol.gt.0) then
            esm=esm+enysol(0)
            teysl = teysl+eysl
!     elseif (itysol.lt.0) then
!     esm=esm+esolan(0)
         else
            eysl=0.d0
         endif
      else 
c     Add excluded volume term and save it in eyvw
         esm=esm+exvlun(0)
         teyvw = teyvw+eyvw
      endif

c The Abagyan entropic corrections depend on the area exposed to the 
c solvent for each residue. So, this term has to be evaluated after the
c solvent term.
      eyab=0.0
      if (ientyp.eq.3) then
         do i = 1,ntlml  
            eyab=eyab+eyabgn(i)
         enddo
      endif
      esm=esm+eyab
c Partial energies for the entire system. If you need the partial
c energies for a single molecule call enyinternal.
      eysm = teysm
      eyel = teyel
      eyvw = teyvw
      eyhb = teyhb
      eyvr = teyvr
      eysl = teysl
      
      if (ientyp.ne.2) then
c     This is temporary. eninteract() does not yet know how to calculate
c     interactions using the Lund potential.
         energy = esm + eninteract()
         return
      endif
      energy=esm
      return
      end

cc Calculates the internal energy for a single molecule.
c  All the partial energies are thus set to their values for molecule
c  nml. 
c
c  @param nml the ID of the molecule
c  @return internal energy of a single molecule
c
c  @author Jan H. Meinke <j.meinke@fz-juelich.de>
      real*8 function enyinternal(nml)

cf2py intent(in) nml
      
      include 'INCL.H'
      esm = 0.d0

      call setvar(nml,vlvr)  

      if (ientyp.eq.0.or.ientyp.eq.3) then
        esm=esm+enyshe(nml)
      else if (ientyp.eq.1) then
        esm=esm+enyflx(nml)
      else if (ientyp.eq.2) then
         esm=esm+enylun(nml)
      endif

      if (ireg.eq.1)  eyrg=enyreg(nml)

      if (ientyp.ne.2) then
         if (itysol.gt.0) then
            esm=esm+enysol(nml)
!     elseif (itysol.lt.0) then
!     esm=esm+esolan(nml)
         else
            eysl=0.d0
         endif
      else 
         esm=esm+exvlun(nml)
      endif
c The Abagyan entropic corrections depend on the area exposed to the 
c solvent for each residue. So, this term has to be evaluated after the
c solvent term.
      eyab=0.0
      if (ientyp.eq.3) then
         do i=1,ntlml
            eyab=eyab + eyabgn(i)
         enddo
      endif
      esm=esm+eyab

      enyinternal = esm
      return
      end
