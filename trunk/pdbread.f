! **************************************************************
!
! This file contains the subroutines: pdbread,pdbvars,atixpdb,getpar
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      subroutine pdbread(pdbfil,ier)

! ....................................................
! PURPOSE: read protein atom coordinates from 'pdbfil'
!          (no Hydrogens, only ATOM records)
!
! RETURNS: 0 = no errors / 1 = error
!
! CALLS: iopfil,iendst
! ......................................................

      include 'INCP.H'
      
      double precision cor

      integer ier, l, iendst, lunpdb, io, iopfil, iat, i

! -------------------------- input
      character*(*) pdbfil
! -------------------------- local
      dimension cor(3)
      character atm*4,rsn*3,rsno*3,chn,chno,
     &          rsid*5,rsido*5,line*132, logString*255

      natp=0
      nchp=0
      nrsp=0

      ier=1

      chno='&'
      rsno='#&#'
      rsido='#&#&#'

      l=iendst(pdbfil)
      if (l.gt.0) then
        lunpdb = 99
      else
        write (logString, '(a)')
     &    ' pdbread> empty file name to read pdb-structure'

        return
      endif

      io=iopfil(lunpdb,pdbfil,'old','formatted')

      if (io.le.0) then
        write (logString, '(a,/,a)')
     &  ' pdbread> ERROR opening file to read pdb-structure: ',
     &  pdbfil(1:iendst(pdbfil))

        return
      endif

    1 read (lunpdb,'(a)',end=3,err=2) line
      l=iendst(line)

      if (l.lt.54.or.index(line(1:4),'ATOM').le.0) goto 1

      if ( line(17:17).ne.' ' )  then
        write (logString, '(a,/,a,/,a,/,2a)')
     &  ' pdbread> found alternate atom location: ',
     &  '                !',
     &  line(:l),' in file: ',pdbfil(1:iendst(pdbfil))

        close(lunpdb)
        return
      endif

      atm=line(13:16)

      if (index(atm(2:2),'H').gt.0) goto 1   ! no H

      read(line,10,err=2) iat,rsn,chn,rsid,(cor(i),i=1,3)

      if ((natp+1).gt.MXATP) then
        write (logString, '(a,i5,a,/,a)')
     &  ' pdbread>  >MXATP (',MXATP,') ATOM lines in PDB file ',
     &  pdbfil(1:iendst(pdbfil))

        close(lunpdb)
        return
      endif

      if (chn.ne.chno) then      ! new chain

        if ((nchp+1).gt.MXCHP) then
          write (logString, '(a,i3,a,/,a)')
     &    ' pdbread>  >MXCHP (',MXCHP,') chains in PDB file ',
     &    pdbfil(1:iendst(pdbfil))

          close(lunpdb)
          return
        endif

        if ((nrsp+1).gt.MXRSP) then
          write (logString, '(a,i3,a,/,a)')
     &    ' pdbread>  >MXRSP (',MXRSP,') residues in PDB file ',
     &    pdbfil(1:iendst(pdbfil))

          close(lunpdb)
          return
        endif

        if (nchp.eq.1) then
          nchrsp(nchp)=nrsp
        elseif (nchp.gt.1) then
          nchrsp(nchp)=nrsp-nchrsp(nchp)
        endif

        nchp=nchp+1
        chno=chn
        chnp(nchp)=chn

        nchrsp(nchp)=nrsp ! -1 1st res.

        if (nrsp.ge.1) then
          nrsatp(nrsp)=natp-irsatp(nrsp)+1
        endif

        rsido=rsid
        rsno=rsn

        nrsp=nrsp+1
        irsatp(nrsp)=natp+1
        rsnmp(nrsp)=rsn
        rsidp(nrsp)=rsid

      elseif (rsid.ne.rsido.or.rsn.ne.rsno) then      ! new residue

        if ((nrsp+1).gt.MXRSP) then
          write (logString, '(a,i3,a,/,a)')
     &    ' pdbread>  >MXRSP (',MXRSP,') residues in PDB file ',
     &    pdbfil(1:iendst(pdbfil))

          close(lunpdb)
          return
        endif

        nrsatp(nrsp)=natp-irsatp(nrsp)+1

        rsido=rsid
        rsno=rsn

        nrsp=nrsp+1
        irsatp(nrsp)=natp+1
        rsnmp(nrsp)=rsn
        rsidp(nrsp)=rsid

      endif

      natp=natp+1

      noatp(natp)=iat
      atnmp(natp)=atm

      xatp(natp)=cor(1)
      yatp(natp)=cor(2)
      zatp(natp)=cor(3)

      goto 1

    2 write (logString, '(a,/,a,/,2a)')
     &  ' pdbread> ERROR reading ATOM line ',
     &  line(:l),
     &  ' from file ',pdbfil(1:iendst(pdbfil))

      close(lunpdb)
      return

    3 close(lunpdb)

      if (natp.gt.0) then

        if (nchp.eq.1) then
          nchrsp(nchp)=nrsp
        elseif (nchp.gt.1) then
          nchrsp(nchp)=nrsp-nchrsp(nchp)
        endif

        nrsatp(nrsp)=natp-irsatp(nrsp)+1
        ier=0

      else

        write (logString, '(a,/,a)')
     &  ' pdbread> NO atom coordinates selected from file ',
     &  pdbfil(1:iendst(pdbfil))

      endif

      return

   10 format(6x,i5,6x,a3,1x,a1,a5,3x,3d8.3)

      end
! **************************************************************

      subroutine pdbvars()

! --------------------------------------------------------------------
! PURPOSE: sequence,indices for selected atoms (data in INCP.H)
!          & torsions from PDB to be used to build SMMP structure
!
!          ixatp(i,)
!          = indices for SMMP atoms pointing to PDB atoms
!            (=0, if atom not selected)
!
! --------------------------------- ref. point & axes
!         ixrfpt(3,),rfpt(3,),xrfax(3,),yrfax(3,),zrfax(3,)
!
! CALLS:  tolost,getmol,bldmol,addend,atixpdb,setmvs,mklist,
!         dihedr,fnd3ba,setsys,getpar,setvar,rmsdopt
! --------------------------------------------------------------------

      include 'INCL.H'
      include 'INCP.H'

      double precision dihedr, rm, av1, av2, rmsd, h, x, y, z

      integer nml, nrs, nc, irb, ire, irs, i, ii, it, i1, i2, i3, i4, j1
      integer j2, j3, j4, inew, j, ix

      character res*3
      dimension rm(3,3),av1(3),av2(3),h(3)


      nml=0
      nrs=0

      do nc=1,nchp  ! PDB chains

! =============================== SMMP molecule
        nml=nml+1
        if (nml.gt.mxml) then
          write (logString, '(a,i4,2a)')' pdbvars> NUMBER of chains > '
     &                          ,mxml,' in ',' ?'
          stop
        endif
        ntlml=nml
! ----------------------------- 'nmml' = ChainID
        nmml(nml)=chnp(nc)

! ======================================== get sequence

        irb=nrs+1
        ire=nrs+nchrsp(nc)
! ----------------------------- # of 1st & last residue
        irsml1(nml)=irb
        irsml2(nml)=ire

        do irs=irb,ire  ! residues of chain 'nc'

          nrs=nrs+1

          if (nrs.gt.mxrs) then
            write (logString, '(a,i4,2a)') 
     &         ' pdbvars> NUMBER of residues > ', mxrs, ' in ',' ?'
            stop
          endif

          res=rsnmp(irs)
          call tolost(res)

          seq(nrs)=res

          if (.not.flex.and.irs.eq.irb.and.seq(nrs)(1:3).eq.'pro')
     &      seq(nrs)='pron'  ! only ECEPP/3

        enddo ! residues

! ======================== get initial coords. for molecule 'nml'
!                          with library values for deg. of freedom

        call getmol(nml)   ! assemble res. data from libraries

        do i=1,6  ! initialize global parameters
          gbpr(i,nml)=zero
        enddo

        call bldmol(nml)                 ! co-ordinates
        call addend(nml,'nh2 ','cooh')   ! modify ends

        call atixpdb(nml)  ! get 'ixatp'

! -------------------------- 'load' SMMP variable information
        call setmvs(nml)   ! moving sets
        call mklist(nml)   ! interaction lists

! ================================= get variables for 'nml'

        ii=ivrml1(nml)
        do i=ii,ii+nvrml(nml)-1       ! SMMP torsions

          isrfvr(i) = .false.
          fxvr(i) = .false.
          idvr(i) = i

          it = ityvr(i)
          i1 = iatvr(i)

          if (it.eq.3) then         ! torsion

            i2=iowat(i1) ! indices of SMMP atoms
            i3=iowat(i2)
            i4=iowat(i3)

            j1=ixatp(i1)  ! inds. for corresp. PDB atoms
            j2=ixatp(i2)
            j3=ixatp(i3)
            j4=ixatp(i4)

            if (j1.le.0.or.j2.le.0.or.j3.le.0.or.j4.le.0) then

              vlvr(i) = toat(i1)  ! default value from library

            else  ! get value from PDB

              xat(i1)=xatp(j1)
              yat(i1)=yatp(j1)
              zat(i1)=zatp(j1)

              xat(i2)=xatp(j2)
              yat(i2)=yatp(j2)
              zat(i2)=zatp(j2)

              xat(i3)=xatp(j3)
              yat(i3)=yatp(j3)
              zat(i3)=zatp(j3)

              xat(i4)=xatp(j4)
              yat(i4)=yatp(j4)
              zat(i4)=zatp(j4)

              vlvr(i) = dihedr(i1,i2,i3,i4)

              isrfvr(i) = .true.

            endif

          elseif (it.eq.2) then  ! b.angle
            vlvr(i)=baat(i1)
          elseif (it.eq.1) then  ! b.length
            vlvr(i)=blat(i1)
          endif

          olvlvr(i) = vlvr(i)

        enddo  ! SMMP vars.

        nvr = ivrml1(ntlml)+nvrml(ntlml)-1

! ================================= global parameters for 'nml'

! +++++++++++
       inew=0

       if (inew.eq.1) then
! ++++++++++++++++++++++++

        call setvar(nml,vlvr)

        nrs = irsml2(nml)-irsml1(nml)+1
        call rmsdopt(nml,1,nrs,ixatp,xatp,yatp,zatp,0,rm,av1,av2,rmsd)

! ---------------------------- retrieve ref. coords.
!                     & transform acc. to opt. rmsd
        do i=1,3
          ii=ixrfpt(i,nml)

          h(1)=xat(ii)-av1(1)
          h(2)=yat(ii)-av1(2)
          h(3)=zat(ii)-av1(3)

          x=av2(1)
          y=av2(2)
          z=av2(3)

          do j=1,3
            x = x + rm(j,1) * h(j)
            y = y + rm(j,2) * h(j)
            z = z + rm(j,3) * h(j)
          enddo

          xat(ii)=x
          yat(ii)=y
          zat(ii)=z
        enddo

        call getpar(nml)
        call bldmol(nml)  ! finally build SMMP molecule

! ++++++++++++++++
       else  ! old
! ++++++++++++++++

        call fnd3ba(nml,i1,i2,i3)  ! three 1st bb atoms in SMMP (e.g. n,ca,c')

        ixrfpt(1,nml)=i1
        ixrfpt(2,nml)=i2
        ixrfpt(3,nml)=i3

! -------------------------------- retrieve ref. coords.
        do i=1,3
          ii=ixrfpt(i,nml)
          ix=ixatp(ii)
          if (ix.gt.0) then
            xat(ii)=xatp(ix)
            yat(ii)=yatp(ix)
            zat(ii)=zatp(ix)
          else
            write (logString, '(3a)') ' pdbvars> missing PDB atom ',
     &         nmat(ii), ' is ref. point for SMMP - cannot proceed !'
          endif
        enddo

        call getpar(nml)
        call setvar(nml,vlvr)  ! finally build SMMP molecule

        nrs = irsml2(nml)-irsml1(nml)+1
        call rmsdopt(nml,1,nrs,ixatp,xatp,yatp,zatp,0,rm,av1,av2,rmsd)

! ++++++++++
       endif
! ++++++++++

       write (logString, *) ' '
       write (logString, *) ' Initial RMSD ',rmsd

      enddo ! chains(molecules)

      return
      end
! ***************************
      subroutine atixpdb(nml)

! --------------------------------------------------------------------
! PURPOSE: get ixatp - pointer of each SMMP atom to corresponding atom
!                      of reference structure loaded in 'INCP.H'
!                      (=0 if no corr. atom in ref. str.)
!
! CALLS:   toupst
! --------------------------------------------------------------------

      include 'INCL.H'
      include 'INCP.H'

      integer irs, nml, i1, i2, iat, ix, i

      character*4 atm


      atm=' '

      do irs=irsml1(nml),irsml2(nml)  ! SMMP residues

        i1=irsatp(irs)
        i2=i1+nrsatp(irs)-1

        do iat=iatrs1(irs),iatrs2(irs)  ! SMMP atoms

          ix=0

          if (nmat(iat)(1:1).ne.'h') then ! ignore hydrogens

            atm(2:4)=nmat(iat)(1:3)
            call toupst(atm)

            do i=i1,i2  ! atoms of PDB residue
              if (index(atnmp(i),atm).gt.0) then
                ix=i
                goto 1
              endif
            enddo

!            write (logString, '(8a)') ' pdbvars> ',atm,' not found in '
!     #       ,chnp(nc),' ',rsidp(irs),' ',rsnmp(irs)

          endif

    1     ixatp(iat)=ix

        enddo   ! SMMP atoms of 'irs'
      enddo   ! residues

      return
      end
! **************************
      subroutine getpar(nml)

      include 'INCL.H'

      double precision tol, h1, h2, h3, x1, x2, x3, d, z1, z2, z3, y1
      double precision y2, y3, yp1, yp2

      integer i1, nml, i2, i3, i

      parameter (TOL = 1.d-12)

! Obtain molecule-fixed system (J,K,L) for 1st 3 bb-atoms,
! -> determine global parameters: shifts dX,dY,dZ
! & angles alpha,beta,gamma [rad], put into 'gbpr'
!
! CALLS: none
!

      i1=ixrfpt(1,nml)  ! from 'INCL.H'
      i2=ixrfpt(2,nml)
      i3=ixrfpt(3,nml)
! -------------------------------------- Shifts
      gbpr(1,nml) = xat(i1)
      gbpr(2,nml) = yat(i1)
      gbpr(3,nml) = zat(i1)

      do i = 4,6
        gbpr(i,nml) = 0.d0
      enddo
! --------------------------------- J
      h1=xat(i2)
      h2=yat(i2)
      h3=zat(i2)

      x1=h1-xat(i1)
      x2=h2-yat(i1)
      x3=h3-zat(i1)

      d=sqrt(x1*x1+x2*x2+x3*x3)

      x1=x1/d
      x2=x2/d
      x3=x3/d
! --------------------------------- L
      h1=xat(i3)-h1
      h2=yat(i3)-h2
      h3=zat(i3)-h3

      z1=x2*h3-x3*h2
      z2=x3*h1-x1*h3
      z3=x1*h2-x2*h1

      d=sqrt(z1*z1+z2*z2+z3*z3)

      z1=z1/d
      z2=z2/d
      z3=z3/d

! ---------------------------------- K
      y1=z2*x3-z3*x2
      y2=z3*x1-z1*x3
      y3=z1*x2-z2*x1

      if ( ( 1.d0 - abs(y3) ) .gt. TOL )  then       ! ============ |beta| < PI/2

! ----------------------------------------------- Y'
        d = sqrt( y1 * y1 + y2 * y2 )
        yp1= y1 / d
        yp2= y2 / d

        gbpr(4,nml)= atan2( -yp1, yp2 )                             ! alpha
        gbpr(5,nml)= atan2( y3, ( y1*yp1+y2*yp2 ) )                 ! beta
        gbpr(6,nml)= atan2( ( z1*yp2-z2*yp1 ),( x1*yp2-x2*yp1 ) )   ! gamma

      else        ! ======================= |beta| = PI/2

        gbpr(4,nml) = atan2( x2, x1 )   ! alpha+gamma

        if ( abs(y3) .lt. 1.d0 ) then  ! beta
          gbpr(5,nml) = asin(y3)
        else if ( y3 .gt. 0.d0 ) then
          gbpr(5,nml) = pi*.5d0
        else
          gbpr(5,nml) = -pi*.5d0
        endif

        gbpr(6,nml) = 0.0

      endif

      return
      end

