! **************************************************************
!
! This file contains the subroutines:  bldmol, fnd3ba,eyring,
!                                      setsys,setgbl
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      subroutine bldmol(nml)

! .................................................
! PURPOSE: calculate coordinates for molecule 'nml'
!
! OUTPUT:  xat,yat,zat,xbaat,ybaat,zbaat,xtoat,ytoat,
!          ztoat (via 'eyring')
!
!          1st backbone atom of 1st residue of 'nml': 
!
!          - it's position: from 'gbpr(1-3,nml)'
!          - it's axes: from 'setgbl' according to
!            global angles 'gbpr(4-5,nml)'
!
! CALLS: eyring, fnd3ba,setgbl,setsys
! .................................................

      include 'INCL.H'
      ! Subroutine arguments
      integer nml
      
      integer i, i1, i2, i3, ifirs, jow, j, jj
      double precision x1, x2, x3, y1, y2, y3, z1, z2, z3
      
      double precision xg, zg
      dimension xg(3),zg(3)


      call fnd3ba(nml,i1,i2,i3)
! ------------------------------ first 3 bb atoms of 'nml'
      ixrfpt(1,nml)=i1
      ixrfpt(2,nml)=i2
      ixrfpt(3,nml)=i3

! ------------------------------ position of 1st bb atom
      xat(i1) = gbpr(1,nml)
      yat(i1) = gbpr(2,nml)
      zat(i1) = gbpr(3,nml)

      rfpt(1,nml)=xat(i1)
      rfpt(2,nml)=yat(i1)
      rfpt(3,nml)=zat(i1)

      call setgbl(nml,i1,i2,i3, xg, zg)

      xbaat(i1) = zg(1)
      ybaat(i1) = zg(2)
      zbaat(i1) = zg(3)

      xtoat(i1) = xg(1)
      ytoat(i1) = xg(2)
      ztoat(i1) = xg(3)

      ifirs=irsml1(nml)

      do i=ifirs,irsml2(nml)
        if (i.eq.ifirs) then  ! not construct
          jj=iatrs1(i)+1      ! 1st bb atom again
        else
          jj=iatrs1(i)
        endif
        do j=jj,iatrs2(i)
          jow=iowat(j)
          call eyring(j,jow)
        enddo
      enddo

      call setsys(i1,i2,i3, x1,x2,x3,y1,y2,y3,z1,z2,z3)

      xrfax(1,nml)=x1      ! J
      xrfax(2,nml)=x2
      xrfax(3,nml)=x3
      yrfax(1,nml)=-y1     ! K
      yrfax(2,nml)=-y2
      yrfax(3,nml)=-y3
      zrfax(1,nml)=-z1     ! L
      zrfax(2,nml)=-z2
      zrfax(3,nml)=-z3

      return
      end
! ***********************************
      subroutine fnd3ba(nml,i1,i2,i3)
 
! .................................................
! PURPOSE: return indices 'i1,i2,i3' of
!          first 3 backbone atoms in molecule 'nml'
!
! CALLS:   fndbrn
! .................................................
 
      include 'INCL.H'
! arguments      
      integer nml, i1, i2, i3
      
      integer ibd, i, irs, ix
      dimension ibd(4)
      logical bb


      irs=irsml1(nml)

! --------------------------- 1st bb atom
      i1=iatrs1(irs)

      call fndbrn(nml,irs,i1,i,ix,i2,bb)

! --------------------------- 2nd bb atom
      i2=i+1

! ------------------------ check for ring

      ibd(1)=iowat(i1)
      ibd(2)=ibdat(1,i1)
      ibd(3)=ibdat(2,i1)
      ibd(4)=ibdat(3,i1)

      ix = 0
      do i=1,nbdat(i1)+1

        if (iowat(ibd(i)).ne.i1) then

          if (ix.ne.0) then
            write (logString, '(2a,i3)') 
     &         ' fnd3ba> Can handle only simple ring at 1st',
     &         ' atom of molecule #',nml
            stop
          endif

          ix=ibd(i)
        endif
      enddo

! --------------------------- 3rd bb atom

      ix=ixatrs(irs)

      if (i2.eq.ix) then

        if (irs.lt.irsml2(nml)) then
          i3 = iatrs1(irs+1)
          return
        endif

      else

        do i = ix,i2+1,-1
          if (iowat(i).eq.i2) then
            i3=i
            return
          endif
        enddo

      endif

      write (logString, '(4a,i4,a,i4)') 
     &   ' fnd3ba> Cannot find backbone atom following ',nmat(i2),
     &   ' of residue ',seq(irs),irs,' in molecule #',nml
      stop

      end
! ***************************
      subroutine eyring(i,ia)

! .........................................................
! PURPOSE:  calculate cartesian coordinates of atom 'i'
!           using EYRING's algorithm
! INPUT:    i - index of atom to be constructed
!               for 'i': blat,csbaat,snbaat,cstoat,sntoat
!           ia- index of atom from which 'i' is to be built
! OUTPUT:   for 'i': xat,yat,zat,xbaat,ybaat,zbaat,xtoat,ytoat,ztoat
!
! CALLS: none
! .........................................................

      include 'INCL.H'
!     arguments
      integer i, ia
      double precision ct, st, ca, sa, bl, x1, x2, x3, z1, z2, z3
      double precision y1, y2, y3, h2, h3, dx, dz
      
      ct=cstoat(i)
      st=sntoat(i)
      ca=csbaat(i)
      sa=snbaat(i)
      bl=blat(i)

      x1=xtoat(ia)
      x2=ytoat(ia)
      x3=ztoat(ia)

      z1=xbaat(ia)
      z2=ybaat(ia)
      z3=zbaat(ia)

      y1=z2*x3-z3*x2
      y2=z3*x1-z1*x3
      y3=z1*x2-z2*x1

      h2=-sa*ct
      h3=-sa*st

      x1=-ca*x1+h2*y1+h3*z1
      x2=-ca*x2+h2*y2+h3*z2
      x3=-ca*x3+h2*y3+h3*z3

      dx=one/sqrt(x1*x1+x2*x2+x3*x3)
      x1=x1*dx
      x2=x2*dx
      x3=x3*dx

      xtoat(i)=x1
      ytoat(i)=x2
      ztoat(i)=x3

      z1=-st*y1+ct*z1
      z2=-st*y2+ct*z2
      z3=-st*y3+ct*z3

      dz=one/sqrt(z1*z1+z2*z2+z3*z3)
      xbaat(i)=z1*dz
      ybaat(i)=z2*dz
      zbaat(i)=z3*dz

      xat(i)=xat(ia)+x1*bl
      yat(i)=yat(ia)+x2*bl
      zat(i)=zat(ia)+x3*bl

      return
      end
! ***********************************************************
      subroutine setsys(i1,i2,i3, x1,x2,x3,y1,y2,y3,z1,z2,z3)

! ..........................................................
!  PURPOSE:  calculate axes X,Y,Z of right-handed orthogonal
!            system given by three atom positions R1, R2, R3
!
!            X = (R2-R1)/ |( )|
!            Z = {X x (R2-R3)} / |{ }|
!            Y = Z x X
!
!  INPUT:    i1, i2, i3 - indices of three atoms
!  OUTPUT:   x1,x2,x3 |
!            y1,y2,y3 | -direction cosines of X,Y,Z
!            z1,z2,z3 |
!
!  CALLS:    none
! ...................................................


      include 'INCL.H'
!     arguments
      integer i1, i2, i3
      double precision x1, x2, x3, y1, y2, y3, z1, z2, z3, h1, h2, h3
      double precision dz, dx
      
      h1=xat(i2)
      h2=yat(i2)
      h3=zat(i2)

      x1=h1-xat(i1)
      x2=h2-yat(i1)
      x3=h3-zat(i1)

      y1=h1-xat(i3)
      y2=h2-yat(i3)
      y3=h3-zat(i3)

      z1=x2*y3-x3*y2
      z2=x3*y1-x1*y3
      z3=x1*y2-x2*y1

      dx=one/sqrt(x1*x1+x2*x2+x3*x3)
      x1=x1*dx
      x2=x2*dx
      x3=x3*dx

      dz=one/sqrt(z1*z1+z2*z2+z3*z3)
      z1=z1*dz
      z2=z2*dz
      z3=z3*dz

      y1=z2*x3-z3*x2
      y2=z3*x1-z1*x3
      y3=z1*x2-z2*x1

      return
      end

! *****************************************
      subroutine setgbl(nml,i1,i2,i3,xg,zg)

! ...................................................
!
! PURPOSE: 1. Obtain global axes (J,K,L)
!             related to x(1 0 0),y(0 1 0),z(0 0 1)
!             by 3 rotations (gbl. parameters #4-#6):
!
!             - round z by angle alpha 
!	      - round x' by a. beta 
!	      - round y" by a. gamma
!
!          2. Return x-axis (xg) & z-axis (zg)
!             for atom #1 in order to orientate J
!             along the bond from backbone atom #1
!             to bb.a. #2 and L according to the
!             cross product [ bond(#1->#2) x
!             bond(#2->#3) ] when using Eyring's
!             algorithm to get the coordinates
!
! CALLS:   none
! ..............................................

      include 'INCL.H'
!     arguments
      integer nml, i1, i2, i3
      double precision xg, zg, ag
      dimension xg(3),zg(3),ag(3,3)

      integer i
      double precision ca, sa, cb,sb, cg, sg, d, ct2, ct3, dx, dz
      double precision x1, x2, x3, y1, y2, y3, z1, z2, z3, st2, sa2, st3

      ca = cos(gbpr(4,nml))  ! alpha
      sa = sin(gbpr(4,nml))
      cb = cos(gbpr(5,nml))  ! beta
      sb = sin(gbpr(5,nml))
      cg = cos(gbpr(6,nml))  ! gamma
      sg = sin(gbpr(6,nml))

! ----------------------------- J
      x1 =  ca*cg - sa*sb*sg
      x2 =  sa*cg + ca*sb*sg
      x3 = -cb*sg

      d = sqrt(x1**2+x2**2+x3**2)
      ag(1,1) = x1/d
      ag(2,1) = x2/d
      ag(3,1) = x3/d
! ----------------------------- K
      y1 = -sa*cb
      y2 =  ca*cb
      y3 =  sb

      d = sqrt(y1**2+y2**2+y3**2)
      ag(1,2) = y1/d
      ag(2,2) = y2/d
      ag(3,2) = y3/d
! ----------------------------- L
      z1 =  ca*sg + sa*sb*cg
      z2 =  sa*sg - ca*sb*cg
      z3 =  cb*cg

      d = sqrt(z1**2+z2**2+z3**2)
      ag(1,3) = z1/d
      ag(2,3) = z2/d
      ag(3,3) = z3/d

! ------------------------------------ X1
      ct2 = cstoat(i2)
      st2 = sntoat(i2)
      sa2 = snbaat(i2)

      x1 = -csbaat(i2)
      x2 = -sa2*ct2
      x3 = -sa2*st2

      dx = sqrt(x1**2+x2**2+x3**2)
      x1 = x1/dx
      x2 = x2/dx
      x3 = x3/dx
! ------------------------------------- Z1
      st3 = sntoat(i3)
      ct3 = cstoat(i3)

      z1 = -st3*(st2*x3+ct2*x2)
      z2 =  st3*ct2*x1 + ct3*st2
      z3 =  st3*st2*x1 - ct3*ct2

      dz = sqrt(z1**2+z2**2+z3**2)
      z1 = z1/dz
      z2 = z2/dz
      z3 = z3/dz
! ------------------------------------- Y1
      y1 = z2 * x3 - z3 * x2
      y3 = z1 * x2 - z2 * x1  ! do not need y2

! ----------------------------- into global system

      xg(1) = ag(1,1)*x1 + ag(1,2)*y1 + ag(1,3)*z1
      xg(2) = ag(2,1)*x1 + ag(2,2)*y1 + ag(2,3)*z1
      xg(3) = ag(3,1)*x1 + ag(3,2)*y1 + ag(3,3)*z1

      dx = sqrt(xg(1)**2+xg(2)**2+xg(3)**2)

      zg(1) = ag(1,1)*x3 + ag(1,2)*y3 + ag(1,3)*z3
      zg(2) = ag(2,1)*x3 + ag(2,2)*y3 + ag(2,3)*z3
      zg(3) = ag(3,1)*x3 + ag(3,2)*y3 + ag(3,3)*z3

      dz = sqrt(zg(1)**2+zg(2)**2+zg(3)**2)

      do i=1,3
        xg(i) = xg(i)/dx
        zg(i) = zg(i)/dz
      enddo

      return
      end

