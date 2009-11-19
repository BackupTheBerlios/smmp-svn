! **************************************************************
!
! This file contains the subroutines: rmsdfun,rmsdopt,fitmol,
!                                     jacobi,rmsinit
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      real*8 function rmsdfun(nml,ir1,ir2,ixat,xrf,yrf,zrf,isl)
!
! --------------------------------------------------------------
! Wrapping function for calculating rmsd
!
! LIMITATION: requires call of rmsinit BEFORE calling this function
!
! CALLS: rmsdopt
!
! ---------------------------------------------------------------
!
      include 'INCL.H'
      include 'INCP.H'
!
! Input
      double precision xrf, yrf, zrf, rm, av1, av2, rssd

      integer nml, ir1, ir2, ixat, isl

      dimension ixat(mxat),xrf(mxatp),yrf(mxatp),zrf(mxatp)
! Local
      dimension rm(3,3),av1(3),av2(3)
      call rmsdopt(nml,ir1,ir2,ixat,xrf,yrf,zrf,isl,rm,av1,av2,rssd)

      rmsdfun = rssd

      return

      end

!*******************************************************************
      subroutine rmsdopt(nml,ir1,ir2,ixat,xrf,yrf,zrf,isl,
     &                   rm,av1,av2,rmsd)

! ---------------------------------------------------------------
! PURPOSE: root mean square deviation (rmsd) between current SMMP
!          structure and reference atom coordinates 'x,y,zrf()'
!          for range of SMMP residues [ir1,ir2] in molecule 'nml'
!
!          ixat(i) - points to the atom in ref. coords., which is
!                    equivalent to atom i of SMMP structure
!                    (=0 if no equivalent in ref. structure exists)
!
!          isl = 0  : select all heavy atoms
!          isl = 1  : backbone atoms n,ca,c
!          isl = 2  : only ca atoms
!
! CALLS:   fitmol [S.K.Kearsley, Acta Cryst. 1989, A45, 208-210]
!
!      NB  uncomment last lines in 'fitmol' to return coordinates
!          in 'x2' after fitting the ref. str. onto SMMP structure
! ----------------------------------------------------------------

      include 'INCL.H'
      include 'INCP.H'
      
      double precision x1, x2, xrf, yrf, zrf, rm, av1, av2, rmsd

      integer nml, nr, na, n, im, ir, ia, ir1, ir2, ix, ixat, isl


!-------------------------------------------------------- input
      dimension ixat(mxat),xrf(mxatp),yrf(mxatp),zrf(mxatp)
!-------------------------------------------------------- output
      dimension rm(3,3),av1(3),av2(3)
!-------------------------------------------------------- local
      dimension x1(3,mxat),x2(3,mxat)
      character*4 atnm


      if (nml.lt.1.or.nml.gt.ntlml) then
        write(*,*) ' rmsdopt>  Sorry, there is no molecule #',nml
        stop
      endif


      nr=0
      na=0
      n=0

      do im=1,ntlml
        do ir=irsml1(im),irsml2(im)
          if (im.eq.nml) nr=nr+1

          do ia=iatrs1(ir),iatrs2(ir)

            na=na+1

            if (im.eq.nml.and.nr.ge.ir1.and.nr.le.ir2) then  ! range of res. for 'nml'

              ix=ixat(na)
              atnm=nmat(ia)

              if (ix.gt.0.and.atnm(1:1).ne.'h') then

                if ( isl.eq.0

     &               .or.

     &              (isl.eq.1.and.(index(atnm,'n   ').gt.0 .or.
     &                             index(atnm,'ca  ').gt.0 .or.
     &                             index(atnm,'c   ').gt.0 ))
     &               .or.

     &              (isl.eq.2.and.index(atnm,'ca  ').gt.0)

     &          ) then

                  n=n+1
                  x1(1,n)=xat(ia)
                  x1(2,n)=yat(ia)
                  x1(3,n)=zat(ia)
                  x2(1,n)=xrf(ix)
                  x2(2,n)=yrf(ix)
                  x2(3,n)=zrf(ix)

                endif  ! atom selection
              endif  ! ix>0 & not 'h'

            endif  ! res. range in mol. 'nml'

          enddo  ! atoms
        enddo  ! residues
      enddo  ! molecules

      if (n.lt.3) then
        write(*,*) ' rmsdopt>  <3 atoms selected !'
        stop
      endif

      call fitmol(n,x1,x2, rm,av1,av2,rmsd)

      return
      end
! *********************************************
      subroutine fitmol(n,x1,x2, rm,a1,a2,rmsd)
!      real*8 function fitmol(n,x1,x2)

! .......................................................
! PURPOSE: compute RMSD of n positions in x1(3,) & x2(3,)
!          [S.K.Kearsley Acta Cryst. 1989,A45,208-210]
!
! CALLS: jacobi
! .......................................................
!f2py intent(out) rmsd

      include 'INCL.H'
!      implicit real*8 (a-h,o-z)
!      implicit integer*4 (i-n)

! ------------------------------------------- input/output
      double precision dn, a1, a2, x1, x2, q, dm, dp, dxm, dym, dzm, dxp
      double precision dyp, dzp, e, v, em, rmsd, rm

      integer n, i, j, ndim4, im

      dimension x1(3,mxat),x2(3,mxat)
! -------------------------------------------------- local
      dimension e(4),q(4,4),v(4,4),dm(3),dp(3),a1(3),a2(3),rm(3,3)

      dn=dble(n)
! ------------------- average of coordinates
      do i=1,3
        a1(i) = 0.d0
        a2(i) = 0.d0
        do j=1,n
          a1(i) = a1(i) + x1(i,j)
          a2(i) = a2(i) + x2(i,j)
        enddo
        a1(i) = a1(i)/dn
        a2(i) = a2(i)/dn
      enddo
! ------------------------- compile quaternion
      do i=1,4
        do j=1,4
          q(i,j)=0.d0
        enddo
      enddo

      do i=1,n

        do j=1,3
          dm(j) = x1(j,i)-a1(j)
          dp(j) = x2(j,i)-a2(j)
        enddo

        dxm = dp(1) - dm(1)
        dym = dp(2) - dm(2)
        dzm = dp(3) - dm(3)
        dxp = dp(1) + dm(1)
        dyp = dp(2) + dm(2)
        dzp = dp(3) + dm(3)

        q(1,1) = q(1,1) + dxm * dxm + dym * dym + dzm * dzm
        q(1,2) = q(1,2) + dyp * dzm - dym * dzp
        q(1,3) = q(1,3) + dxm * dzp - dxp * dzm
        q(1,4) = q(1,4) + dxp * dym - dxm * dyp
        q(2,2) = q(2,2) + dyp * dyp + dzp * dzp + dxm * dxm
        q(2,3) = q(2,3) + dxm * dym - dxp * dyp
        q(2,4) = q(2,4) + dxm * dzm - dxp * dzp
        q(3,3) = q(3,3) + dxp * dxp + dzp * dzp + dym * dym
        q(3,4) = q(3,4) + dym * dzm - dyp * dzp
        q(4,4) = q(4,4) + dxp * dxp + dyp * dyp + dzm * dzm

      enddo

      do i=1,3
        do j=i+1,4
          q(j,i)=q(i,j)
        enddo
      enddo
! ------------------------------ eigenvalues & -vectors
      ndim4=4
      call jacobi(q,ndim4,e,v)
! --------------------------- lowest eigenvalue
      im=1
      em=e(1)
      do i=2,4
         if (e(i).lt.em) then
           em=e(i)
           im=i
         endif
      enddo

      rmsd = sqrt(em/dn)

! ================= uncomment following lines to fit molecule 2 onto 1

! ---------------------------------------------------rotation matrix
      rm(1,1) = v(1,im)**2+v(2,im)**2-v(3,im)**2-v(4,im)**2
      rm(1,2) = 2.d0*( v(2,im)*v(3,im)-v(1,im)*v(4,im) )
      rm(1,3) = 2.d0*( v(2,im)*v(4,im)+v(1,im)*v(3,im) )
      rm(2,1) = 2.d0*( v(2,im)*v(3,im)+v(1,im)*v(4,im) )
      rm(2,2) = v(1,im)**2+v(3,im)**2-v(2,im)**2-v(4,im)**2
      rm(2,3) = 2.d0*( v(3,im)*v(4,im)-v(1,im)*v(2,im) )
      rm(3,1) = 2.d0*( v(2,im)*v(4,im)-v(1,im)*v(3,im) )
      rm(3,2) = 2.d0*( v(3,im)*v(4,im)+v(1,im)*v(2,im) )
      rm(3,3) = v(1,im)**2+v(4,im)**2-v(2,im)**2-v(3,im)**2

!      do i=1,n
!        do j=1,3
!          dm(j) = x2(j,i) - a2(j)
!        enddo
!        do j=1,3
!          dp(j) = a1(j)
!          do k=1,3
!            dp(j) = dp(j) + rm(j,k) * dm(k)
!          enddo
!          x2(j,i) = dp(j)
!        enddo
!      enddo

!      fitmol=rmsd

      return
      end
! ******************************
      subroutine jacobi(a,n,d,v)

! ......................................................
! PURPOSE: for given symmetric matrix 'a(n,n)
!          compute eigenvalues 'd' & eigenvectors 'v(,)'
!
!  [W.H.Press,S.A.Teukolsky,W.T.Vetterling,
!   B.P.Flannery, Numerical Recipes in FORTRAN,
!   Cambridge Univ. Press, 2nd Ed. 1992, 456-462]
!
! CALLS: none
!
! ......................................................

!f2py intent(out) d
!f2py intent(out) v
      integer nmax

      parameter (NMAX=500)

      integer n,nrot,i,ip,iq,j


      real*8 a(n,n),d(n),v(n,n),
     &       c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX),smeps

      smeps=1.0d-6


      do ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
        do iq=1,n
          v(ip,iq)=0.d0
        enddo
        v(ip,ip)=1.d0
      enddo

      nrot=0

      do i=1,500

        sm=0.d0

        do ip=1,n-1
          do iq=ip+1,n
            sm=sm+abs(a(ip,iq))
          enddo
        enddo

        if (sm.le.smeps) return  ! normal end

        if (i.lt.4) then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif

        do ip=1,n-1
          do iq=ip+1,n

            g=100.d0*abs(a(ip,iq))

            if((i.gt.4).and.(abs(d(ip))+

     &g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.d0

            else if(abs(a(ip,iq)).gt.tresh)then

              h=d(iq)-d(ip)

              if (abs(h)+g.eq.abs(h)) then

                t=a(ip,iq)/h

              else

                theta=0.5d0*h/a(ip,iq)
                t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                if (theta.lt.0.d0) t=-t

              endif

              c=1.d0/sqrt(1.d0+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0

              do j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo

              do j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
              enddo
              do j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
              enddo
              nrot=nrot+1

            endif

          enddo
        enddo

        do ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.d0
        enddo

      enddo

      write(*,*) ' jacobi> too many iterations'
      stop

      return
      end

! ***********************************************************

      subroutine rmsinit(nml,string)
!
!------------------------------------------------------------------------------
! Reads in pdb-file 'string' into INCP.H and initalizes
! the files that 'rmdsopt' needs to calculate the rmsd
! of a configuration with the pdb-configuration
!
! CALLS: pdbread,atixpdb
!
! ----------------------------------------------------------------------------
!
      include 'INCL.H'
      include 'INCP.H'

      integer i, nml, ier

      character string*(*)

      if(string.eq.'smmp') then
!
! Compare with a smmp-structure
!
        do i=iatrs1(irsml1(nml)),iatrs2(irsml2(nml))
         if(nmat(i)(1:1).ne.'h') then
            ixatp(i)=i
         else
            ixatp(i) = 0
         end if
        enddo
!
       else
!
! Reference structure is read in from pdb-file
!
         call pdbread(string,ier)
         if(ier.ne.0) stop
         call atixpdb(nml)
!
      end if
      print *,'RMSD initialized with ',string
      return

      end

