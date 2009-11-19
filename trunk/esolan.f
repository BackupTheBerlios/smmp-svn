! **************************************************************
!
! This file contains the subroutines: esolan
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      real*8 function esolan(nmol)

! -----------------------------------------------------------------
!
! Calculates the solvation energy of the protein with 
! solavtion parameters model E=\sum\sigma_iA_i. 
! The solvent accessible surface area per atom
! and its gradients with respect to the Cartesian coordinates of 
! the atoms are calculated using the projection method described in
!
!            J Comp Phys 26, 334-343, 2005
!
! USAGE: eysl=esolan(nmol)-Returns the value of the solvation energy
!
! INPUT: nmol - the order number of the protein chain.
!      The atomic coordinates, radii and solvation parameters  
!      are taken from the global arrays of SMMP
!
! OUTPUT: global array gradan(mxat,3) contains the gradients of 
!         solvation energy per atom
!         local array as(mxat) contains accessible surface area per atom
!
! Output file "solvation.dat' conatains detailed data.
!
! Correctness of this program was last time rechecked with GETAREA
!
! CALLS: none
!
! ----------------------------------------------------------------------
! TODO Test the solvent energy for multiple molecules

      include 'INCL.H'
      
      
      integer nmol, ks0, ks2
Cf2py intent(in) nmol
      
      parameter (ks0=mxat*(mxat+1))
      parameter (ks2=mxat+mxat)
! Functions
      real*8 grnd

      integer neib, neibp, ivrx, iv, ial, i, ij, j, iii, idi,ivr, ifu,
     &          ita2, ita3, ine, kk, j2, jkk, jjj, jj, jk, k, l, ll, 
     &          nb, numat, nyx
      real*8 vertex, ax, pol, as, ayx, ayx1, probe, dd, ddat, dadx, gp,
     &          dalp,dbet, daalp, dabet, vrx, dv, dx, dy, dz, dt, di,
     &          dii, ss, dta, dtb, di1, di2, gs, xold, yold, zold, 
     &          energy, sss, a1, a, ab1, a2, aa, ac2, ac1, ab2, am2,
     &          am, aia, ak, alp, am1, ap1, amibe, amabe, amaal, amial,
     &          ang, D, ddp, cf, b2, arg2, arg1, ap2, b, b1, b2c2, ca,
     &          ba, bcd, bc, bet, c1, c, c2, cc, sfcg, caba, cg, cfsg,
     &          cfcg, daba, ct, cs, d1, d2, dbe, dal, dcotbma, dcbet,
     &          dcalp, dcbetmalp, dcbetpalp, ddp2, div, dsbetmalp, 
     &          dsalp, dsbet, yy, dsbetpalp, enr, p1, p2, r42, r22, 
     &          r22p, r2, pom, prob, r, r42p, ratom, rr, s, sfsg, sf,
     &          sg, vv1, t, ufi, tt, uga, uv, vv, vv2, wv, xx, xv, yv, 
     &          zz, zv
      dimension neib(0:mxat),neibp(0:mxat),ivrx(ks0)
      dimension vertex(ks0,4),ax(ks0,2),
     &          pol(mxat),as(mxat),ayx(ks0,2),
     &          ayx1(ks0),probe(ks0),dd(mxat),ddat(mxat,4)
      dimension dadx(4,3),gp(4),dalp(4),dbet(4),
     &   daalp(4),dabet(4),vrx(ks0,4),dv(4),dx(4),dy(4),dz(4),dt(4),
     &   di(4),dii(4,3),ss(mxat),dta(4),dtb(4),di1(4),di2(4),gs(3)
      dimension xold(-1:mxat),yold(-1:mxat),zold(-1:mxat)
      integer ta2(0:mxat),ta3(0:mxat),fullarc(0:mxat),al(0:ks2)
      real*8 neibor(mxat,4)
      
      real*8, dimension(:, :, :), allocatable :: grad
      
      allocate(grad(mxat, mxat, 3))
      
      
!     open(3,file='solvation.dat')! Output file with solvation data

      energy=0.0d0 ! Total solvation energy
      sss=0.0d0 ! Total solvent accessible surface area
      if (nmol.eq.0) then
         numat=iatrs2(irsml2(ntlml))-iatrs1(irsml1(1))+1 ! Number of atoms
      else 
         numat=iatrs2(irsml2(nmol))-iatrs1(irsml1(nmol))+1 ! Number of atoms
      endif

      do i=1,numat
        do j=1,3
         gradan(i,j)=0.0d0 
      end do
      end do

      do i=1,numat
         xold(i)=xat(i)! Redefine the coordinates just for safety
         yold(i)=yat(i)! 
         zold(i)=zat(i)! 
      pol(i)=rvdw(i)*rvdw(i)!The water radius is already added in 'main' 
      As(i)=pi4*pol(i)! Initially the whole surface of the atom is accessible.
      end do
      
             do iii=1,numat
              do jjj=1,numat
                grad(iii,jjj,1)=0.d0
                grad(iii,jjj,2)=0.d0
                grad(iii,jjj,3)=0.d0
              enddo
             enddo

!***************************
!   Start the computations *
!***************************

      do 1 i=1,numat ! Lop over atoms "i"           
! ----------------------------------------------------
         
       R=rvdw(i) 
       jj=0
        do j=1,numat !Find the neighbours of "i"
       if(i.ne.j) then
       ddp=(xold(j)-xold(i))**2+(yold(j)-yold(i))**2
     &   +(zold(j)-zold(i))**2! Distance between atom i and j
       if(ddp.lt.1e-10) then
         write(*,*)'ERROR in data: centres of two atoms coincide!'
         write(*,*)i,j,xold(i),yold(i),zold(i),rvdw(i),
     &                   xold(j),yold(j),zold(j),rvdw(j) 
         stop 'Centres of atoms coincide!'
       endif
       ddp=dsqrt(ddp)
       if(ddp+R.le.rvdw(j)) then! Atom "i" is enclosed within "j"
         As(i)=0.d0 ! no accessible surface
         goto 1
       endif
       if(.not.(ddp.le.R-rvdw(j).or.ddp.ge.R+rvdw(j))) then
        jj=jj+1 ! The next neighbour of ith atom
        neib(jj)=j ! and its order number
       endif
         end if
         end do!Neighbours
       
       neib(0)=jj   ! The number of neighbors of "i"

       if(neib(0).eq.0) goto 1 !Finish the atom i if it doesn't have neighbors

!----
       R2=R*R                   ! R square
       R22=R2+R2             ! 2 * R square
       R42=R22+R22       ! 4 * R square
       R22p=R22*pi       ! 2 pi R square
       R42p=R22p+R22p       ! 4 pi R square

         do j=1,neib(0)
           k=neib(j)
           ddat(k,2)=xold(k)-xold(i)
           ddat(k,3)=yold(k)-yold(i)
           ddat(k,4)=zold(k)-zold(i)
           ddat(k,1)=rvdw(k)
         enddo
           ddat(i,2)=0d0
           ddat(i,3)=0d0
           ddat(i,4)=0d0
           ddat(i,1)=R
!
!    Verification of the north point of ith sphere
!
       jjj=0
!
!      If jjj=0 then - no transformation
!             else the transformation is necessary
!
       k=neib(1)  ! The order number of the first neighbour
! ddp - the square minimal distance of NP from the neighboring surfaces
       ddp=dabs((ddat(k,2))**2+(ddat(k,3))**2+
     &    (R-ddat(k,4))**2-pol(k))

       do j=2,neib(0)
        k=neib(j)
        ddp2=dabs((ddat(k,2))**2+(ddat(k,3))**2+
     &      (R-ddat(k,4))**2-pol(k))
        if(ddp2.lt.ddp) ddp=ddp2
       enddo
!
!   Check whether the NP of ith sphere is too close to the intersection line
!
       do while (ddp.lt.1.d-6)
        jjj=jjj+1
!
!       generate grndom numbers
!
        uga=grnd() ! Random \gamma angle
        ufi=grnd() ! Random \pi angle 
        uga=pi*uga/2
        ufi=pi*ufi/2
         cf=dcos(ufi)
         sf=dsin(ufi)
         cg=dcos(uga)
         sg=dsin(uga)
         cfsg=cf*sg
         cfcg=cf*cg
        xx=R*cfsg
        yy=R*sf
        zz=R*cfcg
        k=neib(1)
        ddp=dabs((xx-ddat(k,2))**2+(yy-ddat(k,3))**2+
     &     (zz-ddat(k,4))**2-pol(k))
        do j=2,neib(0)
         k=neib(j)
         ddp2=dabs((xx-ddat(k,2))**2+(yy-ddat(k,3))**2+
     &       (zz-ddat(k,4))**2-pol(k))
         if(ddp2.lt.ddp) ddp=ddp2
        end do
      end do
!       
      
       if(jjj.ne.0) then ! Rotation is necessary
         sfsg=sf*sg
         sfcg=sf*cg
         do j=1,neib(0)
           k=neib(j)
           xx=ddat(k,2)
           yy=ddat(k,3)
           zz=ddat(k,4)
           ddat(k,2)=xx*cg-zz*sg            ! (x) Coordinates
           ddat(k,3)=-xx*sfsg+yy*cf-zz*sfcg ! (y) after
           ddat(k,4)=xx*cfsg+yy*sf+zz*cfcg  ! (z) rotation
         enddo
       endif

!  iiiiiiiiiii

       pom=8.d0*pol(i)

! In this loop the parameters a,b,c,d for the equation
!      a*(t^2+s^2)+b*t+c*s+d=0 are calculated (see the reference article)
      do jj=1,neib(0)
      j=neib(jj)
       neibor(j,1)=(ddat(j,2))**2+(ddat(j,3))**2+      
     &            (ddat(j,4)-R)**2-pol(j)           ! a
        neibor(j,2)=-pom*ddat(j,2)                  ! b
        neibor(j,3)=-pom*ddat(j,3)                  ! c
        neibor(j,4)=4d0*pol(i)*((ddat(j,2))**2+
     &      (ddat(j,3))**2+(R+ddat(j,4))**2-pol(j)) ! d
       enddo

!       end of the 1st cleaning

       iv=0

       nb=neib(0)
       k=neib(0)
       do while(k.gt.1)               ! B
        k=k-1
! Analyse mutual disposition of every pair of neighbours                                                      
        do 13 L=neib(0),k+1,-1        ! A

         if(neibor(neib(k),1).gt.0d0.and.neibor(neib(L),1).gt.0d0) then ! 03 a01

          b1=neibor(neib(k),2)/neibor(neib(k),1)
          c1=neibor(neib(k),3)/neibor(neib(k),1)
          d1=neibor(neib(k),4)/neibor(neib(k),1)
          b2=neibor(neib(L),2)/neibor(neib(L),1)
          c2=neibor(neib(L),3)/neibor(neib(L),1)
          d2=neibor(neib(L),4)/neibor(neib(L),1)
          D=4d0*((b1-b2)*(b2*d1-b1*d2)+(c1-c2)*(c2*d1-c1*d2)-
     &      (d1-d2)**2)+(b1*c2-b2*c1)**2
! if D<0 then the circles neib(k) and neib(L) don't intersect
          if(D.le.0d0.and.dsqrt((b1-b2)**2+(c1-c2)**2).le.      ! a01 01
     &     dsqrt(b2*b2+c2*c2-4d0*d2)-dsqrt(b1*b1+c1*c1-4d0*d1)) then ! 04
! Circle neib(L) encloses circle neib(k) and the later is discarded
           neib(0)=neib(0)-1
           do j=k,neib(0)
            neib(j)=neib(j+1)
           enddo
           goto 12
          elseif(D.le.0d0.and.dsqrt((b1-b2)**2+(c1-c2)**2).le.      ! a01 02
     &    -dsqrt(b2*b2+c2*c2-4d0*d2)+dsqrt(b1*b1+c1*c1-4d0*d1)) then
! The circle neib(k) encloses circle neib(L) and the later is discarded 
           neib(0)=neib(0)-1
           do j=L,neib(0)
            neib(j)=neib(j+1)
           enddo
          elseif(D.gt.0d0) then                        ! a01 03
! The circles nieb(L) and neib(k) have two intersection points (IP)
           am=2d0*((b1-b2)**2+(c1-c2)**2)
           t=-2d0*(b1-b2)*(d1-d2)-(b2*c1-b1*c2)*(c1-c2)
           s=-2d0*(c1-c2)*(d1-d2)+(b2*c1-b1*c2)*(b1-b2)
           iv=iv+1
           pom=dsqrt(D)
           p1=(c1-c2)*pom
           p2=(b1-b2)*pom
           vertex(iv,1)=(t+p1)/am ! t coordinate of the first IP
           vertex(iv,2)=(s-p2)/am ! s coordinate of the first IP
           vertex(iv,3)=neib(k)   ! the order number of the first circle
           vertex(iv,4)=neib(L)   ! the order number of the second circle
           iv=iv+1
           vertex(iv,1)=(t-p1)/am ! t coordinate of the second IP
           vertex(iv,2)=(s+p2)/am ! s coordinate of the second IP
           vertex(iv,3)=neib(k)   ! the order number of the first circle
           vertex(iv,4)=neib(L)   ! the order number of the second circle
          endif                     ! 04

         elseif(neibor(neib(k),1).gt.0d0.and.      ! a03
     &              neibor(neib(L),1).lt.0d0) then
          b1=neibor(neib(k),2)/neibor(neib(k),1)
          c1=neibor(neib(k),3)/neibor(neib(k),1)
          d1=neibor(neib(k),4)/neibor(neib(k),1)
          b2=neibor(neib(L),2)/neibor(neib(L),1)
          c2=neibor(neib(L),3)/neibor(neib(L),1)
          d2=neibor(neib(L),4)/neibor(neib(L),1)
          D=4d0*((b1-b2)*(b2*d1-b1*d2)+(c1-c2)*(c2*d1-c1*d2)-
     &        (d1-d2)**2)+(b1*c2-b2*c1)**2
! if D<0 then neib(k) and neib(L) don't intersect
          if(D.le.0d0.and.dsqrt((b1-b2)**2+(c1-c2)**2).le.  ! a03 01
     &     -dsqrt(b2*b2+c2*c2-4d0*d2)+dsqrt(b1*b1+c1*c1-4d0*d1)) then ! 06
! Neighbours neib(k) and neib(L) cover fully the atom "i" and as(i)=0
           As(i)=0.d0
           goto 11
          elseif(D.le.0d0.and.dsqrt((b1-b2)**2+(c1-c2)**2).ge.  ! a03 02
     &      dsqrt(b2*b2+c2*c2-4d0*d2)+dsqrt(b1*b1+c1*c1-4d0*d1)) then
! Don't exclude neib(k) 
           neib(0)=neib(0)-1
           do j=k,neib(0)
            neib(j)=neib(j+1)
           enddo
           goto 12
          elseif(D.gt.0.d0) then                        ! a03 03
! Circles neib(L) and neib(k) have two intersection points.
           am=2d0*((b1-b2)**2+(c1-c2)**2)
           t=-2d0*(b1-b2)*(d1-d2)-(b2*c1-b1*c2)*(c1-c2)
           s=-2d0*(c1-c2)*(d1-d2)+(b2*c1-b1*c2)*(b1-b2)
           iv=iv+1
           pom=dsqrt(D)
           p1=(c1-c2)*pom
           p2=(b1-b2)*pom
           vertex(iv,1)=(t+p1)/am ! t coordinate of the first IP
           vertex(iv,2)=(s-p2)/am ! s coordinate of the first IP
           vertex(iv,3)=neib(k)   ! order number of the first circle
           vertex(iv,4)=neib(L)   ! order number of the second circle
           iv=iv+1
           vertex(iv,1)=(t-p1)/am ! t coordinate of the second IP
           vertex(iv,2)=(s+p2)/am ! s coordinate of the second IP
           vertex(iv,3)=neib(k)   ! order number of the first circle
           vertex(iv,4)=neib(L)   ! order number of the second circle
          endif                                                         ! 06

         elseif(neibor(neib(k),1).lt.0d0.and.neibor(neib(L),1).gt.0d0)! a07
     &           then
          b1=neibor(neib(k),2)/neibor(neib(k),1)
          c1=neibor(neib(k),3)/neibor(neib(k),1)
          d1=neibor(neib(k),4)/neibor(neib(k),1)
          b2=neibor(neib(L),2)/neibor(neib(L),1)
          c2=neibor(neib(L),3)/neibor(neib(L),1)
          d2=neibor(neib(L),4)/neibor(neib(L),1)
          D=4d0*((b1-b2)*(b2*d1-b1*d2)+(c1-c2)*(c2*d1-c1*d2)-
     &        (d1-d2)**2)+(b1*c2-b2*c1)**2
! If D<0 then the circles neib(k) and neib(L) don't intersect
          if(D.le.0d0.and.dsqrt((b1-b2)**2+(c1-c2)**2).le.   ! a07 01
     &     dsqrt(b2*b2+c2*c2-4d0*d2)-dsqrt(b1*b1+c1*c1-4d0*d1)) then    ! 10
! atom "i" is covered fully by neib(k) and neib(L) and as(i)=0.
           As(i)=0.d0
           goto 12
          elseif(D.le.0d0.and.dsqrt((b1-b2)**2+(c1-c2)**2).ge.   ! a07 02
     &      dsqrt(b2*b2+c2*c2-4d0*d2)+dsqrt(b1*b1+c1*c1-4d0*d1)) then
! discard the circle neib(L) 
           neib(0)=neib(0)-1
           do j=L,neib(0)
            neib(j)=neib(j+1)
           enddo
          elseif(D.gt.0d0) then                          ! a07 03
! There are two IP's between neib(k) and neib(L).
           am=2d0*((b1-b2)**2+(c1-c2)**2)
           t=-2d0*(b1-b2)*(d1-d2)-(b2*c1-b1*c2)*(c1-c2)
           s=-2d0*(c1-c2)*(d1-d2)+(b2*c1-b1*c2)*(b1-b2)
           iv=iv+1
           pom=dsqrt(D)
           p1=(c1-c2)*pom
           p2=(b1-b2)*pom
!-- Assign t and s coordinates and the order number of circles
           vertex(iv,1)=(t+p1)/am 
           vertex(iv,2)=(s-p2)/am
           vertex(iv,3)=neib(k)
           vertex(iv,4)=neib(L)
           iv=iv+1
           vertex(iv,1)=(t-p1)/am
           vertex(iv,2)=(s+p2)/am
           vertex(iv,3)=neib(k)
           vertex(iv,4)=neib(L)
!--  
          endif                                                 ! 10

         elseif(neibor(neib(k),1).lt.0d0.and.neibor(neib(L),1).lt.0d0) ! a09
     &           then
          b1=neibor(neib(k),2)/neibor(neib(k),1)
          c1=neibor(neib(k),3)/neibor(neib(k),1)
          d1=neibor(neib(k),4)/neibor(neib(k),1)
          b2=neibor(neib(L),2)/neibor(neib(L),1)
          c2=neibor(neib(L),3)/neibor(neib(L),1)
          d2=neibor(neib(L),4)/neibor(neib(L),1)
          D=4d0*((b1-b2)*(b2*d1-b1*d2)+(c1-c2)*(c2*d1-c1*d2)-
     &        (d1-d2)**2)+(b1*c2-b2*c1)**2
! D<0 - no intersection poit between neib(k) and neib(L)
          if(D.le.0d0.and.dsqrt((b1-b2)**2+(c1-c2)**2).le.
     &   dsqrt(b2*b2+c2*c2-4d0*d2)-dsqrt(b1*b1+c1*c1-4d0*d1)) then!12 ! a09 01
! omit the circle neib(L)
           neib(0)=neib(0)-1
           do j=L,neib(0)
            neib(j)=neib(j+1)
           enddo
          elseif(D.le.0d0.and.dsqrt((b1-b2)**2+(c1-c2)**2).le. ! a09 02
     &     -dsqrt(b2*b2+c2*c2-4d0*d2)+dsqrt(b1*b1+c1*c1-4d0*d1)) then
! Omit the circle neib(k)
           neib(0)=neib(0)-1
           do j=k,neib(0)
            neib(j)=neib(j+1)
           enddo
           goto 12
          elseif(D.le.0.d0) then                            ! a09 03
! The whole surface of atom "i" is covered fully, and as(i)=0.
           As(i)=0.d0
           goto 11
          else                                             ! a09 04
! Two intersection points
           am=2.d0*((b1-b2)**2+(c1-c2)**2)
           t=-2.d0*(b1-b2)*(d1-d2)-(b2*c1-b1*c2)*(c1-c2)
           s=-2.d0*(c1-c2)*(d1-d2)+(b2*c1-b1*c2)*(b1-b2)
           iv=iv+1
           pom=dsqrt(D)
           p1=(c1-c2)*pom
           p2=(b1-b2)*pom
! Assign the t and s coordinates and order numbers of the circles
           vertex(iv,1)=(t+p1)/am
           vertex(iv,2)=(s-p2)/am
           vertex(iv,3)=neib(k)
           vertex(iv,4)=neib(L)
           iv=iv+1
           vertex(iv,1)=(t-p1)/am
           vertex(iv,2)=(s+p2)/am
           vertex(iv,3)=neib(k)
           vertex(iv,4)=neib(L)
          endif                                     ! 12

         endif                                      ! 03
13        continue
12        continue                                  ! A
        enddo                                       ! B

        ita2=0
        ita3=0
        do j=1,neib(0)
        if(neibor(neib(j),1).eq.0.d0) then
! Collect the lines (after rotation it is empty).
         ita2=ita2+1
         ta2(ita2)=j
        endif
        if(neibor(neib(j),1).lt.0.d0) then
! Collect the neighbours with inner part
         ita3=ita3+1
         ta3(ita3)=j
        endif
       enddo

!*** Consider to remove this part.
         if(ita2.gt.0.and.ita3.lt.1) then
           As(i)=R22p
           if(ita2.gt.1) then
            amial=pi
            amaal=-pi
            amibe=pi
            amabe=-pi
            do j=1,ita2
             jj=ta2(j)
             p1=neibor(neib(jj),3)
             p2=neibor(neib(jj),2)
             alp=datan2(p1,p2)
             bet=alp
             if(alp.le.0d0) bet=bet+pi2
             if(amial.gt.alp) amial=alp
             if(amaal.lt.alp) amaal=alp
             if(amibe.gt.bet) amibe=bet
             if(amabe.lt.bet) amabe=bet
            enddo
            dal=amaal-amial
            dbe=amabe-amibe
            if(dal.lt.pi) then
             As(i)=R22*(pi-dal)  
            elseif(dbe.lt.pi) then
             As(i)=R22*(pi-dbe)
            else
             As(i)=0d0
            endif
           endif
         elseif((nb.gt.0.and.neib(0).lt.1).or.ita3.gt.0) then
          As(i)=0d0
         endif
!****** End of could-be-removed part

       neibp(0)=neib(0)
       ifu=0 ! Full arcs (without intersection points).
       ine=0 ! Circles with intersection points.
       do j=1,neib(0)
        neibp(j)=neib(j)
       enddo
       do j=1,neib(0)
! "iv" is the number of intersection point. 
        do k=1,iv
       if(neibp(j).eq.vertex(k,3).or.neibp(j).eq.vertex(k,4)) goto 30
        enddo
         ifu=ifu+1
! Arrays with the ordering numbers of full arcs.
         fullarc(ifu)=neibp(j)
         al(ifu)=neibp(j)
        goto 31
30        continue
         ine=ine+1
         neib(ine)=neibp(j) !Array containing the intersection points.
31        continue

       enddo
       neib(0)=ine ! Number of IP.
       fullarc(0)=ifu ! Number of full arcs.
       al(0)=ifu

! Loop over full arcs.
       do k=1,ifu
         a=neibor(fullarc(k),1)
         if(a.lt.0.d0) then 
          aia=-1.d0
         else
          aia=1.d0
         endif
         b=neibor(fullarc(k),2)
         c=neibor(fullarc(k),3)
         d=neibor(fullarc(k),4)
         if(a.ne.0d0) then !True if rotation has taken place.
! Remove the part of the atomic surface covered cut by this full arc
          As(i)=As(i)-aia*R22p*(1d0+(-d/a-R42)*dabs(a)/dsqrt((R42*a-
     &                 d)**2+R42*(b*b+c*c)))
          dadx(1,1)=2d0*ddat(fullarc(k),2)
          dadx(1,2)=2d0*ddat(fullarc(k),3)
          dadx(1,3)=2d0*(ddat(fullarc(k),4)-R)
          dadx(2,1)=-8d0*pol(i)
          dadx(2,2)=0d0
          dadx(2,3)=0d0
          dadx(3,1)=0d0
          dadx(3,2)=dadx(2,1)
          dadx(3,3)=0d0
          dadx(4,1)=8d0*pol(i)*ddat(fullarc(k),2)
          dadx(4,2)=8d0*pol(i)*ddat(fullarc(k),3)
          dadx(4,3)=8d0*pol(i)*(ddat(fullarc(k),4)+R)
          div=R22*R42p/((R42*a-d)**2+R42*(b*b+c*c))**(1.5)
            gp(1)=(R42*(b*b+c*c-2d0*a*d)+2d0*d*d)*div
          gp(2)=-b*(d+R42*a)*div
          gp(3)=-c*(d+R42*a)*div
          gp(4)=(b*b+c*c-2d0*a*d+(R42+R42)*a*a)*div
          grad(i,fullarc(k),1)=gp(1)*dadx(1,1)+gp(2)*dadx(2,1)+
     &        gp(3)*dadx(3,1)+gp(4)*dadx(4,1)
          grad(i,fullarc(k),2)=gp(1)*dadx(1,2)+gp(2)*dadx(2,2)+
     &        gp(3)*dadx(3,2)+gp(4)*dadx(4,2)
          grad(i,fullarc(k),3)=gp(1)*dadx(1,3)+gp(2)*dadx(2,3)+
     &        gp(3)*dadx(3,3)+gp(4)*dadx(4,3)

         else
          As(i)=As(i)+R22p*d/dsqrt(R42*(b*b+c*c)+d*d)
         endif
       enddo

! Loop over the circles with intersection points.
       do k=1,ine
         jk=neib(k) ! The order number of kth neighbour.
         a=neibor(jk,1) ! a>0 - outer part, a<0-inner part
         ak=a*a
         daba=dabs(a)
         if(a.lt.0d0) then 
          aia=-1d0
          else
          aia=1d0
         endif
         b=neibor(jk,2)
         c=neibor(jk,3)
         d=neibor(jk,4)
         b2c2=b*b+c*c
             vv=dsqrt((R42*a-d)**2+R42*b2c2)
             vv1=1d0/vv
             vv2=vv1*vv1
             wv=dsqrt(b2c2-4d0*a*d)
         nyx=0 ! The number of IP's which lie on this circle.
         ivr=0
         do j=1,iv
          if(vertex(j,3).eq.jk) then
           nyx=nyx+1
           ayx(nyx,1)=vertex(j,1) ! t coordinate of IP
           ayx(nyx,2)=vertex(j,2) ! s coordinate
           ivr=ivr+1
           ivrx(ivr)=j
          endif
          if(vertex(j,4).eq.jk) then
           nyx=nyx+1
           ayx(nyx,1)=vertex(j,1) ! t coordinate
           ayx(nyx,2)=vertex(j,2) ! s coordinate
           aa=vertex(j,3)
           vertex(j,3)=vertex(j,4)! make the fixed ordering number as the third
           vertex(j,4)=aa         ! component in vertex if it was the 4th.
           ivr=ivr+1
           ivrx(ivr)=j
          endif
         enddo

           ial=ifu
           do j=1,ine
           if(neib(j).ne.jk) then
            ial=ial+1
              al(ial)=neib(j) ! Add other neighbours with IP's to this array
           endif
          enddo
          al(0)=ial  ! The number of all "active" circles.

         if(a.ne.0.d0) then ! True after rotation
          aa=0.5d0/a
          ba=b*aa
          ca=c*aa
          bcd=b2c2-2d0*a*d
             xv=(d+R42*a)*vv1
             yv=(bcd+2d0*R42*ak)/daba*vv1
             zv=wv/a*vv1
! Move the (t,s) origo into the centre of this circle.
! Modify t and s coordinates of IP.
          do j=1,nyx
           ayx(j,1)=ayx(j,1)+ba ! t
           ayx(j,2)=ayx(j,2)+ca ! s
          enddo

          ct=-ba
          cs=-ca
          rr=0.5d0*wv/daba
! Calculate the polar angle of IP.
          do j=1,nyx
           ayx1(j)=datan2(ayx(j,2),ayx(j,1))
          enddo
! Sort the angles in ascending order.
          do j=1,nyx-1
           a1=ayx1(j)
           jj=j
           do ij=j+1,nyx
            if(a1.gt.ayx1(ij)) then
             a1=ayx1(ij)
             jj=ij
            endif
           enddo
           if(jj.ne.j) then
            a1=ayx1(jj)
            ayx1(jj)=ayx1(j)
            ayx1(j)=a1
            ivr=ivrx(jj)
            ivrx(jj)=ivrx(j)
            ivrx(j)=ivr
           endif
          enddo

          ayx1(nyx+1)=ayx1(1)+pi2

         do j=1,nyx
           do jkk=1,4
             vrx(j,jkk)=vertex(ivrx(j),jkk) ! vrx-helping array
           enddo
         enddo

           do jkk=1,4
             vrx(nyx+1,jkk)=vrx(1,jkk)
           enddo

          do j=1,nyx

!  
!  Escape 'bad' intersections.
           if(dabs(ayx1(j+1)-ayx1(j)).lt.1d-8) goto 40

           prob=(ayx1(j)+ayx1(j+1))/2.d0
           ap1=ct+rr*dcos(prob) ! The middle point of the arc.
           ap2=cs+rr*dsin(prob)
! Verify if the middle point belongs to a covered part. 
! If yes then omit.
           do ij=1,ial
          if(neibor(al(ij),1)*(ap1*ap1+ap2*ap2)+neibor(al(ij),2)*ap1+
     &         neibor(al(ij),3)*ap2+neibor(al(ij),4).lt.0d0) goto 40
           enddo
             alp=ayx1(j)
             bet=ayx1(j+1)
             dsalp=dsin(alp)
             dcalp=dcos(alp)
             dsbet=dsin(bet)
             dcbet=dcos(bet)
             dcbetmalp=dcos(5d-1*(bet-alp))
             dcbetpalp=dcos(5d-1*(bet+alp))
             dsbetmalp=dsin(5d-1*(bet-alp))
             dsbetpalp=dsin(5d-1*(bet+alp))
             dcotbma=1.d0/dtan(5d-1*(bet-alp))
             uv=daba*(bcd+2d0*R42*ak)*dcbetmalp
     &          -a*dsqrt(b2c2-4d0*a*d)*(b*dcbetpalp+c*dsbetpalp)
             As(i)=As(i)+R2*(aia*(alp-bet)+(d+R42*a)*vv1*(pi-2d0*
     &             datan(uv/(2d0*ak*vv*dsbetmalp))))
             kk=vrx(j,4)
             LL=vrx(j+1,4)
             tt=yv*dcotbma-zv*(b*dcbetpalp+c*dsbetpalp)/dsbetmalp
             a1=neibor(kk,1)
             b1=neibor(kk,2)
             c1=neibor(kk,3)
             d1=neibor(kk,4)
             a2=neibor(LL,1)
             b2=neibor(LL,2)
             c2=neibor(LL,3)
             d2=neibor(LL,4)
             ab1=a*b1-b*a1
             ac1=a*c1-c*a1
             ab2=a*b2-b*a2
             ac2=a*c2-c*a2
             am1=wv*aia*a1*(ab1*dsalp-ac1*dcalp)
             am2=wv*aia*a2*(ab2*dsbet-ac2*dcbet)
             dalp(1)=(b1*ab1+c1*ac1-2d0*d*a1*a1-a*
     &               (b1*b1+c1*c1-4d0*a1*d1)-2d0*d*aia*a1*(ab1*
     &               dcalp+ac1*dsalp)/wv+wv*aia*a1*
     &               (b1*dcalp+c1*dsalp))/am1
             dalp(2)=(-a1*ab1+b*a1*a1+b*aia*a1*(ab1*
     &               dcalp+ac1*dsalp)/wv-wv*aia*a1*a1*dcalp)/am1
             dalp(3)=(-a1*ac1+c*a1*a1+c*aia*a1*(ab1*
     &               dcalp+ac1*dsalp)/wv-wv*aia*a1*a1*dsalp)/am1
             dalp(4)=(-2d0*a*a1*a1-2d0*daba*a1*(ab1*
     &               dcalp+ac1*dsalp)/wv)/am1
             dbet(1)=(b2*ab2+c2*ac2-2d0*d*a2*a2-a*
     &               (b2*b2+c2*c2-4d0*a2*d2)-2d0*d*aia*a2*(ab2*
     &               dcbet+ac2*dsbet)/wv+wv*aia*a2*
     &               (b2*dcbet+c2*dsbet))/am2
             dbet(2)=(-a2*ab2+b*a2*a2+b*aia*a2*(ab2*
     &               dcbet+ac2*dsbet)/wv-wv*aia*a2*a2*dcbet)/am2
             dbet(3)=(-a2*ac2+c*a2*a2+c*aia*a2*(ab2*
     &               dcbet+ac2*dsbet)/wv-wv*aia*a2*a2*dsbet)/am2
             dbet(4)=(-2d0*a*a2*a2-2d0*daba*a2*(ab2*
     &               dcbet+ac2*dsbet)/wv)/am2
             daalp(1)=(-b*ab1-c*ac1+wv*wv*a1+2d0*ak*d1+
     &              wv*aia*((ab1-b*a1)*dcalp+(ac1-c*a1)*dsalp))/am1
             daalp(2)=(a*ab1-ak*b1+wv*daba*a1*dcalp)/am1
             daalp(3)=(a*ac1-ak*c1+wv*daba*a1*dsalp)/am1
             daalp(4)=(2d0*ak*a1)/am1
             dabet(1)=(-b*ab2-c*ac2+wv*wv*a2+2d0*ak*d2+
     &             wv*aia*((ab2-b*a2)*dcbet+(ac2-c*a2)*dsbet))/am2
             dabet(2)=(a*ab2-ak*b2+wv*daba*a2*dcbet)/am2
             dabet(3)=(a*ac2-ak*c2+wv*daba*a2*dsbet)/am2
             dabet(4)=(2d0*ak*a2)/am2
               dv(2)=R42*b*vv1      
               dv(3)=R42*c*vv1      
               dv(4)=(d-R42*a)*vv1      
               dv(1)=-dv(4)*R42      
               dx(1)=R42*vv1-(d+R42*a)*dv(1)*vv2      
             dx(2)=-(d+R42*a)*dv(2)*vv2      
             dx(3)=-(d+R42*a)*dv(3)*vv2      
             dx(4)=1d0*vv1-(d+R42*a)*dv(4)*vv2      
             dy(1)=(-2d0*d+4d0*R42*a)/daba*vv1-(b*b+c*c-2d0*a*d+2d0*
     &             R42*ak)*(aia*vv+daba*dv(1))/ak*vv2
               dy(2)=2d0*b/daba*vv1-(b*b+c*c-2d0*a*d+2d0*
     &             R42*ak)*aia*dv(2)/a*vv2
               dy(3)=2d0*c/daba*vv1-(b*b+c*c-2d0*a*d+2d0*
     &             R42*ak)*aia*dv(3)/a*vv2
               dy(4)=-2.d0*a/daba*vv1-(b*b+c*c-2.d0*a*d+2.d0*
     &             R42*ak)*aia*dv(4)/a*vv2      
             dz(1)=-2d0*d/wv/a*vv1-wv*(vv+a*dv(1))/ak*vv2
             dz(2)=b/wv/a*vv1-wv*dv(2)/a*vv2
             dz(3)=c/wv/a*vv1-wv*dv(3)/a*vv2
             dz(4)=-2d0/wv*vv1-wv*dv(4)/a*vv2
             dt(1)=dcotbma*dy(1)-5d-1*yv*
     &               (dbet(1)-dalp(1))/dsbetmalp**2-((b*dcbetpalp+
     &             c*dsbetpalp)/dsbetmalp)*dz(1)-
     &             5d-1*zv*(((-b*dsbetpalp+c*dcbetpalp)*dsbetmalp)*
     &             (dalp(1)+dbet(1))-((b*dcbetpalp+c*dsbetpalp)*
     &             dcbetmalp)*(dbet(1)-dalp(1)))/dsbetmalp**2
             dt(2)=dcotbma*dy(2)-5d-1*yv*
     &             (dbet(2)-dalp(2))/dsbetmalp**2-(b*dcbetpalp+
     &             c*dsbetpalp)/dsbetmalp*dz(2)-
     &             5d-1*zv*((-b*dsbetpalp+c*dcbetpalp)*dsbetmalp*
     &             (dalp(2)+dbet(2))-(b*dcbetpalp+c*dsbetpalp)*
     &             dcbetmalp*(dbet(2)-dalp(2))+dcbetpalp
     &             *2d0*dsbetmalp)/dsbetmalp**2
             dt(3)=dcotbma*dy(3)-5d-1*yv*
     &             (dbet(3)-dalp(3))/dsbetmalp**2-(b*dcbetpalp+
     &             c*dsbetpalp)/dsbetmalp*dz(3)-
     &             5d-1*zv*((-b*dsbetpalp+c*dcbetpalp)*dsbetmalp*
     &             (dalp(3)+dbet(3))-(b*dcbetpalp+c*dsbetpalp)*
     &             dcbetmalp*(dbet(3)-dalp(3))+dsbetpalp*2d0*
     &             dsbetmalp)/dsbetmalp**2
             dt(4)=dcotbma*dy(4)-5d-1*yv*
     &             (dbet(4)-dalp(4))/dsbetmalp**2-(b*dcbetpalp+
     &             c*dsbetpalp)/dsbetmalp*dz(4)-
     &             5d-1*zv*((-b*dsbetpalp+c*dcbetpalp)*dsbetmalp*
     &             (dalp(4)+dbet(4))-(b*dcbetpalp+c*dsbetpalp)*
     &             dcbetmalp*(dbet(4)-dalp(4)))/dsbetmalp**2
             di(1)=R2*(aia*(dalp(1)-dbet(1))+(pi-2d0*datan(5d-1*tt))*
     &               dx(1)-4d0*xv*dt(1)/(4d0+tt**2))
             di(2)=R2*(aia*(dalp(2)-dbet(2))+(pi-2d0*datan(5d-1*tt))*
     &               dx(2)-4d0*xv*dt(2)/(4d0+tt**2))
             di(3)=R2*(aia*(dalp(3)-dbet(3))+(pi-2d0*datan(5d-1*tt))*
     &               dx(3)-4d0*xv*dt(3)/(4d0+tt**2))
             di(4)=R2*(aia*(dalp(4)-dbet(4))+(pi-2d0*datan(5d-1*tt))*
     &               dx(4)-4d0*xv*dt(4)/(4d0+tt**2))
             dii(1,1)=2d0*ddat(jk,2)
             dii(1,2)=2d0*ddat(jk,3)
             dii(1,3)=2d0*(ddat(jk,4)-R)
             dii(2,1)=-R42-R42
             dii(2,2)=0d0
             dii(2,3)=0d0
             dii(3,1)=0d0
             dii(3,2)=dii(2,1)
             dii(3,3)=0d0
             dii(4,1)=dii(1,1)*R42
             dii(4,2)=dii(1,2)*R42
             dii(4,3)=2d0*(ddat(jk,4)+R)*R42
             
             do idi=1,3
              grad(i,jk,idi)=grad(i,jk,idi)+
     &                   di(1)*dii(1,idi)+di(2)*dii(2,idi)+
     &                   di(3)*dii(3,idi)+di(4)*dii(4,idi)

             enddo             
             do idi=1,4
               dta(idi)=5d-1*daalp(idi)*(((yv-zv*((-b*dsbetpalp
     &                  +c*dcbetpalp)*dsbetmalp+
     &                (b*dcbetpalp+c*dsbetpalp)*
     &                dcbetmalp))/dsbetmalp**2))
               dtb(idi)=5d-1*dabet(idi)*(((-yv-zv*((-b*dsbetpalp
     &                  +c*dcbetpalp)*dsbetmalp-
     &                (b*dcbetpalp+c*dsbetpalp)*
     &                dcbetmalp))/dsbetmalp**2))
             di1(idi)=R2*(aia*daalp(idi)-4d0*xv*dta(idi)/(4d0+tt*tt))
            di2(idi)=R2*(-aia*dabet(idi)-4d0*xv*dtb(idi)/(4d0+tt*tt))
             enddo             

             dii(1,1)=2d0*ddat(kk,2)
             dii(1,2)=2d0*ddat(kk,3)
             dii(1,3)=2d0*(ddat(kk,4)-R)
             dii(4,1)=dii(1,1)*R42
             dii(4,2)=dii(1,2)*R42
             dii(4,3)=2d0*(ddat(kk,4)+R)*R42
             do idi=1,3
              grad(i,kk,idi)=grad(i,kk,idi)+
     &                   di1(1)*dii(1,idi)+di1(2)*dii(2,idi)+
     &                       di1(3)*dii(3,idi)+di1(4)*dii(4,idi)
             enddo

             dii(1,1)=2d0*ddat(LL,2)
             dii(1,2)=2d0*ddat(LL,3)
             dii(1,3)=2d0*(ddat(LL,4)-R)
             dii(4,1)=dii(1,1)*R42
             dii(4,2)=dii(1,2)*R42
             dii(4,3)=2d0*(ddat(LL,4)+R)*R42
             do idi=1,3
              grad(i,LL,idi)=grad(i,LL,idi)+
     &                   di2(1)*dii(1,idi)+di2(2)*dii(2,idi)+
     &                       di2(3)*dii(3,idi)+di2(4)*dii(4,idi)
             enddo

40           continue
          enddo
         
         else ! analyse the case with lines (not necessary after rotation).

           ang=datan2(b,c)
           a1=dsin(ang)
           a2=dcos(ang)
! Rotate the intersection points by the angle "ang".
! After rotation the position of IP will be defined
! by its first coordinate "ax(.,1)".
           do j=1,nyx
            ax(j,1)=ayx(j,1)*a2-ayx(j,2)*a1
            ax(j,2)=ayx(j,1)*a1+ayx(j,2)*a2
           enddo
! Sorting by acsending order.
           do j=1,nyx-1
            a1=ax(j,1)
            jj=j
            do L=j+1,nyx
             if(a1.gt.ax(L,1)) then
              a1=ax(L,1)
              jj=L
             endif
            enddo
              if(jj.ne.j) then
             a1=ax(jj,1)
             a2=ax(jj,2)
             ax(jj,1)=ax(j,1)
             ax(jj,2)=ax(j,2)
             ax(j,1)=a1
             ax(j,2)=a2
             a1=ayx(jj,1)
             a2=ayx(jj,2)
             ayx(jj,1)=ayx(j,1)
             ayx(jj,2)=ayx(j,2)
             ayx(j,1)=a1
             ayx(j,2)=a2
            endif
           enddo

           if(c.ne.0d0) then
            do j=1,nyx
             ayx1(j)=(ayx(j,1)-ayx(1,1))/c
            enddo
           else
            do j=1,nyx
             ayx1(j)=-(ayx(j,2)-ayx(1,2))/b
            enddo
           endif

           probe(1)=-1d0
           do j=1,nyx-1
            probe(j+1)=(ayx1(j)+ayx1(j+1))/2d0
           enddo
           probe(nyx+1)=ayx1(nyx)+1d0
           

           bc=b*b+c*c
           cc=dsqrt(R42*(bc)+d*d)
           caba=c*ayx(1,1)-b*ayx(1,2)

           p1=ayx(1,1)+c*probe(1)
           p2=ayx(1,2)-b*probe(1)
! Verify if the middle point belongs to a covered part.
! Omit, if yes.
           do L=1,ial
            if(neibor(al(L),1)*(p1*p1+p2*p2)+neibor(al(L),2)*p1+
     &         neibor(al(L),3)*p2+neibor(al(L),4).lt.0d0) goto 20
           enddo
           arg2=bc*ayx1(1)+caba
           As(i)=As(i)+R22*d*(datan(arg2/cc)+pi/2d0)/cc
20           continue

           do j=2,nyx
            p1=ayx(1,1)+c*probe(j)
            p2=ayx(1,2)-b*probe(j)
            do L=1,ial
             if(neibor(al(L),1)*(p1*p1+p2*p2)+neibor(al(L),2)*p1+
     &          neibor(al(L),3)*p2+neibor(al(L),4).lt.0d0) goto 21
            enddo
            arg1=bc*ayx1(j-1)+caba
            arg2=bc*ayx1(j)+caba
            As(i)=As(i)+R22*d*(datan(arg2/cc)-datan(arg1/cc))/cc
21            continue
           enddo
           p1=ayx(1,1)+c*probe(nyx+1)
           p2=ayx(1,2)-b*probe(nyx+1)
           do      L=1,ial
            if(neibor(al(L),1)*(p1*p1+p2*p2)+neibor(al(L),2)*p1+
     &         neibor(al(L),3)*p2+neibor(al(L),4).lt.0d0) goto 22
           enddo
           arg1=bc*ayx1(nyx)+caba
           As(i)=As(i)+R22*d*(pi/2d0-datan(arg1/cc))/cc
22           continue
            endif
         enddo       
11      continue

!
!         changed in 2004

      do j=1,ine
      al(ifu+j)=neib(j)
      enddo
      al(0)=ifu+ine

      do k=1,3
      ss(k)=0d0
       do j2=1,al(0)
         j=al(j2)
         ss(k)=ss(k)+grad(i,j,k)
       enddo
      grad(i,i,k)=-ss(k)
      enddo

      if(jjj.ne.0) then
! Restore the original configuration if the molecule has been rotated.
! This is necessary for calculation of gradients.
          xx=grad(i,i,1)
          yy=grad(i,i,2)
          zz=grad(i,i,3)
          grad(i,i,1)=xx*cg-yy*sfsg+zz*cfsg
          grad(i,i,2)=yy*cf+zz*sf
          grad(i,i,3)=-xx*sg-yy*sfcg+zz*cfcg

        do j2=1,al(0)
          j=al(j2)
          xx=grad(i,j,1)
          yy=grad(i,j,2)
          zz=grad(i,j,3)
          grad(i,j,1)=xx*cg-yy*sfsg+zz*cfsg
          grad(i,j,2)=yy*cf+zz*sf
          grad(i,j,3)=-xx*sg-yy*sfcg+zz*cfcg
        enddo

      endif


1     sss=sss+As(i)*sigma(i)

      do j=1,numat
         do k=1,3
          gs(k)=0.d0
            do i=1,numat
               gs(k)=gs(k)+sigma(i)*grad(i,j,k)
            enddo
          gradan(j,k)=gs(k)
         enddo
      enddo

111   continue

!     write(3,*)'   No   Area    sigma   Enrg    gradx   grady',     
!    & '   gradz   Rad    Atom'
!     write(3,*)

      
      do i=1,numat
        enr=As(i)*sigma(i)
        energy=energy+enr
        ratom=rvdw(i)-rwater

        if (nmat(i)(1:1).eq.'h') ratom=0.0

!     write(3,700),i,As(i),sigma(i),
!    &        enr,(gradan(i,k),k=1,3),ratom,nmat(i)

      enddo

!     write(3,*)'Total Area: ',sss
!     write(3,*)'Total solvation energy: ',energy

      esolan=sss
      
      deallocate(grad)
!     close(3)
      return

100   format(20i4)
200   format(2i4,3f16.6)
700   format(i5,7f8.3,5x,a4)
         
      end
