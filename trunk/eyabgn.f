! *********************************************************************
! This file contains eyrccr, init_abgn, eyentr, eyabgn
!
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! Corrections to ECEPP energy terms due to R. A. Abagyan et al.
!
! Two terms are calculated: eyrccr and eyentr, representing respectively
! c a term to slightly shift the backbone dihedral angle preferences in
! the ECEPP potential slightly away from the helix region, and another
! term to estimate the side-chain entropy from a given configuration.
!
!
! *********************************************************************
      real*8 function eyrccr(nml)
      include 'INCL.H'
      double precision et, rsscl

      integer nml, istres, indres, i, ipsi, in, ica, ic, mlvr, iphi

      dimension iN(mxrs),iCa(mxrs),iC(mxrs),mlvr(mxvr)
      dimension iphi(mxrs),ipsi(mxrs)
      common /lundds/iN,iCa,iC,mlvr,iphi,ipsi
      character mynm*4

      if (nml.eq.0) then
         istres=irsml1(1)
         indres=irsml2(ntlml)
      else
         istres=irsml1(nml)
         indres=irsml2(nml)
      endif
      et=0.0
!      print *,'***********'
      do i=istres,indres
         mynm=seq(i)
         call tolost(mynm)
         if ((mynm.eq.'val').or.(mynm.eq.'ile').or.
     &           (mynm.eq.'thr')) then
            rsscl=1.0
         else
            rsscl=0.5
         endif
         et=et+rsscl*(1.0-sin(vlvr(ipsi(i))))
!         print *,'  contribution = ',rsscl*(1.0-sin(vlvr(ipsi(i))))
!         print *,'obtained using scale ',rsscl,' and angle ',
!     #        vlvr(ipsi(i))
      enddo
!      print *,'abagyan dihedral term = ',et
!      print *,'***********'
      eyrccr=et
      return
      end

      subroutine init_abgn
      include 'INCL.H'
!       dimension rsstrg(mxrs)
!       common /abgncor/rsstrg
      double precision xarea, estrg

      integer i, istres, indres, imytyp, j

      dimension xarea(nrsty),estrg(nrsty)
      character mynm*4
!      print *,'Initialization of Abagyan entropic term'
!     Maximum accessible surface areas for different residue types
      data (xarea(i),i=1,nrsty)/
!             1         2         3         4         5
     &     117.417 , 244.686 , 245.582 , 146.467 , 144.485 ,
!             6         7         8         9        10
     &     144.192 , 142.805 , 147.568 , 183.103 , 177.094 ,
!            11       12        13        14        15
     &     186.293 , 83.782 , 187.864 , 187.864 , 187.864 ,
!            16       17        18        19        20
     &     187.864 , 160.887 , 161.741 , 184.644 , 179.334 ,
!            21       22        23        24        25
     &     209.276 , 209.276 , 203.148 , 208.902 , 153.124 ,
!            26       27        28        29        30
     &     153.973 , 153.037 , 158.695 , 157.504 , 157.504 ,
!            31       32        33        34        35
     &     119.786 , 146.488 , 238.641 , 223.299 , 160.283 /
!     Entropic contribution for maximally exposed residue
      data (estrg(i),i=1,nrsty)/
!             1ala      2arg      3arg+     4asn      5asp
     &       0.0   ,   2.13  ,   2.13  ,   0.81  ,   0.61  ,
!             6asp-     7cys      8cyss     9gln     10glu
     &       0.61  ,   1.14  ,   1.14  ,   2.02  ,   1.65  ,
!            11glu-   12gly     13his     14hise    15hisd
     &       1.65  ,  0.0    ,   0.99  ,   0.99  ,   0.99  ,
!            16his+   17hyp     18hypu    19ile     20leu
     &       0.99  ,   0.99  ,   0.99  ,  0.75   ,   0.75  ,
!            21lys    22lys+     23met     24phe     25cpro
     &       2.21  ,   2.21  ,   1.53  ,   0.58  ,   0.0   ,
!            26pro     27cpru    28prou    29pron    30pro+
     &       0.0   ,   0.0   ,   0.0   ,   0.0   ,   0.0   ,
!            31ser     32thr     33trp     34tyr     35val
     &       1.19  ,   1.12  ,   0.97  ,   0.99  ,   0.50  /
      do i=1,mxrs
         rsstrg(i)=0.0
      enddo
      istres=irsml1(1)
      indres=irsml2(ntlml)
      do i=istres,indres
         imytyp=0
         mynm=seq(i)
         call tolost(mynm)
         do j=1,nrsty
            if (rsnmcd(j).eq.mynm) imytyp=j
!            print *,'comparing ',mynm,' with ',rsnmcd(j),imytyp
         enddo
         if (imytyp.eq.0) then
            print *, 'Unknown residue type ',seq(i)
            print *, 'Abagyan term strength set to 0'
            rsstrg(i)=0.0
         else
            rsstrg(i)=estrg(imytyp)/xarea(imytyp)
         endif
!         print *,'residue ',i,seq(i),' type ',imytyp
!         print *, 'strength for residue ',i,seq(i),' is ',rsstrg(i)
      enddo
      print *, 'initialized Abagyan corrections to ECEPP force field'
      end

      real*8 function eyentr(nml)
      include 'INCL.H'
!       dimension rsstrg(mxrs)
!       common /abgncor/rsstrg
!       common/ressurf/surfres(mxrs)
!      common/bet/beta
      double precision eentr, aars, strh

      integer nml, istres, indres, i

      eentr=0
      if (nml.eq.0) then
         istres=irsml1(1)
         indres=irsml2(ntlml)
      else
         istres=irsml1(nml)
         indres=irsml2(nml)
      endif
!      print *,'residue range ',istres,indres
!      print *,'for molecule ',nml
      do i=istres, indres
         aars=surfres(i)
         strh=rsstrg(i)
!        The maximal burial entropies were estimated at temperature 300k
!        The values in the array estrg are k_B * T (=300k) * Entropy
!        Presently we need it at temperature 1/beta, so we need to
!        multiply the strengths in estrg with (1/beta)/(300 kelvin)
!        300 kelvin is approximately 0.59576607 kcal/mol.
         eentr=eentr+aars*strh/(0.59576607*beta)
!         print *,'contribution = ',aars*strh/(0.59576607*beta)
!         print *,'residue, exposed area = ',i,aars
!         print *,'strength = ',strh,' for residue index = ',i
!         print *,'beta = ',beta
      enddo
!      print *,'abagyan entropic term = ',eentr
      eyentr=eentr
      return
      end

      real*8 function eyabgn(nml)
      include 'INCL.H'
      double precision eyrccr, eyentr

      integer nml

      eyabgn=eyrccr(nml)+eyentr(nml)
!      print *,'Abagyan term = ',eyabgn
      return
      end
