! Subroutine initlund: Initializes data structures used frequently 
! in connection with Biased Gaussian Steps and the energy functions
! from Anders Irback's protein folding model. Calls: none.
!
subroutine init_lund
    include 'INCL.H'
    include 'incl_lund.h'
    integer i, npprs, j, k

    logical bgsposs
    do i=1,mxrs
        iN(i)=-1
        iCa(i)=-1
        iC(i)=-1
        iphi(i)=-34
        ipsi(i)=-35
    enddo
!      print *,'total number of variables = ',nvr

    do i=1,ntlml
        npprs=1
        do j=ivrml1(i),ivrml1(i)+nvrml(i)-1
            mlvr(j)=i
            if (nmvr(j).eq.'phi') then 
                iphi(npprs)=j
! Now if the residue is a proline, there is no phi angle in the variable
! list in SMMP, and so iphi(j) will remain at the initial value. 
! So, make sure you never use iphi(i) for proline. 
            endif
            if (nmvr(j).eq.'psi'.or.nmvr(j).eq.'pst') then 
                ipsi(npprs)=j
                npprs=npprs+1
            endif
        enddo
        do j=irsml1(i),irsml2(i)
            iN(j)=iatrs1(j) 
            do k=iatrs1(j),iatrs2(j)
            if (nmat(k)(1:2).eq.'ca') then 
                iCa(j)=k 
            endif
            if (nmat(k)(1:1).eq.'c') then 
                iC(j)=k 
            endif 
            enddo
    !            print *,'determined phi,psi serial for residue',j,' as '
    !     #           ,iphi(j),ipsi(j)
        enddo
    enddo
    abgs=300.0
    bbgs=10.0
    bgsnvar=0
! JHM: Took the following lines out to avoid another dependency.
!     do i=1,nvr
!         if (bgsposs(i)) then
!             bgsnvar=bgsnvar+1
!             bgsvar(bgsnvar)=i
!         endif
!     enddo
end subroutine init_lund

