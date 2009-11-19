! ......................
! contents of a PDB file
! ......................
      integer MXCHP, MXRSP, MXATP, nchp, nchrsp, nrsp, irsatp, nrsatp
      integer natp, noatp
      
      double precision xatp, yatp, zatp

      parameter (MXCHP =100,    ! max. no. of polypeptide chains
     &           MXRSP =1000,   ! max. no. of residues
     &           MXATP =10000)  ! max. no. of atoms

!   nchp      - no. of polypeptide chains
!   nchrsp()  - no. of residues / chain
!   nrsp      - total no. of residues
!   irsatp()  - index of 1st atom / res.
!   nrsatp()  - no. of atoms / res.
!   natp      - total no. of selected atoms
!   noatp()   - atom numbers

      common /pdb_i/ nchp,nchrsp(MXCHP),
     &               nrsp,irsatp(MXRSP),nrsatp(MXRSP),
     &               natp,noatp(MXATP)
      save /pdb_i/

!   chnp()  - chain identifiers
!   rsidp() - residue identifiers (number + insertion code)
!   rsnmp() - residues (sequence, 3-letter code)
!   atnmp() - atom names

      character chnp(MXCHP),
     &          rsidp(MXRSP)*5,rsnmp(MXRSP)*3,
     &          atnmp(MXATP)*4

      common /pdb_c/ chnp,rsnmp,rsidp,atnmp
      save /pdb_c/

!     xatp,yatp,zatp - atom coordinates

      common /pdb_r/ xatp(MXATP),yatp(MXATP),zatp(MXATP)
      save /pdb_r/

! ------------------- code to list all PDB information
!      ir=0
!      do i=1,nchp
!        write(*,*) ' ===== chain |',chnp(i),'|'
!        do j=1,nchrsp(i)
!          ir=ir+1
!          write(*,*) ' ----- ',rsidp(ir),' ',rsnmp(ir),' ',nrsatp(ir)
!          k1=irsatp(ir)
!          k2=k1+nrsatp(ir)-1
!          do k=k1,k2
!            write(*,*) ' ',noatp(k),' ',atnmp(k),' ',(xyzp(l,k),l=1,3)
!          enddo
!        enddo
!      enddo
