! **************************************************************
!
! This file contains the subroutines: enyshe 
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************


      real*8 function enyshe(nml)

! ............................................................................
!
! PURPOSE: Calculate internal energy of molecule 'nml' with ECEPP parameters
!
! CALLS: none
!
! The function loops over all moving sets within the molecule. Within
! this loop it loops over the van-der-Waals domains of each atom in the 
! moving set and finally over the atoms that belong to the 1-4 interaction
! set.
! ............................................................................

      include 'INCL.H'

! If nml == 0 calculate the interaction between all pairs.
      if (nml.eq.0) then
          ntlvr = nvr
      else
        ntlvr=nvrml(nml)
      endif
      
      if (ntlvr.eq.0) then
        write (*,'(a,i4)')
     &           ' enyshe> No variables defined in molecule #',nml
        return
      endif

      enyshe=0.0
      eyel=0.0
      eyvw=0.0
      eyhb=0.0
      eyvr=0.0
      if (nml.eq.0) then
        ifivr = ivrml1(1)
        i1s = imsml1(ntlml) + nmsml(ntlml)
      else 
! Index of first variable in molecule.
        ifivr=ivrml1(nml)
! Index of last moving set in molecule
        i1s=imsml1(nml)+nmsml(nml)
      endif
! Loop over moving sets/variables in reverse order     
      do io=ifivr+ntlvr-1,ifivr,-1  
! The array iorvr contains the variables in an "apropriate" order.
        iv=iorvr(io)       
! Index of the primary moving atom for the variable with index iv
        ia=iatvr(iv)       
! Get the type of variable iv (valence length, valence angle, dihedral angle)
        it=ityvr(iv)       
! Class of variable iv's potential  (Q: What are they)
        ic=iclvr(iv)       
! If iv is a dihedral angle ...
        if (it.eq.3) then      
! Barrier height * 1/2 of the potential of iv.
          e0=e0to(ic)
! Calculate the periodic potential term. sgto is the sign of the barrier, rnto is
! the periodicity and toat is torsion angle(?) associate with atom ia.
          if (e0.ne.0.) 
     &         eyvr=eyvr+e0*(1.0+sgto(ic)*cos(toat(ia)*rnto(ic)))
! else if iv is a valence angle ...
        elseif (it.eq.2) then  
! vr is the valence angle of ia
          vr=baat(ia)
! else if iv is a valence length...
        elseif (it.eq.1) then  
! vr is the length of the valence bond
          vr=blat(ia)
        endif

! ============================================ Energies & Atomic forces
! index of next to last moving set
        i2s=i1s-1
! index of first moving set associated with iv
        i1s=imsvr1(iv) 
! Loop over all moving sets starting from the one associated with vr to the end.
        do ims=i1s,i2s  
! First atom of the current moving set
          i1=latms1(ims)
! Last atom of the current moving set
          i2=latms2(ims)
! Loop over all atoms of the current moving set.
          do i=i1,i2  
! Atom class of current atom
            ity=ityat(i)
! Point charge at current atom
            cqi=conv*cgat(i)
! Cartesian coordinates of current atom
            xi=xat(i)
            yi=yat(i)
            zi=zat(i)
! Loop over the atoms of the van der Waals domain belonging to atom i
            do ivw=ivwat1(i),ivwat2(i)  
! Loop over the atoms of the van der Waals domain of the atoms of the 
! van der Waals domain of atom i
! Q: Which atoms are in these domains?
              do j=lvwat1(ivw),lvwat2(ivw)  
! Atom type of partner
                jty=ityat(j)
! Differences in cartesian coordinates
                xij=xat(j)-xi
                yij=yat(j)-yi
                zij=zat(j)-zi
! Cartesian distance and higher powers
                rij2=xij*xij+yij*yij+zij*zij
                rij4=rij2*rij2
                rij6=rij4*rij2
                rij=sqrt(rij2)
! Are we using a distance dependent dielectric constant?
                if(epsd) then
                 sr=slp*rij
                 ep=plt-(sr*sr+2.0*sr+2.0)*(plt-1.0)*exp(-sr)/2.0
                else
                 ep = 1.0d0
                end if
! Coulomb interaction
                eyel=eyel+cqi*cgat(j)/(rij*ep)
! If the two atoms cannot form a hydrogen bond use 6-12 Lennard-Jones potential
                if (ihbty(ity,jty).eq.0) then
                  eyvw=eyvw+aij(ity,jty)/(rij6*rij6)
     &                     -cij(ity,jty)/rij6
                else
! For hydrogen bonding use 10-12 Lennard-Jones potential
                  eyhb=eyhb+ahb(ity,jty)/(rij6*rij6)
     &                     -chb(ity,jty)/(rij6*rij4)
                endif

              enddo  
            enddo  
            
! Loop over 1-4 interaction partners
! The interactions between atoms that are three bonds apart in the protein are
! dominated by quantum mechanical effects. They are treated separately.
            do i14=i14at1(i),i14at2(i)   
              j=l14at(i14)

              jty=ityat(j)

              xij=xat(j)-xi
              yij=yat(j)-yi
              zij=zat(j)-zi
              rij2=xij*xij+yij*yij+zij*zij
              rij4=rij2*rij2
              rij6=rij4*rij2
              rij = sqrt(rij2)
! Are we using a distance dependent dielectric constant?
              if(epsd) then
               sr=slp*rij
               ep=plt-(sr*sr+2.0*sr+2.0)*(plt-1.0)*exp(-sr)/2.0
              else
               ep=1.0d0
              end if

              eyel=eyel+cqi*cgat(j)/(rij*ep)
! If hydrogen bonding is not possible use 6-12 Lennard-Jones potential.
              if (ihbty(ity,jty).eq.0) then
                eyvw=eyvw+a14(ity,jty)/(rij6*rij6)
     &                   -cij(ity,jty)/rij6
              else
! Use 10-12 Lennard-Jones potential for hydrogen bonds.
                eyhb=eyhb+ahb(ity,jty)/(rij6*rij6)
     &                   -chb(ity,jty)/(rij6*rij4)
              endif

            enddo  ! ... 1-4-partners of i

          enddo  ! ... atoms i
        enddo  ! ... m.s.

      enddo  ! ... variables

      eysm = eyel + eyvw + eyhb + eyvr

      enyshe=eysm
      return
      end

