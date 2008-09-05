! *********************************************************************
! This file contains eninteract
!
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
!
! Description: Calculates the interaction energy between molecules
! The function assumes that all molecules are up-to-date. If in doubt
! call energy first.
! The energy function is based on the ECEPP/3 dataset.
!
! TODO: Intermolecular interaction energy for FLEX and ECEPP/2
      real*8 function eninteract()

        include 'INCL.H'
                
        eysmi=0.0
        eyeli=0.0
        eyvwi=0.0
        eyhbi=0.0

        do iml = 1, ntlml
          do jml = iml + 1, ntlml
            do ires= irsml1(iml), irsml2(iml)
              do jres= irsml1(jml), irsml2(jml)
                do iat = iatrs1(ires), iatrs2(ires)
! Atom class of current atom
                  ity=ityat(iat)
! Point charge at current atom
                  cqi=conv*cgat(iat)
! Cartesian coordinates of current atom
                  xi=xat(iat)
                  yi=yat(iat)
                  zi=zat(iat)

                  do jat = iatrs1(jres), iatrs2(jres)
! Atom type of partner
                    jty=ityat(jat)
! Differences in cartesian coordinates
                    xj=xat(jat)
                    yj=yat(jat)
                    zj=zat(jat)
                    
                    xij=xat(jat)-xi
                    yij=yat(jat)-yi
                    zij=zat(jat)-zi
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
                    eyeli=eyeli+cqi*cgat(jat)/(rij*ep)
! If the two atoms cannot form a hydrogen bond use 6-12 Lennard-Jones potential
                    if (ihbty(ity,jty).eq.0) then
                      eyvwi=eyvwi+aij(ity,jty)/(rij6*rij6)
     &                          -cij(ity,jty)/rij6
                    else
! For hydrogen bonding use 10-12 Lennard-Jones potential
                      eyhbi=eyhbi+ahb(ity,jty)/(rij6*rij6)
     &                          -chb(ity,jty)/(rij6*rij4)
                    endif
                  enddo 
                enddo
              enddo
            enddo
          enddo
        enddo

        eysmi = eyeli + eyvwi + eyhbi
        eninteract = eysmi
        return 
      end   