!**************************************************************
!     
! This file contains the subroutines: partem_p 
! USE WITH main_p, NOT WITH main!!!!!!
!     
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!     
!     **************************************************************

      subroutine  partem_p(num_rep, nequi, nswp, nmes, nsave, newsta,
     &                     switch, rep_id, partem_comm)
!     
!     PURPOSE: SIMULATION OF PROTEINS BY PARALLEL TEMPERING ALGORITHM 
!     ON PARALLEL COMPUTERS USING MPI
!     
!     switch: Choses the starting configuration:
!     -1 - stretched configuration
!     0 - don't change anything
!     1 - random start configuration
!     
!     CALLS:  addang,contacts,energy,hbond,helix,iendst,metropolis,
!     outvar,(rand),rgyr
!     
      include 'INCL.H'
      include 'INCP.H'
      include 'mpif.h'

      logical newsta
      integer switch, partem_comm, rep_id, nsave
!     external rand
      external can_weight

!     nequi:  number of Monte Carlo sweeps for thermalization
!     nswp:   number of Monte Carlo sweeps
!     nmes:   number of Monte Carlo sweeps between measurments
!     newsta: .true. for new simulations, .false. for re-start
      double precision temp, eavm, sph, geavm, gsph, dv, grnd, vr
      double precision addang, dummy, eol, energy, acz, rmsv, rmsdfun
      double precision rgy, ee, tmhb, dham, swp, wij, rd, e_final

      integer ifrrm, nmes, nswp, num_rep, i, j, nresi, iold, inode
      integer intem, iv, jold, idum1, idum2, idum3, mpi_integer
      integer mpi_comm_world, ierr, mpi_double_precision, nsw, nequi
      integer nml, nhel, mhel, nbet, mbet, mhb, imhb, nctot, ncnat
      integer mpi_comm_null, k1, k, nu, no1, in, jn
      
      dimension  eavm(MAX_PROC),sph(MAX_PROC),intem(MAX_PROC),
     &     inode(MAX_PROC), geavm(MAX_PROC), gsph(MAX_PROC)
      double precision    pbe(MAX_PROC),yol(MAX_PROC),acy(MAX_PROC),
     &     acy1(MAX_PROC),acx1(MAX_PROC),
     &     rgyrp(MAX_PROC),rmsdp(MAX_PROC), eol0,acz0
      
      double precision    e_min, e_minp(MAX_PROC), e_minpt(MAX_PROC)
      integer   h_max, h_maxp(MAX_PROC)
!     Order of replica exchange
      integer   odd
!     Counter to keep random number generators in sync
      integer randomCount
      
!     Collect partial energies. Only the root writes to disk. We have to
!     collect the information from the different replicas and provide 
!     arrays to store them.
!     eyslr    storage array for solvent energy
!     eyelp     -      "        - coulomb energy
!     eyvwp     -      "        - van-der-Waals energy
!     eyhbp     -      "        - hydrogen bonding energy
!     eysmi    -      "        - intermolecular interaction energy
!     eyabp     -      "        - Abagyan correction term
      double precision eyslr(MAX_PROC)
      double precision eyelp(MAX_PROC),eyvwp(MAX_PROC),eyhbp(MAX_PROC), 
     &     eyvrp(MAX_PROC),eysmip(MAX_PROC), eyabp(MAX_PROC)
!     Collect information about accessible surface and van-der-Waals volume
!     asap      storage array for solvent accessible surface
!     vdvolp     storage array for van-der-Waals volume
      double precision asa_p(MAX_PROC), vdvolp(MAX_PROC)

      integer nhelp(MAX_PROC),nbetp(MAX_PROC), mhbp(MAX_PROC),
     &     ncnatp(MAX_PROC),nctotp(MAX_PROC)
      integer imhbp(MAX_PROC)
      character*80 filebase, fileNameMP, tbase0,tbase1
!     frame     frame number for writing configurations
!     trackID   configuration that should be tracked and written out
!     dir          direction in random walk
!     -1 - visited highest temperature last
!     1 - visited lowest temperature last
!     0 - haven't visited the boundaries yet.
!     dirp      storage array for directions.
      integer frame, trackID, dir
      integer dirp(MAX_PROC)

      frame = ifrrm
      trackID = 1
      odd = 1
      write (*,*) 'Starting parallel tempering.'
      write (*,*) 'parameters, ',switch,newsta,nmes,nswp,nmes,
     &            rep_id, num_rep, partem_comm, myrank
      call flush(6)
!     
!     
!     File with temperatures 
      open(11,file='temperatures',status='old')
!     File with reference conformation
      tbase0='trj_00000'
      open(18,file=fileNameMP(tbase0,5,9,rep_id),status='unknown')
      if (rep_id.eq.0.and.myrank.eq.0) then
!     File with time series of simulation
         open(14,file='ts.d',status='unknown')
!     Track weights
!      open(16, file='weights.dat', status='unknown')
      endif
      
!     READ IN TEMPERATURES
      do i=1,num_rep
         read(11,*) j,temp
         pbe(j) = 1.0d0/( temp * 1.98773d-3 )
      end do
      close(11)

!     nresi:  number of residues
      nresi=irsml2(1)-irsml1(1)+1
!     
!     Initialize variables
      do i=1,num_rep      
         acx1(i) = 0.0d0
         acy(i) = 0.0d0
         eavm(i) = 0.0d0
         sph(i) = 0.0d0
         geavm(i) =0.0d0
         gsph(i) = 0.0d0
         e_minp(i) = 1.0d15
         h_maxp(i) = 0 
         dirp(i) = 0
      end do
      dirp(1) = 1
      dirp(num_rep) = -1
      e_min = 1.0d15
      h_max = 0
      dir = dirp(rep_id + 1)

!     _________________________________ Initialize Variables
      if(newsta) then
         iold=0
         do i=1,num_rep
            inode(i) = i
            intem(i) = i
         end do
!     _________________________________ initialize starting configuration
         if (switch.ne.0) then
            do i=1,nvr
               iv=idvr(i)       ! provides index of non-fixed variable
               if (switch.gt.0) then
                  dv=axvr(i)*(grnd()-0.5)
                  vr=addang(pi,dv)
               else
                  vr = pi 
               endif 
               vlvr(iv)=vr
            enddo
         endif
      else
         if(rep_id.eq.0.and.myrank.eq.0) then
            open(13,file='par_R.in', status='unknown')
            read(13,*) iold
            do i=1,num_rep
               read(13,*) j,inode(i),intem(i),yol(i),e_minp(i),h_maxp(i)
               write (*,*) "par_R.in:",i,j
            end do
            jold=(iold/nmes)*num_rep
            rewind 14
            do i=1,jold
               read(14,*) idum1,idum2,idum3,dummy, dir
     &              ,dummy, dummy, dummy, dummy, dummy
     &              ,dummy, dummy, dummy, dummy
     &              ,dummy, idum1, idum2, idum3
     &              ,idum1, idum2, idum3, e_min
     &              ,dummy, dummy
               write (*,*) i
               call flush(6)
            end do
            close(13)
         end if
         CALL MPI_BCAST(IOLD,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(INTEM,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(INODE,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YOL,num_rep,MPI_DOUBLE_PRECISION,0,
     &        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(E_MINP, num_rep, MPI_DOUBLE_PRECISION, 0, 
     &        MPI_COMM_WORLD, IERR)
         CALL MPI_BCAST(h_maxp,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,
     &        IERR)
      end if
      
      BETA = pbe(inode(rep_id+1))
      e_min = e_minp(rep_id+1)
      h_max = h_maxp(rep_id+1)
      write (*,*) "E_min=",e_min," for ", rep_id + 1 
      eol=energy()
      if(.not.newsta.and.abs(yol(rep_id + 1) - eol).gt.0.1) then
         write(*,*) rep_id, ' Warning: yol(rep_id).ne.eol:'
         write(*,*) rep_id, yol(rep_id + 1), eol
      endif
!     Start of simulation
      write (*,*) '[',rep_id, myrank, beta, partem_comm,
     &            '] Energy before equilibration:', eol
!     =====================Equilibration by canonical Metropolis
      do nsw=1,nequi
         call metropolis(eol,acz,can_weight)
      end do
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      write (*,*) '[',rep_id,'] Energy after equilibration:', eol
      call flush(6)
!     
!======================Multiple Markov Chains
      acz = 0
      do nsw=1,nswp
!------------First ordinary Metropolis 
         call metropolis(eol,acz,can_weight)
         iold = iold + 1	
         eol0 = eol
         if (myrank.eq.0.and.rep_id.eq.0) then
            write (*,*) "Finished sweep", nsw
            call flush(6)
         endif
         if(mod(iold,nmes).eq.0) then
            if ((rep_id + 1).eq.trackID.and.myrank.eq.0) then
               frame = iold /nmes
               filebase = "frame_00000.pdb"
               call outpdb(0, fileNameMP(filebase, 7, 11, frame))
            endif
            acz0 = acz
!     Evaluate RMSD
            nml = 1
            rmsv = rmsdfun(nml,irsml1(nml),irsml2(nml),ixatp,xatp,yatp, &
     &             zatp,0)
!            print *,myrank,'received RMSD,energy ',rmsv,eyab,beta
!     Measure global radius of gyration
            call rgyr(0,rgy,ee)  
            rgyp = rgy
!     Measure Helicity and Sheetness
            call helix(nhel,mhel,nbet,mbet)
!     Measure Number of hydrogen bonds
            mhb = 0
            do i = 1, ntlml
               call hbond(i,tmhb,-1) 
               mhb = mhb + 1
            enddo
            call interhbond(imhb) 
!     Measure total number of contacts (NCTOT) and number of
!     native contacts (NCNAT)
            call contacts(nctot,ncnat,dham)
!     Add tracking of lowest energy configuration
            if (eol.lt.e_min) then
!     Write out configuration
               i=rep_id+1
               j=inode(i)
               e_min = eol
               filebase = "c_emin_0000.pdb"
               call outpdb(0, fileNameMP(filebase, 8, 11, i))
               filebase = "c_emin_0000.var"
               call outvar(0, fileNameMP(filebase, 8, 11, i))
               filebase = "c_emin_0000.dat"
               open(15, file=fileNameMP(filebase, 8, 11, i),
     &              status="unknown")
!     write(15,'(i8,2i4,f6.2,2f8.2,5i8)') iold,i,j,pbe(i), 
               write(15,*) iold,j,i,beta, 
     &              eol, eyab, eysl, eyel, eyvw, eyhb, eyvr, eysmi,asa,
     &              vdvol, rgy, nhel, nbet, mhb, imhb, nctot,ncnat
               close(15)
            endif
!     Add tracking of configuration with larges hydrogen contents.
            if ((mhb + imhb).gt.h_max) then
!     Write out configuration
               i = rep_id + 1
               j = inode(i)
               h_max = mhb + imhb
               filebase = "c_hmax_0000.pdb"
               call outpdb(0,fileNameMP(filebase,8,11,i))
               filebase = "c_hmax_0000.var"
               call outvar(0,fileNameMP(filebase,8,11,i))
               filebase = "c_hmax_0000.dat"
               open(15, file=fileNameMP(filebase, 8, 11, i),
     &              status="unknown")
!     write(15,'(i8,2i4,f6.2,2f8.2,5i8)') iold,i,j,pbe(i), 
               write(15,*) iold,j,i,beta, 
     &              eol, eyab, eysl, eyel, eyvw, eyhb, eyvr, eysmi,asa,
     &              vdvol, rgy, nhel, nbet, mhb, imhb, nctot,ncnat
               close(15)
            endif

!     
!--------------------Gather measurement data
! I only use the master node of each replica for data collection. The
! variable partem_comm provides the appropriate communicator.
            if (partem_comm.ne.MPI_COMM_NULL) then
               CALL MPI_GATHER(rmsv,1,MPI_DOUBLE_PRECISION,rmsdp,1,
     &              MPI_DOUBLE_PRECISION, 0,partem_comm,IERR)
               CALL MPI_GATHER(eyab,1,MPI_DOUBLE_PRECISION,eyabp,1,
     &              MPI_DOUBLE_PRECISION, 0,partem_comm,IERR)
               CALL MPI_GATHER(RGYP,1,MPI_DOUBLE_PRECISION,RGYRP,1,
     &              MPI_DOUBLE_PRECISION, 0,partem_comm,IERR)
               CALL MPI_GATHER(NHEL,1,MPI_INTEGER,NHELP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(NBET,1,MPI_INTEGER,NBETP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(MHB,1,MPI_INTEGER,MHBP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(iMHB,1,MPI_INTEGER,iMHBP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(NCTOT,1,MPI_INTEGER,NCTOTP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(NCNAT,1,MPI_INTEGER,NCNATP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(dir,1,MPI_INTEGER,dirp,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(acz0,1,MPI_DOUBLE_PRECISION,acy1,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(e_min,1,MPI_DOUBLE_PRECISION,e_minp,1,
     &              MPI_DOUBLE_PRECISION,0, partem_comm,IERR)
               CALL MPI_GATHER(EOL0,1,MPI_DOUBLE_PRECISION,YOL,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eysl,1,MPI_DOUBLE_PRECISION,eyslr,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eyel,1,MPI_DOUBLE_PRECISION,eyelp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eyvw,1,MPI_DOUBLE_PRECISION,eyvwp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eyhb,1,MPI_DOUBLE_PRECISION,eyhbp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eyvr,1,MPI_DOUBLE_PRECISION,eyvrp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eysmi,1,MPI_DOUBLE_PRECISION,eysmip,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(asa,1,MPI_DOUBLE_PRECISION,asa_p,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(vdvol,1,MPI_DOUBLE_PRECISION,vdvolp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)

!     CALL MPI_GATHER(EOL0,1,MPI_DOUBLE_PRECISION,YOL,1,MPI_DOUBLE_PRECISION,
!     &                0,MPI_COMM_WORLD,IERR)
!     CALL MPI_GATHER(E_MIN, 1, MPI_DOUBLE_PRECISION, E_MINP, MPI_DOUBLE_PRECISION,
!     &                0,MPI_COMM_WORLD, IERR)                

!     Write trajectory
               write (18,*) '@@@',iold,inode(rep_id+1)
               call outvbs(0,18)
               write (18,*) '###'
!                call flush(18)
!     Write current configuration
               if ((mod(iold, nsave).eq.0)) then
                  filebase = "conf_0000.var"
                  call outvar(0, fileNameMP(filebase, 6, 9, rep_id+1))
               endif
            endif

            if(rep_id.eq.0.and.myrank.eq.0) then
               randomCount = 0
!  Update acceptance, temperature wise average of E and E^2 used to calculate
!  specific heat. 
               do i=1,num_rep
                  j=intem(i)
                  acy(i)=0.0
!  Above: contents of acy1 are added to acy(i) a few lines down. 
!  acy1(intem(i)) contains information received from the node at temperature
!  i, on how many updates have been accepted in node intem(i). Since acz
!  is not reset to 0 every cycle, acy(i) must be set to 0 here. Else, there 
!  will be serious double counting and the values of acceptance printed 
!  will be simply wrong.
               end do
               do i=1, num_rep
                  j=intem(i)
                  acy(i)=acy(i)+acy1(j)
                  eavm(i)= eavm(i)+yol(j)
                  sph(i) = sph(i)+yol(j)*yol(j)
               enddo


!     Write measurements to the time series file ts.d
               do i=1,num_rep
                  j=intem(i)
                     write(14,*) iold,i,j,pbe(i), dirp(j),
     &                 yol(j),eyslr(j), eyelp(j), eyvwp(j), eyhbp(j), 
     &                    eyvrp(j),eysmip(j), asa_p(j), vdvolp(j),
     &                    rgyrp(j),nhelp(j),nbetp(j),mhbp(j),
     &                 imhbp(j), nctotp(j),ncnatp(j), e_minp(j), 
     &                 eyabp(j),rmsdp(j)
!                      call flush(14)
               end do
!     Write the current parallel tempering information into par_R.in
!               timeLeft = llwrem(2) ! Time left till hard limit
!               if ((mod(iold, nsave).eq.0).or.(timeLeft.lt.minTimeLeft)
               if ((mod(iold, nsave).eq.0))
     &         then
                  open(13,file='par_R.in', status='unknown')
                  write(13,*) iold
                  do i=1,num_rep
                     write(13,*) i,inode(i),intem(i),yol(i),e_minp(i),
     &                    h_maxp(i)
                  end do
!     -------------------------- Various statistics of current run
!               swp=nswp-nequi
                  swp=nsw
                  write(13,*) 'Acceptance rate for change of chains:'
                  do k1=1,num_rep
                     temp=1.0d0/pbe(k1)/0.00198773
                     write(13,*) temp, acx1(k1)*2.0d0*nmes/swp 
!  Above: it's the acceptance rate of exchange of replicas. Since a 
!  replica exchange is attempted only once every nmes sweeps, the 
!  rate should be normalized with (nmes/swp).
                  end do 
                  write(13,*)
                  do k1=1,num_rep
                     k = intem(k1)   
                     temp=1.0d0/pbe(k1)/0.00198773
                     beta = pbe(k1)
                     geavm(k1) = nmes*eavm(k1)/swp
                     gsph(k1)  = (nmes*sph(k1)/swp-geavm(k1)**2)
     &                    *beta*beta/nresi
                     write(13,'(a,2f9.2,i4,f12.3)') 
     &                    'Temperature, Node,local acceptance rate:',
     &                    beta,temp,k,acy(k1)/dble(nsw*nvr)
!  Above: Changed (nswp-nequi) in the denominator of acceptance as 
!  acceptance values are initialized to 0 after equilibration cycles are
!  finished. Note also that since this is being written in the middle of
!  the simulation, it is normalized to nsw instead of nswp.
                     write(13,'(a,3f12.2)') 
     &                    'Last Energy, Average Energy, Spec. Heat:', 
     &                    yol(k),geavm(k1),gsph(k1) 
                     write(13,*) 
                  end do
                  close(13)
! Finally, flush the time series and trajectory files to ensure that we can do 
! a proper restart.
                  call flush(14)
                  call flush(18)
               end if

!--------------------Parallel Tempering  update
!     Swap with right neighbor (odd, even)            
               if(odd.eq.1) then
                  nu=1
                  no1 = num_rep-1
!     Swap with left neighbor (even, odd)
               else
                  nu = 2
                  no1 = num_rep
               end if
               do i=nu,no1,2
                  j=i+1
!     Periodic bc for swaps
                  if(i.eq.num_rep) j=1
                  in=intem(i)
                  jn=intem(j)
                  wij=exp(-pbe(i)*yol(jn)-pbe(j)*yol(in)
     &                 +pbe(i)*yol(in)+pbe(j)*yol(jn))
! The random number generator is getting out of sync here, because
!        the swap is only done on node 0!
! Keep track of number of random numbers used.
                  rd=grnd()
                  randomCount = randomCount + 1
!                  write (16,*) '>', iold, i,j
!     &            ,pbe(i),yol(in), pbe(j), yol(jn), wij, rd
                  if(wij.ge.rd) then
! Next line: Replica exchange only happens after equilibration, 
! which takes place outside this loop over nsw. So, I think nsw.gt.nequi
! is irrelevant for the calculation of acceptance of replica exchanges.
! /Sandipan
!                     if(nsw.gt.nequi) 
                     acx1(i) = acx1(i)+1
                     intem(i) = jn
                     intem(j) = in
                     inode(in)= j
                     inode(jn)= i
                  end if
               end do
!     ---------------- End Loop over nodes which creates a new temperature
!     map for all nodes, at the node with rank 0. 
!     
               odd = 1 - odd
            end if
!     End of "if (myrank.eq.0) ...". The block above includes PT update and 
!     writing of observables into the time series file etc. 
           
!     Below: Communicate new temperature-node map to all nodes
            CALL MPI_BCAST(INTEM,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,
     &           IERR)
            CALL MPI_BCAST(INODE,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,
     &           IERR)
! Synchronize random number generators for replica 0
            if (rep_id.eq.0) then
               CALL MPI_BCAST(randomCount,1,MPI_INTEGER,0,my_mpi_comm,
     &                IERR)
               if (myrank.ne.0) then
!                  write (*,*) '[', myrank,'] Missed', randomCount, 
!     &                            'random numbers.'
                  do i = 1, randomCount
                     rd = grnd()
!                     write (*,*) '[', myrank,'] rd=', rd
                  enddo
               endif
            endif

            BETA=PBE(INODE(rep_id+1))
            if (INODE(rep_id + 1).eq.1) dir = 1
            if (INODE(rep_id + 1).eq.num_rep) dir = -1

         endif
!        End of "if (mod(iold,nmes).eq.0) ..."
      end do
!-----------End Loop over sweeps
!     
!     OUTPUT:
!--------------------For Re-starts:
      nu = rep_id + 1
      filebase = "conf_0000.var"
      call outvar(0, fileNameMP(filebase, 6, 9, nu))
      e_final=energy()
      if (partem_comm.ne.MPI_COMM_NULL) then
         write (*,*) rep_id, ' E_final', e_final
      endif
      eol0 = eol
      acz0 = acz
      if (partem_comm.ne.MPI_COMM_NULL) then
         CALL MPI_GATHER(EOL0,1,MPI_DOUBLE_PRECISION,YOL,1,
     &        MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
         CALL MPI_GATHER(acz0,1,MPI_DOUBLE_PRECISION,acy1,1,
     &     MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
      endif
      
      if(rep_id.eq.0.and.myrank.eq.0) then
         close(14)
         open(13,file='par_R.in', status='unknown')
         write(13,*) iold
         do i=1,num_rep
            write(13,*) i,inode(i),intem(i),yol(i),e_minp(i),h_maxp(i)
         end do
!     -------------------------- Various statistics of current run
         swp=nswp
         write(13,*) 'Acceptance rate for change of chains:'
         do k1=1,num_rep
            temp=1.0d0/pbe(k1)/0.00198773
            write(13,*) temp, acx1(k1)*2.0d0*nmes/swp
         end do 
         write(13,*)
         do k1=1,num_rep
            k = intem(k1)   
            temp=1.0d0/pbe(k1)/0.00198773
            beta = pbe(k1)
            geavm(k1) = nmes*eavm(k1)/swp
            gsph(k1)  = (nmes*sph(k1)/swp-geavm(k1)**2)*beta*beta/nresi
            write(13,'(a,2f9.2,i4,f12.3)') 
     &           'Temperature, Node,local acceptance rate:',
     &           beta,temp,k,acy(k1)/dble((nswp)*nvr)
            write(13,'(a,3f12.2)') 
     &           'Last Energy, Average Energy, Spec. Heat:', 
     &           yol(k),geavm(k1),gsph(k1) 
            write(13,*) 
         end do
         close(13)
!         close(16)
      end if
      close(18)

!     =====================
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      return

      end
