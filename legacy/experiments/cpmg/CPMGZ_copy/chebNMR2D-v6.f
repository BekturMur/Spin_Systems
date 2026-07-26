        program chebspindynParallelWithDDandGates
        use mpi

        implicit none 
        include "chsdpar.h"

        integer narg,ioutstat, ioutGen
        character*32 infile

        integer iargc

        integer rootp, npMax
        parameter (rootp=0)
        parameter (npMax=64)

        integer Nit
        parameter (Nit=180) ! realizations of 2D spin arrangement
c        parameter (Nit=186*2) ! realizations of 2D spin arrangement
c        parameter (Nit=189) ! realizations of 2D spin arrangement
        double precision PIc, MDLT, MDLT2
        parameter (PIc=3.14159265358979323846D0)
        parameter (MDLT=1.0D-15, MDLT2=1.0D-7)

        integer it, itSent,itRecvd,itC,iToGo
        integer ioutArr(npMax-1)
        integer iittag,iouttag,itogotag,vecresouttag
        integer istatus(MPI_STATUS_SIZE), reqiout(npMax-1)
        integer istarray(MPI_STATUS_SIZE,npMax-1)
        integer ierr,lr,np0,np, ip, i
        double precision vecResOut(maxLenVecResOut)
        integer mrkResOut(maxLenVecResOut), mrkresouttag
        integer iptVecResOut, iptvecresouttag
        double precision vecResAv(maxLenVecResOut)
        double precision vecResDev(maxLenVecResOut)

        integer nfiles
        double precision JxT(maxNconn,maxLen), JyT(maxNconn,maxLen)
        double precision JzT(maxNconn,maxLen)
        double precision HxT(maxLtot,maxLen), HyT(maxLtot,maxLen)
        double precision HzT(maxLtot,maxLen), tauT(maxLen)
        integer LsysT(maxLen),LbathT(maxLen),LtotT(maxLen)
        integer nconnT(maxLen)
        integer aspT(maxNconn,maxLen),bspT(maxNconn,maxLen)
        integer nstepsT(maxLen),iRandT(maxLen)
        integer init_stateT(maxLen),autonormT(maxLen),imflagT(maxLen)
        double precision EpsFT(maxLen)
        integer ispecflagT(maxLen)
        integer ifile, nels

        integer imarkT(2,maxLen), nooutputT(maxLen)
        integer iNtypeT(maxLtot,maxLen),intypetag
        logical flagtst

        double precision xnorm,timeC, MaxErr/0.0D0/, xx,yy
        integer irnd,j,k,is
        double precision xrnd,xJz(maxNconn),JzScale/1.D0/,HzScale/1.D2/
        double precision dHz(maxLtot), HzT1(maxLtot,maxLen)
        integer locLtot,locnc,locasp(maxNconn),locbsp(maxNconn),lociRand

c        integer lu
c        character*36 lfname,lfname0,lfnameN,lfname1
        integer ltu
        character*36 ltfname,ltfname0,ltfnameN,ltfname1

c executable statements

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,lr,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np0,ierr)

c       if (lr.ne.rootp) then
c          lu = lr+10
c          lfname0 = 'ptst'
c          lfname1 = '.out'
c          write(lfnameN,'(i0.3)') lr
c          lfname = trim(lfname0)//trim(lfnameN)//trim(lfname1)
c          open(lu,file=trim(lfname))
c          close(lu)
c        endif

        iittag = 1
        iouttag = 2
        itogotag = 4
        vecresouttag = 5
        iptvecresouttag = 6
        mrkresouttag = 7
        intypetag = 8

        np = min(npMax,np0,Nit+1)
        ioutstat=0
        ioutGen = 0
        iToGo = 1
        do i=1,npMax-1
          ioutArr(i)=0
        enddo
        do i=1,maxLenVecResOut
          vecResAv(i) = 0.0D0
          vecResDev(i) = 0.0D0
        enddo

        if (lr.eq.rootp) then 
          open(chlog,file=logfnam)
          write(chlog,'(a)') 
     +      '------------------------------------------------'
          write(chlog,'(a)') 
     +      'START OF Chebyshev-based quantum spin dynamics'
          write(chlog,'(a)') 
     +      '------------------------------------------------'
          close(chlog)
          narg=iargc()
          if (narg.lt.1) then
            open(chlog,file=logfnam,access='append')
            write(chlog,*) 'No name of the input file'
            close(chlog)
            ioutGen = -567
            iToGo=0
            goto 2000 
          endif
          call GetArg(1,infile)
          if (narg.gt.1) then
            open(chlog,file=logfnam,access='append')
            write(chlog,*) 'Several input files, only the 
     +        first one used: ',infile(:len_trim(infile))
            close(chlog)
          endif

          call chparsgenPDDG(ioutstat,infile,nfiles,
     +      aspT,bspT,JxT,JyT,JzT,nconnT,HxT,HyT,HzT,
     +      LsysT,LbathT,LtotT,iRandT,nstepsT,tauT,init_stateT,
     +      autonormT,EpsFT,imflagT,ispecflagT, iNtypeT,
     +      nooutputT,imarkT)

          if (ioutstat.ne.0) then
            open(chlog,file=logfnam,access='append')
            write(chlog,*) 'Parsing error, ioutstat=',ioutstat
            close(chlog)
            ioutGen = ioutstat
            iToGo=0
            goto 2000
          endif

c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
          if (nconnT(1).gt.0.and.dabs(JzT(1,1)).gt.MDLT2) then
            JzScale=JzScale*JzT(1,1)
          endif
          do ifile=2,nfiles
            if (LtotT(ifile).ne.LtotT(1)) then
              ioutstat=-1893
              open(chlog,file=logfnam,access='append')
              write(chlog,*)'LtotT Error in file #',ifile
              write(chlog,*)'LtotT=',LtotT(ifile),' LtotT(1)=',LtotT(1)
              write(chlog,*)'Should be the same, ioutstat=',ioutstat
              close(chlog)
              ioutGen = ioutstat
              iToGo=0
              goto 2000
            endif
          enddo
c END
        endif ! if root

c end parsing and checking input files

2000    continue
        call MPI_BCAST(iToGo,1,MPI_INTEGER,rootp,MPI_COMM_WORLD,ierr)
        if (iToGo.eq.0) goto 1000

        call MPI_BCAST(nfiles,1,MPI_INTEGER,rootp,MPI_COMM_WORLD,ierr)

        nels = maxLen
        call MPI_BCAST(LsysT,nels,MPI_INTEGER,rootp,MPI_COMM_WORLD,
     +    ierr)
        call MPI_BCAST(LbathT,nels,MPI_INTEGER,rootp,MPI_COMM_WORLD,
     +    ierr)
        call MPI_BCAST(LtotT,nels,MPI_INTEGER,rootp,MPI_COMM_WORLD,
     +    ierr)
c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
c        call MPI_BCAST(nconnT,nels,MPI_INTEGER,rootp,MPI_COMM_WORLD,
c     +    ierr)
c END
        call MPI_BCAST(iRandT,nels,MPI_INTEGER,rootp,MPI_COMM_WORLD,
     +    ierr)
        call MPI_BCAST(nstepsT,nels,MPI_INTEGER,rootp,MPI_COMM_WORLD,
     +    ierr)
        call MPI_BCAST(init_stateT,nels,MPI_INTEGER,rootp,
     +    MPI_COMM_WORLD,ierr)
        call MPI_BCAST(autonormT,nels,MPI_INTEGER,rootp,
     +    MPI_COMM_WORLD,ierr)
        call MPI_BCAST(imflagT,nels,MPI_INTEGER,rootp,
     +    MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ispecflagT,nels,MPI_INTEGER,rootp,
     +    MPI_COMM_WORLD,ierr)

        call MPI_BCAST(tauT,nels,MPI_DOUBLE_PRECISION,rootp,
     +    MPI_COMM_WORLD,ierr)
        call MPI_BCAST(EpsFT,nels,MPI_DOUBLE_PRECISION,rootp,
     +    MPI_COMM_WORLD,ierr)

        nels = maxLtot*maxLen
        call MPI_BCAST(HxT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
     +    MPI_COMM_WORLD,ierr)
        call MPI_BCAST(HyT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
     +    MPI_COMM_WORLD,ierr)
        call MPI_BCAST(HzT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
     +    MPI_COMM_WORLD,ierr)
c
c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
c        nels = maxNconn*maxLen
c        call MPI_BCAST(aspT(1,1),nels,MPI_INTEGER,rootp,
c     +    MPI_COMM_WORLD,ierr)
c        call MPI_BCAST(bspT(1,1),nels,MPI_INTEGER,rootp,
c     +    MPI_COMM_WORLD,ierr)
c        call MPI_BCAST(JxT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
c     +    MPI_COMM_WORLD,ierr)
c        call MPI_BCAST(JyT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
c     +    MPI_COMM_WORLD,ierr)
c        call MPI_BCAST(JzT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
c     +    MPI_COMM_WORLD,ierr)
c
c END
        nels=2*maxLen
        call MPI_BCAST(imarkT(1,1),nels,MPI_INTEGER,rootp,  ! imark - number of the file to be executed next
     +    MPI_COMM_WORLD,ierr)                              ! normally: the next file, when ENDCYCLE - file at the start of the cycle

        nels=maxLtot*maxLen
        call MPI_BCAST(iNtypeT(1,1),nels,MPI_INTEGER,rootp,
     +    MPI_COMM_WORLD,ierr)
        nels=maxLen
        call MPI_BCAST(nooutputT,nels,MPI_INTEGER,rootp,MPI_COMM_WORLD,
     +    ierr)

c DATA SENT TO ALL PROCS, START MAIN BODY

        if (lr.eq.rootp) then ! IF ROOT - MAIN BODY

          iToGo=0 ! extra procs - not needed
          do ip=np,np0-1 
          call MPI_SEND(iToGo,1,MPI_INTEGER,ip,itogotag,MPI_COMM_WORLD,
     +      ierr)
          enddo

          iToGo=1 ! procs to work with
          do ip=1,np-1 
            call MPI_SEND(iToGo,1,MPI_INTEGER,ip,itogotag,
     +        MPI_COMM_WORLD,ierr)
            call MPI_IRECV(ioutArr(ip),1,MPI_INTEGER,ip,iouttag,
     +        MPI_COMM_WORLD,reqiout(ip),ierr)
          enddo

          itRecvd=0
          itSent=0
          do while (itRecvd.lt.Nit)
            call MPI_WAITANY(np-1,reqiout,ip,istatus,ierr)
c            if (ip.eq.MPI_UNDEFINED) goto 1000 ! if all procs completed
            if (ip.eq.MPI_UNDEFINED) then
              call MPI_TESTALL(np-1,reqiout,flagtst,istarray,ierr)
              if (flagtst) then
                goto 1000 ! all procs completed
              else
                open(chlog,file=logfnam,access='append')
                write(chlog,*) 'WAITANY undefined'
                close(chlog)
                iToGo=0
              endif
            endif
            call MPI_RECV(itC,1,MPI_INTEGER,ip,iittag,MPI_COMM_WORLD,
     +        istatus,ierr)

            if (ioutArr(ip).ne.0) then
              ioutGen = ioutArr(ip)
              iToGo=0
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Error iout=',ioutArr(ip),' ip=',ip,
     +          ' itC=',itC
              close(chlog)
            endif

            if (itC.ne.-1) then ! if there is output, i.e. not first time
              itRecvd=itRecvd+1 
              call MPI_RECV(iptVecResOut,1,MPI_INTEGER,ip,
     +          iptvecresouttag,MPI_COMM_WORLD,istatus,ierr)
              call MPI_RECV(vecResOut,maxLenVecResOut,
     +          MPI_DOUBLE_PRECISION,ip,vecresouttag,
     +          MPI_COMM_WORLD,istatus,ierr)
              call MPI_RECV(mrkResOut,maxLenVecResOut,
     +          MPI_INTEGER,ip,mrkresouttag,
     +          MPI_COMM_WORLD,istatus,ierr)

              if (ioutArr(ip).eq.0) then ! if no error: process output and put it in proper place
                i=1
                do while (i.le.iptVecResOut)
                  if (dabs(vecResOut(i)).gt.1.0D0/MDLT) then
                    ioutstat=-777
                    iToGo=0
                    open(chlog,file=logfnam,access='append')
                    write(chlog,*)'Error, ip=',ip
                    write(chlog,*)'vecResOut too large, i=',i
                    close(chlog)
                  endif
                  timeC = vecResOut(i)
                  do j=i+1,i+mrkResOut(i)-1
                    if (dabs(vecResOut(j)).gt.1.0D0/MDLT) then
                      ioutstat=-777
                      iToGo=0
                      open(chlog,file=logfnam,access='append')
                      write(chlog,*)'Error, ip=',ip
                      write(chlog,*)'vecResOut too large, i=',i,' j=',j
                      close(chlog)
                    endif
c running average and deviation
                  xx = vecResOut(j) - vecResAv(j)
                  vecResAv(j) = vecResAv(j) + xx/itRecvd ! itRecvd starts from 1
                  vecResDev(j)=vecResDev(j)+(itRecvd-1.D0)*xx*xx/itRecvd
                  enddo
                  i = i + mrkResOut(i)
cc now i marks beginning of the next record, i.e. i=iptVecResOut+1 for the last record
cc here we do not reset iptVecResOut, and do not touch VecResOut or mrkVecResOut
cc all this is done in chebsd
                enddo

c output running averages and deviations every 8th time
                if (mod(itRecvd,8).eq.0) then
                open(71,file='echo2D.cur')
                i=1
                do while (i.le.iptVecResOut)
                  timeC = vecResOut(i)
                  write(71,'(f12.4$)') timeC
                  do j=i+1,i+mrkResOut(i)-1
                   write(71,'(f12.4$)') vecResAv(j)
                   write(71,'(f12.4$)') dsqrt(vecResDev(j)/itRecvd) ! important: we output only every 8th time, so itRecvd-1 is never zero
                  enddo
                  i = i + mrkResOut(i)
                  write(71,'(a)') ' '
                enddo
                close(71)
                endif

c example: test individual vecResOut
c                open(75,file='veclog.out',access='append')
c                write(75,'(a)') '#' ! separator from other instances of it
c                write(75,'(i10)') itC
c                i=1
c                do while (i.le.iptVecResOut)
c                  write(75,'(f10.5$)') vecResOut(i) ! timeC
c                  write(75,'(i5$)') mrkResOut(i) ! timeC
c                  do j=i+1,i+mrkResOut(i)-1
c                    write(75,'(f10.5$)') vecResOut(j) 
c                  enddo
c                  i = i + mrkResOut(i)
c                  write(75,'(a)') ' ' 
c                  write(75,'(a)') ' ' 
c                enddo
c                close(75)
c end of example
c                do i=1,iptVecResOut,4
c                  write(71,'(e18.8$)') vecResOut(i)   ! timeGl
c                  write(71,'(e18.8$)') vecResOut(i+1) ! sx(1)
c                  write(71,'(e18.8$)') vecResOut(i+2) ! sy(1)
c                  write(71,'(e18.8)') vecResOut(i+3)  ! sz(1)
c                enddo

                open(chlog,file=logfnam,access='append')
                write(chlog,*)'Recvd ',itRecvd,' itC=',itC,' ip=',ip
                close(chlog)
              else
c
c               if error - vecResOut contains information about the error
c               analyze this information as needed
c
                open(chlog,file=logfnam,access='append')
                write(chlog,*)'EpsF=',vecResOut(1)
                write(chlog,*)'alp=',vecResOut(2)
                write(chlog,*)'emax=',vecResOut(3)
                write(chlog,*)'nord=',vecResOut(4)
                write(chlog,*)'tau=',vecResOut(5)
                close(chlog)
                iToGo = 0
              endif
            endif

            if (itSent.lt.Nit.and.iToGo.ne.0) then ! if new input to be sent

              call MPI_SEND(iToGo,1,MPI_INTEGER,ip,itogotag,
     +          MPI_COMM_WORLD,ierr)
              call MPI_SEND(itSent,1,MPI_INTEGER,ip,iittag,
     +          MPI_COMM_WORLD,ierr)
c
c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
c The input will be generated locally, see below
c END
              open(chlog,file=logfnam,access='append')
              write(chlog,*)'Sent itSent=',itSent,' ip=',ip
              close(chlog)
              itSent = itSent+1
c
              call MPI_IRECV(ioutArr(ip),1,MPI_INTEGER,ip,iouttag,
     +          MPI_COMM_WORLD,reqiout(ip),ierr)
            else
              iToGo=0
              call MPI_SEND(iToGo,1,MPI_INTEGER,ip,itogotag,
     +          MPI_COMM_WORLD,ierr)
            endif

          enddo ! while itRecvd

          if (ioutGen.ne.0) then
            open(chlog,file=logfnam,access='append')
            write(chlog,*) '!!!=========>'
            write(chlog,*) 'UNSUCCESSFUL, IOUTGEN=',ioutGen
            write(chlog,*) '!!!=========>'
            close(chlog)
          else
            open(71,file='echo2D.out')
            i=1
            do while (i.le.iptVecResOut)
              timeC = vecResOut(i)
              write(71,'(f18.10$)') timeC
c              write(71,'(i15$)') mrkResOut(i) ! used for tests, to check mrkResOut
              do j=i+1,i+mrkResOut(i)-1
                 write(71,'(f18.10$)') vecResAv(j)
                 write(71,'(f18.10$)') dsqrt(vecResDev(j)/Nit)
c                write(71,'(f18.10$)') vecResAv(j)
c                write(71,'(f25.15$)') vecResOut(j) ! used for tests, to check the last vecResOut
              enddo
              i = i + mrkResOut(i)
cc now i marks beginning of the next record, or i=iptVecResOut+1 for the last record
              write(71,'(a)') ' '
            enddo
            close(71)

            open(chlog,file=logfnam,access='append')
            write(chlog,*) 'SUCCESSFUL END'
            close(chlog)
          endif
c          goto 1000

        else ! IF NOT ROOT - MAIN BODY

          call MPI_RECV(iToGo,1,MPI_INTEGER,rootp,itogotag,
     +      MPI_COMM_WORLD,istatus,ierr)
          if (iToGo.eq.0) goto 1000 ! this is an extra proc, not needed

          it=-1 ! first time - no output, just reporting ready
          ioutstat=0
          call MPI_SEND(ioutstat,1,MPI_INTEGER,rootp,iouttag,
     +      MPI_COMM_WORLD,ierr)
          call MPI_SEND(it,1,MPI_INTEGER,rootp,iittag,
     +      MPI_COMM_WORLD,ierr)

          do while (.true.) ! main cycle
            call MPI_RECV(iToGo,1,MPI_INTEGER,rootp,itogotag,
     +        MPI_COMM_WORLD,istatus,ierr)
            if (iToGo.eq.0) goto 1000

            call MPI_RECV(it,1,MPI_INTEGER,rootp,iittag,MPI_COMM_WORLD,
     +        istatus,ierr)

c
c              prepare new input for this process
c
c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
            do ifile=1,nfiles
              iRandT(ifile)=iRandT(ifile)+it*11489 ! 11489 is prime
            enddo
            lociRand=iRandT(1)
            locLtot=LtotT(1)
            call gendipnmr2D(lociRand,locLtot,locnc,xJz,locasp,locbsp)
            do i=1,locLtot
              xx=MDLT+(1.0D0-2*MDLT)*rand()
              yy=rand()
              xrnd=dsqrt(-2.0D0*dlog(xx))*dcos(2*PIc*yy)
              dHz(i) = HzScale*xrnd
            enddo
            do ifile=1,maxLen
              nconnT(ifile)=locnc
              do i=1,locnc
                JzT(i,ifile)=JzScale*xJz(i)
                JxT(i,ifile)=-0.5D0*JzScale*xJz(i)
                JyT(i,ifile)=-0.5D0*JzScale*xJz(i)
                aspT(i,ifile)=locasp(i)
                bspT(i,ifile)=locbsp(i)
              enddo
              do i=locnc+1,maxNconn
                JxT(i,ifile)=0.0D0
                JyT(i,ifile)=0.0D0
                JzT(i,ifile)=0.0D0
                aspT(i,ifile)=0
                bspT(i,ifile)=0
              enddo
              do i=1,locLtot
c                HxT(i,ifile)=0.0D0
c                HyT(i,ifile)=0.0D0
                HzT1(i,ifile)=HzT(i,ifile) + dHz(i)
              enddo
              do i=locLtot+1,maxLtot
                HxT(i,ifile)=0.0D0
                HyT(i,ifile)=0.0D0
                HzT1(i,ifile)=0.0D0
              enddo
            enddo
c END
c
c example of logging the local parameters
c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
c logging the local parameters
            if (mod(4*it,Nit).eq.0) then ! it starts from 0
              ltu = it+1
              ltfname0 = 'plog'
              ltfname1 = '.out'
              write(ltfnameN,'(i0.3)') ltu
              ltfname = trim(ltfname0)//trim(ltfnameN)//trim(ltfname1)
              open(ltu,file=trim(ltfname))
              write(ltu,'(a,i6,a,i6)') 'start steps at lr=',lr,'  it=',it
              write(ltu,'(a,i6)') 'nfiles=',nfiles
              do ifile=1,nfiles
                write(ltu,'(a,i5)') 'ifile=',ifile
                do i=1,LtotT(ifile)
                  write(ltu,'(i5,3f20.5)') 
     +              i,HxT(i,ifile),HyT(i,ifile),HzT1(i,ifile)
                enddo
                write(ltu,'(a)') ' '
c          do i=1,nconnT(1)
                do i=1,locnc
                  write(ltu,'(2i5$)') aspT(i,ifile),bspT(i,ifile)
                  write(ltu,'(3f20.5)') 
     +              JxT(i,ifile),JyT(i,ifile),JzT(i,ifile)
                enddo
                write(ltu,'(a)') ' '
              enddo
              close(ltu)
            endif
c END
c end of example

            call chebsdPDDG(ioutstat,nfiles,
     +        aspT,bspT,JxT,JyT,JzT,nconnT,HxT,HyT,HzT1,
     +        LsysT,LbathT,LtotT,iRandT,nstepsT,tauT,init_stateT,
     +        autonormT,EpsFT,imflagT,ispecflagT, 
     +        iNtypeT,nooutputT,imarkT,
     +        vecResOut,iptVecResOut,mrkResOut)
c
c example of logging the local results
c          open(lu,file=trim(lfname),access='append')
c          write(lu,'(a,i6)') 'local results from lr=',lr
c          i=1
c          do while (i.le.iptVecResOut)
c            timeC = vecResOut(i)
c            write(lu,'(f15.10$)') vecResOut(j)
c            xnorm = vecResOut(i+1)
c            do j=i+2,i+mrkResOut(i)-1
c              write(lu,'(f15.10$)') vecResOut(j)
c            enddo
c            i = i + mrkResOut(i)
c            write(lu,'(a)') ' '
c          enddo
c          write(lu,'(a)') ' '
c          close(lu)
c end of example

            call MPI_SEND(ioutstat,1,MPI_INTEGER,rootp,iouttag,
     +        MPI_COMM_WORLD,ierr)
            call MPI_SEND(it,1,MPI_INTEGER,rootp,iittag,
     +        MPI_COMM_WORLD,ierr)

c            if (ioutstat.eq.0) then
              call MPI_SEND(iptVecResOut,1,MPI_INTEGER,rootp,
     +          iptvecresouttag,MPI_COMM_WORLD,ierr)
              call MPI_SEND(vecResOut,maxLenVecResOut,
     +          MPI_DOUBLE_PRECISION,rootp,vecresouttag,
     +          MPI_COMM_WORLD,ierr)
              call MPI_SEND(mrkResOut,maxLenVecResOut,
     +          MPI_INTEGER,rootp,mrkresouttag,
     +          MPI_COMM_WORLD,ierr)
c            endif
          enddo ! while(.true.) - main cycle

        endif ! IF ROOT OR NOT

1000    continue
        call MPI_FINALIZE(ierr)
        end

