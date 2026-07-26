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
c        parameter (Nit=200)
        integer it, itSent,itRecvd,itC,iToGo
        double precision Dt 
        integer ioutArr(npMax-1)
        integer iittag,iouttag,itogotag,vecresouttag, iDttag
        integer istatus(MPI_STATUS_SIZE), reqiout(npMax-1)
        integer ierr,lr,np0,np, ip, i
        double precision vecResOut(maxLenVecResOut)
        integer mrkResOut(maxLenVecResOut), mrkresouttag
        integer iptVecResOut, iptvecresouttag
        double precision vecResAv(maxLenVecResOut)

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
        integer iNtypeT(maxLtot,maxLen)
        double precision HzT1(maxLtot,maxLen)
        double precision JxT1(maxNconn,maxLen),JyT1(maxNconn,maxLen) 
        integer ihztag,intypetag, ijxtag,ijytag, is1,is2
        logical lt1,lt2,lt3
        double precision PIc
        parameter (PIc=3.14159265358979323846D0)

        double precision xnorm,timeC, MaxErr/0.0D0/
        integer irnd,j,k
        double precision dHres(6)


        integer lu
        character*36 lfname,lfname0,lfnameN,lfname1




c executable statements

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,lr,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np0,ierr)


c        if (lr.ne.rootp) then
c          lu = lr+10
c          lfname0 = 'ptst'
c          lfname1 = '.out'
c          write(lfnameN,'(i0.3)') lr
c          lfname = trim(lfname0)//trim(lfnameN)//trim(lfname1)
c          open(lu,file=trim(lfname))
c          close(lu)
c        endif


        dHres(1) = -2*PIc*0.0D0
        dHres(2) = -2*PIc*0.0D0
        dHres(3) = -2*PIc*0.0D0
        dHres(4) = -2*PIc*6.5D0
        dHres(5) = -2*PIc*0.0D0
        dHres(6) = -2*PIc*0.0D0
c        dHres(1) = -2*PIc*268.6D0
c        dHres(2) = -2*PIc*300.6D0
c        dHres(3) = -2*PIc*390.1D0
c        dHres(4) = -2*PIc*396.6D0
c        dHres(5) = -2*PIc*469.4D0
c        dHres(6) = -2*PIc*493.8D0


#tags
        iittag = 1
        iouttag = 2
        ihztag = 3
        itogotag = 4
        vecresouttag = 5
        iptvecresouttag = 6
        mrkresouttag = 7
        intypetag = 8
        ijxtag = 9
        ijytag = 10


        Nit = 3*(np0-1) #number of iterations
        #np- number of processes


        np = min(npMax,np0,Nit+1)
        ioutstat=0 #it should be zero, status of out
        ioutGen = 0 #it vozvrashat znacheine
        iToGo = 1 # idem , 1 - idem , 0 -ne idem
        do i=1,npMax-1
          ioutArr(i)=0 #current sutuation of process
        enddo
        do i=1,maxLenVecResOut #максимальная длина вектора
          vecResAv(i) = 0.0D0
        enddo
#rank 0 = root,manager, rank neq = 0 working processor
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
            ioutGen = -567 #error, it should be 0
            iToGo=0
            goto 2000 # tushi svet, my zakanchivaem
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
     +      nooutputT,imarkT) # this program should format vector of parameteres of input processes



#тут комменты для того, потому что при компиляции могут быть проблемы -  этот кусок для тестов
c          open(chlog,file=logfnam,access='append')
c          write(chlog,'(a,i6)') 'chk0 nfiles=',nfiles
c          do i=1,LtotT(1)
c            write(chlog,'(i4$)') iNtypeT(i,1)
c          enddo
c          write(chlog,'(a)') ' '
c          do i=1,LtotT(2)
c            write(chlog,'(i4$)') iNtypeT(i,2)
c          enddo
c          write(chlog,'(a)') ' '
c          close(chlog)


#iNtype - index (i - integer)
#input of two types - input-file and data-file
c checking iNtype, modifying iNtypeT array for this specific calculation
          do i=1,maxLtot #size of system
            k=0 ! k=1 if there are files where i <= Ltot
            j=-1 ! j neq -1 if iNtype was set for the spin i
            do ifile=1,nfiles
              if (i.le.LtotT(ifile)) then
                k=1
                #i - number of spin, ifile - number of file
                if (iNtypeT(i,ifile).ne.-1) then
                  j=iNtypeT(i,ifile) ! iNtype is determined by its first value, future changes ignored
                  exit
                endif
              endif
            enddo
            if (k.eq.1) then
              if (j.eq.0) then ! NV
                do ifile=1,nfiles
                iNtypeT(i,ifile)=0 ! line 0 - NV itself
                enddo
              else if (j.eq.1) then ! P1 [111]
                do ifile=1,nfiles
                iNtypeT(i,ifile) = 1 ! line 1 (can also be 3 or 6)
                enddo
              else if (j.ge.2.and.j.le.4) then ! P1 [-111] or other
                do ifile=1,nfiles
                iNtypeT(i,ifile) = 2 ! line 2 (can also be 4 or 5)
                enddo
              else ! error - wrong value of iNtype
                open(chlog,file=logfnam,access='append')
                write(chlog,*)'Wrong value of iNtype !!!'
                write(chlog,*)'i=',i,'  ifile=',ifile,'  iNtype=',j
                write(chlog,*)'Execution interrupted'
                close(chlog)
                ioutstat = -4567
                exit
              endif
            endif
          enddo
                    
          if (ioutstat.ne.0) then
            open(chlog,file=logfnam,access='append')
            write(chlog,*) 'Parsing error, ioutstat=',ioutstat
            close(chlog)
            ioutGen = ioutstat
            iToGo=0
c            goto 2000
          endif
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
        call MPI_BCAST(nconnT,nels,MPI_INTEGER,rootp,MPI_COMM_WORLD,
     +    ierr)
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
c        call MPI_BCAST(HzT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
c     +    MPI_COMM_WORLD,ierr)

        nels = maxNconn*maxLen
        call MPI_BCAST(aspT(1,1),nels,MPI_INTEGER,rootp,
     +    MPI_COMM_WORLD,ierr)
        call MPI_BCAST(bspT(1,1),nels,MPI_INTEGER,rootp,
     +    MPI_COMM_WORLD,ierr)
c        call MPI_BCAST(JxT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
c     +    MPI_COMM_WORLD,ierr)
c        call MPI_BCAST(JyT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
c     +    MPI_COMM_WORLD,ierr)
        call MPI_BCAST(JzT(1,1),nels,MPI_DOUBLE_PRECISION,rootp,
     +    MPI_COMM_WORLD,ierr)

        nels=2*maxLen
        call MPI_BCAST(imarkT(1,1),nels,MPI_INTEGER,rootp,
     +    MPI_COMM_WORLD,ierr)
        nels=maxLtot*maxLen
c        call MPI_BCAST(iNtypeT(1,1),nels,MPI_INTEGER,rootp,
c     +    MPI_COMM_WORLD,ierr)
        nels=maxLen
        call MPI_BCAST(nooutputT,nels,MPI_INTEGER,rootp,MPI_COMM_WORLD,
     +    ierr)

c DATA SENT TO ALL PROCS, START MAIN BODY
#теперь есть два варианта
        if (lr.eq.rootp) then ! IF ROOT - MAIN BODY

          open(71,file='rdm1PDDG.out')
          close(71,status='delete')
          open(71,file='rdm1PDDG.out')
          close(71)
          open(72,file='nrmPDDG.out')
          close(72,status='delete')
          open(72,file='nrmPDDG.out')
          close(72)

          iToGo=0 ! extra procs - not needed
          do ip=np,np0-1 
          call MPI_SEND(iToGo,1,MPI_INTEGER,ip,itogotag,MPI_COMM_WORLD,
     +      ierr)
          enddo
#аналогия  с училкой и школьниками - вызов мастера от рабочего
          iToGo=1 ! procs to work with
          do ip=1,np-1 
            call MPI_SEND(iToGo,1,MPI_INTEGER,ip,itogotag,
     +        MPI_COMM_WORLD,ierr)
            call MPI_IRECV(ioutArr(ip),1,MPI_INTEGER,ip,iouttag,
     +        MPI_COMM_WORLD,reqiout(ip),ierr)
          enddo

          itRecvd=0
          itSent=0
          #ip - number of process which says i'm ready
          #itc - iteration card - текущий номер ситуации
          do while (itRecvd.lt.Nit)
            call MPI_WAITANY(np-1,reqiout,ip,istatus,ierr)
            if (ip.eq.MPI_UNDEFINED) goto 1000 ! if all procs completed

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

              if (ioutArr(ip).eq.0) then ! if no error
c
c               if no error - put the output in proper place and process it
c
                i=1
                do while (i.le.iptVecResOut)
                  timeC = vecResOut(i) #время
                  xnorm = vecResOut(i+1) #нормировака волновой функции, чтобы смотреть насколько точни вычисления
                  MaxErr = max(MaxErr,dabs(xnorm-1.0D0)) #проверка
                  do j=i+2,i+mrkResOut(i)-1 #mark - какой длины запись
                    vecResAv(j) = vecResAv(j) + vecResOut(j)/Nit #усреднение результатов
                  enddo
                  i = i + mrkResOut(i) #увеличили i на длину записи
c now i marks beginning of the next record, or i=iptVecResOut+1 for the last record
                enddo
c
c                write(71,'(a)') '#' ! 
c                write(71,'(a)') ' ' ! separator from other instances of it
c                do i=1,iptVecResOut,4
c                  write(71,'(e18.8$)') vecResOut(i)   ! timeGl
c                  write(71,'(e18.8$)') vecResOut(i+1) ! sx(1)
c                  write(71,'(e18.8$)') vecResOut(i+2) ! sy(1)
c                  write(71,'(e18.8)') vecResOut(i+3)  ! sz(1)
c                enddo
                close(71)
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
                iToGo = 0 #продолжаешь работать
              endif
            endif
#новая порция данных
            if (itSent.lt.Nit.and.iToGo.ne.0) then ! if new input to be sent
              itSent = itSent+1
              open(chlog,file=logfnam,access='append')
              write(chlog,*)'Sent ',itSent,' ip=',ip
              close(chlog)
c
c              prepare new input and send it
c

              do i=1,maxLtot
                j = iNtypeT(i,1)
                if (j.eq.1.or.j.eq.3.or.j.eq.6) then
                  irnd = int(3.0D0*rand())
                  if (irnd.eq.0) j=1
                  if (irnd.eq.1) j=3
                  if (irnd.eq.2) j=6
                else if (j.eq.2.or.j.eq.4.or.j.eq.5) then
                  irnd = int(3.0D0*rand())
                  if (irnd.eq.0) j=2
                  if (irnd.eq.1) j=4
                  if (irnd.eq.2) j=5
                endif
                do ifile=1,nfiles
                  iNtypeT(i,ifile) = j
                  HzT1(i,ifile) = HzT(i,ifile) + dHres(j)
                enddo
              enddo
              do ifile=1,nfiles
                do i=1,nconnT(ifile)
                  is1 = aspT(i,ifile)
                  is2 = bspT(i,ifile)
                  lt1=iNtypeT(is1,ifile).eq.3.or.iNtypeT(is1,ifile).eq.4
                  lt2=iNtypeT(is2,ifile).eq.3.or.iNtypeT(is2,ifile).eq.4
                  lt3=iNtypeT(is1,ifile).eq.iNtypeT(is2,ifile)
                  if ((lt1.and.lt2).or.lt3) then
                    JxT1(i,ifile) = JxT(i,ifile)
                    JyT1(i,ifile) = JyT(i,ifile)
                  else
                    JxT1(i,ifile) = 0.0D0
                    JyT1(i,ifile) = 0.0D0
                  endif
                enddo
              enddo


              call MPI_SEND(iToGo,1,MPI_INTEGER,ip,itogotag,
     +          MPI_COMM_WORLD,ierr)
              call MPI_SEND(itSent,1,MPI_INTEGER,ip,iittag,
     +          MPI_COMM_WORLD,ierr)

              nels = maxLtot*maxLen
              call MPI_SEND(HzT1,nels,MPI_DOUBLE_PRECISION,ip,ihztag,
     +          MPI_COMM_WORLD,ierr)
              call MPI_SEND(iNtypeT,nels,MPI_INTEGER,ip,intypetag,
     +          MPI_COMM_WORLD,ierr)
              nels = maxNconn*maxLen
              call MPI_SEND(JxT1,nels,MPI_DOUBLE_PRECISION,ip,ijxtag,
     +          MPI_COMM_WORLD,ierr)
              call MPI_SEND(JyT1,nels,MPI_DOUBLE_PRECISION,ip,ijytag,
     +          MPI_COMM_WORLD,ierr)

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
            open(72,file='nrmPDDG.out',access='append')
            write(72,'(f25.15)') MaxErr
            close(72)
            open(71,file='rdm1PDDG.out',access='append')
            i=1
            do while (i.le.iptVecResOut)
              timeC = vecResOut(i)
              write(71,'(f15.10$)') timeC
              do j=i+2,i+mrkResOut(i)-1
                write(71,'(f25.15$)') vecResAv(j)
              enddo
              i = i + mrkResOut(i)
c now i marks beginning of the next record, or i=iptVecResOut+1 for the last record
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
#отправили училке задание - я готов
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

            nels = maxLtot*maxLen
            call MPI_RECV(HzT1,nels,MPI_DOUBLE_PRECISION,rootp,ihztag,
     +        MPI_COMM_WORLD,istatus,ierr)
            call MPI_RECV(iNtypeT,nels,MPI_INTEGER,rootp,intypetag,
     +        MPI_COMM_WORLD,istatus,ierr)
            nels = maxNconn*maxLen
            call MPI_RECV(JxT1,nels,MPI_DOUBLE_PRECISION,rootp,ijxtag,
     +        MPI_COMM_WORLD,istatus,ierr)
            call MPI_RECV(JyT1,nels,MPI_DOUBLE_PRECISION,rootp,ijytag,
     +        MPI_COMM_WORLD,istatus,ierr)

c
c prepare input parameters for this specific proc
c here is an example
c



c          open(lu,file=trim(lfname),access='append')
c          write(lu,'(a,i6)') 'start steps at lr=',lr
c          write(lu,'(a,i6)') 'nfiles=',nfiles
c          do i=1,LtotT(1)
c            write(lu,'(i4$)') iNtypeT(i,1)
c          enddo
c          write(lu,'(a)') ' '
c          do i=1,LtotT(2)
c            write(lu,'(i4$)') iNtypeT(i,2)
c          enddo
c          write(lu,'(a)') ' '
c          close(lu)


            call chebsdPDDG(ioutstat,nfiles,
     +        aspT,bspT,JxT1,JyT1,JzT,nconnT,HxT,HyT,HzT1,
     +        LsysT,LbathT,LtotT,iRandT,nstepsT,tauT,init_stateT,
     +        autonormT,EpsFT,imflagT,ispecflagT, 
     +        iNtypeT,nooutputT,imarkT,
     +        vecResOut,iptVecResOut,mrkResOut)
          

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

