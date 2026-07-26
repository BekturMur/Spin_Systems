       program tstsxi
       implicit double precision (a-h,o-z)
       parameter (Nspin=21)
       parameter (Nrho=5)
       parameter (Ns12=6)
       parameter (Ntime=178)
       integer itime, ispin, nSpCur, indDsMax, indTotMax, 
     +       ispTotMax, ifilTotMax
           integer ios
       double precision time, s1, s2, ds, dsMax, totMax
       character*128 nam1, nam2
       character*63 namDir1, namDir2, namsx(8)
       character*60 lj(19)
       logical exOutFile
c
c   character*128 namsx1, namsx2, namsx3, namsx4, namsx5, namsx6, namsx7, namsx8
c
       dimension s1(21), s2(21), dsMax(21), indDsMax(21)

       namsx(1) = 'sxi.out'
       namsx(2) = 'syi.out'
       namsx(3) = 'szi.out'
       namsx(4) = 'sx1sxi.out'
       namsx(5) = 'sy1syi.out'
       namsx(6) = 'sz1szi.out'
       namsx(7) = 's12.out'
       namsx(8) = 'rho.out'

       namDir1 = 'mi/Current/'
       namDir2 = '../mi/Current/'

       if (Nspin.gt.21) then
         write(*,*) 'Nspin > 21, change this code'
         stop
       endif

       totMax = -1.0
       ifilTotMax = -1
       ispTotMax = -1
       indTotMax = -1

       inquire(file='tstsx.dif', exist=exOutFile)
       if (exOutFile) then
         open(23, file='tstsx.dif')
         close(23, status='delete')
       endif

       do ispin=1,21
         dsMax(ispin) = -1.0
         indDsMax(ispin) = -1
       enddo
       do ifil = 1, 3
         do ispin=1,21
           dsMax(ispin) = -1.0
           indDsMax(ispin) = -1
         enddo
         nam1 = (namDir1(:len_trim(namDir1)))//namsx(ifil)
         nam2 = (namDir2(:len_trim(namDir2)))//namsx(ifil)
         write(*,*) nam1
         write(*,*) nam2
         open(21, file=(nam1(:len_trim(nam1))))
         open(22, file=(nam2(:len_trim(nam2))))
         write(*,*) ifil

         nSpCur=Nspin
         if (ifil.eq.7) nSpCur= Ns12
         if (ifil.eq.8) nSpCur = Nrho

         do itime = 1, Ntime
           read(21,*,iostat=ios) time, (s1(ispin), ispin=1,nSpCur)
           if (ios.ne.0) then
             write(*,*) 'error',nam1,'  itime',itime
             stop
           endif
           read(22,*,iostat=ios) time, (s2(ispin), ispin=1,nSpCur)
           if (ios.ne.0) then
             write(*,*) 'error',nam2,'  itime',itime
             stop
           endif
           do ispin = 1, nSpCur
             ds = dabs(s1(ispin)-s2(ispin))
             if (ds.gt.dsMax(ispin)) then
               dsMax(ispin) = ds
               indDsMax(ispin) = itime
               if (ds.gt.totMax) then
                 totMax = ds
                 ifilTotMax = ifil
                 ispTotMax = ispin
                 indTotMax = itime
               endif
             endif
           enddo
         enddo
         close(21)
         close(22)
         open(unit=23, file='tstsx.dif', access='append')
         write(23,*) (namsx(ifil)(:len_trim(namsx(ifil)))),
     +           '    file #',ifil
         write(23,*) (indDsMax(ispin), ispin=1,nSpCur)
         write(23,*) (dsMax(ispin), ispin=1,nSpCur)
         write(23,*) ('-', ispin=1,20)
         close(23)
       enddo
       write(*,*) totMax, ifilTotMax, ispTotMax, indTotMax
       end
