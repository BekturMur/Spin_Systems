        subroutine gendipnmr2D(lociRand,locLtot,locnc,xJz,locasp,locbsp)
        implicit none
        include "chsdpar.h"

c in/out variables

        integer lociRand,locLtot,locnc,locasp(maxNconn),locbsp(maxNconn)
        double precision xJz(maxNconn)

c internal variables

        integer i,isp,jsp,iconn
        double precision Dx,Dy,dMin,x,y,z,xx,xsp(maxLtot),ysp(maxLtot)
        logical xyok

c executable statements

        Dx = 0.5D0*dsqrt(1.0D0*locLtot)
        Dy = 0.5D0*dsqrt(1.0D0*locLtot)
        dMin = 5.0D-2
        do i=1,lociRand
          xx=rand()
        enddo

        do isp=1,locLtot
          xsp(isp) = 0.0D0
          ysp(isp) = 0.0D0
        enddo
        iconn=0
        do isp=2,locLtot
          xyok = .false.
          do while(.not.xyok)
            x = Dx*(2*rand()-1.0D0)
            y = Dy*(2*rand()-1.0D0)
            xyok = .true.
            do jsp=1,isp-1
              if (dabs(x-xsp(jsp)).lt.dMin.and.dabs(y-ysp(jsp))
     +         .lt.dMin) xyok = .false.
            enddo
          enddo
          xsp(isp) = x
          ysp(isp) = y
          do jsp=1,isp-1
            xx=(x-xsp(jsp))*(x-xsp(jsp))
            xx=xx+(y-ysp(jsp))*(y-ysp(jsp))
            iconn=iconn+1
            xJz(iconn)=1.0D0/(xx*dsqrt(xx))
            locbsp(iconn)=isp
            locasp(iconn)=jsp ! asp < bsp
          enddo
        enddo
        locnc=iconn
        end
