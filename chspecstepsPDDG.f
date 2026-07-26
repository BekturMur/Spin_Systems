        subroutine chspecstepsPDDG(asp0,bsp0,Jx,Jy,Jz,nconn0,Hx,Hy,Hz,
     +    Lsys,Lbath,Ltot0,iRand,nsteps,tau,init_state,timeGl,ioutstat,
     +    psiR,psiI,autonorm,firstime,EpsF,phiR,phiI,dpsiR,dpsiI,
     +    ispecflag,nooutput, iNtype,
     +    vecResOut,iptVecResOut,mrkResOut)

        implicit none
        include "chsdpar.h"

c in/out variables

        double precision Jx(maxNconn),Jy(maxNconn),Jz(maxNconn)
        integer asp0(maxNconn), bsp0(maxNconn)
        double precision Hx(maxLtot), Hy(maxLtot), Hz(maxLtot)
        double precision tau, timeGl, EpsF
        integer autonorm, firstime
        integer Lsys,Lbath,Ltot0,nconn0,nsteps,iRand,init_state,ioutstat
        double precision psiI(0:maxStat-1),psiR(0:maxStat-1)

        double precision phiI(0:maxStat-1),phiR(0:maxStat-1)
        double precision dpsiI(0:maxStat-1),dpsiR(0:maxStat-1)

        integer ispecflag

        integer nooutput
        integer iNtype(maxLtot)

        double precision vecResOut(maxLenVecResOut)
        integer iptVecResOut
        integer mrkResOut(maxLenVecResOut)

c internal variables

        integer nconn, Ltot, asp(maxNconn),bsp(maxNconn) 
        double precision jefx(maxNconn),jefy(maxNconn),jefz(maxNconn)
        double precision hefx(maxLtot),hefy(maxLtot),hefz(maxLtot)
        integer nstates, nstat1,nstat1max,i,k,iconn,j,jp, l,i0,ip0,ip1
        integer j0,jp0,jp1,j1,j2,k1,l1
        integer istep,ia,ib, iord 
c        double precision psiI(0:maxStat-1),psiR(0:maxStat-1)
c        double precision phiI(0:maxStat-1),phiR(0:maxStat-1)
c        double precision dpsiI(0:maxStat-1),dpsiR(0:maxStat-1)
c        double precision pbufI(0:maxStat-1),pbufR(0:maxStat-1)
c        double precision dpbI(0:maxStat-1),dpbR(0:maxStat-1)
        double precision xnorm,xnR,xnI, x,y,z
        double precision HmR(2**maxLsys,2**maxLsys)
        double precision HmI(2**maxLsys,2**maxLsys)

        double complex HmT(2**maxLsys,2**maxLsys)
        double complex work((2**maxLsys)*((2**maxLsys)+1))
        double precision w(2**maxLsys),rwork(3*(2**maxLsys)-2)
        integer lwork

        double precision vecR(2**maxLsys,2**maxLsys)
        double precision vecI(2**maxLsys,2**maxLsys)
        double precision tPre, emax,alp,atab(0:nordMax)
        integer nR,nI, nord,ncalc
        logical coeflarge

        double precision xrnd

        double precision dnrm2,ddot
        external dnrm2,ddot
        external daxpy,dscal,dcopy,dswap,ribesl
        external steprnPDDG,addJSxxSyy,addJSzz,addHSx,addHSy,addHSz
        external actSx,actSy,actSz

        double precision n1x,n1y,n1z,n2x,n2y,n2z,th1,th2,ct1,st1,ct2,st2
        double precision ynR,ynI, ar,ai,br,bi,cr,ci,dr,di

c executable statements

        nstates = 2**Ltot0
        Ltot = Ltot0
        nconn = nconn0
c        do i=0,maxStat-1 
c          psiI(i)=psi0I(i)
c          psiR(i)=psi0R(i)
c        enddo
        do i=1,maxNconn
          asp(i) = asp0(i)
          bsp(i) = bsp0(i)
        enddo
        do i=1,2**maxLsys
          do k=1,2**maxLsys
            HmR(k,i) = 0
            HmI(k,i) = 0
          enddo
        enddo

        if (init_state.ne.-10) then ! if initial state is going to be (re)set
          do i=0,nstates-1
            psiR(i) = 0.0D0
            psiI(i) = 0.0D0
          enddo
        endif

        if (init_state.eq.-1) then ! SPIN UP: default value
          psiR(nstates-1) = 1
        else if (init_state.eq.-2) then ! SPIN DOWN
          psiR(0) = 1
        else if (init_state.eq.-3) then ! RANDOM
          do i=0,nstates-1
            xrnd=rand()
            psiR(i) = 2*xrnd-1
          enddo
          do i=0,nstates-1
            xrnd=rand()
            psiI(i) = 2*xrnd-1
          enddo
          xnR = dnrm2(nstates,psiR(0),1)
          xnI = dnrm2(nstates,psiI(0),1)
          xnorm = 1.0D0/sqrt(xnR*xnR + xnI*xnI)
          call dscal(nstates,xnorm,psiR(0),1)
          call dscal(nstates,xnorm,psiI(0),1)
        else if (init_state.eq.-4) then ! GROUND STATE
          nstat1 = 2**Lsys
          nstat1max = 2**maxLsys
          do iconn=1,nconn
            if (asp(iconn).le.Lsys.and.bsp(iconn).le.Lsys) then
              x=-0.25D0*Jx(iconn)
              y=-0.25D0*Jy(iconn)
              z=-0.25D0*Jz(iconn)
              do k=1,nstat1
                psiR(k-1)=1
                call actSx(psiR,psiI,asp(iconn),nstat1,dpsiR,dpsiI)
                call actSx(dpsiR,dpsiI,bsp(iconn),nstat1,phiR,phiI)
                call daxpy(nstat1,x,phiR(0),1,HmR(1,k),1)
                call daxpy(nstat1,x,phiI(0),1,HmI(1,k),1)
                call actSy(psiR,psiI,asp(iconn),nstat1,dpsiR,dpsiI)
                call actSy(dpsiR,dpsiI,bsp(iconn),nstat1,phiR,phiI)
                call daxpy(nstat1,y,phiR(0),1,HmR(1,k),1)
                call daxpy(nstat1,y,phiI(0),1,HmI(1,k),1)
                call actSz(psiR,psiI,asp(iconn),nstat1,dpsiR,dpsiI)
                call actSz(dpsiR,dpsiI,bsp(iconn),nstat1,phiR,phiI)
                call daxpy(nstat1,z,phiR(0),1,HmR(1,k),1)
                call daxpy(nstat1,z,phiI(0),1,HmI(1,k),1)
                psiR(k-1)=0
              enddo
            endif
          enddo
          do i=1,Lsys 
            x=-0.5D0*Hx(i)
            y=-0.5D0*Hy(i)
            z=-0.5D0*Hz(i)
            do k=1,nstat1
              psiR(k-1)=1
              call actSx(psiR,psiI,i,nstat1,phiR,phiI)
              call daxpy(nstat1,x,phiR(0),1,HmR(1,k),1)
              call daxpy(nstat1,x,phiI(0),1,HmI(1,k),1)
              call actSy(psiR,psiI,i,nstat1,phiR,phiI)
              call daxpy(nstat1,y,phiR(0),1,HmR(1,k),1)
              call daxpy(nstat1,y,phiI(0),1,HmI(1,k),1)
              call actSz(psiR,psiI,i,nstat1,phiR,phiI)
              call daxpy(nstat1,z,phiR(0),1,HmR(1,k),1)
              call daxpy(nstat1,z,phiI(0),1,HmI(1,k),1)
              psiR(k-1)=0
            enddo
          enddo

          lwork=((2**maxLsys)*((2**maxLsys)+1))
          do k=1,nstat1
          do i=1,nstat1
            HmT(i,k)=dcmplx(HmR(i,k),HmI(i,k))
          enddo
          enddo
          call zheev('V','U',nstat1,HmT,nstat1max,w,work,lwork,rwork,
     +      ioutstat)
          if (ioutstat.ne.0) then
c            open(chlog,file=logfnam,access='append')
c            write(chlog,*) 'ZHEEV failed, error #', ioutstat
c            close(chlog)
            ioutstat= -201
            return
          endif
          do k=1,nstat1
            vecR(k,1) = dble(HmT(k,1))
            vecI(k,1) = dimag(HmT(k,1))
          enddo

c          jp=max(nstates-1-2**Lsys,0)
          jp=2**Lbath-1
          do k=0,jp
            xrnd=rand()
            dpsiR(k)=2*xrnd-1
          enddo
          do k=0,jp
            xrnd=rand()
            dpsiI(k)=2*xrnd-1
          enddo
          jp = 2**Lbath
          xnR = dnrm2(jp,dpsiR(0),1)
          xnI = dnrm2(jp,dpsiI(0),1)
          xnorm = 1.0D0/sqrt(xnR*xnR + xnI*xnI)
          call dscal(jp,xnorm,dpsiR(0),1)
          call dscal(jp,xnorm,dpsiI(0),1)
          j=1
          jp=0
          do k=0,nstates-1
            psiR(k)= vecR(j,1)*dpsiR(jp) - vecI(j,1)*dpsiI(jp)
            psiI(k)= vecR(j,1)*dpsiI(jp) + vecI(j,1)*dpsiR(jp)
            j=j+1
            if (j.gt.nstat1) then  ! nstat1= 2**Lsys
              j=1
              jp=jp+1
            endif
          enddo
        else if (init_state.ge.0) then
          psiR(init_state)=1
        endif
        if (firstime.ne.0.and.nooutput.eq.0) then
          call steprnPDDG(timeGl,psiR,psiI,Lsys,Ltot,nstates,
     +      phiR,phiI,dpsiR,dpsiI,vecR,vecI,asp,bsp,Jx,Jy,Jz,
     +      Hx,Hy,Hz,nconn,
     +      vecResOut,iptVecResOut,mrkResOut,ioutstat)
          if (ioutstat.ne.0) then
          return
          endif
        endif

c
c psi(t=0) initialized, let's go on
c

c pi-pulse along X-axis
        if (ispecflag.eq.1) then ! pi-pulse along X
          do l=0,nstates-1,2
            xnR=psiR(l)
            xnI=psiI(l)
            psiR(l)=psiI(l+1)
            psiI(l)=-psiR(l+1)
            psiR(l+1)=xnI
            psiI(l+1)=-xnR
          enddo
        else if (ispecflag.eq.2) then ! pi-pulse along Y
          do l=0,nstates-1,2
            xnR=psiR(l)
            xnI=psiI(l)
            psiR(l)=psiR(l+1)
            psiI(l)=psiI(l+1)
            psiR(l+1)=-xnR
            psiI(l+1)=-xnI
          enddo
        else if (ispecflag.eq.3) then ! pi-pulse along Z
          do l=0,nstates-1,2
            xnR=psiR(l)
            psiR(l)=-psiI(l)
            psiI(l)=xnR
            xnR=psiR(l+1)
            psiR(l+1)=psiI(l+1)
            psiI(l+1)=-xnR
          enddo

        else if (ispecflag.eq.101) then ! pi-pulse along X, spin 2
          do l=0,nstates-1,4
            xnR=psiR(l)
            xnI=psiI(l)
            ynR=psiR(l+1)
            ynI=psiI(l+1)
            psiR(l)=psiI(l+2)
            psiI(l)=-psiR(l+2)
            psiR(l+2)=xnI
            psiI(l+2)=-xnR
            psiR(l+1)=psiI(l+3)
            psiI(l+1)=-psiR(l+3)
            psiR(l+3)=ynI
            psiI(l+3)=-ynR
          enddo
        else if (ispecflag.eq.102) then ! pi-pulse along Y, spin 2
          do l=0,nstates-1,4
            xnR=psiR(l)
            xnI=psiI(l)
            ynR=psiR(l+1)
            ynI=psiI(l+1)
            psiR(l)=psiR(l+2)
            psiI(l)=psiI(l+2)
            psiR(l+2)=-xnR
            psiI(l+2)=-xnI
            psiR(l+1)=psiR(l+3)
            psiI(l+1)=psiI(l+3)
            psiR(l+3)=-ynR
            psiI(l+3)=-ynI
          enddo

        else if (ispecflag.eq.11) then ! EPR1
          do l=0,nstates-1,4
            psiR(l)=-psiR(l)
            psiI(l)=-psiI(l)
            psiR(l+3)=-psiR(l+3)
            psiI(l+3)=-psiI(l+3)
            xnR=psiR(l+1)
            xnI=psiI(l+1)
            psiR(l+1)=psiR(l+2)
            psiI(l+1)=psiI(l+2)
            psiR(l+2)=xnR
            psiI(l+2)=xnI
          enddo
        else if (ispecflag.eq.12) then ! EPR2
          do l=0,nstates-1,4
            psiR(l)=-psiR(l)
            psiI(l)=-psiI(l)
            psiR(l+3)=-psiR(l+3)
            psiI(l+3)=-psiI(l+3)
            xnR=psiR(l+1)
            xnI=psiI(l+1)
            psiR(l+1)=-psiR(l+2)
            psiI(l+1)=-psiI(l+2)
            psiR(l+2)=-xnR
            psiI(l+2)=-xnI
          enddo

        else if (ispecflag.eq.13) then ! EPR12, P2 - new
          ar = 0.25D0
          ai = dsqrt(3.0D0)*0.25D0
          br = -0.75D0
          bi = dsqrt(3.0D0)*0.25D0
          cr = 0.25D0
          ci = -dsqrt(3.0D0)*0.25D0
          dr = 0.75D0
          di = dsqrt(3.0D0)*0.25D0
          do l=0,nstates-1,4
            xnR = psiR(l)
            xnI = psiI(l)
            ynR = psiR(l+3)
            ynI = psiI(l+3)
            psiR(l) = xnR*ar + ynR*br - xnI*ai - ynI*bi
            psiI(l) = xnR*ai + ynR*bi + xnI*ar + ynI*br
            psiR(l+3) = xnR*br + ynR*ar - xnI*bi - ynI*ai
            psiI(l+3) = xnR*bi + ynR*ai + xnI*br + ynI*ar

            xnR = psiR(l+1)
            xnI = psiI(l+1)
            ynR = psiR(l+2)
            ynI = psiI(l+2)
            psiR(l+1) = xnR*cr + ynR*dr - xnI*ci - ynI*di
            psiI(l+1) = xnR*ci + ynR*di + xnI*cr + ynI*dr
            psiR(l+2) = xnR*dr + ynR*cr - xnI*di - ynI*ci
            psiI(l+2) = xnR*di + ynR*ci + xnI*dr + ynI*cr
          enddo

        else if (ispecflag.eq.14) then ! EPR2 (singlet), P1 - new
          do l=0,nstates-1,4
            xnR=psiR(l+1)
            xnI=psiI(l+1)
            psiR(l+1)=psiR(l+2)
            psiI(l+1)=psiI(l+2)
            psiR(l+2)=xnR
            psiI(l+2)=xnI
          enddo

        else if (ispecflag.eq.15) then ! All EPRs, P4 - new
          do l=0,nstates-1,4
            xnR = psiR(l)
            xnI = psiI(l)
            ynR = psiR(l+3)
            ynI = psiI(l+3)
            psiR(l) = -0.5D0*(xnR+ynR) - 0.5D0*(xnI-ynI)
            psiI(l) = -0.5D0*(xnI+ynI) + 0.5D0*(xnR-ynR)
            psiR(l+3) = -0.5D0*(xnR+ynR) - 0.5D0*(ynI-xnI)
            psiI(l+3) = -0.5D0*(xnI+ynI) + 0.5D0*(ynR-xnR)

            xnR = psiR(l+1)
            xnI = psiI(l+1)
            ynR = psiR(l+2)
            ynI = psiI(l+2)
            psiR(l+1) = 0.5D0*(xnR+ynR) + 0.5D0*(xnI-ynI)
            psiI(l+1) = 0.5D0*(xnI+ynI) - 0.5D0*(xnR-ynR)
            psiR(l+2) = 0.5D0*(xnR+ynR) + 0.5D0*(ynI-xnI)
            psiI(l+2) = 0.5D0*(xnI+ynI) - 0.5D0*(ynR-xnR)
          enddo

        else if (ispecflag.eq.20.or.ispecflag.eq.21) then ! random rotation of spins 1 and 2
          n1x= -0.826213643123139D0
          n1y= -0.501191209212064D0
          n1z= 0.257251603932297D0
          th1=0.834156688384137D0
          st1=dsin(0.5D0*th1)
          ct1=dcos(0.5D0*th1)
          n2x= -0.507460929675091D0 
          n2y= -0.152638755808748D0
          n2z=-0.848047649061331D0
          th2=4.26670451769502D0
          st2=dsin(0.5D0*th2)
          ct2=dcos(0.5D0*th2)
          if (ispecflag.eq.21) then
            st1=-st1
            st2=-st2
          endif
          do l=0,nstates-1,2 ! rotating spin 1
            xnR=psiR(l)
            xnI=psiI(l)
            ynR=psiR(l+1)
            ynI=psiI(l+1)
            psiR(l) = ct1*xnR - n1z*st1*xnI + n1y*st1*ynR + n1x*st1*ynI
            psiI(l) = ct1*xnI + n1z*st1*xnR + n1y*st1*ynI - n1x*st1*ynR
            psiR(l+1) = ct1*ynR + n1z*st1*ynI - n1y*st1*xnR +n1x*st1*xnI
            psiI(l+1) = ct1*ynI - n1z*st1*ynR - n1y*st1*xnI -n1x*st1*xnR
          enddo
          do k=0,nstates-1,4 ! rotating spin 2
          do l=k,k+1       
            xnR=psiR(l)
            xnI=psiI(l)
            ynR=psiR(l+2)
            ynI=psiI(l+2)
            psiR(l) = ct2*xnR - n2z*st2*xnI + n2y*st2*ynR + n2x*st2*ynI
            psiI(l) = ct2*xnI + n2z*st2*xnR + n2y*st2*ynI - n2x*st2*ynR
            psiR(l+2) = ct2*ynR + n2z*st2*ynI - n2y*st2*xnR +n2x*st2*xnI
            psiI(l+2) = ct2*ynI - n2z*st2*ynR - n2y*st2*xnI -n2x*st2*xnR
          enddo
          enddo
        else
          ioutstat = -4444
c          open(chlog,file=logfnam,access='append')
c          write(chlog,*) 'Wrong value of ispecflag:', ispecflag
c          close(chlog)
          return
        endif

c        i0 = 2
c        j0 = 1
c
c        ip0 = 2**(i0-1)
c        ip1 = 2*ip0  ! 2**i0
c        jp0 = 2**(j0-1)
c        jp1 = 2*jp0
c        do l=0,nstates-1,ip1
c          do k=l,l+ip0-1,jp1
c            do k1=k,k+jp0-1
c              j = k1+ip0
c              j1 = k1+jp0
c              j2 = k1+ip0+jp0
c
c              xnR = psiR(j2)
c              xnI = psiI(j2)
c              psiR(j2) = psiR(k1)
c              psiR(k1) = xnR
c              psiI(j2) = psiI(k1)
c              psiI(k1) = xnI
c            enddo
c          enddo
c        enddo

        if (nooutput.eq.0) then
          call steprnPDDG(timeGl,psiR,psiI,Lsys,Ltot,nstates,
     +      phiR,phiI,dpsiR,dpsiI,vecR,vecI,asp,bsp,Jx,Jy,Jz,
     +      Hx,Hy,Hz,nconn,
     +      vecResOut,iptVecResOut,mrkResOut,ioutstat)
          if (ioutstat.ne.0) then
          return
          endif
        endif

        return
        end
