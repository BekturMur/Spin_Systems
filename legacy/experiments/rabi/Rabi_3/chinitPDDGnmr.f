        subroutine chinitPDDG(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +    Lsys,Lbath,Ltot,iRand,init_state,ioutstat,
     +    psiR,psiI,dpsiR,dpsiI)

        implicit none
        include "chsdpar.h"
        double precision PIc, MDLT, MDLT2
        parameter (PIc=3.14159265358979323846D0)
        parameter (MDLT=1.0D-15, MDLT2=1.0D-7)

c in/out variables

        double precision Jx(maxNconn),Jy(maxNconn),Jz(maxNconn)
        integer asp(maxNconn), bsp(maxNconn)
        double precision Hx(maxLtot), Hy(maxLtot), Hz(maxLtot)
        integer Lsys,Lbath,Ltot,nconn,iRand,init_state,ioutstat
        double precision psiI(0:maxStat-1),psiR(0:maxStat-1)
        double precision dpsiI(0:maxStat-1),dpsiR(0:maxStat-1)

c internal variables

        integer nstates, nstat1,nstat1max,i,k,iconn,j,jp, istat
        integer ia,ib 
        double precision xnorm,xnR,xnI, x,y,z
        double precision HmR(2**maxLsys,2**maxLsys)
        double precision HmI(2**maxLsys,2**maxLsys)

        double complex HmT(2**maxLsys,2**maxLsys)
        double complex work((2**maxLsys)*((2**maxLsys)+1))
        double precision w(2**maxLsys),rwork(3*(2**maxLsys)-2)
        integer lwork
        double precision vecR(2**maxLsys,2**maxLsys)
        double precision vecI(2**maxLsys,2**maxLsys)

        double precision xrnd

        double precision dnrm2,ddot
        external dnrm2,ddot
        external daxpy,dscal,dcopy,dswap,rjbesl
        external addJSxyz,addHSx,addHSy,addHSz

c executable statements

        ioutstat=0

        nstates = 2**Ltot
        do i=1,2**maxLsys
          do k=1,2**maxLsys
            HmR(k,i) = 0
            HmI(k,i) = 0
          enddo
        enddo
        if (init_state.eq.-10) return ! if initial state will not be (re)set

        do i=0,nstates-1
          psiR(i) = 0.0D0
          psiI(i) = 0.0D0
        enddo

        if (init_state.eq.-1) then ! SPIN UP: default value
          psiR(nstates-1) = 1
        else if (init_state.eq.-2) then ! SPIN DOWN
          psiR(0) = 1
        else if (init_state.eq.-3) then ! RANDOM
          xnorm=0.0D0
          do istat=0,nstates-1
            x=MDLT+(1.0D0-2*MDLT)*rand()
            y=rand()
            z=dlog(x)
            xnR=dsqrt(-2.0D0*z)*dcos(2*PIc*y)
            xnI=dsqrt(-2.0D0*z)*dsin(2*PIc*y)
            psiR(istat)=xnR
            psiI(istat)=xnI
            xnorm=xnorm + xnR*xnR + xnI*xnI
          enddo
          xnorm=1.0D0/dsqrt(xnorm)
          do istat=0,nstates-1
            psiR(istat)=psiR(istat)*xnorm
            psiI(istat)=psiI(istat)*xnorm
          enddo
        else if (init_state.eq.-4) then ! GROUND STATE
          nstat1 = 2**Lsys
          nstat1max = 2**maxLsys
          do iconn=1,nconn
            if (asp(iconn).le.Lsys.and.bsp(iconn).le.Lsys) then
              x=-0.25D0*Jx(iconn)
              y=-0.25D0*Jy(iconn)
              z=-0.25D0*Jz(iconn)
              ib=bsp(iconn)
              ia=asp(iconn)
              do k=1,nstat1
                psiR(k-1)=1
           call addJSxyz(x,y,z,psiR,psiI,ib,ia,nstat1,HmR(1,k),HmI(1,k)) ! ib > ia, since asp < bsp
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
              call addHSx(x,psiR,psiI,i,nstat1,HmR(1,k),HmI(1,k))
              call addHSy(y,psiR,psiI,i,nstat1,HmR(1,k),HmI(1,k))
              call addHSz(z,psiR,psiI,i,nstat1,HmR(1,k),HmI(1,k))
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
            ioutstat= -201
            return
          endif
          do k=1,nstat1
            vecR(k,1) = dble(HmT(k,1))
            vecI(k,1) = dimag(HmT(k,1))
          enddo

          jp = 2**Lbath
          xnorm=0.0D0
          do istat=0,jp-1
            x=MDLT+(1.0D0-2*MDLT)*rand()
            y=rand()
            z=dlog(x)
            xnR=dsqrt(-2.0D0*z)*dcos(2*PIc*y)
            xnI=dsqrt(-2.0D0*z)*dsin(2*PIc*y)
            dpsiR(istat)=xnR
            dpsiI(istat)=xnI
            xnorm=xnorm + xnR*xnR + xnI*xnI
          enddo
          xnorm=1.0D0/dsqrt(xnorm)
          do istat=0,jp-1
            dpsiR(istat)=dpsiR(istat)*xnorm
            dpsiI(istat)=dpsiI(istat)*xnorm
          enddo

          do jp=0,nstates-1,nstat1
            do j=1,nstat1    ! because in vec(j,1) index j starts from 1, goes to nstat1
              k=jp+j-1       ! because j is from 1 to nstat1
              psiR(k)= vecR(j,1)*dpsiR(jp) - vecI(j,1)*dpsiI(jp)
              psiI(k)= vecR(j,1)*dpsiI(jp) + vecI(j,1)*dpsiR(jp)
            enddo
          enddo
        else if (init_state.ge.0.and.init_state.lt.nstates) then
          psiR(init_state)=1
        else
          ioutstat=-83 ! bad init_state value, same ioutstat value as in chparsfPDDG.f
        endif
c        if (firstime.ne.0.and.nooutput.eq.0) then
c          call steprnPDDG(timeGl,psiR,psiI,Lsys,Ltot,nstates,
c     +      phiR,phiI,dpsiR,dpsiI,vecR,vecI,asp,bsp,Jx,Jy,Jz,
c     +      Hx,Hy,Hz,nconn,
c     +      vecResOut,iptVecResOut,mrkResOut,ioutstat)
c          if (ioutstat.ne.0) then
c          return
c          endif
        return
        end
