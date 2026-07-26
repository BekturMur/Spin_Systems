        integer LENfnam
        parameter (LENfnam=32) ! max length of all filenames, barring logfnam
        integer LENline
        parameter (LENline=128)! max length of all lines in the input/data files
        character*16 logfnam
        character*15 fsxnam,fsynam,fsznam,fs12nam,frhonam,fnrmnam
	character*16 frdm1nam
        character*18 fs1sxnam,fs1synam,fs1sznam
        parameter (logfnam='Current/chsd.log') ! log file name
        parameter (fsxnam='Current/sxi.out')
        parameter (fsynam='Current/syi.out')
        parameter (fsznam='Current/szi.out')
        parameter (fs1sxnam='Current/sx1sxi.out')
        parameter (fs1synam='Current/sy1syi.out')
        parameter (fs1sznam='Current/sz1szi.out')
        parameter (fs12nam='Current/s12.out')
        parameter (frhonam='Current/rho.out')
        parameter (fnrmnam='Current/nrm.out')
        parameter (frdm1nam='Current/rdm1.out')
        integer chlog, chin, chcur
        integer chsx,chsy,chsz,chs1sx,chs1sy,chs1sz,chs12,chrho,chnrm
	integer chrdm1
        parameter (chlog=8)
        parameter (chin=10)
        parameter (chcur=11)
        parameter (chsx=21)
        parameter (chsy=22)
        parameter (chsz=23)
        parameter (chs1sx=24)
        parameter (chs1sy=25)
        parameter (chs1sz=26)
        parameter (chs12=27)
        parameter (chrho=28)
        parameter (chnrm=29)
        parameter (chrdm1=30)

        integer maxLtot, maxStat,maxLsys, maxNconn, maxLen
        parameter (maxLtot=24)
        parameter (maxStat=2**maxLtot)
        parameter (maxNconn=maxLtot*maxLtot)
        parameter (maxLsys=3)
        parameter (maxLen=128)

        integer maxLenVecResOut
        parameter (maxLenVecResOut=100000)

        double precision Eps
        parameter(Eps=1.0D-7)
        integer extraPrec
        parameter(extraPrec=0)
        integer nordMax ! max order of the Chebyshev expansion
        parameter (nordMax=6000)
