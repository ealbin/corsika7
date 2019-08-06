      implicit none

      integer nnx,nny,nx,ny,ixmin,ixmax,iymin,iymax,ipr,irun,ncol,icol

      real    xmin,xmax,ymin,ymax,xconst,yconst,xxconst,yyconst,
     *        cut,cutem,cutmu,cuthad,radi,higr

      character*30 ich,dsn,dsn2 

      parameter(nnx=500,nny=1000)
      integer ifeld(nnx,nny)

      common /work/ nx,ny,ifeld,ixmin,ixmax,iymin,iymax,
     *              ipr,irun,ncol,icol,
     *              xmin,xmax,ymin,ymax,xconst,yconst,xxconst,yyconst,
     *              cut,cutem,cutmu,cuthad,radi,higr,
     *              ich,dsn,dsn2
