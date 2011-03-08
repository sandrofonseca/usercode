#!/bin/bash 


THETRUEWORKDIR=${WORKDIR}
THETRUEOUTDIR=${OUTDIR}
THETRUEPKGDIR=${PKGDIR}


cd $PKGDIR
eval `scramv1 runtime -sh`

cd $THETRUEWORKDIR
echo " pwd ->"
pwd


#cmsRun $PKGDIR/test/${1}.cfg >& ${1}.log
cmsRun $THETRUEPKGDIR/test/${1}.py >& ${1}.log
RETVAL=$?

#if [ $RETVAL != 0 ]; then
#  tar czf ${1}.log.tgz ${1}.log
#  mv -f *.tgz  $OUTDIR/
#fi

#mv -f ${1}.log $OUTDIR/
#mv -f *.root $OUTDIR/
mv -f ${1}.log $THETRUEOUTDIR/
mv -f *.root $THETRUEOUTDIR/
