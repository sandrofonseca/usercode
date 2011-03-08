#!/bin/csh 

setenv PKGDIR $CMSSW_BASE/src/APAnalyzers/FastMuonAnalyzer/
setenv OUTDIR $PKGDIR/test/RootHisto_`date +%Y%m%d%H%M`

[ -d $OUTDIR ] || mkdir $OUTDIR

#setenv QUEUE 1nd
setenv QUEUE 8nh
#setenv QUEUE 1nh
#setenv QUEUE 8nm
setenv RUNSCRIPT run.sh

  setenv MYJOB1 NewMuonAnalyzer_FASTfromFILE
  bsub -q $QUEUE -oo $OUTDIR/$MYJOB1.out $RUNSCRIPT $MYJOB1
  setenv MYJOB2 NewMuonAnalyzer_FULL_HLT
  bsub -q $QUEUE -oo $OUTDIR/$MYJOB2.out $RUNSCRIPT $MYJOB2
  setenv MYJOB3 NewMuonAnalyzer_FULL_noHLT
  bsub -q $QUEUE -oo $OUTDIR/$MYJOB3.out $RUNSCRIPT $MYJOB3
#  setenv MYJOB4 NewMuonAnalyzer_FULL
#  bsub -q $QUEUE -oo $OUTDIR/$MYJOB4.out $RUNSCRIPT $MYJOB4
#  setenv MYJOB5 NewFastMuonProducerAndAnalyzerFile_cfg
#  bsub -q $QUEUE -oo $OUTDIR/$MYJOB5.out $RUNSCRIPT $MYJOB5
