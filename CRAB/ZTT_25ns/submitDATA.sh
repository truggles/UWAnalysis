#!/bin/sh
voms-proxy-init --voms cms --valid 100:00

#Submit the reprocessed Data
farmoutAnalysisJobs  $1 --skip-existing-output --vsize-limit=8000 --input-files-per-job=1 --input-dbs-path=/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD  SINGLEMUON2015D $CMSSW_BASE $CMSSW_BASE/src/UWAnalysis/CRAB/ZTT_25ns/SUB-JSON.py  

farmoutAnalysisJobs  $1 --skip-existing-output --vsize-limit=8000 --input-files-per-job=1 --input-dbs-path=/SingleElectron/Run2015D-16Dec2015-v1/MINIAOD  SINGLEELE2015D $CMSSW_BASE $CMSSW_BASE/src/UWAnalysis/CRAB/ZTT_25ns/SUB-JSON.py  

