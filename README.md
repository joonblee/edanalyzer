# 0. EDAnalyzer

EDAnalyzer is a CMS framework to analyze MINIAOD foramt (or correspondings).
For more details, one can refer `https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule`.

# 1. Basic environments

Set CMSSW and create a skeleton EDAnalyzer module as follow:
```
export SCRAM_ARCH=slc7_amd64_gcc820
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src/
cmsenv
mkdir demo
cd demo
mkedanlzr DemoAnalyzer
cd DemoAnalyzer/
scram b
```
Note your working directory should be located directly under the `src` directory.
This automatically creates `plugins`, `python`, and `test` directories with `plugins/DemoAnalyzer.cc`.
To run this code, one need to make a python configuration file, e.g. ``python/ConfFile_cfg.py`, where it read like:
```
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source    = cms.Source("PoolSource",
                               fileNames = cms.untracked.vstring('file:/path/to/miniaod/file.root')
                               )
process.demo      = cms.EDAnalyzer('DemoAnalyzer',
                                   tracks = cms.untracked.InputTag('generalTracks')
                                   )

process.p = cms.Path(process.demo)
```
Now it is reaady to run.
```
cmsRun python/ConfFile_cfg.py
```
or equivalently one can submit crab jobs.
See below for crab submission.


# 2. Get a test MINIAOD file

A MINIAOD (.root) file can be downloaded at KNU or LXPLUS as follow:
```
voms-proxy-init --voms cms --valid 168:00
xrdcp root://cms-xrd-global.cern.ch///store/data/Run2016D/SingleMuon/MINIAOD/17Jul2018-v1/00000/2E6ADD2B-C68A-E811-AB71-0CC47AFC3C64.root ./
```
(This MINIAOD file is originated from `https://cmsweb.cern.ch/das/request?input=file%3D%2Fstore%2Fdata%2FRun2016D%2FSingleMuon%2FMINIAOD%2F17Jul2018-v1%2F00000%2F2E6ADD2B-C68A-E811-AB71-0CC47AFC3C64.root&instance=prod/global`.)
Import the file into the current directory.

If you want to study more about the MINIAOD format, one can check the MINIAOD twiki, `https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017`, or simply do `edmDumpEventContent /path/to/.root`, which prints out all contents in the root file.


# 3. Crab submission

One can clone 

or submit crab jobs:
```
cd CRAB

```
Please don't forget you should have voms proxy, if not do `voms-proxy-init --voms cms --valid 168:00`.





# Aggregation of output root files

ROOT_AGGREGATION


