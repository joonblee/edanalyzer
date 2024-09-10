import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options   = cms.untracked.PSet( 
  SkipEvent = cms.untracked.vstring('ProductNotFound'),
  wantSummary = cms.untracked.bool(True) 
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            #fileNames = cms.untracked.vstring('file:/data6/Users/joonblee/MINIAODAnalyzer/CMSSW_10_2_18/src/Demo/DemoAnalyzer/Inputs/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/C24D5677-8BE3-E811-AF3D-001E677927B2.root')
                            #fileNames = cms.untracked.vstring('file:/data6/Users/joonblee/MINIAODAnalyzer/CMSSW_10_2_18/src/Demo/DemoAnalyzer/Inputs/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic.root')
                            #fileNames = cms.untracked.vstring('file:/data6/Users/joonblee/MINIAODAnalyzer/CMSSW_10_2_18/src/Demo/DemoAnalyzer/Inputs/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic.root')
                            #fileNames = cms.untracked.vstring('file:/data6/Users/joonblee/MINIAODAnalyzer/CMSSW_10_2_18/src/Demo/DemoAnalyzer/Inputs/MET_Run2016C-17Jul2018-v1/CCB57FEE-278A-E811-BFAC-0CC47A57CC42.root')
                            fileNames = cms.untracked.vstring('file:/data6/Users/joonblee/MINIAODAnalyzer/CMSSW_10_2_18/src/Demo/DemoAnalyzer/Inputs/SingleMuon_C/00408080-0F98-E811-A7E3-34E6D7BEAF01.root')
                            )

process.TFileService = cms.Service("TFileService", 
                                   #fileName = cms.string("/afs/cern.ch/work/j/joon/private/MINIAODAnalyzer/CMSSW_10_2_18/src/Demo/DemoAnalyzer/Outputs/hist.root"),
                                   fileName = cms.string("./hist.root"),
                                   closeFileFast = cms.untracked.bool(True)
)

process.demo = cms.EDAnalyzer('DemoTrigEff',
															## -- Physics Object -- ##
                              pfCands = cms.untracked.InputTag("packedPFCandidates"),
                              Jet = cms.untracked.InputTag("slimmedJets"),
                              Muon = cms.untracked.InputTag("slimmedMuons"),
                              Electron = cms.untracked.InputTag("slimmedElectrons"),
                              Photon = cms.untracked.InputTag("slimmedPhotons"),
                              MET = cms.untracked.InputTag("slimmedMETs"),
                              PrimaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
															## -- Trigger -- ##
                              Triggers = cms.untracked.InputTag("TriggerResults","","HLT"),
                              L1Muon = cms.untracked.InputTag("gmtStage2Digis"),
                              L3obj = cms.untracked.InputTag("slimmedPatTrigger"),
															## -- Gen Event -- ##
                              #GenEventInfo = cms.untracked.InputTag("generator"),
                             )

process.p = cms.Path(process.demo)
