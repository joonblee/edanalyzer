import CRABClient,os
from shutil import copyfile
from CRABClient.UserUtilities import config 
config = config()

analysis = 'trigeff'
date = '09Sep2024'
version = '1'

config.General.requestName = ''
config.General.workArea = analysis+'_'+date+'_v'+version
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/trigeff_cfg.py'

config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False
config.Data.outputDatasetTag = 'edanalyzer_'+analysis

config.Site.storageSite = 'T3_KR_KNU'

# 'MultiCRAB' part

GoldenJSON = './Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
# MuonPhysJSON = '/u/user/kplee/JSON/Run2016/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'
StartRun = 271036
EndRun = 284044


if __name__ == '__main__':

  from CRABAPI.RawCommand import crabCommand

  # -- MET phi correction for MC -- #
  #src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pfMETmultShiftCorrections_MC_cfi.py')
  #dst= os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','..','..','JetMETCorrections/Type1MET/python/pfMETmultShiftCorrections_cfi.py')
  #copyfile(src,dst)

  # -- QCD Pt-binned -- #
  config.General.requestName = 'QCDMuEnriched_Pt15to20'
  config.Data.inputDataset = '/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config) #
  # 4141251

  config.General.requestName = 'QCDMuEnriched_Pt20to30'
  config.Data.inputDataset = '/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 31878740

  config.General.requestName = 'QCDMuEnriched_Pt30to50'
  config.Data.inputDataset = '/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 29954815

  config.General.requestName = 'QCDMuEnriched_Pt50to80'
  config.Data.inputDataset = '/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 19662175

  config.General.requestName = 'QCDMuEnriched_Pt80to120'
  config.Data.inputDataset = '/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 13895437

  config.General.requestName = 'QCDMuEnriched_Pt80to120_ext1'
  config.Data.inputDataset = '/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 9809949

  config.General.requestName = 'QCDMuEnriched_Pt120to170'
  config.Data.inputDataset = '/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 7897731

  config.General.requestName = 'QCDMuEnriched_Pt170to300'
  config.Data.inputDataset = '/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_backup_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 7947159

  config.General.requestName = 'QCDMuEnriched_Pt170to300_ext1'
  config.Data.inputDataset = '/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 9403072
 
  config.General.requestName = 'QCDMuEnriched_Pt300to470'
  config.Data.inputDataset = '/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 7937590

  config.General.requestName = 'QCDMuEnriched_Pt300to470_ext1'
  config.Data.inputDataset = '/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM'
  crabCommand('submit', config = config)
  # 16462878

  config.General.requestName = 'QCDMuEnriched_Pt300to470_ext2'
  config.Data.inputDataset = '/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 24605508

  config.General.requestName = 'QCDMuEnriched_Pt470to600'
  config.Data.inputDataset = '/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 3972819
  
  config.General.requestName = 'QCDMuEnriched_Pt470to600_ext1'
  config.Data.inputDataset = '/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 5668793
  
  config.General.requestName = 'QCDMuEnriched_Pt470to600_ext2'
  config.Data.inputDataset = '/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 9847664
  
  config.General.requestName = 'QCDMuEnriched_Pt600to800'
  config.Data.inputDataset = '/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_backup_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 4010136

  config.General.requestName = 'QCDMuEnriched_Pt600to800_ext1'
  config.Data.inputDataset = '/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 5971175

  config.General.requestName = 'QCDMuEnriched_Pt800to1000'
  config.Data.inputDataset = '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 3962749

  config.General.requestName = 'QCDMuEnriched_Pt800to1000_ext1'
  config.Data.inputDataset = '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 6011849

  config.General.requestName = 'QCDMuEnriched_Pt800to1000_ext2'
  config.Data.inputDataset = '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 9966149

  config.General.requestName = 'QCDMuEnriched_Pt1000toInf'
  config.Data.inputDataset = '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 3990117

  config.General.requestName = 'QCDMuEnriched_Pt1000toInf_ext1'
  config.Data.inputDataset = '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 9638102

  
  # -- Top -- #
  config.General.requestName = 'TTLJ'
  config.Data.inputDataset = '/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 152669800

  config.General.requestName = 'TTLL'
  config.Data.inputDataset = '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
  crabCommand('submit', config = config)
  # 79140880

  config.General.requestName = 'ST_top_ext1'
  config.Data.inputDataset = '/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
  crabCommand('submit', config = config) # 3,256,650
  
  config.General.requestName = 'ST_top'
  config.Data.inputDataset = '/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
  crabCommand('submit', config = config) # 5,429,452
  
  config.General.requestName = 'ST_antitop_ext1'
  config.Data.inputDataset = '/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
  crabCommand('submit', config = config) # 3,256,407
  
  config.General.requestName = 'ST_antitop'
  config.Data.inputDataset = '/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
  crabCommand('submit', config = config) # 5,428,928

  
  # -- Jet DT -- #
  config.General.requestName = 'JetHT_B-v2'
  config.Data.inputDataset   = '/JetHT/Run2016B-17Jul2018_ver2-v2/MINIAOD'
  config.Data.lumiMask       = GoldenJSON
  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
  crabCommand('submit', config = config)
  
  config.General.requestName = 'JetHT_C'
  config.Data.inputDataset   = '/JetHT/Run2016C-17Jul2018-v1/MINIAOD'
  config.Data.lumiMask       = GoldenJSON
  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
  crabCommand('submit', config = config)
  
  config.General.requestName = 'JetHT_D'
  config.Data.inputDataset   = '/JetHT/Run2016D-17Jul2018-v1/MINIAOD'
  config.Data.lumiMask       = GoldenJSON
  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
  crabCommand('submit', config = config)
  
  config.General.requestName = 'JetHT_E'
  config.Data.inputDataset   = '/JetHT/Run2016E-17Jul2018-v1/MINIAOD'
  config.Data.lumiMask       = GoldenJSON
  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
  crabCommand('submit', config = config)
  
  config.General.requestName = 'JetHT_F'
  config.Data.inputDataset   = '/JetHT/Run2016F-17Jul2018-v1/MINIAOD'
  config.Data.lumiMask       = GoldenJSON
  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
  crabCommand('submit', config = config)
  
  config.General.requestName = 'JetHT_G'
  config.Data.inputDataset   = '/JetHT/Run2016G-17Jul2018-v1/MINIAOD'
  config.Data.lumiMask       = GoldenJSON
  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
  crabCommand('submit', config = config)
  
  config.General.requestName = 'JetHT_H'
  config.Data.inputDataset   = '/JetHT/Run2016H-17Jul2018-v1/MINIAOD'
  config.Data.lumiMask       = GoldenJSON
  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
  crabCommand('submit', config = config)
  

# # -- BJet DT -- #
#  config.General.requestName = 'BTagMu_B-v1'
#  config.Data.inputDataset   = '/BTagMu/Run2016B-17Jul2018_ver1-v1/MINIAOD'
#  config.Data.splitting = 'LumiBased'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  crabCommand('submit', config = config)
#  
#  config.General.requestName = 'BTagMu_B-v2'
#  config.Data.inputDataset   = '/BTagMu/Run2016B-17Jul2018_ver2-v1/MINIAOD'
#  config.Data.splitting = 'LumiBased'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  crabCommand('submit', config = config)
#  
#  config.General.requestName = 'BTagMu_C'
#  config.Data.inputDataset   = '/BTagMu/Run2016C-17Jul2018-v1/MINIAOD'
#  config.Data.splitting = 'LumiBased'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  crabCommand('submit', config = config)
#  
#  config.General.requestName = 'BTagMu_D'
#  config.Data.inputDataset   = '/BTagMu/Run2016D-17Jul2018-v1/MINIAOD'
#  config.Data.splitting = 'LumiBased'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  crabCommand('submit', config = config)
#  
#  config.General.requestName = 'BTagMu_E'
#  config.Data.inputDataset   = '/BTagMu/Run2016E-17Jul2018-v1/MINIAOD'
#  config.Data.splitting = 'LumiBased'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  crabCommand('submit', config = config)
#  
#  config.General.requestName = 'BTagMu_F'
#  config.Data.inputDataset   = '/BTagMu/Run2016F-17Jul2018-v1/MINIAOD'
#  config.Data.splitting = 'LumiBased'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  crabCommand('submit', config = config)
#  
#  config.General.requestName = 'BTagMu_G'
#  config.Data.inputDataset   = '/BTagMu/Run2016G-17Jul2018-v1/MINIAOD'
#  config.Data.splitting = 'LumiBased'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  crabCommand('submit', config = config)
#  
#  config.General.requestName = 'BTagMu_H'
#  config.Data.inputDataset   = '/BTagMu/Run2016H-17Jul2018-v1/MINIAOD'
#  config.Data.splitting = 'LumiBased'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  crabCommand('submit', config = config)




  # -- MET DT -- #
#  config.General.requestName = 'MET_B'
#  config.Data.inputDataset   = '/MET/Run2016B-17Jul2018_ver2-v1/MINIAOD'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  #crabCommand('submit', config = config)
#
#  config.General.requestName = 'MET_C'
#  config.Data.inputDataset   = '/MET/Run2016C-17Jul2018-v1/MINIAOD'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  crabCommand('submit', config = config)
#
#  config.General.requestName = 'MET_D'
#  config.Data.inputDataset   = '/MET/Run2016D-17Jul2018-v1/MINIAOD'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  #crabCommand('submit', config = config)
#
#  config.General.requestName = 'MET_E'
#  config.Data.inputDataset   = '/MET/Run2016E-17Jul2018-v1/MINIAOD'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  #crabCommand('submit', config = config)
#
#  config.General.requestName = 'MET_F'
#  config.Data.inputDataset   = '/MET/Run2016F-17Jul2018-v1/MINIAOD'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  #crabCommand('submit', config = config)
#
#  config.General.requestName = 'MET_G'
#  config.Data.inputDataset   = '/MET/Run2016G-17Jul2018-v1/MINIAOD'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  #crabCommand('submit', config = config)
#
#  config.General.requestName = 'MET_H'
#  config.Data.inputDataset   = '/MET/Run2016H-17Jul2018-v1/MINIAOD'
#  config.Data.lumiMask       = GoldenJSON
#  config.Data.runRange       = '%d-%d' % (StartRun, EndRun)
#  #crabCommand('submit', config = config)

