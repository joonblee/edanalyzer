location=/u/user/joonblee/SE_UserHome/

#ls ${location}/*/edanalyzer_trigeff/240905*

hadd ./storage/JetHT.root ${location}/JetHT/edanalyzer_trigeff/240905_*/*/hist_*.root

hadd ./storage/TTLL.root ${location}/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/edanalyzer_trigeff/240905_*/*/hist_*.root
hadd ./storage/TTLJ.root ${location}/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/edanalyzer_trigeff/240905_*/*/hist_*.root

#hadd ./storage/QCD_Pt-15to20.root ${location}/QCD_Pt-15to*/edanalyzer_trigeff/240905_*/*/*.root
#hadd ./storage/QCD_Pt-20to30.root ${location}/QCD_Pt-20to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-30to50.root ${location}/QCD_Pt-30to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-50to80.root ${location}/QCD_Pt-50to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-80to120.root ${location}/QCD_Pt-80to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-120to170.root ${location}/QCD_Pt-120to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-170to300.root ${location}/QCD_Pt-170to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-300to470.root ${location}/QCD_Pt-300to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-470to600.root ${location}/QCD_Pt-470to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-600to800.root ${location}/QCD_Pt-600to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-800to1000.root ${location}/QCD_Pt-800to*/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/QCD_Pt-1000toInf.root ${location}/QCD_Pt-1000toInf*/edanalyzer_trigeff/240905_*/*/*.root
#hadd ./storage/QCD_Pt-1000toInf.root ${location}/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/edanalyzer_trigeff/240905_3*/*/*.root

hadd ./storage/ST_top.root ${location}/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/edanalyzer_trigeff/240905_*/*/*.root
hadd ./storage/ST_antitop.root ${location}/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/edanalyzer_trigeff/240905_*/*/*.root


