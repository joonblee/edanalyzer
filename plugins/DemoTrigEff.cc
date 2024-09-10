// -*- C++ -*-
//
// Package:    Demo/DemoTrigEff
// Class:      DemoTrigEff
// 
/**\class DemoTrigEff DemoTrigEff.cc Demo/DemoTrigEff/plugins/DemoTrigEff.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jhovanny Andres Mejia Guisao
//         Created:  Tue, 05 Sep 2017 22:37:14 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include <TH2D.h>
#include <TMath.h>
#include <TLorentzVector.h>


/////////////////////
// -- For Muons -- //
/////////////////////
#include "DataFormats/MuonReco/interface/MuonFwd.h" // -- Forward declarations of muon variables
#include "DataFormats/MuonReco/interface/Muon.h" // -- A reconstructed Muon. (tracker alone, muon detector alone, combined muon plus tracker)
#include "DataFormats/PatCandidates/interface/Muon.h" // -- Analysis-level muon class, pat::Muon implements the analysis-level muon class within the 'pat' namespace. PAT(Physics Analysis Toolkit)
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h" // -- Variables shows how muon is cosmic-like

///////////////////////////////////
// -- For Electrons & Photons -- //
///////////////////////////////////
#include "DataFormats/PatCandidates/interface/Photon.h" // -- Analysis-level Photon class, pat::Photon implements the analysis-level photon class within the 'pat' namespace
#include "DataFormats/PatCandidates/interface/Electron.h" // -- Analysis-level electron class,  pat::Electron implements the analysis-level electron class within the 'pat' namespace
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h" // -- It seems as a function refer effective areas of ECAL
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"


///////////////////
// -- For MET -- //
///////////////////
#include "DataFormats/PatCandidates/interface/MET.h" // -- Analysis-level MET class, pat::MET implements an analysis-level missing energy class as a 4-vector within the 'pat' namespace
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h" // --  MET made from Particle Flow candidates, with various fraction functions ex)photonEtFraction()
#include "DataFormats/METReco/interface/PFMETCollection.h"

//////////////////////////
// -- Track & Vertex -- //
//////////////////////////
#include "DataFormats/TrackReco/interface/TrackFwd.h" // -- Forward definitions for tracker variables
#include "DataFormats/TrackReco/interface/Track.h" // -- reconstructed tracks that are stored in the AOD and RECO. also contains a reference to more detailed information(RECO)
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h" // -- A composite Candidate  with error matrix and other vertex fix information.
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h" // -- Foward declaration of VertexCompositeCandidate class's variables
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h" // -- Least-squares vertex fitter implemented in the Kalman Filter formalism
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h" // --
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" // -- Helper class to build TransientTrack from the persistent Track.
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

////////////////////
// -- For Jets -- //
////////////////////
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h" // -- Analysis-level calorimeter jet class, Jet implements the analysis-level calorimeter jet class within the 'pat' namespace.
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

////////////////////////////////
// -- For PackedCandidates -- //
////////////////////////////////
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

//////////////////////////
// -- Track & Vertex -- //
//////////////////////////
#include "DataFormats/TrackReco/interface/TrackFwd.h" // -- Forward definitions for tracker variables
#include "DataFormats/TrackReco/interface/Track.h" // -- reconstructed tracks that are stored in the AOD and RECO. also contains a reference to more detailed information(RECO)
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h" // -- A composite Candidate  with error matrix and other vertex fix information.
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h" // -- Foward declaration of VertexCompositeCandidate class's variables
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h" // -- Least-squares vertex fitter implemented in the Kalman Filter formalism
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h" // --
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" // -- Helper class to build TransientTrack from the persistent Track.
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

////////////////////////////
// -- For GenParticles -- //
////////////////////////////
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h" // -- LHE info like PDF
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h" // -- ??
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h" // -- Pythia weight names
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // -- 
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" // -- ???
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" // -- contains information related to the details of the pileup simulation for a given event, ex) Nvtx, z of PU, sum(pt) of vtxs


///////////////////
// -- Trigger -- //
///////////////////
#include "FWCore/Common/interface/TriggerNames.h" // -- Used to access the names and indices of the triggers corresponding to a particular TriggerResults object
#include "FWCore/Common/interface/TriggerResultsByName.h" // --  Class which provides methods to access trigger results
#include "DataFormats/Common/interface/TriggerResults.h" // -- The trigger path results are maintained here as a sequence of entries, one per trigger path
#include "DataFormats/HLTReco/interface/TriggerEvent.h" // -- The single EDProduct to be saved for each event (AOD case) describing the (HLT) trigger table
#include "DataFormats/HLTReco/interface/TriggerObject.h" // --  A single trigger object (e.g., an isolated muon, or MET) described by its 4-momentum and physics type
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h" // -- Analysis-level trigger object class (stand-alone). (within the 'pat' namespace.)
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" // -- This class provides access routines to get hold of the HLT Configuration
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"

////////////////////////////
// -- MuonBxCollection -- //
////////////////////////////
#include "DataFormats/L1Trigger/interface/Muon.h"

///////////////////////////
// -- For HLT objects -- //
///////////////////////////
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"

////////////////
// -- Else -- //
////////////////
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h" // -- Class to provide lumi weighting for analyzers to weight "flat-to-N" MC samples to data


// -- calculations -- //
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//
class DemoTrigEff : public edm::EDAnalyzer {
   public:
      explicit DemoTrigEff(const edm::ParameterSet&);
      ~DemoTrigEff();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


      // ----------member data ---------------------------
      //edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackTags_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pcToken;
      edm::EDGetTokenT<std::vector<pat::Jet>> JetToken;
      edm::EDGetTokenT<std::vector<pat::Muon>> MuonToken;
      edm::EDGetTokenT<edm::View<pat::Electron>> ElectronToken;
      edm::EDGetTokenT<edm::View<pat::Photon>> PhotonToken;
      edm::EDGetTokenT<std::vector<pat::MET>> MetToken;
      edm::EDGetTokenT<reco::VertexCollection> pvToken;
      //edm::EDGetTokenT< edm::View<reco::Track> > TrackToken;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
      //edm::EDGetTokenT<reco::RecoChargedCandidateCollection> L3ObjToken;
      edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> TriggerObjectToken;
      //edm::EDGetTokenT<GenEventInfoProduct> GenEventInfoToken;

      // ----- histograms ----- //
      TH1D *muonPt, *muonEta, *muonPhi;
      //TH2D *nTrksPerJetPt, *nGoodTrksPerJetPt, *nGoodPVTrksPerJetPt;
      TH1D *Nevt, *CutFlow;
      TH1D *JetPt_den, *MuonPt_den, *MuonEta_den, *dR_den, *mass_den;
      TH1D *JetPt_num, *MuonPt_num, *MuonEta_num, *dR_num, *mass_num;
      TH2D *eff_den, *eff_num;

      // ------ functions ----- //
      std::vector<pat::Jet> jet_selection(edm::Handle<std::vector<pat::Jet>>, std::vector<pat::Muon>);
      std::vector<pat::Muon> GetMuons(const reco::Vertex &vtx, edm::Handle<std::vector<pat::Muon>> MuonHandle, double pt, double eta, double iso, bool niso);
      std::vector<pat::Jet> GetJets(edm::Handle<std::vector<pat::Jet>> JetHandle, double pt, double eta);
      std::vector<pat::Jet> GetBJets(edm::Handle<std::vector<pat::Jet>> JetHandle, double pt, double eta);

      std::vector<TString> SetTrigger(const edm::Event& iEvent, edm::Handle<edm::TriggerResults> triggerBits);
      bool PassTrigger(TString trig, std::vector<TString> TriggerSets);
      bool PassTrigger(std::vector<TString> trigs, std::vector<TString> TriggerSets);
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoTrigEff::DemoTrigEff(const edm::ParameterSet& iConfig) :  

//trackTags_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("tracks")))
pcToken       ( consumes< pat::PackedCandidateCollection > (iConfig.getUntrackedParameter<edm::InputTag>("pfCands"))   ),
JetToken      ( consumes< std::vector<pat::Jet> >          (iConfig.getUntrackedParameter<edm::InputTag>("Jet")) ),
MuonToken     ( consumes< std::vector<pat::Muon> >         (iConfig.getUntrackedParameter<edm::InputTag>("Muon")) ),
ElectronToken ( consumes< edm::View<pat::Electron> >       (iConfig.getUntrackedParameter<edm::InputTag>("Electron")) ),
PhotonToken   ( consumes< edm::View<pat::Photon> >         (iConfig.getUntrackedParameter<edm::InputTag>("Photon")) ),
MetToken      ( consumes< std::vector<pat::MET> >          (iConfig.getUntrackedParameter<edm::InputTag>("MET")) ),
pvToken       ( consumes< reco::VertexCollection >         (iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertex")) ),
triggerResultsToken    ( consumes< edm::TriggerResults >                  (iConfig.getUntrackedParameter<edm::InputTag>("Triggers")) ),
//L3ObjToken            ( consumes< reco::RecoChargedCandidateCollection > (iConfig.getUntrackedParameter<edm::InputTag>("L3Obj")) ),
TriggerObjectToken     ( consumes< std::vector<pat::TriggerObjectStandAlone> > (iConfig.getUntrackedParameter<edm::InputTag>("L3obj")) )
//GenEventInfoToken      ( consumes< GenEventInfoProduct >                       (iConfig.getUntrackedParameter<edm::InputTag>("GenEventInfo")) )
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  Nevt                = fs->make<TH1D>("Nevt", "Nevt", 1,0,1);
  CutFlow             = fs->make<TH1D>("CutFlow", "CutFlow", 4,0,4);
  const Int_t NetaBin = 4; Double_t edge_eta[NetaBin+1] = {0., .8, 1.2, 2.1, 2.4};
  const Int_t NdRBin = 12; Double_t edge_dR[NdRBin+1] = {0., .01, .02, .03, .04, .05, .07, .1, .15, .2, .25, .3, 1.};
  const Int_t Npt = 10; Double_t edge_pt_[Npt+1] = {31., 32., 33., 34., 35., 40., 50.,60., 80., 100., 200.};
  const Int_t Nmass = 17; Double_t edge_mass_[Nmass+1] = {0.,1.,2.,2.9,3.3,4.,5.,6.,8.,10.,15.,20.,25.,30.,35.,40.,50.,100.};

  //MuonPt_den = fs->make<TH1D>("MuonPt_den", "MuonPt_den", Npt, edge_pt_);
  //MuonEta_den = fs->make<TH1D>("MuonEta_den", "MuonEta_den", NetaBin, edge_eta);
  //dR_den = fs->make<TH1D>("dR_den", "dR_den", NdRBin, edge_dR);
  //mass_den = fs->make<TH1D>("mass_den", "mass_den", Nmass, edge_mass_);
  //MuonPt_num = fs->make<TH1D>("MuonPt_num", "MuonPt_num", Npt, edge_pt_);
  //MuonEta_num = fs->make<TH1D>("MuonEta_num", "MuonEta_num", NetaBin, edge_eta);
  //dR_num = fs->make<TH1D>("dR_num", "dR_num", NdRBin, edge_dR);
  //mass_num = fs->make<TH1D>("mass_num", "mass_num", Nmass, edge_mass_);

  JetPt_den = fs->make<TH1D>("JetPt_den", "JetPt_den", 500,0.,1000.);
  MuonPt_den = fs->make<TH1D>("MuonPt_den", "MuonPt_den", 200,0.,200.);
  MuonEta_den = fs->make<TH1D>("MuonEta_den", "MuonEta_den", 48,-2.4,2.4);
  dR_den = fs->make<TH1D>("dR_den", "dR_den", 100,0.,1.);
  mass_den = fs->make<TH1D>("mass_den", "mass_den", 1000,0.,100.);
  JetPt_num = fs->make<TH1D>("JetPt_num", "JetPt_num", 500,0.,1000.);
  MuonPt_num = fs->make<TH1D>("MuonPt_num", "MuonPt_num", 200,0.,200.);
  MuonEta_num = fs->make<TH1D>("MuonEta_num", "MuonEta_num", 48,-2.4,2.4);
  dR_num = fs->make<TH1D>("dR_num", "dR_num", 100,0.,1.);
  mass_num = fs->make<TH1D>("mass_num", "mass_num", 1000,0.,100.);

  eff_den = fs->make<TH2D>("eff_den", "eff_den", Nmass, edge_mass_, Npt, edge_pt_);
  eff_num = fs->make<TH2D>("eff_num", "eff_num", Nmass, edge_mass_, Npt, edge_pt_);
}


DemoTrigEff::~DemoTrigEff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void DemoTrigEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;

  edm::Handle< pat::PackedCandidateCollection > pcHandle;
  iEvent.getByToken(pcToken, pcHandle);
  edm::Handle< std::vector<pat::Jet> > JetHandle;
  iEvent.getByToken(JetToken, JetHandle);
  edm::Handle<std::vector<pat::Muon>> MuonHandle;
  iEvent.getByToken(MuonToken, MuonHandle);
  edm::Handle<edm::View<pat::Electron>> ElectronHandle;
  iEvent.getByToken(ElectronToken, ElectronHandle);
  edm::Handle<edm::View<pat::Photon>> PhotonHandle;
  iEvent.getByToken(PhotonToken, PhotonHandle);
  edm::Handle<std::vector<pat::MET>> MetHandle;
  iEvent.getByToken(MetToken, MetHandle);
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken, pvHandle);

  Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerResultsToken, triggerBits);
  const edm::TriggerNames trigNames = iEvent.triggerNames(*triggerBits);
  //edm::Handle<trigger::TriggerEvent> h_triggerEvent;
  //iEvent.getByToken(t_triggerEvent_,   h_triggerEvent);
  edm::Handle< std::vector<pat::TriggerObjectStandAlone> > triggerObject;
  iEvent.getByToken(TriggerObjectToken, triggerObject);
  //edm::Handle<GenEventInfoProduct> genEvtInfo;
  //iEvent.getByToken(GenEventInfoToken, genEvtInfo);

  // -- options -- //
  bool validation_mode = 1;

  // -- Data MC -- //
  bool b_data = iEvent.id().run() == 1 ? 0 : 1;

  // -- Set weight -- //
  double wgt = 0.;
  Nevt->Fill(0.,1.);

  // -- Print Information -- //
  if(validation_mode) {
    cout<<endl<<endl;
    cout<<"[Event id] run # = "<<iEvent.id().run()<<", evt # = "<<iEvent.id().event()<<", lumi = "<<iEvent.id().luminosityBlock()<<endl;
    //cout<<"weight size = "<<genEvtInfo->weights().size()<<", weight = "<<genEvtInfo->weight()<<endl;
  }

  /*
  edm::Handle<reco::RecoChargedCandidateCollection> L3muon;
  iEvent.getByToken(L3ObjToken, L3muon);

  if( !L3muon.isValid() ) return;
  for(unsigned int iL3=0; iL3<L3muon->size(); iL3++) {
    reco::RecoChargedCandidateRef L3obj(L3muon, iL3);
    cout<<"[L3 muon] ("<<L3obj->pt()<<", "<<L3obj->eta()<<", "<<L3obj->phi()<<") charge = "<<L3obj->charge()<<endl;
    reco::TrackRef trackRef = L3obj->track();
    cout<<"[L3 muon track] pt = "<<trackRef->pt()<<endl;
  }
  */

  // -- Trigger Selection -- //
  if(validation_mode) cout<<"Processing orthogonal trigger selection ... ";

  TString target = "HLT_Mu30_TkMu11_v";
  // vector<TString> orthogonals = {"HLT_PFJet450_v"}; // unprescaled
  vector<TString> orthogonals = {"HLT_PFJet900_v", "HLT_PFJet500_v", "HLT_PFJet450_v", "HLT_PFJet400_v", "HLT_PFJet320_v", "HLT_PFJet260_v", "HLT_PFJet200_v", "HLT_PFJet140_v", "HLT_PFJet80_v", "HLT_PFJet60_v", "HLT_PFJet40_v"}; // all single pf jet triggers
  // vector<TString> orthogonals = {"HLT_BTagMu_Jet300_Mu5_v"};
  TString orthogonal="";

  // check trigger prescales, i.e. (act/eff) lumi, regarding to each year
  // https://cmshltinfo.app.cern.ch/summary
  double lumi = 36.310048842; // fb-1
  map<TString, double> prescales; prescales["HLT_PFJet40_v"]=lumi/0.000266973; prescales["HLT_PFJet60_v"]=lumi/0.000726088; prescales["HLT_PFJet80_v"]=lumi/0.002758456; prescales["HLT_PFJet140_v"]=lumi/0.024188608; prescales["HLT_PFJet200_v"]=lumi/0.103797776; prescales["HLT_PFJet260_v"]=lumi/0.593633969; prescales["HLT_PFJet320_v"]=lumi/1.772108294; prescales["HLT_PFJet400_v"]=lumi/5.193142997; prescales["HLT_PFJet450_v"]=lumi/36.310048842; prescales["HLT_PFJet500_v"]=1.; // prescales["HLT_PFJet900_v"]=1.;

  vector<TString> TriggerSets = SetTrigger(iEvent, triggerBits);
  for(auto trig:orthogonals) {
    if(PassTrigger({trig}, TriggerSets)) {
      orthogonal=trig;
      if(!b_data) wgt=wgt*(1.-1./prescales[trig])+1./prescales[trig];
    }
  }

  if(orthogonal=="") return;
  if(b_data) wgt=1.;
  //if( !PassTrigger(orthogonals, TriggerSets) ) return; // Pass orthogonal trigger
  if(validation_mode) cout<<"Pass "<<orthogonal<<endl;
  //CutFlow->Fill(0.,wgt);

  //cout<<endl<<endl<<endl<<endl;
  //cout<<"============================================================="<<endl;
  //cout<<" ### -- Pass HLT_MET200_v -- ###"<<endl;

  // -- Call objects -- //
  vector<pat::Muon> muons = GetMuons(pvHandle->front(), MuonHandle, 13., 2.4, .3, 1);
  vector<pat::Jet> jets = GetJets(JetHandle, 30., 3.0);
  vector<pat::Jet> bjets = GetBJets(JetHandle, 30., 2.4);

  // -- HLT object matching -- //
  if(validation_mode) cout<<"Check #(triggered jet) > 0 ... ";
  if( !triggerObject.isValid() ) return;

  orthogonal+="*";
  std::string orthogonal_str(orthogonal.Data());

  vector<pat::Jet> trigjets;
  for(auto jet:jets) {
    //if( jet.pt() > 550. && fabs(jet.eta()) < 3. ) {
    //  trigjets.push_back(jet);
    //}
    for(pat::TriggerObjectStandAlone L3obj : *triggerObject) {
      L3obj.unpackPathNames(trigNames);
      L3obj.unpackFilterLabels(iEvent, *triggerBits);
      if(!L3obj.path(orthogonal_str,true,true)) continue;
      if( deltaR(jet, L3obj) < .4 && fabs(jet.pt()-L3obj.pt())/jet.pt() < .5 ) {
        trigjets.push_back(jet);
        break;
      }
    }
  }
  if(trigjets.size()==0) 
    return;
  if(validation_mode) cout<<"Pass"<<endl;

  // -- Event Selection -- //
  if(validation_mode) cout<<"Processing event selection ... ";
  pat::Jet selected_jet, hlt_jet;
  vector<pat::Muon> selected_muons;
  for(auto bjet:bjets) {
    for(unsigned i=0; i<muons.size(); i++) {
      pat::Muon muon_1 = muons[i];
      if( !(muon_1.pt() > 32) ) continue;
      if( !(deltaR(bjet, muon_1) < 0.4) ) continue;
      for(unsigned j=i+1; j<muons.size(); j++) {
        pat::Muon muon_2 = muons[j];
        if( !(deltaR(bjet, muon_2) < 0.4) ) continue;
        if( !(muon_1.charge()+muon_2.charge()==0) ) continue;
        for(auto trigjet:trigjets) {
          if( !(deltaR(trigjet, bjet) > 0.5) ) continue;
          selected_jet = bjet;
          selected_muons = {muon_1, muon_2};
					hlt_jet = trigjet;
          break;
        }
        if(selected_muons.size()==2) break;
      }
      if(selected_muons.size()==2) break;
    }
    if(selected_muons.size()==2) break;
  }
  if(selected_muons.size()!=2) return;
  if(validation_mode) cout<<"Pass"<<endl;

  CutFlow->Fill(1.,wgt); // Pass event selection
  JetPt_den->Fill(selected_jet.pt(),wgt);
  for(auto muon:selected_muons) {
    MuonEta_den->Fill(muon.eta(),wgt);
    MuonPt_den->Fill(muon.pt(),wgt);
  }
  dR_den->Fill(deltaR(selected_muons[0], selected_muons[1]),wgt);

  TLorentzVector muon_1_4vec, muon_2_4vec;
  muon_1_4vec.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), 0.10565837/*selected_muons[0].mass()*/);
  muon_2_4vec.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), 0.10565837/*selected_muons[1].mass()*/);
  double mass = (muon_1_4vec + muon_2_4vec).M();

  mass_den->Fill(mass,wgt);
  eff_den->Fill(mass,selected_muons[0].pt(),wgt);

  // -- Pass HLT -- //
  if(validation_mode) cout<<"Processing target trigger selection ... ";
  if( !PassTrigger({target}, TriggerSets) ) return; // Pass target trigger
  CutFlow->Fill(2.,wgt);
  JetPt_num->Fill(selected_jet.pt(),wgt);
  for(auto muon:selected_muons) {
    MuonEta_num->Fill(muon.eta(),wgt);
    MuonPt_num->Fill(muon.pt(),wgt);
  }
  dR_num->Fill(deltaR(selected_muons[0], selected_muons[1]),wgt);
  mass_num->Fill(mass,wgt);

  eff_num->Fill(mass,selected_muons[0].pt(),wgt);
  if(validation_mode) {
    cout<<"Pass"<<endl;
    cout<<" -- END --"<<endl;
  }
}

std::vector<pat::Jet> DemoTrigEff::jet_selection(edm::Handle<std::vector<pat::Jet>> JetHandle, std::vector<pat::Muon> muons)
{
  std::vector<pat::Jet> jets;
  for(unsigned iJet = 0; iJet < JetHandle->size(); iJet++) {
    const pat::Jet jet = JetHandle->at(iJet);
    bool muonVeto = false;
    for( auto muon : muons ) {
      if( deltaR( jet , muon ) < 0.5 ) { 
        muonVeto = true; break;
      }
    }
    if( muonVeto ) continue;
    if( !(30 < jet.pt()) ) continue;
    if( !(fabs(jet.eta()) < 2.4) ) continue;
    jets.push_back( jet );
  }
  return jets;
}

std::vector<TString> DemoTrigEff::SetTrigger(const edm::Event& iEvent, edm::Handle<edm::TriggerResults> triggerBits){
  vector<TString> out;
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //std::cout << "\n == TRIGGER PATHS= " << std::endl;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    //std::cout << "Trigger " << names.triggerName(i) <<
    //        ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
    //        ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
    //        << std::endl;
    if( triggerBits->accept(i) ) out.push_back(names.triggerName(i));
  }
  return out;
}

bool DemoTrigEff::PassTrigger(TString trig, std::vector<TString> TriggerSets){
  std::vector<TString> tmp_vec;
  tmp_vec.push_back(trig);
  return PassTrigger(tmp_vec, TriggerSets);
}

bool DemoTrigEff::PassTrigger(std::vector<TString> trigs, std::vector<TString> TriggerSets){

  for(unsigned int i=0; i<trigs.size(); i++){
    TString this_check_trig = trigs.at(i);

    //cout << this_check_trig << endl;
    for(unsigned int j=0; j<TriggerSets.size(); j++){

      //cout << TriggerSets.at(j) << endl;
      if( TriggerSets.at(j).Index(this_check_trig)!=-1 ){
        return true;
      }
    }

  }
  return false;
}


std::vector<pat::Muon> DemoTrigEff::GetMuons(const reco::Vertex &vtx, edm::Handle<std::vector<pat::Muon>> MuonHandle, double pt, double eta, double iso, bool niso = 1){
  vector<pat::Muon> out;
  for(unsigned i=0; i<MuonHandle->size(); i++) {
    const pat::Muon muon = MuonHandle->at(i);
    if( !(muon.pt() > pt) ) return out;
    if( !(fabs(muon.eta()) < eta) ) continue;
    double absiso = muon.pfIsolationR04().sumChargedHadronPt+std::max( 0., muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5*muon.pfIsolationR04().sumPUPt );
    if( niso ) {
      if( !(absiso/muon.pt() > iso) ) continue;
    }
    else {
      if( !(absiso/muon.pt() < iso) ) continue;
    }
    if( !(muon.isMediumMuon()) ) continue;
    out.push_back( muon );
  }
  return out;
}

std::vector<pat::Jet> DemoTrigEff::GetJets(edm::Handle<std::vector<pat::Jet>> JetHandle, double pt, double eta){
  vector<pat::Jet> out;
  for(unsigned i=0; i<JetHandle->size(); i++) {
    const pat::Jet jet = JetHandle->at(i);
    if( !(jet.pt() > pt) ) return out;
    double eta_ = jet.eta();
    if( !(fabs(eta_) < eta) ) continue;
    double NHF  = jet.neutralHadronEnergyFraction();
    double NEMF = jet.neutralEmEnergyFraction();
    double CHF  = jet.chargedHadronEnergyFraction();
    double MUF  = jet.muonEnergyFraction();
    double CEMF = jet.chargedEmEnergyFraction();
    double NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    double NumNeutralParticles =jet.neutralMultiplicity();
    double CHM      = jet.chargedMultiplicity();
    bool tightJetID=false;
    //bool tightLepVetoJetID=false;
    if(fabs(eta_)>3.0){
      tightJetID = (NEMF<0.90 && NumNeutralParticles>10 );
      //tightLepVetoJetID = (NEMF<0.90 && NumNeutralParticles>10 );
    }
    else if(fabs(eta_)>2.7){
      tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
      //tightLepVetoJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
    }
    else if(fabs(eta_)>2.4){
      tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1);
    }
    else {
      tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99);
      //tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((fabs(eta_)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || fabs(eta_)>2.4) && fabs(eta_)<=2.7;
    }
    if( !tightJetID ) continue;
    out.push_back( jet );
  }
  return out;
}

std::vector<pat::Jet> DemoTrigEff::GetBJets(edm::Handle<std::vector<pat::Jet>> JetHandle, double pt, double eta){
  vector<pat::Jet> out;
  for(unsigned i=0; i<JetHandle->size(); i++) {
    const pat::Jet jet = JetHandle->at(i);
    if( !(jet.pt() > pt) ) return out;
    double eta_ = jet.eta();
    if( !(fabs(eta_) < eta) ) continue;
    double NHF  = jet.neutralHadronEnergyFraction();
    double NEMF = jet.neutralEmEnergyFraction();
    double CHF  = jet.chargedHadronEnergyFraction();
    double MUF  = jet.muonEnergyFraction();
    double CEMF = jet.chargedEmEnergyFraction();
    double NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    double NumNeutralParticles =jet.neutralMultiplicity();
    double CHM      = jet.chargedMultiplicity();
    bool tightJetID=false;
    //bool tightLepVetoJetID=false;
    if(fabs(eta_)>3.0){
      tightJetID = (NEMF<0.90 && NumNeutralParticles>10 );
      //tightLepVetoJetID = (NEMF<0.90 && NumNeutralParticles>10 );
    }
    else if(fabs(eta_)>2.7){
      tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
      //tightLepVetoJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
    }
    else if(fabs(eta_)>2.4){
      tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1);
    }
    else {
      tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99);
      //tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((fabs(eta_)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || fabs(eta_)>2.4) && fabs(eta_)<=2.7;
    }
    if( !tightJetID ) continue;
    if( !(jet.bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll") > 0.6321) ) continue; // Loose: 0.2217, Medium: 0.6321, Tight: 0.8953

    out.push_back( jet );
  }
  return out;
}


// ------------ method called once each job just before starting event loop  ------------
void DemoTrigEff::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void DemoTrigEff::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
/*
void 
DemoTrigEff::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DemoTrigEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DemoTrigEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DemoTrigEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoTrigEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoTrigEff);
