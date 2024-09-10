// -*- C++ -*-
//
// Package:    Demo/DemoTrigAnls
// Class:      DemoTrigAnls
// 
/**\class DemoTrigAnls DemoTrigAnls.cc Demo/DemoTrigAnls/plugins/DemoTrigAnls.cc

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



// -- calculations -- //
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//
class DemoTrigAnls : public edm::EDAnalyzer {
   public:
      explicit DemoTrigAnls(const edm::ParameterSet&);
      ~DemoTrigAnls();

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
      edm::EDGetTokenT<l1t::MuonBxCollection> L1MuonToken;

     

      TH1D *jetPt, *jetEta, *jetPhi;
      TH1D *muonPt, *muonEta, *muonPhi;
      //TH2D *nTrksPerJetPt, *nGoodTrksPerJetPt, *nGoodPVTrksPerJetPt;
      TH1D *CutFlow, *Den_eta, *PassL1_eta, *PassL3_eta, *PassHLT_eta, *Den_dR, *PassL1_dR, *PassL3_dR, *PassHLT_dR, *Den_pt, *PassL1_pt, *PassL3_pt, *PassHLT_pt;

      std::vector< std::vector<int> > vec_Ntrks, vec_Ngoodtrks, vec_NgoodPVtrks;
      int MaxJetPt;


      // ------ functions ----- //
      std::vector<pat::Muon> W_selection(edm::Handle<std::vector<pat::Muon>>, edm::Handle<std::vector<pat::MET>>, edm::Handle<reco::VertexCollection>);
      std::vector<pat::Muon> Z_selection(edm::Handle<std::vector<pat::Muon>>, edm::Handle<reco::VertexCollection>);
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
DemoTrigAnls::DemoTrigAnls(const edm::ParameterSet& iConfig) :  

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
TriggerObjectToken     ( consumes< std::vector<pat::TriggerObjectStandAlone> > (iConfig.getUntrackedParameter<edm::InputTag>("L3obj")) ),
L1MuonToken            ( consumes< l1t::MuonBxCollection  >               (iConfig.getUntrackedParameter<edm::InputTag>("L1Muon")) )

{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  jetPt               = fs->make<TH1D>("jet pt" , "jet pt" , 200 , 0 , 200 );
  jetEta              = fs->make<TH1D>("jet eta" , "jet eta" , 100 , -5. , 5. );
  jetPhi              = fs->make<TH1D>("jet phi" , "jet phi" , 64 , -3.2 , 3.2 );
  muonPt              = fs->make<TH1D>("muon pt" , "muon pt" , 200 , 0 , 200 );
  muonEta             = fs->make<TH1D>("muon eta" , "muon eta" , 100 , -5. , 5. );
  muonPhi             = fs->make<TH1D>("muon phi" , "muon phi" , 64 , -3.2 , 3.2 );
  CutFlow             = fs->make<TH1D>("CutFlow", "CutFlow", 4,0,4);
  const Int_t NetaBin = 4; Double_t edge_eta[NetaBin+1] = {0., .8, 1.2, 2.1, 2.4};
  const Int_t NdRBin = 12; Double_t edge_dR[NdRBin+1] = {0., .01, .02, .03, .04, .05, .07, .1, .15, .2, .25, .3, 1.};
  const Int_t NptBin = 5; Double_t edge_pt[NptBin+1] = {10., 20., 30., 50., 100., 200.};
  Den_eta = fs->make<TH1D>("Den_eta", "Den_eta", NetaBin, edge_eta);
  PassL1_eta = fs->make<TH1D>("PassL1_eta", "PassL1_eta", NetaBin, edge_eta);
  PassL3_eta = fs->make<TH1D>("PassL3_eta", "PassL3_eta", NetaBin, edge_eta);
  PassHLT_eta = fs->make<TH1D>("PassHLT_eta", "PassHLT_eta", NetaBin, edge_eta);
  Den_dR = fs->make<TH1D>("Den_dR", "Den_dR", NdRBin, edge_dR);
  PassL1_dR = fs->make<TH1D>("PassL1_dR", "PassL1_dR", NdRBin, edge_dR);
  PassL3_dR = fs->make<TH1D>("PassL3_dR", "PassL3_dR", NdRBin, edge_dR);
  PassHLT_dR = fs->make<TH1D>("PassHLT_dR", "PassHLT_dR", NdRBin, edge_dR);
  Den_pt = fs->make<TH1D>("Den_pt", "Den_pt", NptBin, edge_pt);
  PassL1_pt = fs->make<TH1D>("PassL1_pt", "PassL1_pt", NptBin, edge_pt);
  PassL3_pt = fs->make<TH1D>("PassL3_pt", "PassL3_pt", NptBin, edge_pt);
  PassHLT_pt = fs->make<TH1D>("PassHLT_pt", "PassHLT_pt", NptBin, edge_pt);

}


DemoTrigAnls::~DemoTrigAnls()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void DemoTrigAnls::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  edm::Handle<l1t::MuonBxCollection> L1muon;
  iEvent.getByToken(L1MuonToken, L1muon);
  edm::Handle< std::vector<pat::TriggerObjectStandAlone> > triggerObject;
  iEvent.getByToken(TriggerObjectToken, triggerObject);


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
  vector<TString> TriggerSets = SetTrigger(iEvent, triggerBits);
  TString orthogonal = "HLT_MET200_v";
  //TString orthogonal = "HLT_PFJet450_v";
  if( !PassTrigger({orthogonal}, TriggerSets) ) return;

  //cout<<endl<<endl<<endl<<endl;
  //cout<<"============================================================="<<endl;
  //cout<<" ### -- Pass HLT_MET200_v -- ###"<<endl;
  //cout<<"[Event id] run # = "<<iEvent.id().run()<<", evt # = "<<iEvent.id().event()<<endl;

  // -- Event Selection -- //
  vector<pat::Muon> muons = GetMuons(pvHandle->front(), MuonHandle, 13., 2.4, .3, 1);
  //vector<pat::Jet> jets = GetJets(JetHandle, 500., 2.4);
  vector<pat::Jet> bjets = GetBJets(JetHandle, 30., 2.4);

  pat::Jet selected_jet;
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
        //for(auto trig_jet:jets) {
        //  if( !(deltaR(trig_jet, bjet) > 0.4) ) continue;
        //}
        selected_jet = bjet;
        selected_muons = {muon_1, muon_2};
        break;
      }
      if(selected_muons.size()==2) break;
    }
    if(selected_muons.size()==2) break;
  }
  if(selected_muons.size()!=2) return;

  CutFlow->Fill(0);
  for(auto muon:selected_muons) {
    Den_eta->Fill(muon.eta());
    Den_pt->Fill(muon.pt());
  }
  Den_dR->Fill(deltaR(selected_muons[0], selected_muons[1]));


  // -- L1 object matching -- //
  if( !L1muon.isValid() ) return;
  vector<bool> PassL1={0,0};
  for(unsigned i=0; i<selected_muons.size(); i++) {
    pat::Muon muon = selected_muons[i];
    for(int ibx=L1muon->getFirstBX(); ibx<=L1muon->getLastBX(); ++ibx) {
      if(ibx != 0) continue; // -- only take when ibx == 0 -- //
      for(auto it=L1muon->begin(ibx); it!=L1muon->end(ibx); it++) {
        //cout<<"[L1 muon] bx: "<<ibx<<", et: "<<it->et()<<", ("<<it->pt()<<", "<<it->eta()<<", "<<it->phi()<<") charge = "<<it->charge()<<", hwQual = "<<it->hwQual()<<endl;
        if( deltaR(muon, *it) < .3 && muon.charge() == it->charge() ) {
          PassL1[i]=1;
          break;
        }
      }
      if( PassL1[i] ) break;
    }
  }
  if( !(PassL1[0] && PassL1[1]) ) return;

  CutFlow->Fill(1);
  for(auto muon:selected_muons) {
    PassL1_eta->Fill(muon.eta());
    PassL1_pt->Fill(muon.pt());
  }
  PassL1_dR->Fill(deltaR(selected_muons[0], selected_muons[1]));


  // -- L3 object matching -- //
  TString target = "HLT_Mu30_TkMu11_v";
  std::string target_str(target.Data());

  if( !triggerObject.isValid() ) return;
  vector<bool> PassL3={0,0};
  for(unsigned i=0; i<selected_muons.size(); i++) {
    pat::Muon muon = selected_muons[i];
    for(pat::TriggerObjectStandAlone L3obj : *triggerObject) {
      L3obj.unpackPathNames(trigNames);
      L3obj.unpackFilterLabels(iEvent, *triggerBits);
      if(!L3obj.path(target_str+"*",true,true)) continue;
      if( deltaR(muon, L3obj) < .3 && muon.charge() == L3obj.charge() ) {
        PassL3[i]=1;
        break;
      }
    }
  }
  if( !(PassL3[0] && PassL3[1]) ) return;

  CutFlow->Fill(2);
  for(auto muon:selected_muons) {
    PassL3_eta->Fill(muon.eta());
    PassL3_pt->Fill(muon.pt());
  }
  PassL3_dR->Fill(deltaR(selected_muons[0], selected_muons[1]));


  // -- Pass HLT -- //
  if( !PassTrigger({target}, TriggerSets) ) return;
  CutFlow->Fill(3);
  for(auto muon:selected_muons) {
    PassHLT_eta->Fill(muon.eta());
    PassHLT_pt->Fill(muon.pt());
  }
  PassHLT_dR->Fill(deltaR(selected_muons[0], selected_muons[1]));



/*

  vector<pat::Muon> muons = W_selection(MuonHandle, MetHandle, pvHandle);
  //vector<pat::Muon> muons = Z_selection(MuonHandle, pvHandle);
  if( muons.size() == 0 ) return;


  for( auto muon : muons ) { muonPt->Fill(muon.pt()); muonEta->Fill(muon.eta()); muonPhi->Fill(muon.phi()); }
  const pat::MET  met  = MetHandle->front();
  metUncorPt->Fill(met.pt()); metUncorPhi->Fill(met.phi());
  metPt->Fill(met.corPt(pat::MET::Type1XY)); metPhi->Fill(met.corPhi(pat::MET::Type1XY));


  vector<pat::Jet> jets = jet_selection(JetHandle, muons);
  
  for( auto jet : jets ) { jetPt->Fill(jet.pt()); jetEta->Fill(jet.eta()); jetPhi->Fill(jet.phi()); }
  nJets->Fill(jets.size());
  if( jets.size() == 0 ) return;
  
  cout<<endl;
  cout<<"============================================================="<<endl;
  cout<<"[Event id] run # = "<<iEvent.id().run()<<", evt # = "<<iEvent.id().event()<<endl;
  cout<<"[Jet Info] #(AllJets) = "<<JetHandle->size()<<", #(SelectedJets) = "<<jets.size()<<endl;


  MaxJetPt = 500;
  for(int i=0; i<MaxJetPt; i++) { 
    vec_Ntrks.push_back({0,0}); vec_Ngoodtrks.push_back({0,0}); vec_NgoodPVtrks.push_back({0,0});
  }


  for(unsigned iJet=0; iJet < jets.size(); iJet++) {
    const pat::Jet jet = jets[iJet];

    vector<pat::PackedCandidate> jetContents;
    bool isMuon=false;
    for(unsigned iPC=0; iPC < pcHandle->size(); iPC++) {
      const pat::PackedCandidate pc = pcHandle->at(iPC);
      if( deltaR( jet , pc ) < 0.4 ) {
        //if( (abs(pc.pdgId())==11 || abs(pc.pdgId())==13 || abs(pc.pdgId())==15) && pc.pt() > 0.5*jet.pt() ) { isMuon=true; break; }
        jetContents.push_back(pc);
      }
    }

    if( jetContents.size() == 0 ) continue;
    if( isMuon ) continue;

    cout<<" - "<<endl;
    cout<<"[Jet Info] jet"<<iJet<<": ("<<jet.pt()<<", "<<jet.eta()<<", "<<jet.phi()<<")"<<endl;

    double pt=0, trkpt=0, goodtrkpt=0, PVtrkpt=0;
    int nTrks_=0, nGoodTrks_=0, nGoodPVTrks_=0;
    for(unsigned iCont=0; iCont<jetContents.size(); iCont++) {
      pat::PackedCandidate jc = jetContents[iCont];

      pt += jc.pt();

      //bool fromPV = (jc.fromPV()>1 || fabs(jc.dz()) < 0.1);
      //cout<<"[Packed Candidate Info] jc"<<iCont<<": ("<<jc.pt()<<", "<<jc.eta()<<", "<<jc.phi()<<") id: "<<jc.pdgId()<<", track? "<<jc.hasTrackDetails()<<", PV"<<fromPV<<", purity"<<jc.trackHighPurity()<<", nPixHit "<<jc.numberOfPixelHits()<<", nHits "<<jc.numberOfHits()<<", jcLayer "<<jc.trackerLayersWithMeasurement()<<endl;
      //cout<<"[Packed Candidate Info] jc"<<iCont<<": ("<<jc.pt()<<", "<<jc.eta()<<", "<<jc.phi()<<") id: "<<jc.pdgId()<<", track? "<<jc.hasTrackDetails()<<", PV"<<jc.fromPV()<<", PV(dz)"<<fromPV<<", purity"<<jc.trackHighPurity()<<endl;

      if( jc.hasTrackDetails() ) {
        trkpt += jc.pt(); nTrks_++; trkPt->Fill(jc.pt());
        if( jc.trackHighPurity() ) {
          goodtrkpt += jc.pt(); nGoodTrks_++; goodTrkPt->Fill(jc.pt());
          if( jc.fromPV() > 0 || jc.fromPV() == -1 ) {
            PVtrkpt += jc.pt(); nGoodPVTrks_++; goodPVTrkPt->Fill(jc.pt());
            if( 30 < jet.pt() && jet.pt() < 40 ) goodPVTrkPt_pt30->Fill(jc.pt());
            else if( 100 < jet.pt() && jet.pt() < 120 ) goodPVTrkPt_pt100->Fill(jc.pt());
          }
        }
      }

    }

    nTrks->Fill(nTrks_);
    nGoodTrks->Fill(nGoodTrks_);
    nGoodPVTrks->Fill(nGoodPVTrks_);
    if( 30 < jet.pt() && jet.pt() < 40 ) nGoodPVTrks_pt30->Fill(nGoodPVTrks_);
    else if( 100 < jet.pt() && jet.pt() < 120 ) nGoodPVTrks_pt100->Fill(nGoodPVTrks_);

    nTrksPerJetPt->Fill(jet.pt(), nTrks_);
    nGoodTrksPerJetPt->Fill(jet.pt(), nGoodTrks_);
    nGoodPVTrksPerJetPt->Fill(jet.pt(), nGoodPVTrks_);

    trkPtSumRatio->Fill(trkpt/jet.pt());
    goodTrkPtSumRatio->Fill(goodtrkpt/jet.pt());
    goodPVTrkPtSumRatio->Fill(PVtrkpt/jet.pt());

    pcPtRatio->Fill(pt/jet.pt());

    if( jet.pt() < MaxJetPt ) {
      vec_Ntrks[(int)jet.pt()][0]++; vec_Ntrks[(int)jet.pt()][1]+=nTrks_;
      vec_Ngoodtrks[(int)jet.pt()][0]++; vec_Ngoodtrks[(int)jet.pt()][1]+=nGoodTrks_;
      vec_NgoodPVtrks[(int)jet.pt()][0]++; vec_NgoodPVtrks[(int)jet.pt()][1]+=nGoodPVTrks_;
    }

//    cout<<"[Packed Candidate Info] total pt       = "<<pt       <<endl;
//    cout<<"[Packed Candidate Info] track pt       = "<<trkpt    <<", ratio = "<<trkpt/pt<<endl;
//    cout<<"[Packed Candidate Info] good trk pt    = "<<goodtrkpt<<", ratio = "<<goodtrkpt/pt<<endl;
//    cout<<"[Packed Candidate Info] good PV trk pt = "<<PVtrkpt  <<", ratio = "<<PVtrkpt/pt<<endl;
  }


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/   
}

std::vector<pat::Muon> DemoTrigAnls::W_selection(edm::Handle<std::vector<pat::Muon>> MuonHandle, edm::Handle<std::vector<pat::MET>> MetHandle, edm::Handle<reco::VertexCollection> pvHandle)
{
  std::vector<pat::Muon> muons;

  const pat::MET  met  = MetHandle->front();
  if( !(met.pt() > 20) ) return muons;
  TLorentzVector met_4vec;
  //met_4vec.SetPtEtaPhiE(met.pt(), 0., met.phi(), met.pt());
  met_4vec.SetPtEtaPhiE(met.corPt(pat::MET::Type1XY), 0., met.corPhi(pat::MET::Type1XY), met.corPt(pat::MET::Type1XY));

  const reco::Vertex &vtx = pvHandle->front();
  for(unsigned iMuon = 0; iMuon < MuonHandle->size(); iMuon++) {
    const pat::Muon muon = MuonHandle->at(iMuon);

    if( !(muon.pt() > 10 && fabs(muon.eta()) < 2.4) ) continue;
    if( !(muon.isGlobalMuon()) ) continue;

    muons.push_back( muon );
  }

  if( muons.size() != 1 ) { muons.clear(); return muons; }

  pat::Muon muon = muons[0];
  if( muon.pt() > 26 && muon.isolationR03().sumPt/muon.pt() < 0.05 && muon.isTightMuon(vtx) ) {
    TLorentzVector muon_4vec;
    muon_4vec.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
    TLorentzVector W_4vec = muon_4vec + met_4vec;
    double mt = sqrt(2*muon.pt()*met.pt()*(1-cos(muon.phi()-met.phi()))); //W_4vec.Mt();
    //MASS->Fill(W_4vec.M()); M_T->Fill(mt);
    if( 60 < mt && mt < 100 ) return muons;
  }

  muons.clear();
  return muons;
}


std::vector<pat::Muon> DemoTrigAnls::Z_selection(edm::Handle<std::vector<pat::Muon>> MuonHandle, edm::Handle<reco::VertexCollection> pvHandle)
{
  std::vector<pat::Muon> muons;

  const reco::Vertex &vtx = pvHandle->front();
  for(unsigned iMuon = 0; iMuon < MuonHandle->size(); iMuon++) {
    const pat::Muon muon_1 = MuonHandle->at(iMuon);

    if( !(muon_1.pt() > 26) ) return muons;
    if( !(fabs(muon_1.eta()) < 2.4) ) continue;
    if( !(muon_1.isolationR03().sumPt/muon_1.pt() < 0.1) ) continue;
    if( !(muon_1.isTightMuon(vtx)) ) continue;

    for(unsigned jMuon = iMuon+1; jMuon < MuonHandle->size(); jMuon++) {
      const pat::Muon muon_2 = MuonHandle->at(jMuon);

      if( !(muon_2.pt() > 10) ) break;
      if( !(fabs(muon_2.eta()) < 2.4) ) continue;
      if( !(muon_2.isolationR03().sumPt/muon_2.pt() < 0.1) ) continue;
      if( !(muon_2.isTightMuon(vtx)) ) continue;

      TLorentzVector muon_1_4vec, muon_2_4vec;
      muon_1_4vec.SetPtEtaPhiM(muon_1.pt(), muon_1.eta(), muon_1.phi(), muon_1.mass());
      muon_2_4vec.SetPtEtaPhiM(muon_2.pt(), muon_2.eta(), muon_2.phi(), muon_2.mass());
      double mass = (muon_1_4vec + muon_2_4vec).M();
      //MASS->Fill(mass);
      if( 70 < mass && mass < 110 ) { muons.push_back(muon_1); muons.push_back(muon_2); return muons; }
    }
  }

  return muons;
}


std::vector<pat::Jet> DemoTrigAnls::jet_selection(edm::Handle<std::vector<pat::Jet>> JetHandle, std::vector<pat::Muon> muons)
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

std::vector<TString> DemoTrigAnls::SetTrigger(const edm::Event& iEvent, edm::Handle<edm::TriggerResults> triggerBits){
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

bool DemoTrigAnls::PassTrigger(TString trig, std::vector<TString> TriggerSets){
  std::vector<TString> tmp_vec;
  tmp_vec.push_back(trig);
  return PassTrigger(tmp_vec, TriggerSets);
}

bool DemoTrigAnls::PassTrigger(std::vector<TString> trigs, std::vector<TString> TriggerSets){

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


std::vector<pat::Muon> DemoTrigAnls::GetMuons(const reco::Vertex &vtx, edm::Handle<std::vector<pat::Muon>> MuonHandle, double pt, double eta, double iso, bool niso = 1){
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

std::vector<pat::Jet> DemoTrigAnls::GetJets(edm::Handle<std::vector<pat::Jet>> JetHandle, double pt, double eta){
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
    else {
      tightJetID        = (NHF<0.90 && NEMF<0.90 && NumConst>1) &&            ((fabs(eta_)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(eta_)>2.4) && fabs(eta_)<=2.7;
      //tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((fabs(eta_)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || fabs(eta_)>2.4) && fabs(eta_)<=2.7;
    }
    if( !tightJetID ) continue;
    out.push_back( jet );
  }
  return out;
}

std::vector<pat::Jet> DemoTrigAnls::GetBJets(edm::Handle<std::vector<pat::Jet>> JetHandle, double pt, double eta){
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
    else {
      tightJetID        = (NHF<0.90 && NEMF<0.90 && NumConst>1) &&            ((fabs(eta_)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(eta_)>2.4) && fabs(eta_)<=2.7;
      //tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((fabs(eta_)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || fabs(eta_)>2.4) && fabs(eta_)<=2.7;
    }
    if( !tightJetID ) continue;
    if( !(jet.bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll") > 0.6321) ) continue; // Loose: 0.2217, Medium: 0.6321, Tight: 0.8953

    out.push_back( jet );
  }
  return out;
}


// ------------ method called once each job just before starting event loop  ------------
void DemoTrigAnls::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void DemoTrigAnls::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
/*
void 
DemoTrigAnls::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DemoTrigAnls::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DemoTrigAnls::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DemoTrigAnls::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoTrigAnls::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoTrigAnls);
