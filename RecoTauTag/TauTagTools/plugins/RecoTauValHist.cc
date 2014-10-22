#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include <vector>
#include <string>
#include <sstream>

// histograms
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"

typedef edm::View<reco::GenJet> GenJetView;

class RecoTauValHist : public edm::EDFilter {
  public:
  explicit RecoTauValHist(const edm::ParameterSet& pset);
    virtual ~RecoTauValHist() {}
    virtual bool filter(edm::Event& evt, const edm::EventSetup& es);
    virtual void endJob();
  private:
    reco::tau::RecoTauQualityCuts qcuts_;
  std::auto_ptr<reco::tau::RecoTauVertexAssociator> vertexAssociator_;
    edm::InputTag src1_;
    edm::InputTag discDM_;
    edm::InputTag genSrc_;
    edm::InputTag genJetSrc_;
    edm::InputTag genTauSrc_;
    std::vector<edm::InputTag> tauIDSrcs_;
    std::vector<std::string> discNames_;
    double maxDeltaR_;
    size_t eventsExamined_;
    size_t tausExamined_;
    size_t trueTaus_;
    size_t tausMatched_;
    size_t differences_;
    size_t passed1_;
    size_t passed2_;
    size_t allPassed1_;
    size_t allPassed2_;
    bool filter_;
    bool mcMatch_;
    bool useGenTaus_;
    double matchingDistance_;
  bool matchToJets_;
  bool verboseOutput_;
  bool verboseOutputMC_;
  bool background_;
  bool Zmumu_;
  bool Zee_;
  bool onlyHadronic_;
  int requireDecayMode_;
  bool checkMother_;
  edm::InputTag vertexTag_;
  edm::InputTag rhoProducer_;

    TH1D* h_nVx;
    TH1D* h_summary;

    TH1D* h_mc_pt;
  TH1D* h_mc_eta;
  TH1D* h_mc_pt_vis_had;
  TH1D* h_mc_pt_vis_el;
  TH1D* h_mc_pt_vis_mu;
  TH1D* h_mc_eta_vis_had;
  TH1D* h_mc_eta_vis_el;
  TH1D* h_mc_eta_vis_mu;
  TH1D* h_decayMode;
  TH1D* h_rho;
  TH2D* h_nVx_rho;

  std::vector<TH1D*> h_discs_raw;
  std::vector<TProfile*> h_discs_pt;
  std::vector<TProfile*> h_discs_eta;
  std::vector<TProfile*> h_discs_phi;
  std::vector<TProfile*> h_discs_nVx;


};

RecoTauValHist::RecoTauValHist(
						     const edm::ParameterSet& pset): qcuts_(pset.exists("qualityCuts") ? pset.getParameterSet("qualityCuts").getParameterSet("isolationQualityCuts") : pset.getParameterSet("qualityCuts"))
 {
  src1_ = pset.getParameter<edm::InputTag>("src1");
  
  vertexAssociator_.reset(
			  new reco::tau::RecoTauVertexAssociator(pset.getParameterSet("qualityCuts"),consumesCollector()));
  discDM_ = pset.getParameter<edm::InputTag>("discDM");
  eventsExamined_ = 0;
  tausExamined_ = 0;
  tausMatched_ = 0;
  trueTaus_ = 0;
  differences_ = 0;
  passed1_ = 0;
  passed2_ = 0;
  allPassed2_ = 0;
  allPassed1_ = 0;
  discNames_.clear();
  filter_ = pset.exists("filter") ? pset.getParameter<bool>("filter") : false;
  mcMatch_ = pset.exists("mcMatch") ? pset.getParameter<bool>("mcMatch"): true;
  if( mcMatch_ && pset.exists("genSrc")) genSrc_= pset.getParameter<edm::InputTag>("genSrc");
  if( mcMatch_ && pset.exists("genJetSrc")) genJetSrc_= pset.getParameter<edm::InputTag>("genJetSrc");
  if(mcMatch_ && pset.exists("genTauSrc")) genTauSrc_ = pset.getParameter<edm::InputTag>("genTauSrc");
  matchingDistance_ = pset.exists("matchingDistance") ? pset.getParameter<double>("matchingDistance"): 0.1 ;
  matchToJets_ = pset.exists("matchToJets") ? pset.getParameter<bool>("matchToJets"): false;
  useGenTaus_ = pset.exists("useGenTaus") ? pset.getParameter<bool>("useGenTaus"): true;
  vertexTag_ = edm::InputTag("offlinePrimaryVertices", "");
  if (pset.exists("primaryVertexSrc")) vertexTag_ = pset.getParameter<edm::InputTag>("primaryVertexSrc");
  
  verboseOutput_= pset.exists("verboseOutput") ? pset.getParameter<bool>("verboseOutput"): true;
  verboseOutputMC_= pset.exists("verboseOutputMC") ? pset.getParameter<bool>("verboseOutputMC"): true;
  background_=pset.exists("background") ? pset.getParameter<bool>("background"): false;
  onlyHadronic_ = pset.exists("onlyHadronic") ? pset.getParameter<bool>("onlyHadronic"): true;
  rhoProducer_ = pset.getParameter<edm::InputTag>("rhoProducer");
  requireDecayMode_ = pset.exists("requireDecayMode") ? pset.getParameter<int>("requireDecayMode"): 0; // -1 No requirement; 0 = any DM; 1 = 1p; 2 = 1p+X; 3=3p
  checkMother_ = pset.exists("checkMother") ? pset.getParameter<bool>("checkMother"): 1;
  Zmumu_ = pset.exists("Zmumu") ? pset.getParameter<bool>("Zmumu"): false;
  Zee_ = pset.exists("Zee") ? pset.getParameter<bool>("Zee"): false;

  if (pset.existsAs<edm::ParameterSet>("idSrc")) {
    edm::ParameterSet idps = pset.getParameter<edm::ParameterSet>("idSrc");
    discNames_ = idps.getParameterNamesForType<edm::InputTag>();
    for (std::vector<std::string>::const_iterator it = discNames_.begin(), ed = discNames_.end(); it != ed; ++it) {
      std::cout << *it << std::endl;
      tauIDSrcs_.push_back( idps.getParameter<edm::InputTag>(*it));
    }
  }


  edm::Service<TFileService> fs;


  h_nVx = fs->make<TH1D>("h_nVx" , "Number of vertices; Reconstructed vertices" , 40 , -0.5 , 39.5 );
  h_rho = fs->make<TH1D>("h_rho" , "Rho; #rho" , 500 , 0.0 , 50.0 );
  h_nVx_rho = fs->make<TH2D>("h_nVx_rho" , "Number of vertices vs rho; Reconstructed vertices;#rho" , 40 , -0.5 , 39.5, 500, 0.0, 50.0 );
  h_summary = fs->make<TH1D>("h_summary","Summary of results",7,-0.5,6.5);

  h_summary->GetXaxis()->SetBinLabel(1,"Examined");
  h_summary->GetXaxis()->SetBinLabel(2,"Matched");
  h_summary->GetXaxis()->SetBinLabel(3,"Differences");
  h_summary->GetXaxis()->SetBinLabel(4,"Pass 1");
  h_summary->GetXaxis()->SetBinLabel(5,"Exclusive 1");
  h_summary->GetXaxis()->SetBinLabel(6,"Pass 2");
  h_summary->GetXaxis()->SetBinLabel(7,"Exclusive 2");


  TFileDirectory DiscResultsVal = fs->mkdir("DiscResultsVal");
  for (std::vector<std::string>::const_iterator it = discNames_.begin(), ed = discNames_.end(); it != ed; ++it) {
    h_discs_raw.push_back(DiscResultsVal.make<TH1D>(TString((*it)+"_raw"),TString(*it+"_raw; raw output; multiplicity"), 100 , 0.0 , 100.0));
    h_discs_pt.push_back(DiscResultsVal.make<TProfile>(TString((*it)+"_pt"),TString(*it+"_pt; p_{T} [GeV]; efficiency"), 1500 , 0.0 , 1500.0));
    h_discs_eta.push_back(DiscResultsVal.make<TProfile>(TString((*it)+"_eta"),TString(*it+"_eta; #eta; efficiency"), 200 , -5.0 , 5.0));
    h_discs_phi.push_back(DiscResultsVal.make<TProfile>(TString((*it)+"_phi"),TString(*it+"_phi; #phi; efficiency"), 100 , -3.2 , 3.2));
    h_discs_nVx.push_back(DiscResultsVal.make<TProfile>(TString((*it)+"_nVx"),TString(*it+"_nVx; vertices; efficiency"), 200 , -0.5 , 199.5));
  }

  h_mc_pt = fs->make<TH1D>("h_mc_pt", "True tau p_{T}; True #tau p_{T}[GeV]", 150, 0.0, 150.0);
  h_mc_pt_vis_had = fs->make<TH1D>("h_mc_pt_vis_had", "True tau visible p_{T}:hadronic modes; True #tau visible p_{T}[GeV]", 1500, 0.0, 1500.0);
  h_mc_pt_vis_el = fs->make<TH1D>("h_mc_pt_vis_el", "True tau visible p_{T}:electron; True #tau visible p_{T}[GeV]", 1500, 0.0, 1500.0);
  h_mc_pt_vis_mu = fs->make<TH1D>("h_mc_pt_vis_mu", "True tau visible p_{T}:muon; True #tau visible p_{T}[GeV]", 1500, 0.0, 1500.0);

  h_mc_eta = fs->make<TH1D>("h_mc_eta", "True tau #eta; True #tau #eta", 100, -3.0, 3.0);
  h_mc_eta_vis_had = fs->make<TH1D>("h_mc_eta_vis_had", "True tau visible #eta:hadronic modes; True #tau visible #eta", 100, -3.0, 3.0);
  h_mc_eta_vis_el = fs->make<TH1D>("h_mc_eta_vis_el", "True tau visible #eta:electron; True #tau visible #eta", 100, -3.0, 3.0);
  h_mc_eta_vis_mu = fs->make<TH1D>("h_mc_eta_vis_mu", "True tau visible #eta:muon; True #tau visible #eta", 100, -3.0, 3.0);


  h_decayMode = fs->make<TH1D>("h_decayMode","Tau Decay Modes", 4, -1.5, 2.5);

  h_decayMode->GetXaxis()->SetBinLabel(1,"Other");
  h_decayMode->GetXaxis()->SetBinLabel(2,"Hadronic");
  h_decayMode->GetXaxis()->SetBinLabel(3,"Electron");
  h_decayMode->GetXaxis()->SetBinLabel(4,"Muon");

  //MVA histograms
  

}

namespace {
  reco::PFJetRef getJetRef(const reco::PFTau& tau) {
    if (tau.jetRef().isNonnull())
      return tau.jetRef();
    else if (tau.pfTauTagInfoRef()->pfjetRef().isNonnull())
      return tau.pfTauTagInfoRef()->pfjetRef();
    else throw cms::Exception("cant find jet ref");
  }
}

bool RecoTauValHist::filter(
    edm::Event& evt, const edm::EventSetup& es) {
  eventsExamined_++;
  unsigned int nMatched = 0;
  // Get taus
  edm::Handle<reco::PFTauCollection> taus1;
  evt.getByLabel(src1_, taus1);
  edm::Handle<reco::GenParticleCollection> genParticles;
  if(mcMatch_) evt.getByLabel(genSrc_, genParticles);
  edm::Handle<reco::GenJetCollection> genJets;
  if (matchToJets_) evt.getByLabel(genJetSrc_, genJets);
  edm::Handle<std::vector<reco::GenJet> > genTaus;
  if(mcMatch_) evt.getByLabel(genTauSrc_, genTaus);
 // Get discriminators
  std::vector<edm::Handle<reco::PFTauDiscriminator> > tauIDdiscs;
  for ( size_t i = 0; i < tauIDSrcs_.size(); ++i ) {
    edm::Handle<reco::PFTauDiscriminator> discTemp;
    evt.getByLabel(tauIDSrcs_[i],discTemp);
    tauIDdiscs.push_back( discTemp);
  }

  edm::Handle<reco::PFTauDiscriminator> discDM;
  evt.getByLabel(discDM_,discDM);
  edm::Handle<reco::VertexCollection> verticesH_;
  evt.getByLabel(vertexTag_, verticesH_);
  int nVx = verticesH_->size();
   h_nVx->Fill(nVx);
  vertexAssociator_->setEvent(evt);
  
  bool differenceFound = false;
  // Loop over first collection
  std::vector<int> goodCands;
  std::vector<bool> isHadronic;
  std::vector<double> pt_visible;
  std::vector<double> eta_visible;
  std::vector<double> phi_visible;
  std::vector<bool> isMatched;
  
  // edm::Handle<double> hRho;
  // evt.getByLabel(rhoProducer_, hRho);
  // double rho_ = *hRho;
  // h_rho->Fill(rho_);
  // h_nVx_rho->Fill(nVx,rho_);


  if(mcMatch_ && !useGenTaus_){
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const reco::GenParticle & ParticleCand = (*genParticles)[i];
    if(verboseOutputMC_){
    std::cout<<i<<") Particle = "<<ParticleCand.pdgId()<<", Status = "<<ParticleCand.status()<<std::endl;
    if(ParticleCand.pdgId()==15 || ParticleCand.pdgId() == -15 || true){
    const reco::GenParticleRefVector& motherRefs = ParticleCand.motherRefVector();
    for(reco::GenParticleRefVector::const_iterator imr = motherRefs.begin(); imr!= motherRefs.end(); ++imr) {
      std::cout<<"   - Mother "<<(*imr).key()<<" "<<(*imr)->pdgId()<<std::endl;
    }
    }
    }
    if( abs(ParticleCand.pdgId())==15 && ParticleCand.status() == 2 && !Zmumu_ && !Zee_)
      {
	if(verboseOutputMC_) std::cout << "  tau candidate!" << std::endl;
	if(background_ || ! checkMother_ ) goodCands.push_back(i);
	else{
	int mother = 0;
	bool final = false;
	int myKey = i;
	while(mother!=23 && !final)
	  {
	    const reco::GenParticleRefVector& mRefs = ((*genParticles)[myKey]).motherRefVector();
	    if(mRefs.size()==0){final=true; break;}
	    for(reco::GenParticleRefVector::const_iterator imr = mRefs.begin(); imr!= mRefs.end(); ++imr) {
	      if((*imr)->pdgId()==23){ mother = 23; break;}
	      if(abs(ParticleCand.pdgId())==15){  myKey =(*imr).key(); break;}
	    }
	  }
	if(mother==23) goodCands.push_back(i);
	}
      }
    if( abs(ParticleCand.pdgId())==13 && ParticleCand.status() == 1 && Zmumu_ && !Zee_)
	  {
	    if(verboseOutputMC_) std::cout << "  mu candidate!" << std::endl;
	    if(background_ || ! checkMother_ ) goodCands.push_back(i);
	    else{
	      int mother = 0;
	      bool final = false;
	      int myKey = i;
	      while(mother!=23 && !final)
		{
		  const reco::GenParticleRefVector& mRefs = ((*genParticles)[myKey]).motherRefVector();
		  if(mRefs.size()==0){final=true; break;}
		  for(reco::GenParticleRefVector::const_iterator imr = mRefs.begin(); imr!= mRefs.end(); ++imr) {
		    if((*imr)->pdgId()==23){ mother = 23; break;}
		    if(abs(ParticleCand.pdgId())==13){  myKey =(*imr).key(); break;}
		  }
		}
	      if(mother==23) goodCands.push_back(i);
	    }
	  }
    if( abs(ParticleCand.pdgId())==11 && ParticleCand.status() == 1 && Zee_ && !Zmumu_)
      {
	if(verboseOutputMC_) std::cout << "  e candidate!" << std::endl;
	if(background_ || ! checkMother_ ) goodCands.push_back(i);
	else{
	  int mother = 0;
	  bool final = false;
	  int myKey = i;
	  while(mother!=23 && !final)
	    {
	      const reco::GenParticleRefVector& mRefs = ((*genParticles)[myKey]).motherRefVector();
	      if(mRefs.size()==0){final=true; break;}
	      for(reco::GenParticleRefVector::const_iterator imr = mRefs.begin(); imr!= mRefs.end(); ++imr) {
		if((*imr)->pdgId()==23){ mother = 23; break;}
		if(abs(ParticleCand.pdgId())==11){  myKey =(*imr).key(); break;}
	      }
	    }
	  if(mother==23) goodCands.push_back(i);
	}
      }


  }
  }
    if(verboseOutputMC_ && !Zmumu_ && !Zee_) std::cout << "There were " << goodCands.size() << " good taus." << std::endl; 
    if(verboseOutputMC_ && Zmumu_ && !Zee_) std::cout << "There were " << goodCands.size() << " good mus." << std::endl; 
    if(verboseOutputMC_ && Zee_ && !Zmumu_) std::cout << "There were " << goodCands.size() << " good els." << std::endl;
 
if(!background_ && mcMatch_ && !useGenTaus_){
 //calculate visible pt
  for(size_t iMC = 0; iMC < goodCands.size(); ++iMC){
  //daughters                                                                                                                                                                                                                    
  const reco::GenParticleRefVector& mRefs = ((*genParticles)[goodCands.at(iMC)]).daughterRefVector();
  reco::Particle::LorentzVector invisibleP4( 0.0, 0.0, 0.0, 0.0 );
  unsigned int decayMode = -1; // 0 = hadronic, 1=electron, 2=muon                                                                                                                                                               
  unsigned int nNu = 0;
  for(reco::GenParticleRefVector::const_iterator imr = mRefs.begin(); imr!= mRefs.end(); ++imr) {
    if(abs((*imr)->pdgId())==11) decayMode = 1;
    if(abs((*imr)->pdgId())==13) decayMode = 2;
    if(abs((*imr)->pdgId())==16 || abs((*imr)->pdgId())==14 || abs((*imr)->pdgId())==12)
      {
	if(verboseOutputMC_) std::cout << "Found neutrino! " << (*imr)->pdgId() << " momentum is " << (*imr)->p4() << std::endl;
	invisibleP4+=(*imr)->p4();
	nNu++;
      }
  }
  if(verboseOutputMC_) std::cout << "Total invisible momentum is " << invisibleP4 << std::endl;
  if(nNu==1) decayMode = 0;
  h_decayMode->Fill(decayMode);
  h_mc_pt->Fill(((*genParticles)[goodCands.at(iMC)]).pt());
  h_mc_eta->Fill(((*genParticles)[goodCands.at(iMC)]).eta());

 reco::Particle::LorentzVector visibleP4 = ((*genParticles)[goodCands.at(iMC)]).p4()-invisibleP4;
  switch(decayMode){
   case 0:
     h_mc_pt_vis_had->Fill(visibleP4.pt());
     h_mc_eta_vis_had->Fill(visibleP4.eta());
     isHadronic.push_back(true);
     break;
   case 1:
     h_mc_pt_vis_el->Fill(visibleP4.pt());
     h_mc_eta_vis_el->Fill(visibleP4.eta());
     isHadronic.push_back(false);
     break;
   case 2: 
     h_mc_pt_vis_mu->Fill(visibleP4.pt());
     h_mc_eta_vis_mu->Fill(visibleP4.eta());
     isHadronic.push_back(false);
     break;
   default:
     isHadronic.push_back(false);
     break;
  }
  pt_visible.push_back(visibleP4.pt());
  eta_visible.push_back(visibleP4.eta());
  phi_visible.push_back(visibleP4.phi());
  isMatched.push_back(false);
  }
 }
 if(useGenTaus_){
   for(size_t i = 0; i < genTaus->size(); ++ i) {
     const reco::GenJet & TauCand = (*genTaus)[i];
     const std::vector <const reco::GenParticle*> mRefs = TauCand.getGenConstituents();
     unsigned int decayMode = 0; // 0 = hadronic, 1=electron, 2=muon 
     for(size_t igTauD =0; igTauD < mRefs.size(); igTauD++) {
       if(abs(mRefs[igTauD]->pdgId())==11) decayMode = 1;
       if(abs(mRefs[igTauD]->pdgId())==13) decayMode = 2;
     }
     if(decayMode==0) trueTaus_++;
     h_decayMode->Fill(decayMode);
     h_mc_pt->Fill(((*genTaus)[i]).pt());
     h_mc_eta->Fill(((*genTaus)[i]).eta());

     reco::Particle::LorentzVector visibleP4 = ((*genTaus)[i]).p4();
     switch(decayMode){
     case 0:
       h_mc_pt_vis_had->Fill(visibleP4.pt());
       h_mc_eta_vis_had->Fill(visibleP4.eta());
       isHadronic.push_back(true);
       break;
     case 1:
       h_mc_pt_vis_el->Fill(visibleP4.pt());
       h_mc_eta_vis_el->Fill(visibleP4.eta());
       isHadronic.push_back(false);
       break;
     case 2: 
       h_mc_pt_vis_mu->Fill(visibleP4.pt());
       h_mc_eta_vis_mu->Fill(visibleP4.eta());
       isHadronic.push_back(false);
       break;
     default:
       isHadronic.push_back(false);
       break;
     }
     pt_visible.push_back(visibleP4.pt());
     eta_visible.push_back(visibleP4.eta());
     phi_visible.push_back(visibleP4.phi());
     isMatched.push_back(false);
   }
 }
  for (size_t iTau1 = 0; iTau1 < taus1->size(); ++iTau1) {
    tausExamined_++;
    reco::PFTauRef tau1(taus1, iTau1);
    // Find the best match in the other collection
    if(verboseOutputMC_ && mcMatch_) std::cout << "Finding true match for tau with pt " << tau1->pt() << " !" << std::endl;
    bool matched = false;
    if(background_) matched = true;
    if(!background_ && !useGenTaus_){
      if(isHadronic.size()!=goodCands.size() || pt_visible.size()!=goodCands.size() || 
       isMatched.size()!=goodCands.size() || eta_visible.size()!=goodCands.size()){  edm::LogWarning("tauDifferenceAnalyzer") << "Wrong array size! "; break;}
    }
    double pt_vis = -999.;
    double eta_vis = -999.;
    double phi_vis = -999.;
    int matchIndex = -1;
    for(size_t iMC = 0; iMC < goodCands.size(); ++iMC){
      if(!Zmumu_ && !Zee_ && !background_ && onlyHadronic_ && !isHadronic.at(iMC))
	{
	  if(verboseOutputMC_) std::cout << "Not hadronic tau!" << std::endl;
	    continue;
	}
      if(!background_ && isMatched.at(iMC)) continue;
        double dR_MC = deltaR(tau1->p4(),((*genParticles)[goodCands.at(iMC)]).p4());
      if(verboseOutputMC_) std::cout << " Tau distance is " << dR_MC << std::endl;
      if(dR_MC < matchingDistance_) 
	{ 
	  if(!background_){
	    matched = true;
	    pt_vis=pt_visible.at(iMC);
	    eta_vis = eta_visible.at(iMC);
	    phi_vis = phi_visible.at(iMC);
	  }
	  else matched = false;
	  matchIndex=iMC;
       	}
    }
    if(useGenTaus_)
      {
	for(size_t iMC = 0; iMC < genTaus->size() && matched==false; ++ iMC) {
	  if(!Zmumu_ && !Zee_ && !background_ && onlyHadronic_ && !isHadronic.at(iMC))
	    {
	      if(verboseOutputMC_) std::cout << "Not hadronic tau!" << std::endl;
	      continue;
	    }
	  if(!background_ && isMatched.at(iMC)) continue;
	  double dR_MC = deltaR(tau1->p4(),((*genTaus)[iMC]).p4());
	  if(verboseOutputMC_) std::cout << " Tau distance is " << dR_MC << std::endl;
	  if(dR_MC < matchingDistance_) 
	    { 
	      if(!background_){
		matched = true;
		pt_vis=pt_visible.at(iMC);
		eta_vis = eta_visible.at(iMC);
		phi_vis = phi_visible.at(iMC);
		if(verboseOutputMC_) std::cout << "This reco tau is matched to a true tau with pt " << pt_vis << std::endl;
	      }
	      else matched = false;
	      matchIndex=iMC;
	    }
	}
      }
    if(!background_ && matched) isMatched.at(matchIndex)=true;
    if(matchToJets_){
      matched = false;
    for(size_t iJet=0; iJet < genJets->size(); ++iJet)
    {
	const reco::GenJet & JetCand = (*genJets)[iJet];
	double dR_jet = deltaR(tau1->p4(),JetCand.p4());
	if(dR_jet < matchingDistance_) matched = true;
    }
}


    if((mcMatch_ && !matched)||(matchToJets_&& !matched)) continue;
    tausMatched_++; nMatched++;
    // now we have a reco tau, that is matched to truth
    // See what's up with the discriminators
    double resultDM = (*discDM)[tau1];
    double pt1 = tau1->pt();
    double pt_mon, eta_mon, phi_mon;
 
      pt_mon=pt_vis;
      eta_mon=eta_vis;
      phi_mon=phi_vis;
 


    if(pt_mon < 5.0) continue;
    
    if(requireDecayMode_>=0 && resultDM){

      if(pt_mon>5.0){
	for ( size_t i = 0; i < tauIDdiscs.size(); ++i ) {
	double result = (*(tauIDdiscs[i]))[tau1];
	h_discs_raw[i]->Fill(result);
	if(pt1>20.0){
	  h_discs_pt[i]->Fill(pt_mon,result);
	  h_discs_eta[i]->Fill(eta_mon,result);
	  h_discs_phi[i]->Fill(phi_mon,result);
	  h_discs_nVx[i]->Fill(nVx,result);
	}
	else h_discs_pt[i]->Fill(pt_mon,0);
	}
      }


    }
  }

  
 return (filter_ ? differenceFound : true);
}


void RecoTauValHist::endJob() {
  if(verboseOutput_){
  std::cout <<  " RECO TAU DIFFERENCE SUMMARY: " << std::endl;
  std::cout <<  " Examined " << tausExamined_ << " taus in "
    << eventsExamined_ << " events." << std::endl;
  std::cout << " There were " << trueTaus_ << " true hadronic taus " << std::endl;
  std::cout << " There were " << tausMatched_ << " taus matched to MC." << std::endl; 
  }
  h_summary->SetBinContent(1,tausExamined_);
  h_summary->SetBinContent(2,tausMatched_);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoTauValHist);
