/*
 * "Unembed" the pizeros in a reco::PFTau.
 *
 * This converts a collection of PFTaus which have their PiZeros stored
 * as std::vector<RecoTauPiZeros>s to an output collection which has
 * the PiZeros stored in a separate product, with the PiZeros stored as Refs
 * within the tau.  This will improve the de-serialization speed of the taus.
 *
 * Author: Evan K. Friis, UW Madison
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/TauReco/interface/RecoTauPiZeroFwd.h"


class RecoTauPiZeroUnembedder : public edm::stream::EDProducer<> {
  public:
    RecoTauPiZeroUnembedder(const edm::ParameterSet& pset);
    virtual ~RecoTauPiZeroUnembedder(){}
    void produce(edm::Event& evt, const edm::EventSetup& es) override;
  private:
    edm::InputTag src_;
    edm::EDGetTokenT<reco::CandidateView> token; 
};

RecoTauPiZeroUnembedder::RecoTauPiZeroUnembedder(const edm::ParameterSet& pset) {
  src_ = pset.getParameter<edm::InputTag>("src");
  token = consumes<reco::CandidateView>(src_);
  produces<reco::RecoTauPiZeroCollection>("pizeros");
  produces<reco::PFRecoTauChargedHadronCollection>("recoTauChHadrons");
  produces<reco::PFTauCollection>();
}
void RecoTauPiZeroUnembedder::produce(edm::Event& evt, const edm::EventSetup& es) {
  std::auto_ptr<reco::RecoTauPiZeroCollection> piZerosOut(
      new reco::RecoTauPiZeroCollection);
  std::auto_ptr<reco::PFRecoTauChargedHadronCollection> chHadronsOut( new reco::PFRecoTauChargedHadronCollection);
  std::auto_ptr<reco::PFTauCollection> tausOut(new reco::PFTauCollection);

  edm::Handle<reco::CandidateView> tauView;
  evt.getByToken(token, tauView);

  reco::PFTauRefVector taus =
      reco::tau::castView<reco::PFTauRefVector>(tauView);

  // Get the reference to the products of where the final pizeros and chHadrons will end up
  reco::RecoTauPiZeroRefProd piZeroProd =
    evt.getRefBeforePut<reco::RecoTauPiZeroCollection>("pizeros");

  reco::PFRecoTauChargedHadronRefProd chHadronProd = evt.getRefBeforePut<reco::PFRecoTauChargedHadronCollection>("recoTauChHadrons");

  for (size_t iTau = 0; iTau < taus.size(); ++iTau) {
    // Make a copy
    reco::PFTau myTau = *taus[iTau];
    // The ref vectors that will be filled
    reco::RecoTauPiZeroRefVector signalPiZeroRefs;
    reco::RecoTauPiZeroRefVector isolationPiZeroRefs;

    reco::PFRecoTauChargedHadronRefVector signalTauChargedHadronRefs;
    reco::PFRecoTauChargedHadronRefVector isolationTauChargedHadronRefs;

    // Copy the PiZeros and chHadrons into the new vectors, while updating what refs they will
    // have
    const reco::RecoTauPiZeroCollection& signalPiZeros =
      myTau.signalPiZeroCandidates();

    for (size_t iPiZero = 0; iPiZero < signalPiZeros.size(); ++iPiZero) {
      piZerosOut->push_back(signalPiZeros[iPiZero]);
      // Figure out what the ref for this pizero will be in the new coll.
      signalPiZeroRefs.push_back(
          reco::RecoTauPiZeroRef(piZeroProd, piZerosOut->size()-1));
    }

    const reco::PFRecoTauChargedHadronCollection& signalChHadrons = 
      myTau.signalTauChargedHadronCandidates();

    for (size_t iChHadron = 0; iChHadron < signalChHadrons.size(); ++iChHadron) {
      std::cout << "HEH, UNEMBEDDING DUDE..." << iChHadron << " pt= " << signalChHadrons[iChHadron].pt() << std::endl;
      chHadronsOut->push_back(signalChHadrons[iChHadron]);
      // Figure out what the ref for this chHadron will be in the new coll.
      signalTauChargedHadronRefs.push_back(
				 reco::PFRecoTauChargedHadronRef(chHadronProd, chHadronsOut->size()-1));
    }

    const reco::RecoTauPiZeroCollection& isolationPiZeroCandidates =
      myTau.isolationPiZeroCandidates();
    for (size_t iPiZero = 0; iPiZero < isolationPiZeroCandidates.size(); ++iPiZero) {
      piZerosOut->push_back(isolationPiZeroCandidates[iPiZero]);
      // Figure out what the ref for this pizero will be in the new coll.
      isolationPiZeroRefs.push_back(
          reco::RecoTauPiZeroRef(piZeroProd, piZerosOut->size()-1));
    }

    const reco::PFRecoTauChargedHadronCollection& isolationChHadrons = 
      myTau.isolationTauChargedHadronCandidates();

    for (size_t iChHadron = 0; iChHadron < isolationChHadrons.size(); ++iChHadron) {
      std::cout << "HEH, UNEMBEDDING isoDUDE..." << iChHadron << " pt= " << isolationChHadrons[iChHadron].pt() << std::endl;
      chHadronsOut->push_back(isolationChHadrons[iChHadron]);
      // Figure out what the ref for this chHadron will be in the new coll.
      isolationTauChargedHadronRefs.push_back(
					reco::PFRecoTauChargedHadronRef(chHadronProd, chHadronsOut->size()-1));
    }

    myTau.setSignalPiZeroCandidatesRefs(signalPiZeroRefs);
    myTau.setIsolationPiZeroCandidatesRefs(isolationPiZeroRefs);
    myTau.setSignalTauChargedHadronCandidatesRefs(signalTauChargedHadronRefs);
    myTau.setIsolationTauChargedHadronCandidatesRefs(isolationTauChargedHadronRefs);

    tausOut->push_back(myTau);
  }

  evt.put(piZerosOut, "pizeros");
  evt.put(chHadronsOut, "recoTauChHadrons");
  evt.put(tausOut);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoTauPiZeroUnembedder);
