/**
  \class    pat::PATTauSlimmer PATTauSlimmer.h "PhysicsTools/PatAlgos/interface/PATTauSlimmer.h"
  \brief    Slimmer of PAT Taus 
*/


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

namespace pat {

  class PATTauSlimmer : public edm::EDProducer {
    public:
      explicit PATTauSlimmer(const edm::ParameterSet & iConfig);
      virtual ~PATTauSlimmer() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      edm::EDGetTokenT<edm::View<pat::Tau> > src_;
      bool linkToPackedPF_;
      edm::EDGetTokenT<edm::Association<pat::PackedCandidateCollection>> pf2pc_;

  };

} // namespace

pat::PATTauSlimmer::PATTauSlimmer(const edm::ParameterSet & iConfig) :
    src_(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("src"))),
    linkToPackedPF_(iConfig.getParameter<bool>("linkToPackedPFCandidates"))
{
    produces<std::vector<pat::Tau> >();
    if (linkToPackedPF_) pf2pc_ = consumes<edm::Association<pat::PackedCandidateCollection>>(iConfig.getParameter<edm::InputTag>("packedPFCandidates"));
}

void 
pat::PATTauSlimmer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<pat::Tau> >      src;
    iEvent.getByToken(src_, src);

    Handle<edm::Association<pat::PackedCandidateCollection>> pf2pc;
    if (linkToPackedPF_) iEvent.getByToken(pf2pc_, pf2pc);

    auto_ptr<vector<pat::Tau> >  out(new vector<pat::Tau>());
    out->reserve(src->size());
    std::cout << "Now making slimmed tau collection" << std::endl;
    for (View<pat::Tau>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
      std::cout << " pushing tau to slimmedTauCollection" << *it << std::endl;
        out->push_back(*it);
        if (linkToPackedPF_) {
            pat::Tau & tau = out->back();
	    std::cout << "number of source signal (PF) cands are " << tau.signalPFCands() << " / " << tau.signalCands().size() << std::endl;
	    std::cout << "number of source signalChargedHadr (PF) cands are " << tau.signalChargedHadrPFCands() << " / " << tau.signalChargedHadrCands().size() << std::endl;
	    std::cout << "number of source signalNeutrHadr (PF) cands are " << tau.signalNeutrHadrPFCands() << " / " << tau.signalNeutrHadrCands().size() << std::endl;
	    std::cout << "number of source signalGamma (PF) cands are " << tau.signalGammaPFCands() << " / " << tau.signalGammaCands().size() << std::endl;

	    std::cout << "number of source isolation (PF) cands are " << tau.isolationPFCands() << " / " << tau.isolationCands().size() << std::endl;
	    std::cout << "number of source isolationChargedHadr (PF) cands are " << tau.isolationChargedHadrPFCands() << " / " << tau.isolationChargedHadrCands().size() << std::endl;
	    std::cout << "number of source isolationNeutrHadr (PF) cands are " << tau.isolationNeutrHadrPFCands() << " / " << tau.isolationNeutrHadrCands().size() << std::endl;
	    std::cout << "number of source isolationGamma (PF) cands are " << tau.isolationGammaPFCands() << " / " << tau.isolationGammaCands().size() << std::endl;

	    std::cout << "linking to packedPF " << tau << std::endl;
            reco::CandidatePtrVector signalPtrs, signalChHPtrs, signalNHPtrs, signalGammaPtrs, isolationPtrs, isolationChHPtrs, isolationNHPtrs, isolationGammaPtrs;
            for (const reco::PFCandidatePtr &p : tau.signalPFCands()) {
                signalPtrs.push_back(edm::refToPtr((*pf2pc)[p]));
            }
            tau.setSignalCands(signalPtrs);

	    for (const reco::PFCandidatePtr &p : tau.signalChargedHadrPFCands()) {
	      signalChHPtrs.push_back(edm::refToPtr((*pf2pc)[p]));
            }
            tau.setSignalChargedHadronCands(signalChHPtrs);

	    for (const reco::PFCandidatePtr &p : tau.signalNeutralHadrPFCands()) {
              signalNHPtrs.push_back(edm::refToPtr((*pf2pc)[p]));
            }
            tau.setSignalNeutralHadronCands(signalNHPtrs);

	    for (const reco::PFCandidatePtr &p : tau.signalGammaPFCands()) {
              signalGammaPtrs.push_back(edm::refToPtr((*pf2pc)[p]));
            }
            tau.setSignalGammaCands(signalGammaPtrs);

	    for (const reco::PFCandidatePtr &p : tau.isolationPFCands()) {
                isolationPtrs.push_back(edm::refToPtr((*pf2pc)[p]));
            }
            tau.setIsolationCands(isolationPtrs);

	    for (const reco::PFCandidatePtr &p : tau.isolationChargedHadrPFCands()) {
              isolationChHPtrs.push_back(edm::refToPtr((*pf2pc)[p]));
            }
            tau.setIsolationChargedHadronCands(isolationChHPtrs);

            for (const reco::PFCandidatePtr &p : tau.isolationNeutralHadrPFCands()) {
              isolationNHPtrs.push_back(edm::refToPtr((*pf2pc)[p]));
            }
            tau.setIsolationNeutralHadronCands(isolationNHPtrs);

            for (const reco::PFCandidatePtr &p : tau.isolationGammaPFCands()) {
              isolationGammaPtrs.push_back(edm::refToPtr((*pf2pc)[p]));
            }
            tau.setIsolationGammaCands(isolationGammaPtrs);

	    std::cout << "  After linking: " << std::endl;
	    std::cout << "   number of source signal cands are " << tau.signalCands().size() << std::endl;
	    std::cout << "   number of source signalChargedHadr cands are " << tau.signalChargedHadrCands().size() << std::endl;
	    std::cout << "   number of source signalNeutrHadr cands are " << tau.signalNeutrHadrCands().size() << std::endl;
	    std::cout << "   number of source signalGamma cands are " << tau.signalGammaCands().size() << std::endl;
	    
	    std::cout << "number of source isolation cands are " << tau.isolationCands().size() << std::endl;
	    std::cout << "   number of source isolationChargedHadr cands are " << tau.isolationChargedHadrCands().size() << std::endl;
	    std::cout << "   number of source isolationNeutrHadr cands are " << tau.isolationNeutrHadrCands().size() << std::endl;
	    std::cout << "   number of source isolationGamma cands are " << tau.isolationGammaCands().size() << std::endl;
        }
    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace pat;
DEFINE_FWK_MODULE(PATTauSlimmer);
