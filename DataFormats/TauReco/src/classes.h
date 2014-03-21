#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/TauReco/interface/BaseTauTagInfo.h"
#include "DataFormats/TauReco/interface/CaloTauTagInfo.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"
#include "DataFormats/TauReco/interface/BaseTau.h"
#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/TauReco/interface/RecoTauPiZeroFwd.h"
#include "DataFormats/TauReco/interface/CaloTauDiscriminatorByIsolation.h"
#include "DataFormats/TauReco/interface/CaloTauDiscriminator.h"
#include "DataFormats/TauReco/interface/CaloTauDiscriminatorAgainstElectron.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminatorByIsolation.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/TauReco/interface/JetPiZeroAssociation.h"
#include "DataFormats/TauReco/interface/PFTauDecayModeAssociation.h"
#include "DataFormats/TauReco/interface/L2TauInfoAssociation.h"
#include "DataFormats/TauReco/interface/HLTTau.h"
#include "DataFormats/TauReco/interface/PFRecoTauChargedHadron.h"
#include "DataFormats/TauReco/interface/PFRecoTauChargedHadronFwd.h"
#include "DataFormats/TauReco/interface/PFJetChargedHadronAssociation.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterAssociation.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterFwd.h"
#include "DataFormats/TauReco/interface/PFTau3ProngSummaryFwd.h"
#include "DataFormats/TauReco/interface/PFTau3ProngSummaryAssociation.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include <vector>
#include <map>
#include "TLorentzVector.h"

namespace {
  struct dictionary {
    reco::L2TauIsolationInfo                                    l2iI;
    reco::L2TauInfoAssociation                                  l2ts;
    edm::Wrapper<reco::L2TauInfoAssociation>                    wl2ts;

    std::vector<reco::BaseTauTagInfo>                           btti_v;
    edm::Wrapper<std::vector<reco::BaseTauTagInfo> >            btti_w;
    edm::Ref<std::vector<reco::BaseTauTagInfo> >                btti_r;
    edm::RefProd<std::vector<reco::BaseTauTagInfo> >            btti_rp;
    edm::RefVector<std::vector<reco::BaseTauTagInfo> >          btti_rv;

    std::vector<reco::CaloTauTagInfo>                           calotti_v;
    edm::Wrapper<std::vector<reco::CaloTauTagInfo> >            calotti_w;
    edm::Ref<std::vector<reco::CaloTauTagInfo> >                calotti_r;
    edm::RefProd<std::vector<reco::CaloTauTagInfo> >            calotti_rp;
    edm::RefVector<std::vector<reco::CaloTauTagInfo> >          calotti_rv;

    std::vector<reco::PFTauTagInfo>                             pftti_v;
    edm::Wrapper<std::vector<reco::PFTauTagInfo> >              pftti_w;
    edm::Ref<std::vector<reco::PFTauTagInfo> >                  pftti_r;
    edm::RefProd<std::vector<reco::PFTauTagInfo> >              pftti_rp;
    edm::RefVector<std::vector<reco::PFTauTagInfo> >            pftti_rv;

    std::vector<reco::BaseTau>                                  bt_v;
    edm::Wrapper<std::vector<reco::BaseTau> >                   bt_w;
    edm::Ref<std::vector<reco::BaseTau> >                       bt_r;
    edm::RefProd<std::vector<reco::BaseTau> >                   bt_rp;
    edm::RefVector<std::vector<reco::BaseTau> >                 bt_rv;
    edm::reftobase::Holder<reco::Candidate,reco::BaseTauRef>    bt_rb;

    std::vector<reco::CaloTau>                                  calot_v;
    edm::Wrapper<std::vector<reco::CaloTau> >                   calot_w;
    edm::Ref<std::vector<reco::CaloTau> >                       calot_r;
    edm::RefProd<std::vector<reco::CaloTau> >                   calot_rp;
    edm::RefVector<std::vector<reco::CaloTau> >                 calot_rv;
    edm::reftobase::Holder<reco::BaseTau,reco::CaloTauRef>      calot_rb;
    edm::RefToBaseVector<reco::CaloTau> calot_rftbv;
    edm::Wrapper<edm::RefToBaseVector<reco::CaloTau> > calot_rftbv_w;

    std::vector<reco::PFTau>                                    pft_v;
    edm::Wrapper<std::vector<reco::PFTau> >                     pft_w;
    edm::Ref<std::vector<reco::PFTau> >                         pft_r;
    edm::RefProd<std::vector<reco::PFTau> >                     pft_rp;
    edm::RefVector<std::vector<reco::PFTau> >                   pft_rv;
    edm::Wrapper<edm::RefVector<std::vector<reco::PFTau> > >         pft_rvw;
    edm::reftobase::Holder<reco::BaseTau,reco::PFTauRef>        pft_rb;

    edm::View<reco::PFTau>  jv;
    edm::RefToBaseProd<reco::PFTau> jrtbp;

    edm::RefToBase<reco::PFTau>  rtbj;
    edm::reftobase::IndirectHolder<reco::PFTau> ihj;
    edm::reftobase::Holder<reco::PFTau, reco::PFTauRef> hcj;
    edm::reftobase::RefHolder<reco::PFTauRef> rhtj;
    edm::RefToBaseVector<reco::PFTau> jrtbv;
    edm::Wrapper<edm::RefToBaseVector<reco::PFTau> > jrtbv_w;
    edm::reftobase::BaseVectorHolder<reco::PFTau> * bvhj_p;    // pointer since it's pure virtual

    std::vector<reco::PFTauDecayMode>                                           pftdm_v;
    edm::Wrapper<std::vector<reco::PFTauDecayMode> >                            pftdm_w;
    edm::Ref<std::vector<reco::PFTauDecayMode> >                                pftdm_r;
    edm::RefProd<std::vector<reco::PFTauDecayMode> >                            pftdm_rp;
    edm::RefVector<std::vector<reco::PFTauDecayMode> >                          pftdm_rv;
    edm::reftobase::Holder<reco::CompositeCandidate,reco::PFTauDecayModeRef>    pftdm_rb;
    edm::Association<std::vector<reco::PFTauDecayMode> >                        pftdm_assoc_v;
    edm::Association<std::vector<reco::PFTau> >                                 pftau_assoc_v; // used for matching
    edm::Wrapper<edm::Association<std::vector<reco::PFTauDecayMode> > >         pftdm_assoc_v_wrapper;
    edm::Wrapper<edm::Association<std::vector<reco::PFTau> > >                  pftau_assoc_v_wrapper;

    std::vector<reco::RecoTauPiZero>                                           recoTauPiZero_v;
    edm::Wrapper<std::vector<reco::RecoTauPiZero> >                            recoTauPiZero_w;
    edm::Ref<std::vector<reco::RecoTauPiZero> >                                recoTauPiZero_r;
    edm::RefProd<std::vector<reco::RecoTauPiZero> >                            recoTauPiZero_rp;
    edm::RefVector<std::vector<reco::RecoTauPiZero> >                          recoTauPiZero_rv;
    edm::reftobase::Holder<reco::CompositePtrCandidate, reco::RecoTauPiZeroRef>    recoTauPiZero_rb;

    std::vector<reco::PFRecoTauChargedHadron>                                           pfrecoTauChH_v;
    edm::Wrapper<std::vector<reco::PFRecoTauChargedHadron> >                            pfrecoTauChH_w;
    edm::Ref<std::vector<reco::PFRecoTauChargedHadron> >                                pfrecoTauChH_r;
    edm::RefProd<std::vector<reco::PFRecoTauChargedHadron> >                            pfrecoTauChH_rp;
    edm::RefVector<std::vector<reco::PFRecoTauChargedHadron> >                          pfrecoTauChH_rv;
    edm::reftobase::Holder<reco::CompositePtrCandidate, reco::PFRecoTauChargedHadronRef>        pfrecoTauChH_rb;

    reco::CaloTauDiscriminatorByIsolationBase                   calotdi_b;
    reco::CaloTauDiscriminatorByIsolation                       calotdi_o;
    reco::CaloTauDiscriminatorByIsolationRef                    calotdi_r;
    reco::CaloTauDiscriminatorByIsolationRefProd                calotdi_rp;
    reco::CaloTauDiscriminatorByIsolationRefVector              calotdi_rv;
    edm::Wrapper<reco::CaloTauDiscriminatorByIsolation>         calotdi_w;

    reco::CaloTauDiscriminatorBase                   calodi_b;
    reco::CaloTauDiscriminator                       calodi_o;
    reco::CaloTauDiscriminatorRef                    calodi_r;
    reco::CaloTauDiscriminatorRefProd                calodi_rp;
    reco::CaloTauDiscriminatorRefVector              calodi_rv;
    edm::Wrapper<reco::CaloTauDiscriminator>         calodi_w;

    reco::CaloTauDiscriminatorAgainstElectronBase               calotde_b;
    reco::CaloTauDiscriminatorAgainstElectron                   calotde_o;
    reco::CaloTauDiscriminatorAgainstElectronRef                calotde_r;
    reco::CaloTauDiscriminatorAgainstElectronRefProd            calotde_rp;
    reco::CaloTauDiscriminatorAgainstElectronRefVector          calotde_rv;
    edm::Wrapper<reco::CaloTauDiscriminatorAgainstElectron>     calotde_w;

    std::pair<reco::CaloTauRef, int>                            calotd_p;
    std::vector<std::pair<reco::CaloTauRef, int> >              calotd_v;

    reco::PFTauDiscriminatorByIsolationBase                     pftdi_b;
    reco::PFTauDiscriminatorByIsolation                         pftdi_o;
    reco::PFTauDiscriminatorByIsolationRef                      pftdi_r;
    reco::PFTauDiscriminatorByIsolationRefProd                  pftdi_rp;
    reco::PFTauDiscriminatorByIsolationRefVector                pftdi_rv;
    edm::Wrapper<reco::PFTauDiscriminatorByIsolation>           pftdi_w;

    std::pair<reco::PFTauRef, int>                              pftd_p;
    std::vector<std::pair<reco::PFTauRef, int> >                pftd_v;

    reco::PFTauDiscriminatorBase                     pftdiscr_b;
    reco::PFTauDiscriminator                         pftdiscr_o;
    reco::PFTauDiscriminatorRef                      pftdiscr_r;
    reco::PFTauDiscriminatorRefProd                  pftdiscr_rp;
    reco::PFTauDiscriminatorRefVector                pftdiscr_rv;
    edm::Wrapper<reco::PFTauDiscriminator>           pftdiscr_w;

    std::pair<reco::PFTauRef, float>                              pftdiscr_p;
    std::vector<std::pair<reco::PFTauRef, float> >                pftdiscr_v;

    reco::JetPiZeroAssociationBase                     jetPiZeroAssoc_b;
    reco::JetPiZeroAssociation                         jetPiZeroAssoc_o;
    reco::JetPiZeroAssociationRef                      jetPiZeroAssoc_r;
    reco::JetPiZeroAssociationRefProd                  jetPiZeroAssoc_rp;
    reco::JetPiZeroAssociationRefVector                jetPiZeroAssoc_rv;
    edm::Wrapper<reco::JetPiZeroAssociation>           jetPiZeroAssoc_w;

    std::pair<reco::PFJetRef, std::vector<reco::RecoTauPiZero> >                              jetPiZeroAssoc_p;
    std::vector<std::pair<reco::PFJetRef, std::vector<reco::RecoTauPiZero> > >                jetPiZeroAssoc_v;

    std::vector<std::vector<reco::RecoTauPiZero> >                jetPiZeroAssoc_v_v;
    
    reco::PFJetChargedHadronAssociationBase                     jetChHAssoc_b;
    reco::PFJetChargedHadronAssociation                         jetChHAssoc_o;
    reco::PFJetChargedHadronAssociationRef                      jetChHAssoc_r;
    reco::PFJetChargedHadronAssociationRefProd                  jetChHAssoc_rp;
    reco::PFJetChargedHadronAssociationRefVector                jetChHAssoc_rv;
    edm::Wrapper<reco::PFJetChargedHadronAssociation>           jetChHAssoc_w;

    std::pair<reco::PFJetRef, std::vector<reco::PFRecoTauChargedHadron> >                              jetChHAssoc_p;
    std::vector<std::pair<reco::PFJetRef, std::vector<reco::PFRecoTauChargedHadron> > >                jetChHAssoc_v;

    std::vector<std::vector<reco::PFRecoTauChargedHadron> >                jetChHAssoc_v_v;

    reco::PFTauDecayModeAssociation                         pftdecaymodeass_o;
    reco::PFTauDecayModeAssociationRef                      pftdecaymodeass_r;
    reco::PFTauDecayModeAssociationRefProd                  pftdecaymodeass_rp;
    reco::PFTauDecayModeAssociationRefVector                pftdecaymodeass_rv;
    edm::Wrapper<reco::PFTauDecayModeAssociation>           pftdecaymodeass_w;

    std::pair<reco::PFTauRef, reco::PFTauDecayMode>                              pftdecaymodeass_p;
    std::vector<std::pair<reco::PFTauRef, reco::PFTauDecayMode> >                pftdecaymodeass_v;
    std::pair<reco::CaloTauRef, float>                              calodiscr_p;
    std::vector<std::pair<reco::CaloTauRef, float> >                calodiscr_v;

    //Needed only in HLT-Open
    std::vector<reco::HLTTau>                                  ht_v;
    edm::Wrapper<std::vector<reco::HLTTau> >                   ht_w;
    edm::Ref<std::vector<reco::HLTTau> >                       ht_r;
    edm::RefProd<std::vector<reco::HLTTau> >                   ht_rp;
    edm::RefVector<std::vector<reco::HLTTau> >                 ht_rv;

    edm::Ptr<reco::BaseTau>	 ptr_t;
    edm::PtrVector<reco::BaseTau>	 ptrv_t;

    reco::PFTauTransverseImpactParameter                                                        pftautip_o;
    std::vector<reco::PFTauTransverseImpactParameter>                                           pftautip_ov;
    edm::Wrapper<std::vector<reco::PFTauTransverseImpactParameter> >                            pftautip_vw;
    edm::Ref<std::vector<reco::PFTauTransverseImpactParameter> >                                pftautip_vr;
    edm::RefProd<std::vector<reco::PFTauTransverseImpactParameter> >                            pftautip_vrp;
    edm::RefVector<std::vector<reco::PFTauTransverseImpactParameter> >                          pftautip_vrv;
    edm::reftobase::Holder<reco::PFTauTransverseImpactParameter,edm::Ref<std::vector<reco::PFTauTransverseImpactParameter> > >    pftautip_rb;
    edm::Association<std::vector<reco::PFTauTransverseImpactParameter> >                        pftautip_assoc_v;
    edm::Wrapper<edm::Association<std::vector<reco::PFTauTransverseImpactParameter> > >         pftautip_assoc_v_wrapper;

    reco::PFTauTransverseImpactParameterCollection                                                        pftautip_oc;
    std::vector<reco::PFTauTransverseImpactParameterCollection>                                           pftautip_vc;
    edm::Wrapper<std::vector<reco::PFTauTransverseImpactParameterCollection> >                            pftautip_vwc;
    edm::Ref<std::vector<reco::PFTauTransverseImpactParameterCollection> >                                pftautip_vrc;
    edm::RefProd<std::vector<reco::PFTauTransverseImpactParameterCollection> >                            pftautip_vrpc;
    edm::RefVector<std::vector<reco::PFTauTransverseImpactParameterCollection> >                          pftautip_vrvc;
    edm::reftobase::Holder<reco::PFTauTransverseImpactParameterCollection,edm::Ref<std::vector<reco::PFTauTransverseImpactParameterCollection> > >    pftautip_rbc;
    edm::Association<std::vector<reco::PFTauTransverseImpactParameterCollection> >                        pftautip_assoc_vc;
    edm::Wrapper<edm::Association<std::vector<reco::PFTauTransverseImpactParameterCollection> > >         pftautip_assoc_vc_wrapper;

    reco::PFTauTransverseImpactParameterRef                                                        pftautip_or;
    std::vector<reco::PFTauTransverseImpactParameterRef>                                           pftautip_vrb;
    edm::Wrapper<std::vector<reco::PFTauTransverseImpactParameterRef> >                            pftautip_vwr;
    edm::Ref<std::vector<reco::PFTauTransverseImpactParameterRef> >                                pftautip_vrr;
    edm::RefProd<std::vector<reco::PFTauTransverseImpactParameterRef> >                            pftautip_vrpr;
    edm::RefVector<std::vector<reco::PFTauTransverseImpactParameterRef> >                          pftautip_vrvr;
    edm::reftobase::Holder<reco::PFTauTransverseImpactParameterRef,edm::Ref<std::vector<reco::PFTauTransverseImpactParameterRef> > >    pftautip_rbr;
    edm::Association<std::vector<reco::PFTauTransverseImpactParameterRef> >                        pftautip_assoc_vr;
    edm::Wrapper<edm::Association<std::vector<reco::PFTauTransverseImpactParameterRef> > >         pftautip_assoc_vr_wrapper;

    std::vector<reco::VertexRef>                                                          pftauvertex_o;
    edm::Wrapper<std::vector<reco::VertexRef> >                                           pftauvertex_w;
    edm::Ref<std::vector<reco::VertexRef> >                                               pftauvertex_r;
    edm::RefProd<std::vector<reco::VertexRef> >                                           pftauvertex_rp;
    edm::RefVector<std::vector<reco::VertexRef> >                                         pftauvertex_rv;
    edm::reftobase::Holder<reco::VertexRef,edm::Ref<std::vector<reco::VertexRef> > >      pftauvertex_rb;
    edm::Association<std::vector<reco::VertexRef> >                                       pftauvertex_assoc_v;
    edm::Wrapper<edm::Association<std::vector<reco::VertexRef> > >                        pftauvertex_assoc_v_wrapper;

    std::vector<std::vector<reco::VertexRef> >                                                                      pftauvertexv_v;
    edm::Wrapper<std::vector<std::vector<reco::VertexRef> > >                                                       pftauvertexv_w;
    edm::Ref<std::vector<std::vector<reco::VertexRef> > >                                                           pftauvertexv_r;
    edm::RefProd<std::vector<std::vector<reco::VertexRef> > >                                                       pftauvertexv_rp;
    edm::RefVector<std::vector<std::vector<reco::VertexRef> > >                                                     pftauvertexv_rv;
    edm::reftobase::Holder<std::vector<reco::VertexRef>, edm::Ref<std::vector<std::vector<reco::VertexRef> > > >    pftauvertexv_rb;
    edm::Association<std::vector<std::vector<reco::VertexRef> > >                                                   pftauvertexv_assoc_v;
    edm::Wrapper<edm::Association<std::vector<std::vector<reco::VertexRef> > > >                                    pftauvertexv_assoc_v_wrapper;

    reco::PFTauTIPAssociation                         pftautipass_o;
    reco::PFTauTIPAssociationRef                      pftautipass_r;
    reco::PFTauTIPAssociationRefProd                  pftautipass_rp;
    reco::PFTauTIPAssociationRefVector                pftautipass_rv;
    edm::Wrapper<reco::PFTauTIPAssociation>           pftautipass_w;

    std::pair<reco::PFTauRef, std::vector<reco::PFTauTransverseImpactParameterRef> >                pftaupairtip_o;
    std::vector<std::pair<reco::PFTauRef, std::vector<reco::PFTauTransverseImpactParameterRef> > >  pftaupairtip_v;

    reco::PFTauVertexAssociation                      pftauvertexass_o;
    reco::PFTauVertexAssociationRef                   pftauvertexass_r;
    reco::PFTauVertexAssociationRefProd               pftauvertexass_rp;
    reco::PFTauVertexAssociationRefVector             pftauvertexass_rv;
    edm::Wrapper<reco::PFTauVertexAssociation>        pftauvertexass_w;
    std::pair<reco::PFTauRef, std::vector<reco::VertexRef> >                 pftaupairvertex_o;
    std::vector<std::pair<reco::PFTauRef, std::vector<reco::VertexRef> > >   pftaupairvertex_v;

    reco::PFTauVertexVAssociation                     pftauvertexvass_o;
    reco::PFTauVertexVAssociationRef                  pftauvertexvass_r;
    reco::PFTauVertexVAssociationRefProd              pftauvertexvass_rp;
    reco::PFTauVertexVAssociationRefVector            pftauvertexvass_rv;
    edm::Wrapper<reco::PFTauVertexVAssociation>       pftauvertexvass_w;
    std::pair<reco::PFTauRef, std::vector<std::vector<reco::VertexRef> > >                pftaupairvertexv_o;
    std::vector<std::pair<reco::PFTauRef, std::vector<std::vector<reco::VertexRef> > > >  pftaupairvertexv_v;

    reco::PFTau3ProngSummary                                                        pftau3prong_o;
    std::vector<reco::PFTau3ProngSummary>                                           pftau3prong_v;
    edm::Wrapper<std::vector<reco::PFTau3ProngSummary> >                            pftau3prong_w;
    edm::Ref<std::vector<reco::PFTau3ProngSummary> >                                pftau3prong_r;
    edm::RefProd<std::vector<reco::PFTau3ProngSummary> >                            pftau3prong_rp;
    edm::RefVector<std::vector<reco::PFTau3ProngSummary> >                          pftau3prong_rv;
    edm::reftobase::Holder<reco::PFTau3ProngSummary,edm::Ref<std::vector<reco::PFTau3ProngSummary> > >    pftau3prong_rb;
    edm::Association<std::vector<reco::PFTau3ProngSummary> >                        pftau3prong_assoc_v;
    edm::Wrapper<edm::Association<std::vector<reco::PFTau3ProngSummary> > >         pftau3prong_assoc_v_wrapper;

    reco::PFTau3ProngSummaryCollection                                                        pftau3prong_oc;
    std::vector<reco::PFTau3ProngSummaryCollection>                                           pftau3prong_vc;
    edm::Wrapper<std::vector<reco::PFTau3ProngSummaryCollection> >                            pftau3prong_wc;
    edm::Ref<std::vector<reco::PFTau3ProngSummaryCollection> >                                pftau3prong_rc;
    edm::RefProd<std::vector<reco::PFTau3ProngSummaryCollection> >                            pftau3prong_rpc;
    edm::RefVector<std::vector<reco::PFTau3ProngSummaryCollection> >                          pftau3prong_rvc;
    edm::reftobase::Holder<reco::PFTau3ProngSummaryCollection,edm::Ref<std::vector<reco::PFTau3ProngSummaryCollection> > >    pftau3prong_rbc;
    edm::Association<std::vector<reco::PFTau3ProngSummaryCollection> >                        pftau3prong_assoc_vc;
    edm::Wrapper<edm::Association<std::vector<reco::PFTau3ProngSummaryCollection> > >         pftau3prong_assoc_vc_wrapper;

    reco::PFTau3ProngSummaryRef                                                        pftau3prong_or;
    std::vector<reco::PFTau3ProngSummaryRef>                                           pftau3prong_vr;
    edm::Wrapper<std::vector<reco::PFTau3ProngSummaryRef> >                            pftau3prong_wr;
    edm::Ref<std::vector<reco::PFTau3ProngSummaryRef> >                                pftau3prong_rr;
    edm::RefProd<std::vector<reco::PFTau3ProngSummaryRef> >                            pftau3prong_rpr;
    edm::RefVector<std::vector<reco::PFTau3ProngSummaryRef> >                          pftau3prong_rvr;
    edm::reftobase::Holder<reco::PFTau3ProngSummaryRef,edm::Ref<std::vector<reco::PFTau3ProngSummaryRef> > >    pftau3prong_rbr;
    edm::Association<std::vector<reco::PFTau3ProngSummaryRef> >                        pftau3prong_assoc_vr;
    edm::Wrapper<edm::Association<std::vector<reco::PFTau3ProngSummaryRef> > >         pftau3prong_assoc_vr_wrapper;

    reco::PFTau3ProngSumAssociation                         pftau3prongass_o;
    reco::PFTau3ProngSumAssociationRef                      pftau3prongass_r;
    reco::PFTau3ProngSumAssociationRefProd                  pftau3prongass_rp;
    reco::PFTau3ProngSumAssociationRefVector                pftau3prongass_rv;
    edm::Wrapper<reco::PFTau3ProngSumAssociation>           pftau3prongass_w;

    std::pair<reco::PFTauRef, std::vector<reco::PFTau3ProngSummaryRef> >                pftaupair3prong_o;
    std::vector<std::pair<reco::PFTauRef, std::vector<reco::PFTau3ProngSummaryRef> > >  pftaupair3prong_v;

    std::vector<std::vector<TLorentzVector> > vlv_o;

    // jet -> PFCandidate associations, needed for boosted tau reconstruction
    edm::AssociationMap<edm::OneToMany<std::vector<reco::PFJet>,std::vector<reco::PFCandidate>,unsigned int> > jetPFCandidateAssociation_o;
    edm::Wrapper<edm::AssociationMap<edm::OneToMany<std::vector<reco::PFJet>,std::vector<reco::PFCandidate>,unsigned int> > > jetPFCandidateAssociation_w;
    edm::helpers::KeyVal<edm::RefProd<std::vector<reco::PFJet> >,edm::RefProd<std::vector<reco::PFCandidate> > > jetPFCandidateAssociation_kv;
    edm::helpers::KeyVal<edm::Ref<std::vector<reco::PFJet>,reco::PFJet,edm::refhelper::FindUsingAdvance<std::vector<reco::PFJet>,reco::PFJet> >,edm::RefVector<std::vector<reco::PFCandidate>,reco::PFCandidate,edm::refhelper::FindUsingAdvance<std::vector<reco::PFCandidate>,reco::PFCandidate> > > jetPFCandidateAssociation_kv2;
    std::map<unsigned int,edm::helpers::KeyVal<edm::Ref<std::vector<reco::PFJet>,reco::PFJet,edm::refhelper::FindUsingAdvance<std::vector<reco::PFJet>,reco::PFJet> >,edm::RefVector<std::vector<reco::PFCandidate>,reco::PFCandidate,edm::refhelper::FindUsingAdvance<std::vector<reco::PFCandidate>,reco::PFCandidate> > > > jetPFCandidateAssociation_mkv;
  };
}
