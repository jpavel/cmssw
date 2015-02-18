),
#'file:/afs/cern.ch/work/j/jez/ntuples/tauID/relval/CMSSW_7_1_0_pre4_AK4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v2/00000/08F63CD1-30B5-E311-B53E-003048B835A2.root'),
#'/store/relval/CMSSW_7_1_0_pre4/RelValZTT_13/GEN-SIM-RECO/POSTLS171_V1-v2/00000/D4AF85E0-85AA-E311-95E5-02163E00E782.root'),
    skipEvents = cms.untracked.uint32(sEvents)
)
## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(eMax) )

#### included from standard pat template
## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.load("Configuration.StandardSequences.MagneticField_cff")



process.options.allowUnscheduled = cms.untracked.bool(False)
#process.Tracer = cms.Service("Tracer")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.load("RecoJets.JetProducers.ak5GenJets_cfi")
process.load("RecoJets.Configuration.GenJetParticles_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *

switchToPFTauHPS(process)
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts

tauIDSources = cms.PSet(
        # configure many IDs as InputTag <someName> = <someTag>
        # you can comment out those you don't want to save some
        # disk space
        decayModeFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding","",procName),
        decayModeFindingNewDMs =cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs","",procName),
        chargedIsoPtSum = cms.InputTag("hpsPFTauMVA3IsolationChargedIsoPtSum","",procName),
        neutralIsoPtSum = cms.InputTag("hpsPFTauMVA3IsolationNeutralIsoPtSum","",procName),
        puCorrPtSum = cms.InputTag("hpsPFTauMVA3IsolationPUcorrPtSum","",procName),
        byIsolationMVA3oldDMwoLTraw = cms.InputTag('hpsPFTauDiscriminationByIsolationMVA3oldDMwoLTraw',"",procName),
        byVLooseIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwoLT',"",procName),
        byLooseIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwoLT',"",procName),
        byMediumIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwoLT',"",procName),
        byTightIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByTightIsolationMVA3oldDMwoLT',"",procName),
        byVTightIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwoLT',"",procName),
        byVVTightIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwoLT',"",procName),                     
        byIsolationMVA3oldDMwLTraw = cms.InputTag('hpsPFTauDiscriminationByIsolationMVA3oldDMwLTraw',"",procName),
        byVLooseIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwLT',"",procName),
        byLooseIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwLT',"",procName),
        byMediumIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwLT',"",procName),
        byTightIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByTightIsolationMVA3oldDMwLT',"",procName),
        byVTightIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwLT',"",procName),
        byVVTightIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwLT',"",procName),                             
        byIsolationMVA3newDMwoLTraw = cms.InputTag('hpsPFTauDiscriminationByIsolationMVA3newDMwoLTraw',"",procName),
        byVLooseIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwoLT',"",procName),
        byLooseIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByLooseIsolationMVA3newDMwoLT',"",procName),
        byMediumIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByMediumIsolationMVA3newDMwoLT',"",procName),
        byTightIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByTightIsolationMVA3newDMwoLT',"",procName),
        byVTightIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVTightIsolationMVA3newDMwoLT',"",procName),
        byVVTightIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwoLT',"",procName),                             
        byIsolationMVA3newDMwLTraw = cms.InputTag('hpsPFTauDiscriminationByIsolationMVA3newDMwLTraw',"",procName),
        byVLooseIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwLT',"",procName),
        byLooseIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByLooseIsolationMVA3newDMwLT',"",procName),
        byMediumIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByMediumIsolationMVA3newDMwLT',"",procName),
        byTightIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByTightIsolationMVA3newDMwLT',"",procName),
        byVTightIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByVTightIsolationMVA3newDMwLT',"",procName),
        byVVTightIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwLT',"",procName),                             
        againstElectronLoose = cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection","",procName),
        againstElectronMedium = cms.InputTag("hpsPFTauDiscriminationByMediumElectronRejection","",procName),
        againstElectronTight = cms.InputTag("hpsPFTauDiscriminationByTightElectronRejection","",procName),
        againstMuonLoose = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection","",procName),
        againstMuonMedium = cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection","",procName),
        againstMuonTight = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection","",procName),
        againstMuonLoose2 = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection2","",procName),
        againstMuonMedium2 = cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection2","",procName),
        againstMuonTight2 = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection2","",procName),
        againstMuonLoose3 = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3","",procName),
        againstMuonTight3 = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3","",procName),
        againstMuonMVAraw = cms.InputTag('hpsPFTauDiscriminationByMVArawMuonRejection',"",procName),                                                            
        againstMuonLooseMVA = cms.InputTag('hpsPFTauDiscriminationByMVALooseMuonRejection',"",procName),
        againstMuonMediumMVA = cms.InputTag('hpsPFTauDiscriminationByMVAMediumMuonRejection',"",procName),
        againstMuonTightMVA = cms.InputTag('hpsPFTauDiscriminationByMVATightMuonRejection',"",procName),           
        byCombinedIsolationDeltaBetaCorrRaw3Hits = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits","",procName),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits","",procName),
        byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits","",procName),
        byTightCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits","",procName),
        againstElectronMVA5raw = cms.InputTag("hpsPFTauDiscriminationByMVA5rawElectronRejection","",procName),
        againstElectronMVA5category = cms.InputTag("hpsPFTauDiscriminationByMVA5rawElectronRejection:category","",procName),
        againstElectronVLooseMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5VLooseElectronRejection","",procName),
        againstElectronLooseMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5LooseElectronRejection","",procName),
        againstElectronMediumMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5MediumElectronRejection","",procName),
        againstElectronTightMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5TightElectronRejection","",procName),
        againstElectronVTightMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5VTightElectronRejection","",procName)
)



process.tauDifferenceAnalyzer = cms.EDFilter("RecoTauValHist",
                                             src1 = cms.InputTag("hpsPFTauProducer","",procName),
                                             src2 = cms.InputTag("hpsPFTauProducer","","HLT"),
                                             disc1 = cms.InputTag(discriminator, "", procName),
                                             disc2 = cms.InputTag(discriminator2, "", "HLT"),
                                             discLoose = cms.InputTag(discriminatorLoose, "", procName),
                                             discLoose_2 = cms.InputTag(discriminatorLoose2, "", procName),
                                             discLoose_3 = cms.InputTag(discriminatorLoose3, "", procName),
                                             discMedium = cms.InputTag(discriminatorMedium, "", procName),
                                             discMedium_2 = cms.InputTag(discriminatorMedium2, "", procName),
                                             discMedium_3 = cms.InputTag(discriminatorMedium3, "", procName),
                                             discTight = cms.InputTag(discriminatorTight, "", procName),
                                             discTight_2 = cms.InputTag(discriminatorTight2, "", procName),
                                             discTight_3 = cms.InputTag(discriminatorTight3, "", procName),
                                             chIso1 = cms.InputTag("hpsPFTauMVA3IsolationChargedIsoPtSum","",procName),
                                             chIso2 = cms.InputTag("hpsPFTauMVA3IsolationChargedIsoPtSum","","HLT"),
                                             nIso1 = cms.InputTag("hpsPFTauMVA3IsolationNeutralIsoPtSum","",procName),
                                             nIso2 = cms.InputTag("hpsPFTauMVA3IsolationNeutralIsoPtSum","","HLT"),
                                             PUIso1 = cms.InputTag("hpsPFTauMVA3IsolationPUcorrPtSum","",procName),
                                             PUIso2 = cms.InputTag("hpsPFTauMVA3IsolationPUcorrPtSum","","HLT"),
                                             cmbIso1 = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits","",procName),
                                             cmbIso2 = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits","","HLT"),
                                             idSrc = tauIDSources,
                                             genSrc = cms.InputTag("genParticles"),
                                             genTauSrc = cms.InputTag("tauGenJets"),
                                             mcMatch = cms.bool(True),
                                             primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
                                             verboseOutput = cms.bool(True),
                                             verboseOutputMC = cms.bool(False),
                                             matchingDistance = cms.double(0.1),
                                             background = cms.bool(False),
                                             rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
                                             requireDecayMode = cms.int32(decayMode),
                                             checkMother = cms.bool(True),
                                             Zmumu=cms.bool(False),
                                             onlyHadronic=cms.bool(True),
                                             qualityCuts=PFTauQualityCuts
                                                                            )
process.tauDifferenceAnalyzer.qualityCuts.isolationQualityCuts.minTrackHits = cms.uint32(3)
process.tauDifferenceAnalyzer.qualityCuts.isolationQualityCuts.minTrackPt = cms.double(0.5)
process.tauDifferenceAnalyzer.qualityCuts.isolationQualityCuts.minGammaEt = cms.double(0.5)


## let it run

process.p = cms.Path(
        process.tauGenJets*
        process.PFTau*
        process.tauDifferenceAnalyzer
#            process.patDefaultSequence
            )


process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root"))

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               ## save only events passing the full path
                               #SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               ## save PAT output; you need a '*' to unpack the list of commands
                               ## 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *', *patEventContentNoCleaning )
                               )

#process.outpath = cms.EndPath(process.out)

process.options.wantSummary = False
