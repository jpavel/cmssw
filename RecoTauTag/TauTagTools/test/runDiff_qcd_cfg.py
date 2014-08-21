import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register ('skipEvents',
                  0,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "events to skip")

#parsing arguments
options.parseArguments()

procName="DMtest"

eMax = 100
if options.maxEvents:
    eMax = options.maxEvents


sEvents=0
if options.skipEvents:
    sEvents=options.skipEvents

discName2='LooseCombinedIsolationDBSumPtCorr3Hits'
discriminator2="hpsPFTauDiscriminationBy"+discName2

discName='LooseCombinedIsolationDBSumPtCorr3Hits'
discriminator="hpsPFTauDiscriminationBy"+discName


discNameLoose3='LooseMuonRejection'
discriminatorLoose3="hpsPFTauDiscriminationBy"+discNameLoose3

discNameLoose2='MVA5LooseElectronRejection'
discriminatorLoose2="hpsPFTauDiscriminationBy"+discNameLoose2

discNameLoose='LooseCombinedIsolationDBSumPtCorr3Hits'
discriminatorLoose="hpsPFTauDiscriminationBy"+discNameLoose

#medium
discNameMedium3='MediumMuonRejection'
discriminatorMedium3="hpsPFTauDiscriminationBy"+discNameMedium3

discNameMedium2='MVA5MediumElectronRejection'
discriminatorMedium2="hpsPFTauDiscriminationBy"+discNameMedium2

discNameMedium='MediumCombinedIsolationDBSumPtCorr3Hits'
discriminatorMedium="hpsPFTauDiscriminationBy"+discNameMedium

#tight
discNameTight3='TightMuonRejection'
discriminatorTight3="hpsPFTauDiscriminationBy"+discNameTight3

discNameTight2='MVA5TightElectronRejection'
discriminatorTight2="hpsPFTauDiscriminationBy"+discNameTight2

discNameTight='TightCombinedIsolationDBSumPtCorr3Hits'
discriminatorTight="hpsPFTauDiscriminationBy"+discNameTight


decayMode = 0

DMname=""
if decayMode == 0:
    DMname="DMall"
elif decayMode == 1:
    DMname="DM1p"
elif decayMode == 2:
    DMname="DM1pX"
elif decayMode == 3:
    DMname="DM3p"

#process definition'
procName=procName+DMname
process = cms.Process(procName)

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# process.MessageLogger = cms.Service("MessageLogger",
#        destinations   = cms.untracked.vstring(
#                                              'detailedInfo'
#                                                ),

#        detailedInfo   = cms.untracked.PSet(
#                       threshold  = cms.untracked.string('DEBUG') 
#        ),
#    debugModules = cms.untracked.vstring(
#                                        'PFTau'
#   )
# )


## Options and Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )


## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         '/store/relval/CMSSW_7_1_0_pre6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/44142FEC-20C7-E311-9576-003048678B38.root',
       '/store/relval/CMSSW_7_1_0_pre6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/44C31F40-23C7-E311-AE6D-002618943857.root',
       '/store/relval/CMSSW_7_1_0_pre6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/6A1AB6F1-1CC7-E311-B9EC-002590596498.root',
       '/store/relval/CMSSW_7_1_0_pre6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/78FD5859-21C7-E311-A8BC-0025905A60D2.root',
       '/store/relval/CMSSW_7_1_0_pre6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/7E3742E4-21C7-E311-BA12-0026189438CC.root',
       '/store/relval/CMSSW_7_1_0_pre6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/A42FAB95-25C7-E311-8058-002618943829.root',
       '/store/relval/CMSSW_7_1_0_pre6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/CE038E3C-20C7-E311-9D33-0030486792AC.root',
       '/store/relval/CMSSW_7_1_0_pre6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/D40C037A-1DC7-E311-A321-0025905A60B8.root',
       '/store/relval/CMSSW_7_1_0_pre6/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/F4624392-21C7-E311-86DD-0026189438AE.root'
       #   '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/065CA478-17B6-E311-83AE-0025905A611C.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/08492532-1BB6-E311-90F7-0025905938AA.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/0ACD0B5C-E3B6-E311-B784-003048679076.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/242FEF51-0CB6-E311-8DD3-0025905938B4.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/247306DB-0EB6-E311-8073-0025905A611C.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/28945BC0-22B6-E311-A3CF-0026189438A0.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/86A79188-2BB6-E311-AC16-003048678B7C.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/8A45DEC0-18B6-E311-B9BB-0025905A60B6.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/B4BA04DD-EAB6-E311-B255-0025905A6136.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/C03E8E1C-15B6-E311-A7FE-0025905A48F2.root',
       # '/store/relval/CMSSW_7_1_0_pre5/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/ECE6D91F-15B6-E311-80B2-0025905A6082.root' 
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

process.tauDifferenceAnalyzer = cms.EDFilter("RecoTauDifferenceAnalyzer",
                                             src1 = cms.InputTag("hpsPFTauProducer","",procName),
                                             src2 = cms.InputTag("hpsPFTauProducer","","RECO"),
                                             disc1 = cms.InputTag(discriminator, "", procName),
                                             disc2 = cms.InputTag(discriminator2, "", "RECO"),
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
                                             chIso2 = cms.InputTag("hpsPFTauMVA3IsolationChargedIsoPtSum","","RECO"),
                                             nIso1 = cms.InputTag("hpsPFTauMVA3IsolationNeutralIsoPtSum","",procName),
                                             nIso2 = cms.InputTag("hpsPFTauMVA3IsolationNeutralIsoPtSum","","RECO"),
                                             PUIso1 = cms.InputTag("hpsPFTauMVA3IsolationPUcorrPtSum","",procName),
                                             PUIso2 = cms.InputTag("hpsPFTauMVA3IsolationPUcorrPtSum","","RECO"),
                                             cmbIso1 = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits","",procName),
                                             cmbIso2 = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits","","RECO"),
                                             genSrc = cms.InputTag("genParticles"),
                                             genTauSrc = cms.InputTag("ak5GenJets"),
                                             mcMatch = cms.bool(True),
                                             primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
                                             verboseOutput = cms.bool(True),
                                             verboseOutputMC = cms.bool(False),
                                             matchingDistance = cms.double(0.5),
                                             background = cms.bool(False),
                                             rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
                                             requireDecayMode = cms.int32(decayMode),
                                             checkMother = cms.bool(False),
                                             Zmumu=cms.bool(False),
                                             onlyHadronic=cms.bool(False),
                                             qualityCuts=PFTauQualityCuts
                                                                            )
process.tauDifferenceAnalyzer.qualityCuts.isolationQualityCuts.minTrackHits = cms.uint32(3)
process.tauDifferenceAnalyzer.qualityCuts.isolationQualityCuts.minTrackPt = cms.double(0.5)
process.tauDifferenceAnalyzer.qualityCuts.isolationQualityCuts.minGammaEt = cms.double(0.5)

## let it run

process.p = cms.Path(
#        process.tauGenJets*
        process.genParticlesForJets*
        process.ak5GenJets*
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

process.outpath = cms.EndPath(process.out)

process.options.wantSummary = False
