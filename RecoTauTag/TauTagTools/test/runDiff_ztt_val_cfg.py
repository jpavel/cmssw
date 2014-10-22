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


discNameLoose3='MVA5LooseElectronRejection'
discriminatorLoose3="hpsPFTauDiscriminationBy"+discNameLoose3

discNameLoose2='LooseElectronRejection'
discriminatorLoose2="hpsPFTauDiscriminationBy"+discNameLoose2

discNameLoose='LooseCombinedIsolationDBSumPtCorr3Hits'
discriminatorLoose="hpsPFTauDiscriminationBy"+discNameLoose

#medium
discNameMedium3='MVA5MediumElectronRejection'
discriminatorMedium3="hpsPFTauDiscriminationBy"+discNameMedium3

discNameMedium2='MediumElectronRejection'
discriminatorMedium2="hpsPFTauDiscriminationBy"+discNameMedium2

discNameMedium='MediumCombinedIsolationDBSumPtCorr3Hits'
discriminatorMedium="hpsPFTauDiscriminationBy"+discNameMedium

#tight
discNameTight3='MVA5TightElectronRejection'
discriminatorTight3="hpsPFTauDiscriminationBy"+discNameTight3

discNameTight2='TightElectronRejection'
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
'file:ref.root'
#        'root://xrootd.unl.edu//store/mc/Spring14dr/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00B3D5C6-C2CB-E311-9026-002481E0D480.root'
#        '/store/relval/CMSSW_7_1_0_pre4_AK4/RelValZTT_13/GEN-SIM-RECO/POSTLS171_V1-v2/00000/2491D051-51B5-E311-BFA1-003048679168.root',
#       '/store/relval/CMSSW_7_1_0_pre4_AK4/RelValZTT_13/GEN-SIM-RECO/POSTLS171_V1-v2/00000/BEC396FA-4CB5-E311-B68C-002354EF3BDD.root',
#       '/store/relval/CMSSW_7_1_0_pre4_AK4/RelValZTT_13/GEN-SIM-RECO/POSTLS171_V1-v2/00000/E8DA9320-52B5-E311-AB3B-0025905A6088.root'
# '/store/relval/CMSSW_7_1_0_pre6/RelValZTT_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/7AC27528-32C7-E311-95C9-003048FFD770.root',
# '/store/relval/CMSSW_7_1_0_pre6/RelValZTT_13/GEN-SIM-RECO/PRE_LS171_V3-v1/00000/C42C13E7-30C7-E311-B583-003048678DD6.root'
# '/store/relval/CMSSW_7_1_0_pre5/RelValZTT_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/8E8B2309-E6B6-E311-A8C6-002618943896.root',
# '/store/relval/CMSSW_7_1_0_pre5/RelValZTT_13/GEN-SIM-RECO/POSTLS171_V1-v1/00000/E84A5F0C-E3B6-E311-8515-003048D3C010.root'
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

tauIDSourcesRECO = cms.PSet(
        # configure many IDs as InputTag <someName> = <someTag>                                                                                                                                                                       
        # you can comment out those you don't want to save some                                                                                                                                                                       
        # disk space                                                                                                                                                                                                                  
        decayModeFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding","","RECO"),
        decayModeFindingNewDMs =cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs","","RECO"),
        chargedIsoPtSum = cms.InputTag("hpsPFTauMVA3IsolationChargedIsoPtSum","","RECO"),
        neutralIsoPtSum = cms.InputTag("hpsPFTauMVA3IsolationNeutralIsoPtSum","","RECO"),
        puCorrPtSum = cms.InputTag("hpsPFTauMVA3IsolationPUcorrPtSum","","RECO"),
        byIsolationMVA3oldDMwoLTraw = cms.InputTag('hpsPFTauDiscriminationByIsolationMVA3oldDMwoLTraw',"","RECO"),
        byVLooseIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwoLT',"","RECO"),
        byLooseIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwoLT',"","RECO"),
        byMediumIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwoLT',"","RECO"),
        byTightIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByTightIsolationMVA3oldDMwoLT',"","RECO"),
        byVTightIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwoLT',"","RECO"),
        byVVTightIsolationMVA3oldDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwoLT',"","RECO"),                     
        byIsolationMVA3oldDMwLTraw = cms.InputTag('hpsPFTauDiscriminationByIsolationMVA3oldDMwLTraw',"","RECO"),
        byVLooseIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwLT',"","RECO"),
        byLooseIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwLT',"","RECO"),
        byMediumIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwLT',"","RECO"),
        byTightIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByTightIsolationMVA3oldDMwLT',"","RECO"),
        byVTightIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwLT',"","RECO"),
        byVVTightIsolationMVA3oldDMwLT = cms.InputTag('hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwLT',"","RECO"),                             
        byIsolationMVA3newDMwoLTraw = cms.InputTag('hpsPFTauDiscriminationByIsolationMVA3newDMwoLTraw',"","RECO"),
        byVLooseIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwoLT',"","RECO"),
        byLooseIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByLooseIsolationMVA3newDMwoLT',"","RECO"),
        byMediumIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByMediumIsolationMVA3newDMwoLT',"","RECO"),
        byTightIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByTightIsolationMVA3newDMwoLT',"","RECO"),
        byVTightIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVTightIsolationMVA3newDMwoLT',"","RECO"),
        byVVTightIsolationMVA3newDMwoLT = cms.InputTag('hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwoLT',"","RECO"),                             
        byIsolationMVA3newDMwLTraw = cms.InputTag('hpsPFTauDiscriminationByIsolationMVA3newDMwLTraw',"","RECO"),
        byVLooseIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwLT',"","RECO"),
        byLooseIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByLooseIsolationMVA3newDMwLT',"","RECO"),
        byMediumIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByMediumIsolationMVA3newDMwLT',"","RECO"),
        byTightIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByTightIsolationMVA3newDMwLT',"","RECO"),
        byVTightIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByVTightIsolationMVA3newDMwLT',"","RECO"),
        byVVTightIsolationMVA3newDMwLT = cms.InputTag('hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwLT',"","RECO"),                             
        againstElectronLoose = cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection","","RECO"),
        againstElectronMedium = cms.InputTag("hpsPFTauDiscriminationByMediumElectronRejection","","RECO"),
        againstElectronTight = cms.InputTag("hpsPFTauDiscriminationByTightElectronRejection","","RECO"),
        againstMuonLoose = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection","","RECO"),
        againstMuonMedium = cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection","","RECO"),
        againstMuonTight = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection","","RECO"),
        againstMuonLoose2 = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection2","","RECO"),
        againstMuonMedium2 = cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection2","","RECO"),
        againstMuonTight2 = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection2","","RECO"),
        againstMuonLoose3 = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3","","RECO"),
        againstMuonTight3 = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3","","RECO"),
        againstMuonMVAraw = cms.InputTag('hpsPFTauDiscriminationByMVArawMuonRejection',"","RECO"),                                                            
        againstMuonLooseMVA = cms.InputTag('hpsPFTauDiscriminationByMVALooseMuonRejection',"","RECO"),
        againstMuonMediumMVA = cms.InputTag('hpsPFTauDiscriminationByMVAMediumMuonRejection',"","RECO"),
        againstMuonTightMVA = cms.InputTag('hpsPFTauDiscriminationByMVATightMuonRejection',"","RECO"),           
        byCombinedIsolationDeltaBetaCorrRaw3Hits = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits","","RECO"),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits","","RECO"),
        byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits","","RECO"),
        byTightCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits","","RECO"),
        againstElectronMVA5raw = cms.InputTag("hpsPFTauDiscriminationByMVA5rawElectronRejection","","RECO"),
        againstElectronMVA5category = cms.InputTag("hpsPFTauDiscriminationByMVA5rawElectronRejection:category","","RECO"),
        againstElectronVLooseMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5VLooseElectronRejection","","RECO"),
        againstElectronLooseMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5LooseElectronRejection","","RECO"),
        againstElectronMediumMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5MediumElectronRejection","","RECO"),
        againstElectronTightMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5TightElectronRejection","","RECO"),
        againstElectronVTightMVA5 = cms.InputTag("hpsPFTauDiscriminationByMVA5VTightElectronRejection","","RECO")
)

 


process.tauDifferenceAnalyzer = cms.EDFilter("RecoTauValHist",
                                             src1 = cms.InputTag("hpsPFTauProducer","",procName),
                                             src2 = cms.InputTag("hpsPFTauProducer","","RECO"),
                                             discDM = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs","",procName),
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

process.tauDifferenceAnalyzerRECO = process.tauDifferenceAnalyzer.clone(
    src1=cms.InputTag("hpsPFTauProducer","","RECO"),
    discDM = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding","","RECO"),
    # disc1 = cms.InputTag(discriminator, "", "RECO"),
    #  discLoose = cms.InputTag(discriminatorLoose, "", "RECO"),
    #  discLoose_2 = cms.InputTag(discriminatorLoose2, "", "RECO"),
    #  discLoose_3 = cms.InputTag(discriminatorLoose3, "", "RECO"),
    #  discMedium = cms.InputTag(discriminatorMedium, "", "RECO"),
    #  discMedium_2 = cms.InputTag(discriminatorMedium2, "", "RECO"),
    #  discMedium_3 = cms.InputTag(discriminatorMedium3, "", "RECO"),
    #  discTight = cms.InputTag(discriminatorTight, "", "RECO"),
    #  discTight_2 = cms.InputTag(discriminatorTight2, "", "RECO"),
    #  discTight_3 = cms.InputTag(discriminatorTight3, "", "RECO"),
    #  chIso1 = cms.InputTag("hpsPFTauMVA3IsolationChargedIsoPtSum","","RECO"),
    #  nIso1 = cms.InputTag("hpsPFTauMVA3IsolationNeutralIsoPtSum","","RECO"),
    #  PUIso1 = cms.InputTag("hpsPFTauMVA3IsolationPUcorrPtSum","","RECO"),
    #  cmbIso1 = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits","","RECO"),
                                           

    idSrc = tauIDSourcesRECO
)

from RecoTauTag.Configuration.switchMVAtoDB_cfi import switchMVAtoDB
process = switchMVAtoDB(process)

process.hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwoLT.verbosity = cms.int32(1)
process.hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwoLT.mapping[0].cut = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff90')
process.hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwoLT.loadMVAfromDB = cms.bool(True)
process.hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwoLT.mvaOutput_normalization = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_mvaOutput_normalization')

process.hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwoLT.verbosity = cms.int32(1)
# process.hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwoLT.mapping[0].cut = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff50')
# process.hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwoLT.loadMVAfromDB = cms.bool(True)
# process.hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwoLT.mvaOutput_normalization = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_mvaOutput_normalization')

## let it run

process.p = cms.Path(
        process.tauGenJets*
        process.PFTau*
        process.tauDifferenceAnalyzer*
        process.tauDifferenceAnalyzerRECO
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
