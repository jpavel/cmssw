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
