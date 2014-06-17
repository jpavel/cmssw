import FWCore.ParameterSet.Config as cms

tauPFIsoValueCharged05 = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("tauIsoDepositPFChargedHadrons"),
             deltaR = cms.double(0.5),
             weight = cms.string('1'),
             vetos = cms.vstring(),
             skipDefaultVeto = cms.bool(True),
             mode = cms.string('sum'),
             PivotCoordinatesForEBEE = cms.bool(False)
             )
         )
     )


tauPFIsoValueGamma05 = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("tauIsoDepositPFGammas"),
             deltaR = cms.double(0.5),
             weight = cms.string('1'),
             vetos = cms.vstring(),
             skipDefaultVeto = cms.bool(True),
             mode = cms.string('sum'),
             PivotCoordinatesForEBEE = cms.bool(False)
             )
         )
     )

tauPFIsoValueNeutral05 = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("tauIsoDepositPFNeutralHadrons"),
             deltaR = cms.double(0.5),
             weight = cms.string('1'),
             vetos = cms.vstring(),
             skipDefaultVeto = cms.bool(True),
             mode = cms.string('sum'),
             PivotCoordinatesForEBEE = cms.bool(False)
             )
         )
     )

tauPFIsoValuePU05 = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(
            cms.PSet(
            src = cms.InputTag("tauIsoDepositPFCandidates"),
             deltaR = cms.double(0.5),
             weight = cms.string('1'),
             vetos = cms.vstring(),
             skipDefaultVeto = cms.bool(True),
             mode = cms.string('sum'),
             PivotCoordinatesForEBEE = cms.bool(False)
             )
         )
     )


tauPFIsolationValuesSequence = cms.Sequence (
    tauPFIsoValueCharged05+
    tauPFIsoValueGamma05+
    tauPFIsoValueNeutral05+
    tauPFIsoValuePU05
)
