import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'

#added in etau and mutau triggers
from UWAnalysis.Configuration.tools.analysisToolsHTauTauSync import *
defaultReconstructionMC(process,'HLT',
                      [
			'HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1',#etau
			'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2'#muTau
                      ])

                      

#EventSelection
process.load("UWAnalysis.Configuration.hTauTauSync_cff")

process.metCalibration.applyCalibration = cms.bool(False)

process.eventSelectionET = cms.Path(process.selectionSequenceET)


createGeneratedParticles(process,
                         'genDaughters',
                          [
                           "keep++ pdgId = {Z0}",
                           "keep pdgId = {tau+}",
                           "keep pdgId = {tau-}",
                           "keep pdgId = {mu+}",
                           "keep pdgId = {mu-}",
                           "keep pdgId = 6",
                           "keep pdgId = -6",
                           "keep pdgId = 11",
                           "keep pdgId = -11",
                          ]
)


createGeneratedParticles(process,
                         'genTauCands',
                          [
                           "keep pdgId = {tau+} & mother.pdgId()= {Z0}",
                           "keep pdgId = {tau-} & mother.pdgId() = {Z0}"
                          ]
)



from UWAnalysis.Configuration.tools.ntupleToolsSync import addEleTauEventTree
addEleTauEventTree(process,'eleTauEventTree')

addEventSummary(process,True,'ET','eventSelectionET')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       $inputFileNames
		),
		inputCommands=cms.untracked.vstring(
						'keep *',
		)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.TFileService.fileName=cms.string("$outputFileName")
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("$outputFileName")
)
