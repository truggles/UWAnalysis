import FWCore.ParameterSet.Config as cms


#Import tool that creates the cut sequence
from UWAnalysis.Configuration.tools.CutSequenceProducer import *

###############################			Ele-Tau 		###################################
ETanalysisConfigurator = CutSequenceProducer(initialCounter  = 'initialEventsET',
                                  pyModuleName = __name__,
                                  pyNameSpace  = locals())

#ESTausID
ETanalysisConfigurator.addSmearing('patOverloadedTaus','triggeredPatMuons','triggeredPatElectrons','patOverloadedJets','slimmedMETs','ET')

#create dielectrons
ETanalysisConfigurator.addDiCandidateModule('diElectrons','PATElePairProducer','smearedElectronsET','smearedElectronsET','smearedMETET','','smearedJetsET',0,9999,text = '',leadingObjectsOnly = False,dR = 0.15,recoMode = "")
ETanalysisConfigurator.addSelector('osDiElectrons','PATElePairSelector','abs(leg1.eta())<2.5&&abs(leg2.eta())<2.5&&abs(leg1.userFloat("dXY"))<0.045&&abs(leg2.userFloat("dXY"))<0.045&&abs(leg1.userFloat("dZ"))<0.2&&abs(leg2.userFloat("dZ"))<0.2&&leg1.pt()>15&&leg2.pt()>15&&charge==0&&leg2.userFloat("CBIDVeto")>0&&leg1.userFloat("CBIDVeto")>0&&leg1.userFloat("dBRelIso03")<0.3&&leg2.userFloat("dBRelIso03")<0.3','DiElectronCreation',0,100)

#Make DiTaus
ETanalysisConfigurator.addDiCandidateModule('eleTaus','PATEleTauPairProducer','smearedElectronsET','smearedTausET','smearedMETET','smearedTausET','smearedJetsET',1,9999,text = 'AtLeastOneEleTau',leadingObjectsOnly = False,dR = 0.5,recoMode = "",genParticles='genDaughters')

ETanalysisConfigurator.addSelector('eleTausElePtEta','PATEleTauPairSelector','leg1.pt()>23&&abs(leg1.eta())<2.1','electronPtEta',1)
ETanalysisConfigurator.addSelector('eleTausEleID','PATEleTauPairSelector','leg1.userFloat("eleMVAIDnonTrig80")>0','ElectronID',1)
ETanalysisConfigurator.addSelector('eleTausEleVertices','PATEleTauPairSelector','abs(leg1.userFloat("dZ"))<0.2&&abs(leg1.userFloat("dXY"))<0.045','electronVertices',1)
ETanalysisConfigurator.addSelector('eleTausEleConvRej','PATEleTauPairSelector','leg1.userInt("eleConversion")==0','electronConvRej',1)
ETanalysisConfigurator.addSelector('eleTausTauPtEta','PATEleTauPairSelector','leg2.pt()>20&&abs(leg2.eta())<2.3','ETTauPtEta',1)
ETanalysisConfigurator.addSelector('eleTausDecayFound','PATEleTauPairSelector','abs(leg2.userFloat("taudZ"))<0.2&&leg2.tauID("decayModeFindingNewDMs")>0.5','ETTauDecayFound',1)
ETanalysisConfigurator.addSelector('eleTausChargeNot2','PATEleTauPairSelector','abs(leg2.charge())==1','ETChargeIsABS1',1)
#ETanalysisConfigurator.addEleTauSVFitSA('eleTausSVFit')
ETanalysisConfigurator.addSelector('eleTausEleTrigMatch','PATEleTauPairSelector','(leg1.pt()>33&&leg1.userFloat("hltEle32WP75GsfTrackIsoFilter")>0)||(leg2.userFloat("hltPFTau20TrackLooseIso")>0&&leg1.userFloat("hltEle22WP75L1IsoEG20erTau20erGsfTrackIsoFilter")>0)  ','ETEleTrigMatch',1)
ETanalysisConfigurator.addSorter('eleTausSorted','PATEleTauPairSorterByIso')
#ETanalysisConfigurator.addSelector('eleTausOS','PATEleTauPairSelector','charge==0','ETOS',1)

#create the sequence
selectionSequenceET =ETanalysisConfigurator.returnSequence()

################################		Mu-Tau			###################################

MTanalysisConfigurator = CutSequenceProducer(initialCounter  = 'initialEventsMT',
                                  pyModuleName = __name__,
                                  pyNameSpace  = locals())

#Add smearing
MTanalysisConfigurator.addSmearing('patOverloadedTaus','triggeredPatMuons','triggeredPatElectrons','patOverloadedJets','slimmedMETs','MT')

#Create di muon pairs for veto purposes
MTanalysisConfigurator.addDiCandidateModule('diMuons','PATMuPairProducer','smearedMuonsMT','smearedMuonsMT','smearedMETMT','','smearedJetsMT',0,9999,text = '',leadingObjectsOnly = False,dR = 0.15,recoMode = "")
#veto possible second zpeak
MTanalysisConfigurator.addSelector('osDiMuons','PATMuPairSelector','leg1.isPFMuon&&leg2.isPFMuon&&abs(leg1.eta())<2.4&&abs(leg2.eta())<2.4&&abs(leg1.userFloat("dZ"))<0.2&&abs(leg2.userFloat("dZ"))<0.2&&abs(leg2.userFloat("dXY"))<0.045&&abs(leg2.userFloat("dXY"))<0.045&&charge==0&&leg1.isTrackerMuon&&leg2.isTrackerMuon&&leg1.isGlobalMuon&&leg2.isGlobalMuon&&leg1.pt()>15&&leg2.pt()>15&&leg1.userFloat("dBRelIso03")<0.3 &&leg2.userFloat("dBRelIso03")<0.3','DiMuonCreation',0,100)

#Make DiTaus   
MTanalysisConfigurator.addDiCandidateModule('diTaus','PATMuTauPairProducer','smearedMuonsMT','smearedTausMT','smearedMETMT','smearedTausMT','smearedJetsMT',1,9999,text='AtLeastOneDiTau',leadingObjectsOnly = False,dR = 0.5,recoMode ="",genParticles='genDaughters')
MTanalysisConfigurator.addSelector('diTausMuonVertices','PATMuTauPairSelector','abs(leg1.userFloat("dZ"))<0.2&&abs(leg1.userFloat("dXY"))<0.045','MuonVertices',1)
#tightMuonSelectionForNow
MTanalysisConfigurator.addSelector('diTausMuonMediumID','PATMuTauPairSelector','leg1.isMediumMuon()','MuonMediumID',1)
#MTanalysisConfigurator.addSelector('diTausMuonMediumID','PATMuTauPairSelector','leg1.userInt("mediumID")','MuonMediumID',1)
MTanalysisConfigurator.addSelector('diTausMuonPtEta','PATMuTauPairSelector','leg1.pt()>18&&abs(leg1.eta())<2.1','MuonPtEta',1)
MTanalysisConfigurator.addSelector('diTausTauPtEta','PATMuTauPairSelector','leg2.pt()>20&&abs(leg2.eta())<2.3','MTTauPtEta',1)
MTanalysisConfigurator.addSelector('diTausDecayFound','PATMuTauPairSelector','abs(leg2.userFloat("taudZ"))<0.2&&leg2.tauID("decayModeFindingNewDMs")>0.5','MTTauDecayFound',1)
MTanalysisConfigurator.addSelector('diTausChargeNot2','PATMuTauPairSelector','abs(leg2.charge())==1','MTChargeIsABS1',1)
#MTanalysisConfigurator.addMuTauSVFitSA('diTausSVFit')
MTanalysisConfigurator.addSelector('diTausTausMuTrigMatch','PATMuTauPairSelector','(leg1.userFloat("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09")>0&&leg2.userFloat("hltPFTau20TrackLooseIsoAgainstMuon")>0)||(leg1.userFloat("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09")>0&&leg1.pt()>25)','MTMuTrigMatch',1)
#MTanalysisConfigurator.addSelector('diTausTausMuTrigMatch','PATMuTauPairSelector','(leg1.userFloat("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09")>0)||(leg1.userFloat("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09")>0&&leg1.pt()>25)','MTMuTrigMatch',1)
MTanalysisConfigurator.addSorter('diTausSorted','PATMuTauPairSorterByIso')

#create the sequence
selectionSequenceMT =MTanalysisConfigurator.returnSequence()