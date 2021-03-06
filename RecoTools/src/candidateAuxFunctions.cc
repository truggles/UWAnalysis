#include "UWAnalysis/RecoTools/interface/candidateAuxFunctions.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TMath.h"
#include "TLorentzVector.h"

const reco::GenParticle* findGenParticle(const reco::Candidate::LorentzVector& direction,
					 const reco::GenParticleCollection& genParticles, double dRmax, int status,
					 const std::vector<int>* pdgIds, bool pdgIdStrict)
{
  bool bestMatchMatchesPdgId = false;
  const reco::GenParticle* bestMatch = 0;
  
  float lowDiff = 9999.0;
  
  for ( reco::GenParticleCollection::const_iterator genParticle = genParticles.begin();
	genParticle != genParticles.end(); ++genParticle ) {
    bool matchesPdgId = false;
    if ( pdgIds ) {
      for ( std::vector<int>::const_iterator pdgId = pdgIds->begin(); pdgId != pdgIds->end(); ++pdgId ) {
	if ( genParticle->pdgId() == (*pdgId) ) {
	  matchesPdgId = true;
	  break;
        }
      }
    }
    
   
    if ( abs(genParticle->pdgId()) == 12 || abs(genParticle->pdgId()) == 14 || abs(genParticle->pdgId()) == 16 || abs(genParticle->pdgId()) == 22 || abs(genParticle->pdgId()) == 23 ) continue;
    
    // If we require strict PDG id checking, skip it if it doesn't match
    if ( pdgIds && !matchesPdgId && pdgIdStrict ) continue;
    
    // Check if status matches - if not, skip it.
    bool statusMatch = (status == -1 || genParticle->status() == status);
    if ( !statusMatch ) continue;

    double dR = reco::deltaR(direction, genParticle->p4());
    if ( dR > dRmax ) continue;

	float diff = fabs( direction.pt() - genParticle->pt() );

// 	std::cout << "Attempted genMatch: PDGID = " << genParticle->pdgId() << " , status = " << genParticle->status() << " , pT = " << 
// 	genParticle->pt() << " , eta = " << genParticle->eta() << std::endl;
 	
    // Check if old bestMatch was not a prefered ID and the new one is.
    if ( bestMatchMatchesPdgId ) {
      // If the old one matches, only use the new one if it is a better
      // energy match
      if ( matchesPdgId && (diff < lowDiff) ){
        bestMatch = &(*genParticle);
        lowDiff = diff;
      }
    } 
    else {
      // The old one doesn't match the prefferred list if it is either
      // a better energy match or better pdgId match
      if ( (diff < lowDiff) || matchesPdgId ) {
          bestMatch = &(*genParticle);
          lowDiff = diff;

          if ( matchesPdgId ) bestMatchMatchesPdgId = true;
      }
    }
  }
//   if( bestMatch != 0 ){
// 	std::cout << "Best Match: PDGID = " << bestMatch->pdgId() << " , status = " << bestMatch->status() << " , dR = " << 
// 	reco::deltaR(direction, bestMatch->p4()) << " , pT Diff = " << fabs( direction.pt() - bestMatch->pt() ) << std::endl;
//   }
//   else{
//   	std::cout << "No Fucking Match For " << direction.pt() << " , " << direction.eta() << std::endl;
//   }
  return bestMatch;
  
}

float getGenBosonMass(const reco::GenParticleCollection& genParticles) {

    float mass = 0.0;
    std::vector<int> ids;

	ids.push_back(23);
	ids.push_back(22);
	ids.push_back(25);
	ids.push_back(35);
	ids.push_back(36);   

    //get the boson with the highest mass
    if(genParticles.size()>0){
    	for(reco::GenParticleCollection::const_iterator i = genParticles.begin();i!=genParticles.end();++i){ 
      		for(unsigned int j=0;j<ids.size();++j){
				if(abs(i->pdgId())==ids.at(j) && i->mass()>mass){
					for(unsigned int k=0;k<i->numberOfDaughters();++k){
						if(abs(i->daughter(k)->pdgId())==15)
					  		mass = i->mass();
					}
				}
			}
		}
	}
    return mass;
}


reco::Candidate::Point getDecayVertex(const reco::GenParticle* mother)
{
  std::vector<const reco::GenParticle*> stableDaughters;
  findDaughters(mother, stableDaughters, 1);
  // If no daughters, return production vertex
  if ( !stableDaughters.size() ) {
    return mother->vertex();
  } else {
    // Try to sure we aren't going through any intermediate decays
    for ( std::vector<const reco::GenParticle*>::const_iterator dau = stableDaughters.begin();
	  dau != stableDaughters.end(); ++dau ) {
      if ( (*dau)->mother() == mother ) return (*dau)->vertex();
    }
    // If everything had an intermediate decay, just take the first one. :/
    return stableDaughters[0]->vertex();
  }
}

void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters, int status)
{
  unsigned numDaughters = mother->numberOfDaughters();
  for ( unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();

    if ( status == -1 || daughter->status() == status ) daughters.push_back(daughter);

    findDaughters(daughter, daughters, status);
  }
}

bool isNeutrino(const reco::GenParticle* daughter)
{
  return ( TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 );
}

reco::Candidate::LorentzVector getVisMomentum(const std::vector<const reco::GenParticle*>& daughters, int status)
{
  reco::Candidate::LorentzVector p4Vis(0,0,0,0);

  for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
	daughter != daughters.end(); ++daughter ) {
    if ( (status == -1 || (*daughter)->status() == status) && !isNeutrino(*daughter) ) {
      //std::cout << "adding daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","
      //	  << " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
      p4Vis += (*daughter)->p4();
    }
  }

  //std::cout << "--> vis. Momentum: Pt = " << p4Vis.pt() << ", eta = " << p4Vis.eta() << ", phi = " << p4Vis.phi() << std::endl;

  return p4Vis;
}

reco::Candidate::LorentzVector getVisMomentum(const reco::GenParticle* genLeg, const reco::GenParticleCollection* genParticles)
{
  std::vector<const reco::GenParticle*> stableDaughters;
  findDaughters(genLeg, stableDaughters, 1);

  reco::Candidate::LorentzVector p4Vis = getVisMomentum(stableDaughters);

  return p4Vis;
}

reco::Candidate::LorentzVector getInvisMomentum(const std::vector<const reco::GenParticle*>& daughters, int status)
{
  reco::Candidate::LorentzVector p4Invis(0,0,0,0);

  for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
	daughter != daughters.end(); ++daughter ) {
    if ( (status == -1 || (*daughter)->status() == status) && isNeutrino(*daughter) ) {
      //std::cout << "adding daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","
      //	  << " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
      p4Invis += (*daughter)->p4();
    }
  }

  //std::cout << "--> invis. Momentum: Pt = " << p4Invis.pt() << ", phi = " << p4Invis.phi() << std::endl;

  return p4Invis;
}

reco::Candidate::LorentzVector getInvisMomentum(const reco::GenParticle* genLeg, const reco::GenParticleCollection* genParticles)
{
  std::vector<const reco::GenParticle*> stableDaughters;
  findDaughters(genLeg, stableDaughters, 1);

  reco::Candidate::LorentzVector p4Invis = getInvisMomentum(stableDaughters);

  return p4Invis;
}

void compX1X2byCollinearApprox(double& x1, double& x2, double pxLeg1, double pyLeg1, double pxLeg2, double pyLeg2, double pxMEt, double pyMEt)
{
  double x1_numerator = pxLeg1*pyLeg2 - pxLeg2*pyLeg1;
  double x1_denominator = pyLeg2*(pxLeg1 + pxMEt) - pxLeg2*(pyLeg1 + pyMEt);
  x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
  //std::cout << "x1 = " << x1 << std::endl;

  double x2_numerator = x1_numerator;
  double x2_denominator = pxLeg1*(pyLeg2 + pyMEt) - pyLeg1*(pxLeg2 + pxMEt);
  x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
  //std::cout << "x2 = " << x2 << std::endl;
}

double getPhysX(double x, bool& isWithinPhysRange)
{
  double physX = x;

  isWithinPhysRange = true;

  if ( x < 0. ) {
    physX = 0.;
    isWithinPhysRange = false;
  }

  if ( x > 1. ) {
    physX = 1.;
    isWithinPhysRange = false;
  }

  return physX;
}

reco::Candidate::LorentzVector boostToRestFrame(const reco::Candidate::LorentzVector& p4daughter,
						const reco::Candidate::LorentzVector& p4mother)
{
  TLorentzVector p4mother_tlorentzvector(p4mother.px(), p4mother.py(), p4mother.pz(), p4mother.energy());

  TVector3 boostVector = -(p4mother_tlorentzvector.BoostVector());

  TLorentzVector p4daughter_tlorentzvector(p4daughter.px(), p4daughter.py(), p4daughter.pz(), p4daughter.energy());
  p4daughter_tlorentzvector.Boost(boostVector);

  return reco::Candidate::LorentzVector(p4daughter_tlorentzvector.Px(), p4daughter_tlorentzvector.Py(), p4daughter_tlorentzvector.Pz(),
					p4daughter_tlorentzvector.Energy());
}

//
//-------------------------------------------------------------------------------
//

TVector2 getDiTauBisectorDirection(const reco::Candidate::LorentzVector& leg1P4, const reco::Candidate::LorentzVector& leg2P4)
{
  double leg1CosPhi = TMath::Cos(leg1P4.phi());
  double leg1SinPhi = TMath::Sin(leg1P4.phi());

  double leg2CosPhi = TMath::Cos(leg2P4.phi());
  double leg2SinPhi = TMath::Sin(leg2P4.phi());

  double diTauPx_unnormalized = leg1CosPhi + leg2CosPhi;
  double diTauPy_unnormalized = leg1SinPhi + leg2SinPhi;
  double diTauP = TMath::Sqrt(diTauPx_unnormalized*diTauPx_unnormalized + diTauPy_unnormalized*diTauPy_unnormalized);
  double diTauPx_normalized = ( diTauP > 0. ) ? (diTauPx_unnormalized/diTauP) : diTauPx_unnormalized;
  double diTauPy_normalized = ( diTauP > 0. ) ? (diTauPy_unnormalized/diTauP) : diTauPy_unnormalized;

  return TVector2(diTauPx_normalized, diTauPy_normalized);
}

void computeMEtProjection(const reco::PFCandidateCollection& pfCandidates, const TVector2& direction,
			  double& sumEt, double& sumP_par, double& sumP_perp)
{
  sumEt = 0.;

  sumP_par = 0.;
  sumP_perp = 0.;

  double dirCosPhi = direction.X();
  double dirSinPhi = direction.Y();

  for ( reco::PFCandidateCollection::const_iterator pfCandidate = pfCandidates.begin();
	pfCandidate != pfCandidates.end(); ++pfCandidate ) {
//--- skip particle candidate of unknown type,
//    as it is done in the "official" (pf)MET reconstruction code
//    implemented in RecoMET/METAlgorithms/src/PFSpecificAlgo.cc
    if ( !(pfCandidate->particleId() == reco::PFCandidate::h        ||
	   pfCandidate->particleId() == reco::PFCandidate::e        ||
	   pfCandidate->particleId() == reco::PFCandidate::mu       ||
	   pfCandidate->particleId() == reco::PFCandidate::gamma    ||
	   pfCandidate->particleId() == reco::PFCandidate::h0       ||
	   pfCandidate->particleId() == reco::PFCandidate::h_HF     ||
	   pfCandidate->particleId() == reco::PFCandidate::egamma_HF) ) continue;

    sumEt += pfCandidate->et();

    sumP_par += TMath::Abs(pfCandidate->px()*dirCosPhi + pfCandidate->py()*dirSinPhi);
    sumP_perp += TMath::Abs(pfCandidate->px()*dirSinPhi - pfCandidate->py()*dirCosPhi);
  }
}

//
//-------------------------------------------------------------------------------
//

std::string getTauDecayModeName(int tauDecayMode)
{
  // "translate" from enum defined in DataFormats/TauReco/interface/PFTauDecayMode.h
  // to string ( in format defined in PhysicsTools/JetMCUtils/src/JetMCTag.cc )
  //
  // NOTE: the "unphysical" 2-prong modes take into account
  //       track reconstruction inefficiencies (migration from 3-prong decays),
  //       fake tracks and tracks from the underlying event (migration from 1-prong decays)
  //
  if      ( tauDecayMode == reco::PFTauDecayMode::tauDecaysElectron           ) return std::string("electron");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecayMuon                ) return std::string("muon");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion0PiZero ) return std::string("oneProng0Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero ) return std::string("oneProng1Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion2PiZero ) return std::string("oneProng2Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion3PiZero ) return std::string("oneProng3Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion4PiZero ) return std::string("oneProng4Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion0PiZero ) return std::string("twoProng0Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion1PiZero ) return std::string("twoProng1Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion2PiZero ) return std::string("twoProng2Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion3PiZero ) return std::string("twoProng3Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion4PiZero ) return std::string("twoProng4Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero ) return std::string("threeProng0Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion1PiZero ) return std::string("threeProng1Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion2PiZero ) return std::string("threeProng2Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion3PiZero ) return std::string("threeProng3Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion4PiZero ) return std::string("threeProng4Pi0");
  else if ( tauDecayMode == reco::PFTauDecayMode::tauDecayOther               ) return std::string("other");
  else {
    edm::LogError ("getTauDecayModeName")
      << " Invalid tau decay mode = " << tauDecayMode << " !!";
    return std::string("unknown");
  }
}

//
//-------------------------------------------------------------------------------
//

const reco::Candidate* getDistPion(const pat::Tau& recTauJet)
{
  if ( !recTauJet.isPFTau() ) {
    //edm::LogWarning ("getDistPion")
    //  << " Cannot identify 'distinguishable' pion for CaloTaus/TCTaus --> returning NULL pointer !!";
    return 0;
  }

  //const reco::PFCandidateRefVector& recTauJetChargedConstituents = recTauJet.signalPFChargedHadrCands();
  const std::vector<edm::Ptr<reco::PFCandidate> >& recTauJetChargedConstituents = recTauJet.signalPFChargedHadrCands();

  if ( recTauJetChargedConstituents.size() == 1 ) {

//--- tau- --> one-prong case (in particular rho- --> pi- pi0 or tau- --> a1- --> pi- pi0 pi0);
//    the "distinguishable" pion is the leading charged hadron
    return recTauJet.leadPFChargedHadrCand().get();
  } else if ( recTauJetChargedConstituents.size() == 3 ) {
    double recTauJetCharge = recTauJet.charge();
    /*
    for ( reco::PFCandidateRefVector::const_iterator recTauJetChargedConstituent = recTauJetChargedConstituents.begin();
         recTauJetChargedConstituent != recTauJetChargedConstituents.end(); ++recTauJetChargedConstituent ) {
//--- tau- --> three-prong case (in particular a1- --> pi- pi+ pi-);
//    the "distinguishable" pion is the pion of charge opposite to the tau-jet charge
      if ( (*recTauJetChargedConstituent)->charge()*recTauJetCharge < 0 ) {
        return recTauJetChargedConstituent->get();
      }
      }*/
  } else {
    //edm::LogWarning ("getDistPion")
    //  << " Unsupported rec. tau decay mode = " << recTauJet.decayMode() << " --> returning NULL pointer !!";
    return 0;
  }
    
  //edm::LogWarning ("getDistPion")
  //  << " Failed to identify 'distinguishable' rec. pion --> returning NULL pointer !!";
  return 0;
}

const reco::Candidate* getDistPion(const reco::GenJet& genTauJet)
{
  std::string genTauDecayMode = JetMCTagUtils::genTauDecayMode(genTauJet);

  std::vector<const reco::GenParticle*> genTauJetConstituents = genTauJet.getGenConstituents();

  if ( genTauDecayMode == "oneProng1Pi0" ||
       genTauDecayMode == "oneProng2Pi0" ) {
    for ( std::vector<const reco::GenParticle*>::const_iterator genTauJetConstituent = genTauJetConstituents.begin();
	  genTauJetConstituent != genTauJetConstituents.end(); ++genTauJetConstituent ) {

//--- tau- --> rho- --> pi- pi0 or tau- --> a1- --> pi- pi0 pi0 case;
//    the "distinguishable" pion is the charged one
      if ((*genTauJetConstituent)->charge()) return (*genTauJetConstituent);
    }
  } else if ( genTauDecayMode == "threeProng0Pi0" ) {
    double genTauJetCharge = genTauJet.charge();

    for ( std::vector<const reco::GenParticle*>::const_iterator genTauJetConstituent = genTauJetConstituents.begin();
	  genTauJetConstituent != genTauJetConstituents.end(); ++genTauJetConstituent ) {

//--- tau- --> a1- --> pi- pi+ pi- case
//    the "distinguishable" pion is the pion of charge opposite to the tau-jet charge
      if ( (*genTauJetConstituent)->charge() != 0. &&
	   ((*genTauJetConstituent)->charge()*genTauJetCharge) < 0. ) {
	return *genTauJetConstituent;
      }
    }
  } else {
    //edm::LogWarning ("getDistPion")
    //  << " Unsupported gen. tau decay mode = " << genTauDecayMode << " --> returning NULL pointer !!";
    return 0;
  }
  //edm::LogWarning ("getDistPion")
  //  << " Failed to identify 'distinguishable' gen. pion --> returning NULL pointer !!";
  return 0;
}

//
//-------------------------------------------------------------------------------
//

std::pair<double, double> compMEtProjU(const reco::Candidate::LorentzVector& zP4, double metPx, double metPy, int& errorFlag)
{
  if ( zP4.pt() == 0. ) {
    edm::LogWarning ("compMEtProjU")
      << " Failed to compute projection, because Z0 candidate has zero Pt --> returning dummy solution !!";
    errorFlag = 1;
    return std::pair<double, double>(0., 0.);
  }
  
  double qX = zP4.px();
  double qY = zP4.py();
  double qT = TMath::Sqrt(qX*qX + qY*qY);
  
  double uX = -(qX + metPx);
  double uY = -(qY + metPy);
  
  double u1 = (uX*qX + uY*qY)/qT;
  double u2 = (uX*qY - uY*qX)/qT;
  
  return std::pair<double, double>(u1, u2);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

std::vector<double> compTrackPtSums(const reco::VertexCollection& vertices)
{
  size_t numVertices = vertices.size();
  std::vector<double> trackPtSums(numVertices);
  for ( size_t iVertex = 0; iVertex < numVertices; ++iVertex ) {
    const reco::Vertex& vertex = vertices[iVertex];

    double trackPtSum = 0.;
    for ( reco::Vertex::trackRef_iterator track = vertex.tracks_begin();
	  track != vertex.tracks_end(); ++track ) {
      trackPtSum += (*track)->pt();
    }
    
    trackPtSums[iVertex] = trackPtSum;
  }

  return trackPtSums;
}

size_t getNumVerticesPtGtThreshold(const std::vector<double>& trackPtSums, double ptThreshold)
{
  size_t numVertices = 0;
  for ( std::vector<double>::const_iterator trackPtSum = trackPtSums.begin();
	trackPtSum != trackPtSums.end(); ++trackPtSum ) {
    if ( (*trackPtSum) > ptThreshold ) ++numVertices;
  }
  
  return numVertices;
}


