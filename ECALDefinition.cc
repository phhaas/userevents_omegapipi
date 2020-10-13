#include "ECALDefinition.h"


ECALDefinition::ECALDefinition(TTree& tree)
	: _phast(Phast::Ref()),
	  _ECAL1               (PaSetup::Ref().Calorimeter(0)),
	  _ECAL2               (PaSetup::Ref().Calorimeter(1)),
	  _neutralClusterNumber(0)
{
	_phast.h_file->cd();

	// ECAL 1

	// ECAL 2

	// Both ECALs
	tree.Branch("ECAL_neutralClusterNumber",      &_neutralClusterNumber, "neutralClusterNumber/i");
	tree.Branch("ECAL_neutralClusterIndex",       &_neutralClusterIndex);
	tree.Branch("ECAL_neutralClusterX",           &_neutralClusterX);
	tree.Branch("ECAL_neutralClusterY",           &_neutralClusterY);
	tree.Branch("ECAL_neutralClusterZ",           &_neutralClusterZ);
	tree.Branch("ECAL_neutralClusterXError",      &_neutralClusterXError);
	tree.Branch("ECAL_neutralClusterYError",      &_neutralClusterYError);
	tree.Branch("ECAL_neutralClusterZError",      &_neutralClusterZError);
	tree.Branch("ECAL_neutralClusterEnergy",      &_neutralClusterEnergy);
	tree.Branch("ECAL_neutralClusterEnergyError", &_neutralClusterEnergyError);
	tree.Branch("ECAL_neutralClusterTime",        &_neutralClusterTime);
	tree.Branch("ECAL_neutralClusterSize",        &_neutralClusterSize);
	tree.Branch("ECAL_IndexCell",                 &_neutralClusterIndexCell);
	tree.Branch("ECAL_clusterXInCell",            &_neutralClusterXInCell);
	tree.Branch("ECAL_clusterYInCell",            &_neutralClusterYInCell);
	tree.Branch("ECAL_isolation",                 &_neutralClusterIsolation);
}


void
ECALDefinition::fill(const PaEvent& event)
{
	getIsolation(event);

	// Both ECALs
	_neutralClusterNumber = getNumberNeutralClusters(event);
	_neutralClusterIndex.clear();
	_neutralClusterX.clear();
	_neutralClusterY.clear();
	_neutralClusterZ.clear();
	_neutralClusterXError.clear();
	_neutralClusterYError.clear();
	_neutralClusterZError.clear();
	_neutralClusterEnergy.clear();
	_neutralClusterEnergyError.clear();
	_neutralClusterTime.clear();
	_neutralClusterSize.clear();
	_neutralClusterIndexCell.clear();
	_neutralClusterXInCell.clear();
	_neutralClusterYInCell.clear();
	_neutralClusterIsolation.clear();
	_neutralClusterIndex.reserve      (_neutralClusterNumber);
	_neutralClusterX.reserve          (_neutralClusterNumber);
	_neutralClusterY.reserve          (_neutralClusterNumber);
	_neutralClusterZ.reserve          (_neutralClusterNumber);
	_neutralClusterXError.reserve     (_neutralClusterNumber);
	_neutralClusterYError.reserve     (_neutralClusterNumber);
	_neutralClusterZError.reserve     (_neutralClusterNumber);
	_neutralClusterEnergy.reserve     (_neutralClusterNumber);
	_neutralClusterEnergyError.reserve(_neutralClusterNumber);
	_neutralClusterTime.reserve       (_neutralClusterNumber);
	_neutralClusterSize.reserve       (_neutralClusterNumber);
	_neutralClusterIndexCell.reserve  (_neutralClusterNumber);
	_neutralClusterXInCell.reserve    (_neutralClusterNumber);
	_neutralClusterYInCell.reserve    (_neutralClusterNumber);
	_neutralClusterIsolation.reserve    (_neutralClusterNumber);
	
	for (size_t i = 0; i < _vectorNeutrals.size(); ++i) {
		const PaParticle& particle = event.vParticle(_vectorNeutrals[i]);
		const PaCaloClus& cluster  = event.vCaloClus(particle.iCalorim(0));
		double cellX;
		double cellY;
		double dX;
		double dY;
		double isolation = 50000;
		if (_ECAL1.IsMyCluster(cluster)) {
			// ECAL 1
			_neutralClusterIndex.push_back(1);
			_neutralClusterIndexCell.push_back(_ECAL1.iCell(cluster.X(), cluster.Y(), cellX, cellY, 0));
			dX = cluster.X() - cellX;
			dY = cluster.Y() - cellY;
			_neutralClusterXInCell.push_back(dX);
			_neutralClusterYInCell.push_back(dY);
			// calculate isolation for ECAL1
			for (size_t j = 0; j < _chargedCoordECAL1.size(); ++j) {
				double tmpDist = std::sqrt((_chargedCoordECAL1[j].X - cluster.X())*(_chargedCoordECAL1[j].X - cluster.X()) 
					+ (_chargedCoordECAL1[j].Y - cluster.Y())*(_chargedCoordECAL1[j].Y - cluster.Y()));
				if (tmpDist < isolation) isolation = tmpDist;
			}
		} else if (_ECAL2.IsMyCluster(cluster)) {
			// ECAL 2
			_neutralClusterIndex.push_back(2);
			_neutralClusterIndexCell.push_back(_ECAL2.iCell(cluster.X(), cluster.Y(), cellX, cellY, 0));
			dX = cluster.X() - cellX;
			dY = cluster.Y() - cellY;
			_neutralClusterXInCell.push_back(dX);
			_neutralClusterYInCell.push_back(dY);
			// calculate isolation for ECAL2
			for (size_t j = 0; j < _chargedCoordECAL2.size(); ++j) {
				double tmpDist = std::sqrt((_chargedCoordECAL2[j].X - cluster.X())*(_chargedCoordECAL2[j].X - cluster.X()) 
					+ (_chargedCoordECAL2[j].Y - cluster.Y())*(_chargedCoordECAL2[j].Y - cluster.Y()));
				if (tmpDist < isolation) isolation = tmpDist;
			}
		}
		// Both ECALs
		_neutralClusterX.push_back          (cluster.X()     );
		_neutralClusterY.push_back          (cluster.Y()     );
		_neutralClusterZ.push_back          (cluster.Z()     );
		_neutralClusterEnergy.push_back     (cluster.E()     );
		_neutralClusterXError.push_back     (cluster.Cov()[0]);
		_neutralClusterYError.push_back     (cluster.Cov()[2]);
		_neutralClusterZError.push_back     (cluster.Cov()[5]);
		_neutralClusterEnergyError.push_back(cluster.Eerr()  );
		_neutralClusterTime.push_back       (cluster.Time()  );
		_neutralClusterSize.push_back       (cluster.Size()  );
		_neutralClusterIsolation.push_back  ( isolation);
		
	}  // loop over _vectorNeutrals
}

const std::vector<int>&
ECALDefinition::getVectorNeutrals(const PaEvent& event)
{
	_vectorNeutrals.clear();
	_vectorNeutrals.reserve((unsigned long)event.NCaloClus());
	for (int i = 0; i < event.NParticle(); ++i) {
		PaParticle particle = event.vParticle(i);
		// Neutral particle
		if (particle.Q() != 0) {
			continue;
		}

		// No track associated to neutral particle
		if (particle.iTrack() != -1) {
			continue;
		}

		// Only one cluster allowed for photons
		if (particle.NCalorim() != 1) {
			continue;
		}

		// First calorimeter cluster -> photons only have one!
		const PaCaloClus& cluster = event.vCaloClus(particle.iCalorim(0));
		// Ensure "E"CAL cluster
		if (cluster.CalorimName()[0] != 'E') {
			continue;
		}

		_vectorNeutrals.push_back(i);
	}

	return _vectorNeutrals;
}

// function to calculate distance to the next charged track
void
ECALDefinition::getIsolation(const PaEvent& event)
{
	const PaVertex& vertex = event.vVertex(event.iBestPrimaryVertex());
	_chargedCoordECAL1.clear();
	_chargedCoordECAL1.reserve(vertex.NOutParticles());
	_chargedCoordECAL2.clear();
	_chargedCoordECAL2.reserve(vertex.NOutParticles());
	// iterate over all charged particles with an track and vertex at the primary vertex to find closest hit
	for (int i = 0; i < vertex.NOutParticles(); ++i) {
		PaParticle particle = event.vParticle(vertex.iOutParticle(i));
		// charged particle
		if (particle.Q() == 0) {
			continue;
		}
		// particle has track
		if (particle.iTrack() == -1) {
			continue;
		}
		const PaTrack& track = event.vTrack(particle.iTrack());
		PaTPar trackPar;
		//check if tracks can be extrapolated to same z as ECAL1
		if (track.Extrapolate(1385, trackPar, false)) 
		{
			_chargedCoordECAL1.push_back({trackPar.X(), trackPar.Y()});
		}
		// same for ECAL2 
		if (track.Extrapolate(3315, trackPar, false)) 
		{
			_chargedCoordECAL2.push_back({trackPar.X(), trackPar.Y()});
		}
	}
}