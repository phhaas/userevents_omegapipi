#include "ECALDefinition.h"


ECALDefinition::ECALDefinition(TTree& tree)
	: _phast(Phast::Ref()),
	  _ECAL1               (PaSetup::Ref().Calorimeter(0)),
	  _ECAL2               (PaSetup::Ref().Calorimeter(1)),
	  _clusterNumber(0)
{
	_phast.h_file->cd();

	// ECAL 1

	// ECAL 2

	// Both ECALs
	tree.Branch("ECAL_clusterNumber",      &_clusterNumber, "clusterNumber/i");
	tree.Branch("ECAL_clusterIndex",       &_clusterIndex);
	tree.Branch("ECAL_clusterX",           &_clusterX);
	tree.Branch("ECAL_clusterY",           &_clusterY);
	tree.Branch("ECAL_clusterZ",           &_clusterZ);
	tree.Branch("ECAL_clusterXError",      &_clusterXError);
	tree.Branch("ECAL_clusterYError",      &_clusterYError);
	tree.Branch("ECAL_clusterZError",      &_clusterZError);
	tree.Branch("ECAL_clusterEnergy",      &_clusterEnergy);
	tree.Branch("ECAL_clusterEnergyError", &_clusterEnergyError);
	tree.Branch("ECAL_clusterTime",        &_clusterTime);
	tree.Branch("ECAL_clusterSize",        &_clusterSize);
	tree.Branch("ECAL_IndexCell",                 &_clusterIndexCell);
	tree.Branch("ECAL_clusterXInCell",            &_clusterXInCell);
	tree.Branch("ECAL_clusterYInCell",            &_clusterYInCell);
	tree.Branch("ECAL_clusterChargeDistance",     &_clusterChargeDistance);
	tree.Branch("ECAL_closestChargeDistance",     &_closestChargeDistance, "closestChargeDistance/d");
}


void
ECALDefinition::fill(const PaEvent& event)
{
	getChargeDistance(event);

	// Both ECALs
	_clusterNumber = getNumberNeutralClusters(event);
	_clusterIndex.clear();
	_clusterX.clear();
	_clusterY.clear();
	_clusterZ.clear();
	_clusterXError.clear();
	_clusterYError.clear();
	_clusterZError.clear();
	_clusterEnergy.clear();
	_clusterEnergyError.clear();
	_clusterTime.clear();
	_clusterSize.clear();
	_clusterIndexCell.clear();
	_clusterXInCell.clear();
	_clusterYInCell.clear();
	_clusterChargeDistance.clear();
	_clusterIndex.reserve          (_clusterNumber);
	_clusterX.reserve              (_clusterNumber);
	_clusterY.reserve              (_clusterNumber);
	_clusterZ.reserve              (_clusterNumber);
	_clusterXError.reserve         (_clusterNumber);
	_clusterYError.reserve         (_clusterNumber);
	_clusterZError.reserve         (_clusterNumber);
	_clusterEnergy.reserve         (_clusterNumber);
	_clusterEnergyError.reserve    (_clusterNumber);
	_clusterTime.reserve           (_clusterNumber);
	_clusterSize.reserve           (_clusterNumber);
	_clusterIndexCell.reserve      (_clusterNumber);
	_clusterXInCell.reserve        (_clusterNumber);
	_clusterYInCell.reserve        (_clusterNumber);
	_clusterChargeDistance.reserve (_clusterNumber);
	_closestChargeDistance = 50000;

	for (size_t i = 0; i < _vectorNeutrals.size(); ++i) {
		const PaParticle& particle = event.vParticle(_vectorNeutrals[i]);
		const PaCaloClus& cluster  = event.vCaloClus(particle.iCalorim(0));
		double cellX;
		double cellY;
		double dX;
		double dY;
		double chargeDistance = 50000;
		if (_ECAL1.IsMyCluster(cluster)) {
			// ECAL 1
			_clusterIndex.push_back(1);
			_clusterIndexCell.push_back(_ECAL1.iCell(cluster.X(), cluster.Y(), cellX, cellY, 0));
			dX = cluster.X() - cellX;
			dY = cluster.Y() - cellY;
			_clusterXInCell.push_back(dX);
			_clusterYInCell.push_back(dY);
			// calculate closest distance to a charge for ECAL2
			for (size_t j = 0; j < _chargedCoordECAL1.size(); ++j) {
				double tmpDist = std::sqrt((_chargedCoordECAL1[j].X - cluster.X())*(_chargedCoordECAL1[j].X - cluster.X())
					+ (_chargedCoordECAL1[j].Y - cluster.Y())*(_chargedCoordECAL1[j].Y - cluster.Y()));
				if (tmpDist < chargeDistance) chargeDistance = tmpDist;
			}
		} else if (_ECAL2.IsMyCluster(cluster)) {
			// ECAL 2
			_clusterIndex.push_back(2);
			_clusterIndexCell.push_back(_ECAL2.iCell(cluster.X(), cluster.Y(), cellX, cellY, 0));
			dX = cluster.X() - cellX;
			dY = cluster.Y() - cellY;
			_clusterXInCell.push_back(dX);
			_clusterYInCell.push_back(dY);
			// calculate closest distance to a charge for ECAL2
			for (size_t j = 0; j < _chargedCoordECAL2.size(); ++j) {
				double tmpDist = std::sqrt((_chargedCoordECAL2[j].X - cluster.X())*(_chargedCoordECAL2[j].X - cluster.X())
					+ (_chargedCoordECAL2[j].Y - cluster.Y())*(_chargedCoordECAL2[j].Y - cluster.Y()));
				if (tmpDist < chargeDistance) chargeDistance = tmpDist;
			}
		}
		// check if the distance to the next charge is smaller than for previous checked hits
		if (chargeDistance < _closestChargeDistance) _closestChargeDistance = chargeDistance;
		// Both ECALs
		_clusterX.push_back                  (cluster.X()     );
		_clusterY.push_back                  (cluster.Y()     );
		_clusterZ.push_back                  (cluster.Z()     );
		_clusterEnergy.push_back             (cluster.E()     );
		_clusterXError.push_back             (cluster.Cov()[0]);
		_clusterYError.push_back             (cluster.Cov()[2]);
		_clusterZError.push_back             (cluster.Cov()[5]);
		_clusterEnergyError.push_back        (cluster.Eerr()  );
		_clusterTime.push_back               (cluster.Time()  );
		_clusterSize.push_back               (cluster.Size()  );
		_clusterChargeDistance.push_back     (chargeDistance);


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
ECALDefinition::getChargeDistance(const PaEvent& event)
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
		// TODO check Z coord
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
