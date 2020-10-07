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
}


void
ECALDefinition::fill(const PaEvent& event)
{
	// ECAL 1
	
	// ECAL 2

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

	for (size_t i = 0; i < _vectorNeutrals.size(); ++i) {
		const PaParticle& particle = event.vParticle(_vectorNeutrals[i]);
		const PaCaloClus& cluster  = event.vCaloClus(particle.iCalorim(0));
		double cellX;
		double cellY;
		double dX;
		double dY;
		if (_ECAL1.IsMyCluster(cluster)) {
			// ECAL 1
			_neutralClusterIndex.push_back(1);
			_neutralClusterIndexCell.push_back(_ECAL1.iCell(cluster.X(), cluster.Y(), cellX, cellY, 0));
			dX = cluster.X() - cellX;
			dY = cluster.Y() - cellY;
			_neutralClusterXInCell.push_back(dX);
			_neutralClusterYInCell.push_back(dY);
			
		} else if (_ECAL2.IsMyCluster(cluster)) {
			// ECAL 2
			_neutralClusterIndex.push_back(2);
			_neutralClusterIndexCell.push_back(_ECAL2.iCell(cluster.X(), cluster.Y(), cellX, cellY, 0));
			dX = cluster.X() - cellX;
			dY = cluster.Y() - cellY;
			_neutralClusterXInCell.push_back(dX);
			_neutralClusterYInCell.push_back(dY);
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

		//TODO add information for isolation cut, i.e. closest distance of
		//     any of the charged tracks to ECAL cluster
		_vectorNeutrals.push_back(i);
	}

	return _vectorNeutrals;
}
