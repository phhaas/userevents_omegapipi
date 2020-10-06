#include "ECALDefinition.h"


ECALDefinition::ECALDefinition(TTree& tree)
	: _phast(Phast::Ref()),
	  _ECAL1               (PaSetup::Ref().Calorimeter(0)),
	  _ECAL2               (PaSetup::Ref().Calorimeter(1)),
	  _ECAL1ClusterNumber  (0),
	  _ECAL2ClusterNumber  (0),
	  _neutralClusterNumber(0)
{
	_phast.h_file->cd();

	// ECAL 1
	tree.Branch("ECAL1_clusterNumber",      &_ECAL1ClusterNumber, "ECAL1ClusterNumber/i");
	tree.Branch("ECAL1_clusterX",           &_ECAL1ClusterX);
	tree.Branch("ECAL1_clusterY",           &_ECAL1ClusterY);
	tree.Branch("ECAL1_clusterZ",           &_ECAL1ClusterZ);
	tree.Branch("ECAL1_clusterXError",      &_ECAL1ClusterXError);
	tree.Branch("ECAL1_clusterYError",      &_ECAL1ClusterYError);
	tree.Branch("ECAL1_clusterZError",      &_ECAL1ClusterZError);
	tree.Branch("ECAL1_clusterEnergy",      &_ECAL1ClusterEnergy);
	tree.Branch("ECAL1_clusterEnergyError", &_ECAL1ClusterEnergyError);
	tree.Branch("ECAL1_clusterTime",        &_ECAL1ClusterTime);
	tree.Branch("ECAL1_clusterSize",        &_ECAL1ClusterSize);
	// ECAL 2
	tree.Branch("ECAL2_clusterNumber",      &_ECAL2ClusterNumber, "ECAL2ClusterNumber/i");
	tree.Branch("ECAL2_clusterX",           &_ECAL2ClusterX);
	tree.Branch("ECAL2_clusterY",           &_ECAL2ClusterY);
	tree.Branch("ECAL2_clusterZ",           &_ECAL2ClusterZ);
	tree.Branch("ECAL2_clusterXError",      &_ECAL2ClusterXError);
	tree.Branch("ECAL2_clusterYError",      &_ECAL2ClusterYError);
	tree.Branch("ECAL2_clusterZError",      &_ECAL2ClusterZError);
	tree.Branch("ECAL2_clusterEnergy",      &_ECAL2ClusterEnergy);
	tree.Branch("ECAL2_clusterEnergyError", &_ECAL2ClusterEnergyError);
	tree.Branch("ECAL2_clusterTime",        &_ECAL2ClusterTime);
	tree.Branch("ECAL2_clusterSize",        &_ECAL2ClusterSize);
	tree.Branch("ECAL2_IndexCell",          &_ECAL2IndexCell);
	tree.Branch("ECAL2_clusterXInCell",     &_ECAL2ClusterXInCell);
	tree.Branch("ECAL2_clusterYInCell",     &_ECAL2ClusterYInCell);
	// Both ECALs
	tree.Branch("ECAL_neutralClusterNumber",      &_neutralClusterNumber, "neutralClusterNumber/i");
	tree.Branch("ECAL_neutralClusterX",           &_neutralClusterX);
	tree.Branch("ECAL_neutralClusterY",           &_neutralClusterY);
	tree.Branch("ECAL_neutralClusterZ",           &_neutralClusterZ);
	tree.Branch("ECAL_neutralClusterXError",      &_neutralClusterXError);
	tree.Branch("ECAL_neutralClusterYError",      &_neutralClusterYError);
	tree.Branch("ECAL_neutralClusterZError",      &_neutralClusterZError);
	tree.Branch("ECAL_neutralClusterEnergy",      &_neutralClusterEnergy);
	tree.Branch("ECAL_neutralClusterEnergyError", &_neutralClusterEnergyError);
	tree.Branch("ECAL_neutralClusterTime",        &_neutralClusterTime);
}


void
ECALDefinition::fill(const PaEvent& event)
{
	// ECAL 1
	_ECAL1ClusterNumber = getECAL1ClusterNumber(event);
	_ECAL1ClusterX.clear();
	_ECAL1ClusterY.clear();
	_ECAL1ClusterZ.clear();
	_ECAL1ClusterXError.clear();
	_ECAL1ClusterYError.clear();
	_ECAL1ClusterZError.clear();
	_ECAL1ClusterEnergy.clear();
	_ECAL1ClusterEnergyError.clear();
	_ECAL1ClusterTime.clear();
	_ECAL1ClusterSize.clear();
	_ECAL1ClusterX.reserve          (_ECAL1ClusterNumber);
	_ECAL1ClusterY.reserve          (_ECAL1ClusterNumber);
	_ECAL1ClusterZ.reserve          (_ECAL1ClusterNumber);
	_ECAL1ClusterXError.reserve     (_ECAL1ClusterNumber);
	_ECAL1ClusterYError.reserve     (_ECAL1ClusterNumber);
	_ECAL1ClusterZError.reserve     (_ECAL1ClusterNumber);
	_ECAL1ClusterEnergy.reserve     (_ECAL1ClusterNumber);
	_ECAL1ClusterEnergyError.reserve(_ECAL1ClusterNumber);
	_ECAL1ClusterTime.reserve       (_ECAL1ClusterNumber);
	_ECAL1ClusterSize.reserve       (_ECAL1ClusterNumber);
	// ECAL 2
	_ECAL2ClusterNumber = getECAL2ClusterNumber(event);
	_ECAL2ClusterX.clear();
	_ECAL2ClusterY.clear();
	_ECAL2ClusterZ.clear();
	_ECAL2ClusterXError.clear();
	_ECAL2ClusterYError.clear();
	_ECAL2ClusterZError.clear();
	_ECAL2ClusterEnergy.clear();
	_ECAL2ClusterEnergyError.clear();
	_ECAL2ClusterTime.clear();
	_ECAL2ClusterSize.clear();
	_ECAL2IndexCell.clear();
	_ECAL2ClusterXInCell.clear();
	_ECAL2ClusterYInCell.clear();
	_ECAL2ClusterX.reserve          (_ECAL2ClusterNumber);
	_ECAL2ClusterY.reserve          (_ECAL2ClusterNumber);
	_ECAL2ClusterZ.reserve          (_ECAL2ClusterNumber);
	_ECAL2ClusterXError.reserve     (_ECAL2ClusterNumber);
	_ECAL2ClusterYError.reserve     (_ECAL2ClusterNumber);
	_ECAL2ClusterZError.reserve     (_ECAL2ClusterNumber);
	_ECAL2ClusterEnergy.reserve     (_ECAL2ClusterNumber);
	_ECAL2ClusterEnergyError.reserve(_ECAL2ClusterNumber);
	_ECAL2ClusterTime.reserve       (_ECAL2ClusterNumber);
	_ECAL2ClusterSize.reserve       (_ECAL2ClusterNumber);
	_ECAL2IndexCell.reserve         (_ECAL2ClusterNumber);
	_ECAL2ClusterXInCell.reserve    (_ECAL2ClusterNumber);
	_ECAL2ClusterYInCell.reserve    (_ECAL2ClusterNumber);
	// Both ECALs
	_neutralClusterNumber = getNumberNeutralClusters(event);
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

	for (size_t i = 0; i < _vectorNeutrals.size(); ++i) {
		const PaParticle& particle = event.vParticle(_vectorNeutrals[i]);
		const PaCaloClus& cluster  = event.vCaloClus(particle.iCalorim(0));
		if (_ECAL1.IsMyCluster(cluster)) {
			// ECAL 1
			_ECAL1ClusterX.push_back          (cluster.X()     );
			_ECAL1ClusterY.push_back          (cluster.Y()     );
			_ECAL1ClusterZ.push_back          (cluster.Z()     );
			_ECAL1ClusterXError.push_back     (cluster.Cov()[0]);
			_ECAL1ClusterYError.push_back     (cluster.Cov()[2]);
			_ECAL1ClusterZError.push_back     (cluster.Cov()[5]);
			_ECAL1ClusterEnergy.push_back     (cluster.E()     );
			_ECAL1ClusterEnergyError.push_back(cluster.Eerr()  );
			_ECAL1ClusterTime.push_back       (cluster.Time()  );
			_ECAL1ClusterSize.push_back       (cluster.Size()  );
		} else if (_ECAL2.IsMyCluster(cluster)) {
			// ECAL 2
			_ECAL2ClusterX.push_back          (cluster.X()     );
			_ECAL2ClusterY.push_back          (cluster.Y()     );
			_ECAL2ClusterZ.push_back          (cluster.Z()     );
			_ECAL2ClusterXError.push_back     (cluster.Cov()[0]);
			_ECAL2ClusterYError.push_back     (cluster.Cov()[2]);
			_ECAL2ClusterZError.push_back     (cluster.Cov()[5]);
			_ECAL2ClusterEnergy.push_back     (cluster.E()     );
			_ECAL2ClusterEnergyError.push_back(cluster.Eerr()  );
			_ECAL2ClusterTime.push_back       (cluster.Time()  );
			_ECAL2ClusterSize.push_back       (cluster.Size()  );
			double cellX;
			double cellY;
			_ECAL2IndexCell.push_back(_ECAL2.iCell(cluster.X(), cluster.Y(), cellX, cellY, 0));
			double dX = cluster.X() - cellX;
			double dY = cluster.Y() - cellY;
			_ECAL2ClusterXInCell.push_back(dX);
			_ECAL2ClusterYInCell.push_back(dY);
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


unsigned long
ECALDefinition::getECAL1ClusterNumber(const PaEvent& event)
{
	_ECAL1ClusterNumber = 0;
	for (size_t i = 0; i < _vectorNeutrals.size(); ++i) {
		const PaParticle& particle = event.vParticle(_vectorNeutrals[i]);
		const PaCaloClus& cluster  = event.vCaloClus(particle.iCalorim(0));
		if (_ECAL1.IsMyCluster(cluster)) {
			_ECAL1ClusterNumber++;
		}
	}
	return _ECAL1ClusterNumber;
}


unsigned long
ECALDefinition::getECAL2ClusterNumber(const PaEvent& event)
{
	_ECAL2ClusterNumber = 0;
	for (size_t i = 0; i < _vectorNeutrals.size(); ++i) {
		const PaParticle& particle = event.vParticle(_vectorNeutrals[i]);
		const PaCaloClus& cluster  = event.vCaloClus(particle.iCalorim(0));
		if (_ECAL2.IsMyCluster(cluster)) {
			_ECAL2ClusterNumber++;
		}
	}
	return _ECAL2ClusterNumber;
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
