#include "PhotonDefinition.h"

PhotonDefinition::PhotonDefinition(TTree& tree)
	: _phast(Phast::Ref()),
	  _numberPhotonsMC(0),
	  _momentumMCX(MAX_CLUSTERS),
	  _momentumMCY(MAX_CLUSTERS),
	  _momentumMCZ(MAX_CLUSTERS)
{
	_phast.h_file->cd();

	tree.Branch("Photon_numberPhotonsMC", &_numberPhotonsMC, "numberPhotonsMC/I");
	tree.Branch("Photon_momentumMCX",     &_momentumMCX);
	tree.Branch("Photon_momentumMCY",     &_momentumMCY);
	tree.Branch("Photon_momentumMCZ",     &_momentumMCZ);
}


void
PhotonDefinition::fillMC(const PaEvent&    event,
                         const PaMCvertex& vertex)
{
	fillMCPairsOmegaPi0(event, vertex);
}


void
PhotonDefinition::fillMCPairsOmegaPi0(const PaEvent&    event,
                                      const PaMCvertex& vertex)
{

	std::vector<std::pair<int, int> > decayVerticesAndTracks;
	for (int i = 0; i < vertex.NMCtrack(); ++i) {
		const PaMCtrack& track = event.vMCtrack(vertex.iMCtrack(i));
		if (not track.IsBeam()
		    and (track.Pid() == MC_PID_PI0 or track.Pid() == MC_PID_OMEGA)
		    and track.vMCvertex().size() == 1) {
			decayVerticesAndTracks.push_back(std::pair<int, int>(track.vMCvertex()[0], i));
		}
	}

	_numberPhotonsMC = 0;
	_momentumMCX.clear();
	_momentumMCY.clear();
	_momentumMCZ.clear();
	for (std::vector<std::pair<int, int> >::const_iterator iterator = decayVerticesAndTracks.begin(); iterator != decayVerticesAndTracks.end(); ++iterator)
	for (size_t i = 0; i < decayVerticesAndTracks.size(); ++i) {
		const PaMCvertex& decayVertex = event.vMCvertex(decayVerticesAndTracks[i].first);
		if (decayVertex.IsPileup() or decayVertex.NMCtrack() != 2) {
			continue;
		}
		const PaMCtrack& track1 = event.vMCtrack(decayVertex.iMCtrack(0));
		const PaMCtrack& track2 = event.vMCtrack(decayVertex.iMCtrack(1));
		if (track1.Pid() != MC_PID_PHOTON or track2.Pid() != MC_PID_PHOTON) {
			continue;
		}

		if (_numberPhotonsMC > MAX_CLUSTERS) {
			continue;
		}
		const TVector3 momentumVectorMC1 = track1.Mom3();
		_momentumMCX.push_back(momentumVectorMC1.X());
		_momentumMCY.push_back(momentumVectorMC1.Y());
		_momentumMCZ.push_back(momentumVectorMC1.Z());
		_numberPhotonsMC++;

		if (_numberPhotonsMC > MAX_CLUSTERS) {
			continue;
		}
		const TVector3 momentumVectorMC2 = track2.Mom3();
		_momentumMCX.push_back(momentumVectorMC2.X());
		_momentumMCY.push_back(momentumVectorMC2.Y());
		_momentumMCZ.push_back(momentumVectorMC2.Z());
		_numberPhotonsMC++;
	}
}
