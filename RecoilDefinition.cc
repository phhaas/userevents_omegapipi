#include "RecoilDefinition.h"

RecoilDefinition::RecoilDefinition(TTree& tree)
	: _phast(Phast::Ref()),
	  _momentumMCX(0),
	  _momentumMCY(0),
	  _momentumMCZ(0)
{
	_phast.h_file->cd();

	tree.Branch("Recoil_momentumMCX", &_momentumMCX, "momentumMCX/D");
	tree.Branch("Recoil_momentumMCY", &_momentumMCY, "momentumMCY/D");
	tree.Branch("Recoil_momentumMCZ", &_momentumMCZ, "momentumMCZ/D");
}


void
RecoilDefinition::fillMC(const PaEvent&    event,
                         const PaMCvertex& vertex)
{
	const PaMCtrack& track = getRecoilTrackMC(event, vertex);
	_momentumMCX = track.P(0);
	_momentumMCY = track.P(1);
	_momentumMCZ = track.P(2);
}


const PaMCtrack&
RecoilDefinition::getRecoilTrackMC(const PaEvent&    event,
                                   const PaMCvertex& vertex)
{
	for (int i = 0; i < vertex.NMCtrack(); ++i) {
		const PaMCtrack& track = event.vMCtrack(vertex.iMCtrack(i));
		if (track.IsPileup() or track.IsBeam()) {
			continue;
		}
		if (track.Pid() == MC_PID_PROTON) {
			return track;
		}
	}
	std::cerr << "No recoil MC track found! Aborting..." << std::endl;
	throw;
}
