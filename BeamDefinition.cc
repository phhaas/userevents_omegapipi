#include "BeamDefinition.h"

BeamDefinition::BeamDefinition(TTree& tree)
	: _phast(Phast::Ref()),
	  _gradX      (0),
	  _gradY      (0),
	  _time       (0),
	  _firstZ     (0),
	  _lastZ      (0),
	  _minimalZ   (0),
	  _maximumZ   (0),
	  _chi2       (0),
	  _charge     (0),
	  _pid        (0),
	  _momentumMCX(0),
	  _momentumMCY(0),
	  _momentumMCZ(0)
{
	_phast.h_file->cd();

	tree.Branch("Beam_graddXdZ",    &_gradX,       "gradX/D");
	tree.Branch("Beam_graddYdZ",    &_gradY,       "gradY/D");
	tree.Branch("Beam_time",        &_time,        "time/D");
	tree.Branch("Beam_firstZ",      &_firstZ,      "firstZ/D");
	tree.Branch("Beam_lastZ",       &_lastZ,       "lastZ/D");
	tree.Branch("Beam_minimalZ",    &_minimalZ,    "minimalZ/D");
	tree.Branch("Beam_maximumZ",    &_maximumZ,    "maximumZ/D");
	tree.Branch("Beam_chi2",        &_chi2,        "chi2/D");
	tree.Branch("Beam_charge",      &_charge,      "charge/I");
	tree.Branch("Beam_pid",         &_pid,         "pid/I");
	tree.Branch("Beam_momentumMCX", &_momentumMCX, "momentumMCX/D");
	tree.Branch("Beam_momentumMCY", &_momentumMCY, "momentumMCY/D");
	tree.Branch("Beam_momentumMCZ", &_momentumMCZ, "momentumMCZ/D");
}


void
BeamDefinition::fill(const PaEvent&    event)
{
	const PaVertex& vertex = event.vVertex(event.iBestPrimaryVertex());
	const PaParticle& beamParticle = getBeamParticle(event, vertex);
	
	const int      vertexIndex = vertex.MyIndex();
	const PaTrack& track       = event.vTrack(beamParticle.iTrack());
	_gradX    = beamParticle.ParInVtx(vertexIndex).dXdZ();
	_gradY    = beamParticle.ParInVtx(vertexIndex).dYdZ();
	_time     = track.MeanTime();
	_firstZ   = track.ZFirst();
	_lastZ    = track.ZLast();
	_minimalZ = track.Zmin();
	_maximumZ = track.Zmax();
	_chi2     = track.Chi2tot() / track.Ndf();
	_charge   = beamParticle.Q();
	_pid      = beamParticle.PID();
}

void
BeamDefinition::fillMC(const PaMCtrack& track)
{
	// in MC beam particle goes the "wrong" direction
	_momentumMCX = (-1) * track.P(0);
	_momentumMCY = (-1) * track.P(1);
	_momentumMCZ = (-1) * track.P(2);
}


const PaMCtrack&
BeamDefinition::getBeamTrackMC(const PaEvent&    event,
                               const PaMCvertex& vertex)
{
	for (int i = 0; i < vertex.NMCtrack(); ++i) {
		const PaMCtrack& track = event.vMCtrack(vertex.iMCtrack(i));
		if (track.IsBeam() and not track.IsPileup() and track.Pid() == MC_PID_BEAM) {
			return track;
		}
	}
	std::cerr << "No MC beam track found! Aborting..." << std::endl;
	throw;
}
