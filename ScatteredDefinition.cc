#include "ScatteredDefinition.h"

namespace {

	template<typename T>
	bool
	setParticleBranches(TTree&             tree,
	                    std::vector<T>&    data,
	                    const std::string& branchName,
	                    const std::string& leafName,
	                    const std::string& leafType)
	{
		bool success = true;
		for (size_t i = 0; i < data.size(); ++i) {
			std::stringstream name;
			name << branchName << i + 1;
			std::stringstream leaf;
			leaf << leafName << i + 1 << "/" << leafType;
			if (not tree.Branch(name.str().c_str(), &(data[i]), leaf.str().c_str())) {
				success = false;
			}
		}
		return success;
	}

}


ScatteredDefinition::ScatteredDefinition(TTree&    tree,
                                         const int numberChargedParticles,
                                         const int numberNeutralParticles)
	: _phast(Phast::Ref()),
	  _numberChargedParticles(numberChargedParticles),
	  _numberNeutralParticles(numberNeutralParticles)
{
	_phast.h_file->cd();

	_numberTotalParticles = _numberChargedParticles + _numberNeutralParticles;
	_momentumX.resize  (_numberChargedParticles);
	_momentumY.resize  (_numberChargedParticles);
	_momentumZ.resize  (_numberChargedParticles);
	_charge.resize     (_numberChargedParticles);
	_pid.resize        (_numberChargedParticles);
	_minimalZ.resize   (_numberChargedParticles);
	_maximalZ.resize   (_numberChargedParticles);
	_time.resize       (_numberChargedParticles);
	_chi2.resize       (_numberChargedParticles);
	_XX0.resize        (_numberChargedParticles);
	_momentumMCX.resize(_numberTotalParticles);
	_momentumMCY.resize(_numberTotalParticles);
	_momentumMCZ.resize(_numberTotalParticles);
	_chargeMC.resize   (_numberTotalParticles);
	_pidMC.resize      (_numberTotalParticles);

	setParticleBranches<double>(tree, _momentumX,   "Scattered_Momentum_X",    "momentumX",   "D");
	setParticleBranches<double>(tree, _momentumY,   "Scattered_Momentum_Y",    "momentumY",   "D");
	setParticleBranches<double>(tree, _momentumZ,   "Scattered_Momentum_Z",    "momentumZ",   "D");
	setParticleBranches<int>   (tree, _charge,      "Scattered_Charge",        "charge",      "I");
	setParticleBranches<int>   (tree, _pid,         "Scattered_Pid",           "pid",         "I");
	setParticleBranches<double>(tree, _minimalZ,    "Scattered_Minimal_Z",     "minimalZ",    "D");
	setParticleBranches<double>(tree, _maximalZ,    "Scattered_Maximal_Z",     "maximalZ",    "D");
	setParticleBranches<double>(tree, _time,        "Scattered_Time",          "time",        "D");
	setParticleBranches<double>(tree, _chi2,        "Scattered_Chi2",          "chi2",        "D");
	setParticleBranches<double>(tree, _XX0,         "Scattered_XX0",           "XX0",         "D");
	setParticleBranches<double>(tree, _momentumMCX, "Scattered_Momentum_MC_X", "momentumMCX", "D");
	setParticleBranches<double>(tree, _momentumMCY, "Scattered_Momentum_MC_Y", "momentumMCY", "D");
	setParticleBranches<double>(tree, _momentumMCZ, "Scattered_Momentum_MC_Z", "momentumMCZ", "D");
	setParticleBranches<int>   (tree, _chargeMC,    "Scattered_Charge_MC",     "chargeMC",    "I");
	setParticleBranches<int>   (tree, _pidMC,       "Scattered_Pid_MC",        "pidMC",       "I");
}


void
ScatteredDefinition::fill(const PaEvent&            event,
                          const vector<PaParticle>& particles,
                          const PaVertex&           vertex)
{
	const size_t nmbParticles = particles.size();
	_momentumX.resize(nmbParticles);
	_momentumY.resize(nmbParticles);
	_momentumZ.resize(nmbParticles);
	_charge.resize   (nmbParticles);
	_pid.resize      (nmbParticles);
	_minimalZ.resize (nmbParticles);
	_maximalZ.resize (nmbParticles);
	_time.resize     (nmbParticles);
	_chi2.resize     (nmbParticles);
	_XX0.resize      (nmbParticles);
	for (size_t i = 0; i < nmbParticles; ++i) {
		const PaParticle& particle       = particles[i];
		const PaTrack&    track          = event.vTrack(particle.iTrack());
		const PaTPar&     parameter      = particle.ParInVtx(vertex.MyIndex());
		const TVector3    momentumVector = parameter.Mom3();

		_momentumX[i] = momentumVector.X();
		_momentumY[i] = momentumVector.Y();
		_momentumZ[i] = momentumVector.Z();
		_charge[i]    = particle.Q();
		_pid[i]       = particle.PID();
		_minimalZ[i]  = track.Zmin();
		_maximalZ[i]  = track.Zmax();
		_time[i]      = track.MeanTime();
		_chi2[i]      = track.Chi2tot() / track.Ndf();
		_XX0[i]       = track.XX0();
	}
}


void
ScatteredDefinition::fillMC(const vector<PaMCtrack>& tracks,
                            const PaMCvertex&        vertex)
{
	const size_t nmbTracks = tracks.size();
	_momentumMCX.resize(nmbTracks);
	_momentumMCY.resize(nmbTracks);
	_momentumMCZ.resize(nmbTracks);
	_chargeMC.resize   (nmbTracks);
	_pidMC.resize      (nmbTracks);
	for (size_t i = 0; i < nmbTracks; ++i) {
		const PaMCtrack& track            = tracks[i];
		const TVector3   momentumVectorMC = track.Mom3();

		_momentumMCX[i] = momentumVectorMC.X();
		_momentumMCY[i] = momentumVectorMC.Y();
		_momentumMCZ[i] = momentumVectorMC.Z();
		_chargeMC[i]    = track.Q();
		_pidMC[i]       = track.Pid();
	}
}


vector<PaParticle>
ScatteredDefinition::getScatteredParticles(const PaEvent&  event,
                                           const PaVertex& vertex)
{
	vector<PaParticle> scatteredParticles;
	for (int i = 0; i < vertex.NOutParticles(); ++i) {
		scatteredParticles.push_back(event.vParticle(vertex.iOutParticle(i)));
	}
	return scatteredParticles;
}


vector<PaTrack>
ScatteredDefinition::getScatteredParticleTracks(const PaEvent&  event,
                                                const PaVertex& vertex)
{
	vector<PaTrack>    tracks;
	vector<PaParticle> scatteredParticles = getScatteredParticles(event, vertex);
	for (size_t i = 0; i < scatteredParticles.size(); ++i) {
		const PaParticle& particle = scatteredParticles[i];
		tracks.push_back(event.vTrack(particle.iTrack()));
	}
	return tracks;
}


vector<PaMCtrack>
ScatteredDefinition::getScatteredParticleMCtracks(const PaEvent&    event,
                                                  const PaMCvertex& vertex)
{
	vector<PaMCtrack> tracks;
	for (int i = 0; i < vertex.NMCtrack(); ++i) {
		const PaMCtrack& track = event.vMCtrack(vertex.iMCtrack(i));
		if (track.IsBeam() or track.Pid() == MC_PID_PROTON) {
			continue;
		} else {
			tracks.push_back(track);
		}
	}
	return tracks;
}
