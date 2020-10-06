#include "RPD_Helper.h"
#include "SelectionHelper.h"


SelectionHelper::SelectionHelper(TTree& tree)
	: _phast(Phast::Ref()),
	  _runNumber    (0),
	  _spillNumber  (0),
	  _eventInSpill (0),
	  _triggerMask  (0),
	  _masterTrigger(0),
	  _DT0          (0),
	  _tcsPhase     (0),
	  _timeInSpill  (0),
	  _isAcceptedMC (0)
{
	_phast.h_file->cd();
	_phast.HistFileDir("UserEvent");

	_hStatistics = new TH1D("statistics", "statistics", 20, -0.5, 19.5);

	tree.Branch("Selection_runNumber",     &_runNumber,     "runNumber/I");
	tree.Branch("Selection_spillNumber",   &_spillNumber,   "spillNumber/I");
	tree.Branch("Selection_eventInSpill",  &_eventInSpill,  "eventInSpill/I");
	tree.Branch("Selection_triggerMask",   &_triggerMask,   "triggerMask/I");
	tree.Branch("Selection_masterTrigger", &_masterTrigger, "masterTrigger/I");
	tree.Branch("Selection_DT0",           &_DT0,           "DT0/I");
	tree.Branch("Selection_tcsPhase",      &_tcsPhase,      "tcsPhase/D");
	tree.Branch("Selection_timeInSpill",   &_timeInSpill,   "timeInSpill/D");
	tree.Branch("Selection_isAcceptedMC",  &_isAcceptedMC,  "isAcceptedMC/I");
}


void
SelectionHelper::fill(const PaEvent& event)
{
	_runNumber     = event.RunNum();
	_spillNumber   = event.SpillNum();
	_eventInSpill  = event.EvInSpill();
	if (event.IsMC()) {
		if (isDT0()) {
			_triggerMask   = 1;
			_masterTrigger = 1;
		} else {
			_triggerMask   = 0;
			_masterTrigger = 0;
		}
	} else {
		_triggerMask = (event.TrigMask() & 0xfff);
	}
	_masterTrigger = event.MasterTriggerMask();
	//TODO _DT0 is never set
	_tcsPhase    = event.TCSphase();
	_timeInSpill = event.TimeInSpill();
}


void
SelectionHelper::handleMCEvent(TTree&         tree,
                               const bool     accepted,
                               const PaEvent& event)
{
	fill(event);
	setIsAcceptedMC(accepted);
	tree.Fill();
	Phast::Ref().h_file = tree.GetCurrentFile();
}


bool
SelectionHelper::isCorrectCharge(const PaParticle&         beamParticle,
                                 const vector<PaParticle>& particles) const
{
	int qSum = 0;
	for (vector<PaParticle>::const_iterator particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator) {
		const PaParticle& particle = *particleIterator;
		qSum += particle.Q();
	}
	return qSum == beamParticle.Q();
}


bool
SelectionHelper::isDT0() const
{
	RPD& rpd = RPD::Instance();
	return rpd.IsDT0();
}
