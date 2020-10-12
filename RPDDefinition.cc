#include "TString.h"

#include "RPDDefinition.h"

// Make this non-static file name
const TString PROTON_RANGES_FILE = "protonRanges.root";
static TLorentzVector* proton = new TLorentzVector(-10, -10, -10, -1);
enum method {noCorrection, oldCorrection, newCorrection};
enum targetType {hydrogenTarget, diskTarget};
enum material {hydrogen, aluminium, mylar, scintillator, lead, tungsten};

RPDDefinition::RPDDefinition(TTree& tree)
	: _phast(Phast::Ref()),
	  _rpd(RPD::Instance()),
	  _numberTracks        (0),
	  _hasTracks           (0),
	  _indexBestProton     (0),
	  _momentumX           (0),
	  _momentumY           (0),
	  _momentumZ           (0),
	  _beta                (0),
	  _phi                 (0),
	  _theta               (0),
	  _energyNoCorrection  (0),
	  _energyOldCorrection (0),
	  _energyNewCorrection (0),
	  _hitTime             (0),
	  _hitsRingA           (0),
	  _hitsRingB           (0),
	  _zHitPosition        (0),
	  _zHitPositionRingA   (0),
	  _zHitPositionRingB   (0),
	  _energyLossRingA     (0),
	  _energyLossRingB     (0)
{
	_phast.h_file->cd();
	
	tree.Branch("RPD_numberTracks",        &_numberTracks,        "numberTracks/I");
	tree.Branch("RPD_hasTracks",           &_hasTracks,           "hasTracks/I");
	tree.Branch("RPD_indexBestProton",     &_indexBestProton,     "indexBestProton/I");
	tree.Branch("RPD_momentumX",           &_momentumX,           "momentumX/D");
	tree.Branch("RPD_momentumY",           &_momentumY,           "momentumY/D");
	tree.Branch("RPD_momentumZ",           &_momentumZ,           "momentumZ/D");
	tree.Branch("RPD_beta",                &_beta,                "beta/D");
	tree.Branch("RPD_phi",                 &_phi,                 "phi/D");
	tree.Branch("RPD_theta",               &_theta,               "theta/D");
	tree.Branch("RPD_energyNoCorrection",  &_energyNoCorrection,  "energyNoCorrection/D");
	tree.Branch("RPD_energyOldCorrection", &_energyOldCorrection, "energyOldCorrection/D");
	tree.Branch("RPD_energyNewCorrection", &_energyNewCorrection, "energyNewCorrection/D");
	tree.Branch("RPD_hitTime",             &_hitTime,             "hitTime/D");
	tree.Branch("RPD_hitsRingA",           &_hitsRingA,           "hitsRingA/I");
	tree.Branch("RPD_hitsRingB",           &_hitsRingB,           "hitsRingB/I");
	tree.Branch("RPD_zHitPosition",        &_zHitPosition,        "zHitPosition/D");
	tree.Branch("RPD_zHitPositionRingA",   &_zHitPositionRingA,   "zHitPositionRingA/D");
	tree.Branch("RPD_zHitPositionRingB",   &_zHitPositionRingB,   "zHitPositionRingB/D");
	tree.Branch("RPD_energyLossRingA",     &_energyLossRingA,     "energyLossRingA/D");
	tree.Branch("RPD_energyLossRingB",     &_energyLossRingB,     "energyLossRingB/D");
}

RPDDefinition::~RPDDefinition()
{
  delete _momentumVector;
}

// needed to initialize RPD so _rpd.Search doesnt need to be called by other members
void RPDDefinition::initRPD(const PaEvent& event)
{
	const PaVertex& vertex = event.vVertex(event.iBestPrimaryVertex());
	_rpd.Search(event, vertex, RPD_CORRECT_MOMENTUM, RPD_CORRECT_ANGLES);
}

void
RPDDefinition::fill(const PaEvent&  event,
                    const PaVertex& vertex)
{
  
	_protons             = getProtons();
	_numberTracks        = _rpd.nTrack();
	_hasTracks           = _rpd.HasTracks() ? 1 : 0;
	_indexBestProton     = _rpd.iBestProtonTrack();
	_bestProton          = getBestProton();
	_momentumX           = _bestProton.X();
	_momentumY           = _bestProton.Y();
	_momentumZ           = _bestProton.Z();
	_momentumVector      = new TVector3(_momentumX, _momentumY, _momentumY);
	_beta                = _bestProton.Beta();
	_phi                 = getPhi();
	_theta               = getTheta();
	_energyNoCorrection  = _bestProton.E();
	_energyOldCorrection = getCorrectedEnergy(oldCorrection, event, vertex);
	_energyNewCorrection = getCorrectedEnergy(newCorrection, event, vertex);
	_hitTime             = _rpd.PV_T()[_indexBestProton];
	if (_hasTracks) {
		for (int i = 0; i < _numberTracks; ++i) {
			_hitsRingA += _rpd.Hits()[i].first;
			_hitsRingB += _rpd.Hits()[i].second;
		}
	}
	_zHitPosition      = _rpd.PV_Z()[_indexBestProton];
	_zHitPositionRingA = _rpd.ZA()[_indexBestProton];
	_zHitPositionRingB = _rpd.ZB()[_indexBestProton];
	_energyLossRingA   = _rpd.dE_A()[_indexBestProton];
	_energyLossRingB   = _rpd.dE_B()[_indexBestProton];
}


bool
RPDDefinition::cutHasTracks() const
{
	return _rpd.HasTracks();
}


bool
RPDDefinition::cutHasBestProton() const
{
	return _rpd.iBestProtonTrack() > -1;
}


const TLorentzVector&
RPDDefinition::getBestProton() const
{
	if (_indexBestProton < 0 or _indexBestProton > (int)_protons.size()) {
		return* proton;
	}
	return _protons[_indexBestProton];
}


double
RPDDefinition::getPhi() const
{
	if (_indexBestProton < 0 or  _indexBestProton > (int)_protons.size()) {
		return -5;
	}
	return _bestProton.Phi();
}


double
RPDDefinition::getTheta() const
{
	if (_indexBestProton < 0 or _indexBestProton > (int) _protons.size()) {
		return -5;
	}
	return _bestProton.Theta();
}


double
RPDDefinition::getEnergy(const int useMethod) const // Method: 0=noCorrection, 1=oldCorrection, 2=newCorrection
{
	switch (useMethod) {
	// Method: 0=noCorrection, 1=oldCorrection, 2=newCorrection
		case noCorrection: // No correction
			return _energyNoCorrection;
		case oldCorrection:  // Old correction
			return _energyOldCorrection;
		case newCorrection:  // New correction
			return _energyNewCorrection;
		default:
			std::cerr << "Undefined energy type in function:" << __PRETTY_FUNCTION__ << std::endl;
			return -1;
	}
}


// assumes that fill() has been called already
double
RPDDefinition::getCorrectedEnergy(const int       useMethod, // Method: 0=noCorrection, 1=oldCorrection, 2=newCorrection
                                  const PaEvent&  event,
                                  const PaVertex& vertex) const
{
	const double rpdRingAThickness       = 5.0;    // Thickness of scintillator [mm]
	const double rpdChaussetteThickness  = 0.2;    // Thickness of Mylar        [mm]
	// Hydrogen target
	const double targetAluminumThickness = 1.8;    // Thickness of Aluminium    [mm]
	const double targetCellThickness     = 0.125;  // Thickness of Mylar        [mm]
	const double targetCellRadius        = 17.5;   // radius of H_2 target cell [mm]
	// Disk target
	const double targetDiskThicknesses[16] = {
		0.125, 0.125,  0.125,  0.125,
		0.125, 0.125,  0.0625, 0.0625,
		0.025, 0.0125, 0.025,  0.0125,
		0.025, 0.0125, 0.025,  0.0125
	};  // [mm]
	const double targetDiskCenterZPos[18] = {
		-100.0, 0.0, 24.182, 47.582,
		72.182, 96.182, 120.882, 144.1195,
		168.1195, 192.082, 216.0695, 240.0820,
		264.0695, 288.082, 312.0695, 336.082,
		360.0695, 400.0
	};  // [cm]
	
	 // 0=hydrogen, 1=aluminium, 2=mylar, 3=scintillator, 4=lead, 5=tungesten
	const int targetMaterial[16] = {
		tungsten, tungsten, tungsten, tungsten,
		tungsten, tungsten, tungsten, tungsten,
		tungsten, tungsten, 6, 6,
		tungsten, tungsten, 6, 6
	};

	// get uncorrected momentum
	double       pRPD = _energyNoCorrection;
	const double theta = getTheta();
	const double phi   = getPhi();

	const int runNumber = event.RunNum();
	int useDiskTarget = diskTarget; // diskTarget=1
	if (runNumber < 80608 or runNumber > 81093) {
		useDiskTarget = hydrogenTarget; // hydrogenTarget=0
	}
	switch (useMethod) {
	// Method: 0=noCorrection, 1=oldCorrection, 2=newCorrection
	  
		// No correction
		case noCorrection:
			break;
		// Old correction
		case oldCorrection:
			if (pRPD >= 0.596) {
				break;
			}
			// RPD
			pRPD = getInterpolatedEnergy(pRPD, scintillator, rpdRingAThickness      / TMath::Sin(theta));
			pRPD = getInterpolatedEnergy(pRPD, mylar, rpdChaussetteThickness / TMath::Sin(theta));
			switch (useDiskTarget) {
				case hydrogenTarget: {  // 1=hydrogenTarget
					pRPD = getInterpolatedEnergy(pRPD, aluminium, targetAluminumThickness / TMath::Sin(theta));
					pRPD = getInterpolatedEnergy(pRPD, mylar, targetCellThickness     / TMath::Sin(theta));
					double traversedHydrogen = TMath::Sqrt(pow(targetCellRadius * TMath::Cos(phi) - vertex.X(), 2) + pow(targetCellRadius * TMath::Sin(phi) - vertex.Y(), 2));
					traversedHydrogen        = traversedHydrogen / TMath::Sin(theta);
					pRPD = getInterpolatedEnergy(pRPD, hydrogen, traversedHydrogen);
					break;
				}
				case diskTarget:  // 2=diskTarget
					for (size_t i = 1; i < 17; ++i) {
						if (    (targetDiskCenterZPos[i - 1] + targetDiskCenterZPos[i]    ) / 2. < (vertex.Z() + 66.7)
						    and (targetDiskCenterZPos[i]     + targetDiskCenterZPos[i + 1]) / 2. > (vertex.Z() + 66.7)) {
							const int targetDiskIndex = i - 1;
							pRPD = getInterpolatedEnergy(pRPD, targetMaterial[targetDiskIndex], targetDiskThicknesses[targetDiskIndex] / (2. * TMath::Cos(theta)));
							break;
						}
					}
					break;
				default:
					std::cerr << "Undefined run type in function:" << __PRETTY_FUNCTION__ << std::endl;
					return -1;
			}
			break;
		// New correction
		case newCorrection:
			// RPD
			pRPD = getCorrectedEnergyNew(pRPD, scintillator, rpdRingAThickness      / TMath::Sin(theta));
			pRPD = getCorrectedEnergyNew(pRPD, mylar, rpdChaussetteThickness / TMath::Sin(theta));
			switch (useDiskTarget) {
				case hydrogenTarget: {  // 1=hydrogenTarget
					pRPD = getCorrectedEnergyNew(pRPD, aluminium, targetAluminumThickness / TMath::Sin(theta));
					pRPD = getCorrectedEnergyNew(pRPD, mylar, targetCellThickness     / TMath::Sin(theta));
					double traversedHydrogen = TMath::Sqrt(pow(targetCellRadius * TMath::Cos(phi) - vertex.X(), 2) + pow(targetCellRadius * TMath::Sin(phi) - vertex.Y(), 2));
					traversedHydrogen        = traversedHydrogen / TMath::Sin(theta);
					pRPD = getCorrectedEnergyNew(pRPD, hydrogen, traversedHydrogen);
					break;
				}
				case diskTarget:  // Disk target
					for (size_t i = 1; i < 17; ++i) {
						if (    (targetDiskCenterZPos[i - 1] + targetDiskCenterZPos[i]    ) / 2. < (vertex.Z() + 66.7)
						    and (targetDiskCenterZPos[i]     + targetDiskCenterZPos[i + 1]) / 2. > (vertex.Z() + 66.7)) {
							const int targetDiskIndex = i - 1;
							pRPD = getCorrectedEnergyNew(pRPD, targetMaterial[targetDiskIndex], targetDiskThicknesses[targetDiskIndex] / (2.0 * TMath::Cos(theta)));
							break;
						}
					}
					break;
				default:
					std::cerr << "Undefined run type in function:" << __PRETTY_FUNCTION__ << std::endl;
					return -1;
			}
			break;
		default:
			std::cerr << "Undefined m ethod type in function:" << __PRETTY_FUNCTION__ << std::endl;
			return -1;
	}  // switch (useMethod)

	return pRPD;
}


double
RPDDefinition::getCorrectedEnergyNew(const double pRPD,                      // proton momentum before correction [GeV]
                                     const int    materialType,              // 0=hydrogen, 1=aluminium, 2=mylar, 3=scintillator, 4=lead, 5=tungesten
                                     const double traversedThickness) const  // [mm]
{
	// do not correct if outside range if correction tables
	if (pRPD < 2e-3 or pRPD >= 8e3) {
		return pRPD;
	}

	static TFile* protonRanges = new TFile(PaUtils::PhastHome() + "/user/" + PROTON_RANGES_FILE, "READ");
	if (not protonRanges or protonRanges->IsZombie()) {
		std::cerr << "Error in function " << __PRETTY_FUNCTION__ << ". "
		          << "Could not open file with proton range table at '" << PROTON_RANGES_FILE << "'. Aborting..." << std::endl;
		throw;
	}

	static TSpline3* rangeFunction  = NULL;
	static TSpline3* energyFunction = NULL;

	static TSpline3* rangeFunctionAluminium     = (TSpline3*) protonRanges->Get("ALUMINIUM/splineEtoR");
	static TSpline3* energyFunctionAluminium    = (TSpline3*) protonRanges->Get("ALUMINIUM/splineRtoE");
	static TSpline3* rangeFunctionHydrogen      = (TSpline3*) protonRanges->Get("HYDROGEN/splineEtoR");
	static TSpline3* energyFunctionHydrogen     = (TSpline3*) protonRanges->Get("HYDROGEN/splineRtoE");
	static TSpline3* rangeFunctionMylar         = (TSpline3*) protonRanges->Get("POLYETHYLENE TEREPHTHALATE (MYLAR)/splineEtoR");
	static TSpline3* energyFunctionMylar        = (TSpline3*) protonRanges->Get("POLYETHYLENE TEREPHTHALATE (MYLAR)/splineRtoE");
	static TSpline3* rangeFunctionScintillator  = (TSpline3*) protonRanges->Get("PLASTIC SCINTILLATOR/splineEtoR");
	static TSpline3* energyFunctionScintillator = (TSpline3*) protonRanges->Get("PLASTIC SCINTILLATOR/splineRtoE");
	static TSpline3* rangeFunctionLead          = (TSpline3*) protonRanges->Get("LEAD/splineEtoR");
	static TSpline3* energyFunctionLead         = (TSpline3*) protonRanges->Get("LEAD/splineRtoE");
	static TSpline3* rangeFunctionTungsten      = (TSpline3*) protonRanges->Get("TUNGSTEN/splineEtoR");
	static TSpline3* energyFunctionTungsten     = (TSpline3*) protonRanges->Get("TUNGSTEN/splineRtoE");

	switch (materialType) { // 0=hydrogen, 1=aluminium, 2=mylar, 3=scintillator, 4=lead, 5=tungesten
		case hydrogen: // Hydogen
			rangeFunction  = rangeFunctionHydrogen;
			energyFunction = energyFunctionHydrogen;
			break;
		case aluminium: // Aluminium
			rangeFunction  = rangeFunctionAluminium;
			energyFunction = energyFunctionAluminium;
			break;
		case mylar: // Mylar
			rangeFunction  = rangeFunctionMylar;
			energyFunction = energyFunctionMylar;
			break;
		case scintillator: // Scintillator
			rangeFunction  = rangeFunctionScintillator;
			energyFunction = energyFunctionScintillator;
			break;
		case lead: // Lead
			rangeFunction  = rangeFunctionLead;
			energyFunction = energyFunctionLead;
			break;
		case tungsten: // Tungsten
			rangeFunction  = rangeFunctionTungsten;
			energyFunction = energyFunctionTungsten;
			break;
		default:
			std::cerr << "Undefined material type in function:" << __PRETTY_FUNCTION__ << std::endl;
			return -1;
	}

	const double EkinBefore = (std::sqrt(pRPD * pRPD + PROTON_MASS * PROTON_MASS) - PROTON_MASS) * 1000;  // kinetic energy before material in [MeV]
	const double range      = rangeFunction->Eval(EkinBefore) + traversedThickness / 10;
	const double EkinAfter  = energyFunction->Eval(range) / 1000;  // kinetic energy after passing through material [GeV]

	return std::sqrt(EkinAfter * EkinAfter + 2 * EkinAfter * PROTON_MASS);
}


double
RPDDefinition::getInterpolatedEnergy(const double pRPD,
                                     const int    materialType,  // 0=hydrogen, 1=aluminium, 2=mylar, 3=scintillator, 4=lead, 5=tungesten
                                     const double traversedThickness) const
{
	if (pRPD < 0) {
		return pRPD;
	}

	const size_t nData = 151;
	static const double Ek[nData] = {
		4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0,
		40.0, 44.0, 48.0, 52.0, 56.0, 60.0, 64.0, 68.0, 72.0, 76.0,
		80.0, 84.0, 88.0, 92.0, 96.0, 100.0, 104.0, 108.0, 112.0,
		116.0, 120.0, 124.0, 128.0, 132.0, 136.0, 140.0, 144.0, 148.0,
		152.0, 156.0, 160.0, 164.0, 168.0, 172.0, 176.0, 180.0, 184.0,
		188.0, 192.0, 196.0, 200.0, 204.0, 208.0, 212.0, 216.0, 220.0,
		224.0, 228.0, 232.0, 236.0, 240.0, 244.0, 248.0, 252.0, 256.0,
		260.0, 264.0, 268.0, 272.0, 276.0, 280.0, 284.0, 288.0, 292.0,
		296.0, 300.0, 304.0, 308.0, 312.0, 316.0, 320.0, 324.0, 328.0,
		332.0, 336.0, 340.0, 344.0, 348.0, 352.0, 356.0, 360.0, 364.0,
		368.0, 372.0, 376.0, 380.0, 384.0, 388.0, 392.0, 396.0, 400.0,
		404.0, 408.0, 412.0, 416.0, 420.0, 424.0, 428.0, 432.0, 436.0,
		440.0, 444.0, 448.0, 452.0, 456.0, 460.0, 464.0, 468.0, 472.0,
		476.0, 480.0, 484.0, 488.0, 492.0, 496.0, 500.0, 504.0, 508.0,
		512.0, 516.0, 520.0, 524.0, 528.0, 532.0, 536.0, 540.0, 544.0,
		548.0, 552.0, 556.0, 560.0, 564.0, 568.0, 572.0, 576.0, 580.0,
		584.0, 588.0, 592.0, 596.0, 600.0};

	static const double range_H2[nData] = {
		0.1373E+01, 0.4879E+01, 0.1030E+02,
		0.1751E+02, 0.2638E+02, 0.3700E+02, 0.4917E+02, 0.6288E+02,
		0.7799E+02, 0.9461E+02, 0.1130E+03, 0.1323E+03, 0.1534E+03,
		0.1756E+03, 0.1994E+03, 0.2239E+03, 0.2507E+03, 0.2774E+03,
		0.3067E+03, 0.3368E+03, 0.3668E+03, 0.4001E+03, 0.4338E+03,
		0.4674E+03, 0.5035E+03, 0.5411E+03, 0.5787E+03, 0.6164E+03,
		0.6574E+03, 0.6993E+03, 0.7413E+03, 0.7832E+03, 0.8269E+03,
		0.8735E+03, 0.9202E+03, 0.9669E+03, 0.1014E+04, 0.1061E+04,
		0.1113E+04, 0.1165E+04, 0.1217E+04, 0.1268E+04, 0.1320E+04,
		0.1373E+04, 0.1431E+04, 0.1488E+04, 0.1545E+04, 0.1602E+04,
		0.1659E+04, 0.1716E+04, 0.1776E+04, 0.1839E+04, 0.1902E+04,
		0.1965E+04, 0.2028E+04, 0.2091E+04, 0.2154E+04, 0.2216E+04,
		0.2282E+04, 0.2351E+04, 0.2420E+04, 0.2489E+04, 0.2558E+04,
		0.2627E+04, 0.2696E+04, 0.2765E+04, 0.2834E+04, 0.2905E+04,
		0.2980E+04, 0.3055E+04, 0.3130E+04, 0.3205E+04, 0.3281E+04,
		0.3356E+04, 0.3431E+04, 0.3506E+04, 0.3581E+04, 0.3657E+04,
		0.3736E+04, 0.3818E+04, 0.3899E+04, 0.3981E+04, 0.4062E+04,
		0.4144E+04, 0.4226E+04, 0.4307E+04, 0.4389E+04, 0.4470E+04,
		0.4552E+04, 0.4633E+04, 0.4718E+04, 0.4806E+04, 0.4894E+04,
		0.4982E+04, 0.5070E+04, 0.5158E+04, 0.5246E+04, 0.5334E+04,
		0.5422E+04, 0.5510E+04, 0.5598E+04, 0.5686E+04, 0.5774E+04,
		0.5862E+04, 0.5951E+04, 0.6045E+04, 0.6140E+04, 0.6234E+04,
		0.6329E+04, 0.6423E+04, 0.6517E+04, 0.6612E+04, 0.6706E+04,
		0.6801E+04, 0.6895E+04, 0.6989E+04, 0.7084E+04, 0.7178E+04,
		0.7273E+04, 0.7367E+04, 0.7462E+04, 0.7560E+04, 0.7661E+04,
		0.7762E+04, 0.7862E+04, 0.7963E+04, 0.8064E+04, 0.8164E+04,
		0.8265E+04, 0.8366E+04, 0.8466E+04, 0.8567E+04, 0.8668E+04,
		0.8768E+04, 0.8869E+04, 0.8970E+04, 0.9070E+04, 0.9171E+04,
		0.9272E+04, 0.9372E+04, 0.9477E+04, 0.9584E+04, 0.9690E+04,
		0.9797E+04, 0.9904E+04, 0.1001E+05, 0.1012E+05, 0.1022E+05,
		0.1033E+05, 0.1044E+05};

	static const double range_Aluminium[nData] = {
		0.1290E+00, 0.4173E+00, 0.8474E+00,
		0.1408E+01, 0.2087E+01, 0.2893E+01, 0.3809E+01, 0.4833E+01,
		0.5958E+01, 0.7188E+01, 0.8542E+01, 0.9960E+01, 0.1151E+02,
		0.1312E+02, 0.1486E+02, 0.1664E+02, 0.1857E+02, 0.2051E+02,
		0.2262E+02, 0.2479E+02, 0.2695E+02, 0.2934E+02, 0.3175E+02,
		0.3416E+02, 0.3674E+02, 0.3943E+02, 0.4211E+02, 0.4480E+02,
		0.4771E+02, 0.5069E+02, 0.5367E+02, 0.5665E+02, 0.5975E+02,
		0.6305E+02, 0.6635E+02, 0.6965E+02, 0.7295E+02, 0.7633E+02,
		0.7998E+02, 0.8362E+02, 0.8727E+02, 0.9092E+02, 0.9456E+02,
		0.9832E+02, 0.1023E+03, 0.1063E+03, 0.1104E+03, 0.1144E+03,
		0.1184E+03, 0.1224E+03, 0.1266E+03, 0.1310E+03, 0.1354E+03,
		0.1398E+03, 0.1442E+03, 0.1486E+03, 0.1530E+03, 0.1574E+03,
		0.1620E+03, 0.1668E+03, 0.1716E+03, 0.1764E+03, 0.1812E+03,
		0.1860E+03, 0.1908E+03, 0.1956E+03, 0.2004E+03, 0.2053E+03,
		0.2105E+03, 0.2158E+03, 0.2210E+03, 0.2262E+03, 0.2314E+03,
		0.2366E+03, 0.2419E+03, 0.2471E+03, 0.2523E+03, 0.2575E+03,
		0.2630E+03, 0.2687E+03, 0.2743E+03, 0.2800E+03, 0.2856E+03,
		0.2912E+03, 0.2969E+03, 0.3025E+03, 0.3082E+03, 0.3138E+03,
		0.3195E+03, 0.3251E+03, 0.3309E+03, 0.3370E+03, 0.3431E+03,
		0.3491E+03, 0.3552E+03, 0.3613E+03, 0.3673E+03, 0.3734E+03,
		0.3795E+03, 0.3855E+03, 0.3916E+03, 0.3977E+03, 0.4037E+03,
		0.4098E+03, 0.4159E+03, 0.4224E+03, 0.4289E+03, 0.4354E+03,
		0.4419E+03, 0.4484E+03, 0.4549E+03, 0.4614E+03, 0.4678E+03,
		0.4743E+03, 0.4808E+03, 0.4873E+03, 0.4938E+03, 0.5003E+03,
		0.5068E+03, 0.5133E+03, 0.5197E+03, 0.5265E+03, 0.5334E+03,
		0.5403E+03, 0.5472E+03, 0.5541E+03, 0.5610E+03, 0.5679E+03,
		0.5748E+03, 0.5817E+03, 0.5886E+03, 0.5955E+03, 0.6024E+03,
		0.6093E+03, 0.6161E+03, 0.6230E+03, 0.6299E+03, 0.6368E+03,
		0.6437E+03, 0.6506E+03, 0.6578E+03, 0.6651E+03, 0.6723E+03,
		0.6796E+03, 0.6869E+03, 0.6942E+03, 0.7014E+03, 0.7087E+03,
		0.7160E+03, 0.7233E+03};

	static const double range_Mylar[nData] = {
		0.1851E+00, 0.6286E+00, 0.1300E+01,
		0.2182E+01, 0.3259E+01, 0.4540E+01, 0.6001E+01, 0.7640E+01,
		0.9442E+01, 0.1142E+02, 0.1360E+02, 0.1588E+02, 0.1837E+02,
		0.2099E+02, 0.2379E+02, 0.2667E+02, 0.2981E+02, 0.3294E+02,
		0.3637E+02, 0.3988E+02, 0.4339E+02, 0.4728E+02, 0.5120E+02,
		0.5512E+02, 0.5932E+02, 0.6369E+02, 0.6806E+02, 0.7244E+02,
		0.7719E+02, 0.8206E+02, 0.8692E+02, 0.9178E+02, 0.9683E+02,
		0.1022E+03, 0.1076E+03, 0.1130E+03, 0.1184E+03, 0.1239E+03,
		0.1299E+03, 0.1359E+03, 0.1418E+03, 0.1478E+03, 0.1538E+03,
		0.1599E+03, 0.1665E+03, 0.1731E+03, 0.1797E+03, 0.1862E+03,
		0.1928E+03, 0.1994E+03, 0.2062E+03, 0.2135E+03, 0.2207E+03,
		0.2279E+03, 0.2351E+03, 0.2424E+03, 0.2496E+03, 0.2568E+03,
		0.2644E+03, 0.2723E+03, 0.2801E+03, 0.2880E+03, 0.2959E+03,
		0.3038E+03, 0.3117E+03, 0.3196E+03, 0.3275E+03, 0.3356E+03,
		0.3442E+03, 0.3528E+03, 0.3614E+03, 0.3700E+03, 0.3786E+03,
		0.3872E+03, 0.3958E+03, 0.4044E+03, 0.4130E+03, 0.4215E+03,
		0.4306E+03, 0.4399E+03, 0.4492E+03, 0.4585E+03, 0.4678E+03,
		0.4771E+03, 0.4864E+03, 0.4957E+03, 0.5050E+03, 0.5143E+03,
		0.5236E+03, 0.5329E+03, 0.5425E+03, 0.5525E+03, 0.5626E+03,
		0.5726E+03, 0.5826E+03, 0.5926E+03, 0.6026E+03, 0.6126E+03,
		0.6226E+03, 0.6327E+03, 0.6427E+03, 0.6527E+03, 0.6627E+03,
		0.6727E+03, 0.6829E+03, 0.6936E+03, 0.7043E+03, 0.7150E+03,
		0.7257E+03, 0.7365E+03, 0.7472E+03, 0.7579E+03, 0.7686E+03,
		0.7794E+03, 0.7901E+03, 0.8008E+03, 0.8115E+03, 0.8222E+03,
		0.8330E+03, 0.8437E+03, 0.8544E+03, 0.8656E+03, 0.8770E+03,
		0.8884E+03, 0.8998E+03, 0.9112E+03, 0.9226E+03, 0.9341E+03,
		0.9455E+03, 0.9569E+03, 0.9683E+03, 0.9797E+03, 0.9911E+03,
		0.1002E+04, 0.1014E+04, 0.1025E+04, 0.1037E+04, 0.1048E+04,
		0.1060E+04, 0.1071E+04, 0.1083E+04, 0.1095E+04, 0.1107E+04,
		0.1119E+04, 0.1131E+04, 0.1143E+04, 0.1155E+04, 0.1167E+04,
		0.1179E+04, 0.1191E+04};

	static const double range_Scintillator[nData] = {
		0.2301E+00, 0.7883E+00, 0.1636E+01,
		0.2752E+01, 0.4114E+01, 0.5737E+01, 0.7590E+01, 0.9670E+01,
		0.1196E+02, 0.1446E+02, 0.1723E+02, 0.2014E+02, 0.2331E+02,
		0.2663E+02, 0.3019E+02, 0.3385E+02, 0.3784E+02, 0.4183E+02,
		0.4620E+02, 0.5066E+02, 0.5513E+02, 0.6008E+02, 0.6508E+02,
		0.7007E+02, 0.7541E+02, 0.8098E+02, 0.8655E+02, 0.9212E+02,
		0.9818E+02, 0.1044E+03, 0.1106E+03, 0.1168E+03, 0.1232E+03,
		0.1301E+03, 0.1370E+03, 0.1438E+03, 0.1507E+03, 0.1578E+03,
		0.1654E+03, 0.1730E+03, 0.1806E+03, 0.1882E+03, 0.1958E+03,
		0.2036E+03, 0.2120E+03, 0.2204E+03, 0.2288E+03, 0.2372E+03,
		0.2456E+03, 0.2540E+03, 0.2627E+03, 0.2719E+03, 0.2812E+03,
		0.2904E+03, 0.2996E+03, 0.3088E+03, 0.3180E+03, 0.3272E+03,
		0.3369E+03, 0.3470E+03, 0.3570E+03, 0.3671E+03, 0.3772E+03,
		0.3873E+03, 0.3973E+03, 0.4074E+03, 0.4175E+03, 0.4278E+03,
		0.4388E+03, 0.4498E+03, 0.4608E+03, 0.4717E+03, 0.4827E+03,
		0.4937E+03, 0.5046E+03, 0.5156E+03, 0.5266E+03, 0.5375E+03,
		0.5491E+03, 0.5610E+03, 0.5729E+03, 0.5848E+03, 0.5966E+03,
		0.6085E+03, 0.6204E+03, 0.6323E+03, 0.6442E+03, 0.6560E+03,
		0.6679E+03, 0.6798E+03, 0.6920E+03, 0.7048E+03, 0.7176E+03,
		0.7304E+03, 0.7432E+03, 0.7560E+03, 0.7688E+03, 0.7816E+03,
		0.7944E+03, 0.8072E+03, 0.8200E+03, 0.8328E+03, 0.8456E+03,
		0.8584E+03, 0.8713E+03, 0.8850E+03, 0.8987E+03, 0.9124E+03,
		0.9261E+03, 0.9398E+03, 0.9535E+03, 0.9672E+03, 0.9809E+03,
		0.9946E+03, 0.1008E+04, 0.1022E+04, 0.1036E+04, 0.1049E+04,
		0.1063E+04, 0.1077E+04, 0.1090E+04, 0.1105E+04, 0.1119E+04,
		0.1134E+04, 0.1149E+04, 0.1163E+04, 0.1178E+04, 0.1192E+04,
		0.1207E+04, 0.1221E+04, 0.1236E+04, 0.1251E+04, 0.1265E+04,
		0.1280E+04, 0.1294E+04, 0.1309E+04, 0.1324E+04, 0.1338E+04,
		0.1353E+04, 0.1367E+04, 0.1382E+04, 0.1398E+04, 0.1413E+04,
		0.1429E+04, 0.1444E+04, 0.1459E+04, 0.1475E+04, 0.1490E+04,
		0.1506E+04, 0.1521E+04};

	static const double range_Lead[nData] = {
		0.0313921, 0.102026, 0.206344, 0.341674, 0.506432, 0.699295,
		0.918943, 1.16564, 1.437, 1.73392, 2.05463, 2.39912, 2.76652,
		3.15595, 3.56828, 4.00176, 4.45639, 4.93216, 5.42819, 5.94449,
		6.48106, 7.03612, 7.61057, 8.20352, 8.81938, 9.44493, 10.0881,
		10.7577, 11.4361, 12.141, 12.8546, 13.5859, 14.3348, 15.1013,
		15.8855, 16.6784, 17.489, 18.3172, 19.163, 20.0176, 20.8899,
		21.7709, 22.6696, 23.5859, 24.511, 25.4537, 26.4053, 27.3744,
		28.3524, 29.348, 30.3524, 31.3744, 32.4053, 33.4449, 34.5022,
		35.5683, 36.6432, 37.7357, 38.837, 39.9471, 41.0749, 42.2026,
		43.348, 44.511, 45.674, 46.8546, 48.0441, 49.2423, 50.4493,
		51.6652, 52.8899, 54.1322, 55.3833, 56.6344, 57.9031, 59.1806,
		60.467, 61.7533, 63.0573, 64.37, 65.6916, 67.022, 68.3612,
		69.7093, 71.0573, 72.4229, 73.7974, 75.1718, 76.5639, 77.9559,
		79.3568, 80.7665, 82.185, 83.6123, 85.0485, 86.4846, 87.9383,
		89.4273, 90.837, 92.3348, 93.8326, 95.2423, 96.7401, 98.2379,
		99.7357, 101.322, 102.819, 104.317, 105.815, 107.401, 108.899,
		110.485, 112.07, 113.568, 115.154, 116.74, 118.326, 119.912,
		121.498, 123.084, 124.67, 126.256, 127.93, 129.515, 131.189,
		132.775, 134.449, 136.035, 137.709, 139.383, 140.969, 142.643,
		144.317, 145.991, 147.665, 149.339, 151.013, 152.687, 154.449,
		156.123, 157.797, 159.559, 161.233, 162.996, 164.67, 166.432,
		168.106, 169.868, 171.63, 173.392};

	static const double range_Tungsten[nData] = {
		0.0743877, 0.214537, 0.408106, 0.65022, 0.937445,
		1.26784, 1.63877, 2.04934, 2.4978, 2.98238, 3.5022,
		4.05639, 4.64405, 5.26432, 5.9163, 6.6, 7.31366,
		8.05815, 8.82819, 9.62996, 10.467, 11.3216, 12.2115,
		13.1189, 14.0617, 15.0308, 16.0176, 17.0396, 18.0793,
		19.1454, 20.2379, 21.348, 22.4846, 23.6476, 24.8282,
		26.0352, 27.2687, 28.511, 29.7885, 31.0749, 32.3877,
		33.7269, 35.0749, 36.4493, 37.8414, 39.2511, 40.6872,
		42.1322, 43.6035, 45.0837, 46.5903, 48.1145, 49.6564,
		51.207, 52.7841, 54.37, 55.9824, 57.6035, 59.2423,
		60.8987, 62.5727, 64.2555, 65.9648, 67.6828, 69.4097,
		71.1542, 72.9163, 74.696, 76.4846, 78.2907, 80.1057,
		81.9383, 83.7885, 85.6476, 87.5154, 89.4273, 91.2775,
		93.2159, 95.1542, 97.0925, 99.0308, 100.969, 102.907,
		104.934, 106.96, 108.899, 110.925, 112.952, 115.066,
		117.093, 119.119, 121.233, 123.26, 125.374, 127.489,
		129.604, 131.718, 133.833, 136.035, 138.15, 140.352,
		142.555, 144.67, 146.872, 149.075, 151.278, 153.568,
		155.771, 157.974, 160.264, 162.467, 164.758, 167.048,
		169.339, 171.63, 173.921, 176.211, 178.502, 180.881,
		183.172, 185.551, 187.841, 190.22, 192.599, 194.978,
		197.357, 199.736, 202.115, 204.493, 206.872, 209.339,
		211.718, 214.185, 216.564, 219.031, 221.498, 223.877,
		226.344, 228.811, 231.278, 233.744, 236.3, 238.767,
		241.233, 243.789, 246.256, 248.722, 251.278, 253.833, 256.3};

	const double EkinBefore = (pRPD - PROTON_MASS) * 1000.;
	size_t i;
	for (i = 1; i < nData - 1; ++i) {
		if (EkinBefore < Ek[i]) {
			break;
		}
	}
	// reached end of table; do not correct
	if (i == nData - 1) {
		return pRPD / 1000.;
	}

	double a, b, c;
	 // 0=hydrogen, 1=aluminium, 2=mylar, 3=scintillator, 4=lead, 5=tungesten
	switch (materialType) {
		case hydrogen:
			a = range_H2[i - 1];
			b = range_H2[i];
			c = range_H2[i + 1];
			break;
		case aluminium:
			a = range_Aluminium[i - 1];
			b = range_Aluminium[i];
			c = range_Aluminium[i + 1];
			break;
		case mylar:
			a = range_Mylar[i - 1];
			b = range_Mylar[i];
			c = range_Mylar[i + 1];
			break;
		case scintillator: 
			a = range_Scintillator[i - 1];
			b = range_Scintillator[i];
			c = range_Scintillator[i + 1];
			break;
		case lead:
			a = range_Lead[i - 1];
			b = range_Lead[i];
			c = range_Lead[i + 1];
			break;
		case tungsten:
			a = range_Tungsten[i - 1];
			b = range_Tungsten[i];
			c = range_Tungsten[i + 1];
			break;
		default:
			std::cerr << "Undefined material type in function:" << __PRETTY_FUNCTION__ << std::endl;
			return -1;
	}

	const double range     = interpolate3(Ek[i - 1], Ek[i], Ek[i + 1], a, b, c, EkinBefore) + traversedThickness;
	const double EkinAfter = interpolate3(a, b, c, Ek[i - 1], Ek[i], Ek[i + 1], range) / 1000;
	if (std::isnan(EkinAfter) or std::isinf(EkinAfter) or EkinAfter < 0) {
		return 0;
	} else {
		return std::sqrt(EkinAfter * EkinAfter + 2 * EkinAfter * PROTON_MASS);
	}
}


double
RPDDefinition::interpolate3(const double xx1,
                            const double xx2,
                            const double xx3,
                            const double yy1,
                            const double yy2,
                            const double yy3,
                            const double x) const
{
	return   yy1 * (x - xx2) * (x - xx3) / ((xx1 - xx2) * (xx1 - xx3))
	       + yy2 * (x - xx1) * (x - xx3) / ((xx2 - xx1) * (xx2 - xx3))
	       + yy3 * (x - xx1) * (x - xx2) / ((xx3 - xx1) * (xx3 - xx2));
}
