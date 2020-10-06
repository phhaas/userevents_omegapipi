#include "RPD_Helper.h"

#include <fstream>
#include <sstream>

#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TTree.h>

#include <Phast.h>


namespace {
	const bool debug(false); // creating some histograms
	const bool draw(false); // drawing the histograms during reconstruction
	const bool verbose(false); // printing some cout info

	const int NELEMENTS[2] = {12,24}; // ring A and B's quantization
	float mP;
	float tofOffset[12][24]; // weg
	float positionOffsetForA[12], positionOffsetForB[24]; //< calibrate
	float lightSpeedInA[12], lightSpeedInB[24]; //< calibrate
	float dEBcusp[24], dEAcusp[12], zVertexOffset[12][24], calTA[12]; //< calibrate

	const float rA = 12.5;
	const float rB = 75.;
	const float deltaR = rB - rA;
	const float lightSpeed = 30.; // cm/ns
	const float deltaTmin = deltaR / lightSpeed; // minimal possible time for a
	                                             // particle to travel from
	                                             // ring A to ring B
	const float dEAcuspMeV = 8.5;
	const float dEBcuspMeV = 24.;
	const float tofGlobalOffset = 2.2; //< calibrate
	const float T0_ALL = 970.;
	// better individual t0 calibration values:
	// [ring:<A|B>][element:<0..11|0..24>][pm:<up|down>]
	// float t0_const[2][24][2];
	const float zTargetCenter = -48.; //< calibrate
	// see http://wwwcompass.cern.ch/compass/gpd/meetings/2008RPD/NdH_forRPD2008.ppt [1]
	// for details of reconstruction where the element offset of
	// the upstream end is used (page 6)
	const float ELEMENTOFFSET_UP[2] = {13., 8.}; // eeg
	// ******* MC simulation part *******
	// the length of the light guides from private communication with Pawel Schneider
	// ring B has active light guides
	const double LIGHTGUIDELENGTH[2] = {21.5, 24.5};
	// last round waveguide leading to PM (not active for both rings)
	const double LIGHTGUIDELENGTHPM[2] = {6.8+3.5, 3.5};
	// http://wwwcompass.cern.ch/compass/collaboration/2008/co_0811/pdf/burtin_2008-11-20.pdf [2] page 5
	const double ELEMENTLENGTH[2] = {50. +LIGHTGUIDELENGTH[0]+LIGHTGUIDELENGTHPM[0],
								  106.+LIGHTGUIDELENGTH[0]+LIGHTGUIDELENGTHPM[0]}; // ring A and B
	// offset of the middle of the element in respect to the middle of the target
	//const double ELEMENTOFFSET[2] = {2.,12.}; // [cm] as in reference [2]
	const double ELEMENTOFFSET[2] = {2.-15.,12.-20.}; // [cm] set for correct MC simulation
	//const double tdcresolution = 0.058856; // [ns/ch] already decoded
	// attenuation length [cm] {140.,140.} from j.bernhards measurements (diploma thesis), for Ring B possibly wrong
	const double ATTENUATIONLENGTH[2] = {64., 68.}; // from [1] page 9
	// attenuation length correction factor for ring A or B
	const double CORRECTION_FACTOR[2] = {(1./TMath::Exp(-ELEMENTLENGTH[0]/2./ATTENUATIONLENGTH[0])),
			(1./TMath::Exp(-ELEMENTLENGTH[1]/2./ATTENUATIONLENGTH[1]))};
	// measured time resolution from [1] page 7 and 8
	const double TIME_RESOLUTION[2] = {0.209, 0.400}; // ns sigma
	//const double TIME_RESOLUTION[2] = {0., 0.}; // no smearing
	// so we obtain the time smearing factor for 1 PM
	// under the ussumption sigma² = sigma_pmup²+sigma_pmdown²
	// for sigma_pmup = sigma_pmdown
	const double TIME_SMEARING[2] = {TIME_RESOLUTION[0]/TMath::Sqrt(2.), TIME_RESOLUTION[1]/TMath::Sqrt(2.)};

	// database for RPD calibration values
	//std::map<int, Calib_const> Calib_constants_db;
	// current T0's
	Calib_const rpd_constants;
	// last checked runumber
	int last_DB_runnumb;

	// DB as text file for t0's (getting replaced by const_db)
	const float const_t0[] = {
		#include "RPD_DB/t0_calib.db"
	};

	// DB as text file for all calibration constants
	const float const_db[] = {
		#include "RPD_DB/calib.db"
	};
}


// initialization of singleton pointer
RPD* RPD::instance = NULL;


// usage: RPD& rpd = RPD::Instance();
RPD& RPD::Instance(){
	if (!instance){
		// create new instance if no existent
		instance = new RPD();
	}
	return *instance;
}


std::string
RPD::ChannelToPMT(const std::string& RpdPlanes, const int chId, int &thisPMT) const
{
	std::stringstream nameIs;
	thisPMT = -1;

	if (RpdPlanes == "RP01TA__") {
		if (chId >= 0 && chId < 12) { // Au0 -Au11
			thisPMT = chId;
			nameIs << "TAu" << thisPMT;
		} else if (chId >= 32 && chId < 44) { // Ad0 -Ad11
			thisPMT = chId - 32;
			nameIs << "TAd" << thisPMT;
		} else {
			// unknown A channel
			thisPMT = -1;
			nameIs << "Unkn";
		}
	} else if (RpdPlanes == "RP01TBl_") {
		if (chId >= 0 && chId < 12) { // Bu0 -Bu11
			thisPMT = chId;
			nameIs << "TBu" << thisPMT;
		} else if (chId >= 16 && chId < 28) { // Bu12 -Bu23
			thisPMT = 12 + chId - 16;
			nameIs << "TBu" << thisPMT;
		} else if (chId >= 32 && chId < 44) { // Bd0 -Bd11
			thisPMT = chId - 32;
			nameIs << "TBd" << thisPMT;
		} else if (chId >= 48 && chId < 60) { // Bd12 -Bd23
			thisPMT = 12 + chId - 48;
			nameIs << "TBd" << thisPMT;
		} else {
			// unknown B channel
			thisPMT = -1;
			nameIs << "Unkn";
		}
	} else if (RpdPlanes == "RP01Ql__") {
		if (chId >= 0 && chId < 12) { // Au0- Au11
			thisPMT = chId;
			nameIs << "QAu" << thisPMT;
		} else if (chId >= 12 && chId < 24) { // Ad0- Ad11
			thisPMT = chId - 12;
			nameIs << "QAd" << thisPMT;
		} else if (chId >= 24 && chId < 48) { // Bu0- Bu23
			thisPMT = chId - 24;
			nameIs << "QBu" << thisPMT;
		} else if (chId >= 48 && chId < 72) { // Bd0- Bd23
			thisPMT = chId - 24 - 24;
			nameIs << "QBd" << thisPMT;
		} else {
			// unknown SADC channel
			thisPMT = -1;
			nameIs << "Unkn";
		}
	} else if (RpdPlanes == "RP01Qh__") {
		if (chId >= 0 && chId < 12) { // Au0- Au11
			thisPMT = chId;
			nameIs << "QAu" << thisPMT;
		} else if (chId >= 12 && chId < 24) { // Ad0- Ad11
			thisPMT = chId - 12;
			nameIs << "QAd" << thisPMT;
		} else if (chId >= 24 && chId < 48) { // Bu0- Bu23
			thisPMT = chId - 24;
			nameIs << "QBu" << thisPMT;
		} else if (chId >= 48 && chId < 72) { // Bd0- Bd23
			thisPMT = chId - 24 - 24;
			nameIs << "QBd" << thisPMT;
		} else {
			// unknown SADC channel
			thisPMT = -1;
			nameIs << "Unkn"; // unknown SADC channel
		}
	} else {
		thisPMT = -1;
		nameIs << "Unkn";
	}

	return nameIs.str();
}


std::string
RPD::GetChannel(const int ring, const int element, const bool tdc, int &channel_pmup, int &channel_pmdo) const
{
	std::string result;
	channel_pmup = -1;
	channel_pmdo = -1;
	// ************** tdc part *************
	if (tdc){
		if (ring == 0){
			result = "RP01TA__";
			channel_pmup = element;	// channel 0..11
			channel_pmdo = element +32;// channel 32..43
		}
		else {
			result = "RP01TBl_";
			// channel 0..11 for element 0..11 and an offset of 4 for element 12..23
			// thus channels 16..27
			channel_pmup = (element < 12) ? element 	 : element+4;
			// channel 32..43 for element 0..11 and channel 48..59 for element 12..23
			channel_pmdo = (element < 12) ? element+32 : element+36;
		}
	}
	// ************** adc part *************
	else {
		result = "RP01Ql__";
		if (ring == 0){
			channel_pmup = element;		// channel 0..11
			channel_pmdo = element + 12;  	// channel 12..24
		}
		else {
			channel_pmup = element + 24;	// channel 24..48
			channel_pmdo = element + 48;	// channel 48..72
		}
	}
	return result;
}


double
RPD::correctEnergyLoss(double pRpd, const double sinTheta, const double cosTheta, const double phi, const double xVertex, const double yVertex, const double zVertex, const bool debug) const
{
		const Phast& phast = Phast::Ref();
		int runnb = phast.event.RunNum();
		const double rTarget = 17.5;
		const double chaussetteThickness = 0.2;
		const double aluminiumThickness = 1.8;
		const double ringAThickness = 5.;
		const double cellThickness = 0.125;
		double lengthInH2;
		pRpd = correction_energyPol(pRpd, 4, ringAThickness / sinTheta, debug); // scint A
		pRpd = correction_energyPol(pRpd, 3, chaussetteThickness / sinTheta, debug); // mylar
		if (runnb < 80608 || runnb > 81093) { //begin of lead data run 80608 (16.10.2009 13:39), ended with run 81093 (28.10.2009 07:56)
		pRpd = correction_energyPol(pRpd, 2, aluminiumThickness / sinTheta, debug); // aluminium
		pRpd = correction_energyPol(pRpd, 3, cellThickness / sinTheta, debug); // mylar
		lengthInH2 = TMath::Sqrt(pow(rTarget * TMath::Cos(phi) - xVertex, 2) + pow(rTarget * TMath::Sin(phi) - yVertex, 2));
		lengthInH2 = lengthInH2 / sinTheta;
		pRpd = correction_energyPol(pRpd, 1, lengthInH2, debug); // hydrogen
		}
		else {
		// materials: 5...lead, 6...tungsten
		int target_material[16]={5,5,5,5,5,5,5,5,5,5,6,6,5,5,6,6};
		double target_thickness[16]={0.125,0.125,0.125,0.125,0.125,0.125,0.0625,0.0625,0.025,0.0125,0.025,0.0125,0.025,0.0125,0.025,0.0125};
		double target_poscenter[18]={-100.,0.,24.182,47.582,72.182,96.182,120.882,144.1195,168.1195,192.082,216.0695,240.0820,264.0695,288.082,312.0695,336.082,360.0695,400.};
		int target_nb = -1;
		for (int i=1; i<17; i++) {
			if ((target_poscenter[i-1]+target_poscenter[i])/2. < (zVertex + 66.7) && (target_poscenter[i]+target_poscenter[i+1])/2. > (zVertex + 66.7)) {
			target_nb=i-1;
			break;
			} else {
			return (std::isnan(pRpd) || std::isinf(pRpd) || pRpd < 0) ? 0 : pRpd;
			}
		} // a lot of hard-coded stuff here, no good documentation of lead target, measured with vertex distributions 
		//double length = TMath::Sqrt(pow(1.25 * TMath::Cos(phi) - xVertex, 2) + pow(1.25 * TMath::Sin(phi) - yVertex, 2));
		//length = length / sinTheta; //<- first try, forcing the vertex in the target center and using half the target thickness gives better results
		double length = target_thickness[target_nb] / (2. * cosTheta);
		pRpd = correction_energyPol(pRpd, target_material[target_nb], length, debug);
		}
		return (std::isnan(pRpd) || std::isinf(pRpd) || pRpd < 0) ? 0 : pRpd;
	}


double
RPD::correction_energyPol(const double p, const int mat, const double thickness, const bool /*debug*/) const
{
	static TF1* rangeFunction = new TF1("rangeFunction", "(((x<=50)*pol3(0))+((x>50)*((x<500)*pol4(4))))+((x>=500)*pol2(9))");
	static double par_Aluminium[12] = {
	              0.0, 0.0265061853053, 0.0038016393958, -1.24938207462e-06,
	             -2.15740719589, 0.109089884982, 0.00309037606361, -2.86225813051e-06,
	              1.27840792335e-09, -120.238187882, 0.977375559516, 0.000714600752871};
	static double par_H2[12] = {
	              0.0, 0.244101692132, 0.0526487276028, -1.3402701485e-05,
	              -25.9373985022, 1.21875716039, 0.0448317358735, -4.03128421537e-05,
	              1.77395193779e-08, -1738.88278436, 13.7138693978, 0.0109717517097};
	static double par_Lead[12] = {
	              0.0, 0.00666252555628, 0.000908999694029, -3.03111361196e-07,
	              -0.519629317233, 0.0265914384531, 0.000736149016099, -6.7564156876e-07,
	              3.00262905778e-10, -28.6755305789, 0.233059709245, 0.000173081021626};
	static double par_Mylar[12] = {
	              0.0, 0.0369143441703, 0.00614414918572, -1.41754876221e-06,
	              -3.40696750414, 0.165662144104, 0.00508770150398, -4.63404128079e-06,
	              2.04132425803e-09, -199.913829222, 1.59986884667, 0.00119862596984};
	static double par_Scintillator[12] = {
	              0.0, 0.0448832760685, 0.00786238109716, -2.36168110019e-06,
	              -4.21884944032, 0.205344172788, 0.00651313737621, -5.94182623541e-06,
	              2.62912777891e-09, -253.985670322, 2.03258607044, 0.00154408969146};
	static double par_Tungsten[12] = {
	              0.0, 0.0173610957949, 0.00145011007172, -1.35044853985e-06,
	              -0.870415711979, 0.0531812877676, 0.00106290615147, -9.93219313867e-07,
	              4.41961714127e-10, -42.1550085467, 0.356134059042, 0.000236019751861};
	static double* parArray[6] = {par_H2, par_Aluminium, par_Mylar, par_Scintillator, par_Lead, par_Tungsten};
	rangeFunction->SetParameters(parArray[mat-1]);

	double E = (TMath::Sqrt(p * p + mP * mP) - mP) * 1000.;

	double range = rangeFunction->Eval(E) + thickness;
	double energy = rangeFunction->GetX(range, 0, 1e6) / 1000.;

	return  TMath::Sqrt(energy * energy + 2 * energy * mP);
}


double
RPD::correction_energy(const double p, const int mat, const double thickness, const bool debug) const
{
		const int nData = 151;

		// --- kinetic energy in MeV

		double Ek[nData] = { 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0,
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
				584.0, 588.0, 592.0, 596.0, 600.0 };

		// --- protons range in H2 mm

		double range_H2[nData] = { 0.1373E+01, 0.4879E+01, 0.1030E+02,
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
				0.1033E+05, 0.1044E+05 };

		double range_Aluminium[nData] = { 0.1290E+00, 0.4173E+00, 0.8474E+00,
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
				0.7160E+03, 0.7233E+03 };

		// --- protons range in Mylar mm

		double range_Mylar[nData] = { 0.1851E+00, 0.6286E+00, 0.1300E+01,
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
				0.1179E+04, 0.1191E+04 };

		// --- protons range in Scintillator mm

		double range_Scintillator[nData] = { 0.2301E+00, 0.7883E+00, 0.1636E+01,
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
				0.1506E+04, 0.1521E+04 };
				
		double range_Lead[nData] = {0.0313921,0.102026,0.206344,0.341674,0.506432,0.699295,0.918943,1.16564,1.437,1.73392,2.05463,2.39912,2.76652,3.15595,3.56828,4.00176,4.45639,4.93216,5.42819,5.94449,6.48106,7.03612,7.61057,8.20352,8.81938,9.44493,10.0881,10.7577,11.4361,12.141,12.8546,13.5859,14.3348,15.1013,15.8855,16.6784,17.489,18.3172,19.163,20.0176,20.8899,21.7709,22.6696,23.5859,24.511,25.4537,26.4053,27.3744,28.3524,29.348,30.3524,31.3744,32.4053,33.4449,34.5022,35.5683,36.6432,37.7357,38.837,39.9471,41.0749,42.2026,43.348,44.511,45.674,46.8546,48.0441,49.2423,50.4493,51.6652,52.8899,54.1322,55.3833,56.6344,57.9031,59.1806,60.467,61.7533,63.0573,64.37,65.6916,67.022,68.3612,69.7093,71.0573,72.4229,73.7974,75.1718,76.5639,77.9559,79.3568,80.7665,82.185,83.6123,85.0485,86.4846,87.9383,89.4273,90.837,92.3348,93.8326,95.2423,96.7401,98.2379,99.7357,101.322,102.819,104.317,105.815,107.401,108.899,110.485,112.07,113.568,115.154,116.74,118.326,119.912,121.498,123.084,124.67,126.256,127.93,129.515,131.189,132.775,134.449,136.035,137.709,139.383,140.969,142.643,144.317,145.991,147.665,149.339,151.013,152.687,154.449,156.123,157.797,159.559,161.233,162.996,164.67,166.432,168.106,169.868,171.63,173.392};
		
		double range_Tungsten[nData] = {0.0743877,0.214537,0.408106,0.65022,0.937445,1.26784,1.63877,2.04934,2.4978,2.98238,3.5022,4.05639,4.64405,5.26432,5.9163,6.6,7.31366,8.05815,8.82819,9.62996,10.467,11.3216,12.2115,13.1189,14.0617,15.0308,16.0176,17.0396,18.0793,19.1454,20.2379,21.348,22.4846,23.6476,24.8282,26.0352,27.2687,28.511,29.7885,31.0749,32.3877,33.7269,35.0749,36.4493,37.8414,39.2511,40.6872,42.1322,43.6035,45.0837,46.5903,48.1145,49.6564,51.207,52.7841,54.37,55.9824,57.6035,59.2423,60.8987,62.5727,64.2555,65.9648,67.6828,69.4097,71.1542,72.9163,74.696,76.4846,78.2907,80.1057,81.9383,83.7885,85.6476,87.5154,89.4273,91.2775,93.2159,95.1542,97.0925,99.0308,100.969,102.907,104.934,106.96,108.899,110.925,112.952,115.066,117.093,119.119,121.233,123.26,125.374,127.489,129.604,131.718,133.833,136.035,138.15,140.352,142.555,144.67,146.872,149.075,151.278,153.568,155.771,157.974,160.264,162.467,164.758,167.048,169.339,171.63,173.921,176.211,178.502,180.881,183.172,185.551,187.841,190.22,192.599,194.978,197.357,199.736,202.115,204.493,206.872,209.339,211.718,214.185,216.564,219.031,221.498,223.877,226.344,228.811,231.278,233.744,236.3,238.767,241.233,243.789,246.256,248.722,251.278,253.833,256.3};

		double E = (TMath::Sqrt(p * p + mP * mP) - mP) * 1000.;
		double beta = p / TMath::Sqrt(p * p + mP * mP);

		double range = 0.;
		double energy = 0.;
		double a = 0., b = 0., c = 0.;
		//int mat=3;

		// Search for first tabulated energy above the proton's energy.
		int i;
		for (i = 1; i < nData - 1; i++)
			if (E < (double) Ek[i])
				break;

		// If the energy is no longer tabulated, assume no energy loss.
		if (i == nData - 1)
			return p;

		// E is inside table range, calculate correction.
		if (mat == 1) {
			if (debug)
				printf("***********************\n in hydrogen\n");
			a = range_H2[i - 1];
			b = range_H2[i];
			c = range_H2[i + 1];
		} else if (mat == 2) {
			if (debug)
				printf("***********************\n in Aluminium\n");
			a = range_Aluminium[i - 1];
			b = range_Aluminium[i];
			c = range_Aluminium[i + 1];
		} else if (mat == 3) {
			if (debug)
				printf("***********************\n in Mylar\n");
			a = range_Mylar[i - 1];
			b = range_Mylar[i];
			c = range_Mylar[i + 1];
		} else if (mat == 4) {
			if (debug)
				printf("***********************\n in Scintillator\n");
			a = range_Scintillator[i - 1];
			b = range_Scintillator[i];
			c = range_Scintillator[i + 1];
		} else if (mat == 5) {
			if (debug)
				printf("***********************\n in Lead\n");
			a = range_Lead[i - 1];
			b = range_Lead[i];
			c = range_Lead[i + 1];
		} else if (mat == 6) {
			if (debug)
				printf("***********************\n in Tungsten\n");
			a = range_Tungsten[i - 1];
			b = range_Tungsten[i];
			c = range_Tungsten[i + 1];
		}
		if (debug)
			printf( "p=%f beta=%f E=%3.2f Ekin=%3.2f,%3.2f,%3.2f range=%4.2f,%4.2f,%4.2f\n", p * 1000., beta, E, Ek[i - 1], Ek[i], Ek[i + 1], a, b, c);
		range = terpol3(Ek[i - 1], Ek[i], Ek[i + 1], a, b, c, E);
		if (debug)
			printf("interpolated range = %4.2f \n", range);
		range += thickness;
		energy = terpol3(a, b, c, Ek[i - 1], Ek[i], Ek[i + 1], range);
		if (debug) printf("energy:%f p=%f\n", energy, TMath::Sqrt(energy * energy + 2 * energy * mP * 1000.));

		return (std::isnan(energy) || std::isinf(energy) || energy < 0) ? 0 : TMath::Sqrt(energy * energy + 2 * energy * mP * 1000.) / 1000.;
	}


double
RPD::terpol3(const double xx1, const double xx2, const double xx3, const double yy1, const double yy2, const double yy3, const double x) const
{
		double val = yy1 * (x - xx2) * (x - xx3) / ((xx1 - xx2) * (xx1 - xx3))
			  + yy2 * (x - xx1) * (x - xx3) / ((xx2 - xx1) * (xx2 - xx3))
			  + yy3 * (x - xx1) * (x - xx2) / ((xx3 - xx1) * (xx3 - xx2));
		return val;
	}


float
RPD::treatADC(const PaDigit& d) const
{
		//calculate pedestal
		float fPedestal = 0.;
		float maxAmplitude = 0.;
		float sumAmplitude = 0.;
		int nSamples = 32;
		int samples = 0;
		for (int i = 2; i < 6; i++) {
			fPedestal += d.DigInfo(i);
		}
		fPedestal = fPedestal / 4.;
		//analyze waveform
		for (int i = 6; i < nSamples + 2; i++) {
			if (d.DigInfo(i) > maxAmplitude)
				maxAmplitude = d.DigInfo(i);
			sumAmplitude += d.DigInfo(i);
			samples++;
		}
		return sumAmplitude - samples * fPedestal;
	}


RPD::RPD()
{
		Load_DB();  // move this behind the next line if map will become global

		mP = 0.938272;

		tofOffset[0][23] = 4.44;
		positionOffsetForA[0] = 5.59;
		tofOffset[0][0] = -0.64;
		positionOffsetForB[0] = 26.52;
		tofOffset[0][1] = -0.78;
		positionOffsetForB[1] = 13.06;

		tofOffset[1][1] = -1.97;
		positionOffsetForA[1] = 17.62;
		tofOffset[1][2] = -2.55;
		positionOffsetForB[2] = 24.81;
		tofOffset[1][3] = -2.37;
		positionOffsetForB[3] = 28.24;

		tofOffset[2][3] = -0.10;
		positionOffsetForA[2] = 48.34;
		tofOffset[2][4] = -0.86;
		positionOffsetForB[4] = 39.47;
		tofOffset[2][5] = 1.35;
		positionOffsetForB[5] = 35.72;

		tofOffset[3][5] = -3.74;
		positionOffsetForA[3] = 30.46;
		tofOffset[3][6] = -4.97;
		positionOffsetForB[6] = 11.37;
		tofOffset[3][7] = -5.70;
		positionOffsetForB[7] = 15.20;

		tofOffset[4][7] = -4.92;
		positionOffsetForA[4] = 41.52;
		tofOffset[4][8] = -4.09;
		positionOffsetForB[8] = 49.74;
		tofOffset[4][9] = -0.68;
		positionOffsetForB[9] = 14.60;

		tofOffset[5][9] = -0.44;
		positionOffsetForA[5] = 52.19;
		tofOffset[5][10] = -2.27;
		positionOffsetForB[10] = 20.41;
		tofOffset[5][11] = 0.40;
		positionOffsetForB[11] = 3.78;

		tofOffset[6][11] = 0.72;
		positionOffsetForA[6] = 26.58;
		tofOffset[6][12] = -4.99;
		positionOffsetForB[12] = 25.08;
		tofOffset[6][13] = -0.82;
		positionOffsetForB[13] = 42.70;

		tofOffset[7][13] = -2.75;
		positionOffsetForA[7] = 24.33;
		tofOffset[7][14] = -4.92;
		positionOffsetForB[14] = -0.18;
		tofOffset[7][15] = -3.33;
		positionOffsetForB[15] = 17.47;

		tofOffset[8][15] = 2.39;
		positionOffsetForA[8] = 21.59;
		tofOffset[8][16] = -0.23;
		positionOffsetForB[16] = -10.47;
		tofOffset[8][17] = 5.51;
		positionOffsetForB[17] = 27.07;

		tofOffset[9][17] = 0.92;
		positionOffsetForA[9] = 4.85;
		tofOffset[9][18] = -3.50;
		positionOffsetForB[18] = 17.68;
		tofOffset[9][19] = -2.11;
		positionOffsetForB[19] = 18.48;

		tofOffset[10][19] = -2.16;
		positionOffsetForA[10] = 13.08;
		tofOffset[10][20] = -2.62;
		positionOffsetForB[20] = 1.65;
		tofOffset[10][21] = 0.32;
		positionOffsetForB[21] = 9.61;

		//tofOffset[11][21] = -1.3;  positionOffsetForA[11]=7.31;
		tofOffset[11][21] = -1.01;
		positionOffsetForA[11] = 7.31;
		tofOffset[11][22] = -5.85;
		positionOffsetForB[22] = -0.17;
		tofOffset[11][23] = -0.78;
		positionOffsetForB[23] = -17.07;

		// effective lightspeed in the material [cm/ns]
		lightSpeedInA[0] = 12.;
		lightSpeedInA[1] = 12.;
		lightSpeedInA[2] = 11.5;
		lightSpeedInA[3] = 11.5;
		lightSpeedInA[4] = 12.;
		lightSpeedInA[5] = 12.;
		lightSpeedInA[6] = 12.;
		lightSpeedInA[7] = 11.;
		lightSpeedInA[8] = 11.;
		lightSpeedInA[9] = 12.;
		lightSpeedInA[10] = 12.;
		lightSpeedInA[11] = 12.;

		for (unsigned int i = 0; i < (unsigned int) (sizeof(lightSpeedInB) / sizeof(float)); i++) {lightSpeedInB[i] = 13.;}

		// dEBcusp[ 0]= 2280.; //old calibration
		// dEBcusp[ 1]=  2150.;
		// dEBcusp[ 2]= 1830.;
		// dEBcusp[ 3]=  2700.;
		// dEBcusp[ 4]= 2180.;
		// dEBcusp[ 5]=  2100.;
		// dEBcusp[ 6]= 2340.;
		// dEBcusp[ 7]=  2000.;
		// dEBcusp[ 8]= 2270.;
		// dEBcusp[ 9]=  2300.;
		// dEBcusp[10]= 2280.;
		// dEBcusp[11]=  2050.;
		// dEBcusp[12]= 2165.;
		// dEBcusp[13]=  2200.;
		// dEBcusp[14]= 2250.;
		// dEBcusp[15]=  2100.;
		// dEBcusp[16]= 2430.;
		// dEBcusp[17]=  2000.;
		// dEBcusp[18]= 2440.;
		// dEBcusp[19]=  2000.;
		// dEBcusp[20]= 2360.;
		// dEBcusp[21]=  2400.;
		// dEBcusp[22]= 2340.;
		// dEBcusp[23]=  2000.;

		dEBcusp[0] = 2403.;
		dEBcusp[1] = 2117.;
		dEBcusp[2] = 1886.;
		dEBcusp[3] = 2690.2;
		dEBcusp[4] = 2244.;
		dEBcusp[5] = 2072.;
		dEBcusp[6] = 2409.;
		dEBcusp[7] = 2022.;
		dEBcusp[8] = 2318.;
		dEBcusp[9] = 2330.;
		dEBcusp[10] = 2317.;
		dEBcusp[11] = 2047.;
		dEBcusp[12] = 2238.;
		dEBcusp[13] = 2247.;
		dEBcusp[14] = 2332.;
		dEBcusp[15] = 2242.;
		dEBcusp[16] = 2490.;
		dEBcusp[17] = 2052.;
		dEBcusp[18] = 2549.;
		dEBcusp[19] = 2050.;
		dEBcusp[20] = 2495.;
		dEBcusp[21] = 2368.;
		dEBcusp[22] = 2371.;
		dEBcusp[23] = 2048.;

		dEAcusp[0] = 1470.;
		dEAcusp[1] = 1280.;
		dEAcusp[2] = 1100.;
		dEAcusp[3] = 1280.;
		dEAcusp[4] = 1190.;
		dEAcusp[5] = 1260.;
		dEAcusp[6] = 1260.;
		dEAcusp[7] = 1480.;
		dEAcusp[8] = 800.;
		dEAcusp[9] = 1370.;
		dEAcusp[10] = 1230.;
		dEAcusp[11] = 1240.;

		zVertexOffset[0][23] = -1.80;
		zVertexOffset[0][0] = -2.40;
		zVertexOffset[0][1] = -3.12;
		zVertexOffset[1][1] = -3.31;
		zVertexOffset[1][2] = -3.59;
		zVertexOffset[1][3] = -3.84;
		zVertexOffset[2][3] = -4.26;
		zVertexOffset[2][4] = -4.43;
		zVertexOffset[2][5] = -4.71;
		zVertexOffset[3][5] = -3.83;
		zVertexOffset[3][6] = -3.75;
		zVertexOffset[3][7] = -3.36;
		zVertexOffset[4][7] = -1.88;
		zVertexOffset[4][8] = -1.56;
		zVertexOffset[4][9] = -1.04;
		zVertexOffset[5][9] = 0.53;
		zVertexOffset[5][10] = 1.07;
		zVertexOffset[5][11] = 2.13;
		zVertexOffset[6][11] = 2.49;
		zVertexOffset[6][12] = 2.94;
		zVertexOffset[6][13] = 3.44;
		zVertexOffset[7][13] = 4.78;
		zVertexOffset[7][14] = 4.69;
		zVertexOffset[7][15] = 4.67;
		zVertexOffset[8][15] = 4.62;
		zVertexOffset[8][16] = 4.73;
		zVertexOffset[8][17] = 4.47;
		zVertexOffset[9][17] = 5.00;
		zVertexOffset[9][18] = 4.58;
		zVertexOffset[9][19] = 4.35;
		zVertexOffset[10][19] = 3.63;
		zVertexOffset[10][20] = 3.07;
		zVertexOffset[10][21] = 2.56;
		zVertexOffset[11][21] = 0.77;
		zVertexOffset[11][22] = 0.10;
		zVertexOffset[11][23] = -0.62;

		calTA[ 0] = 5.;
		calTA[ 1] = 4.;
		calTA[ 2] = 6.;
		calTA[ 3] = 1.;
		calTA[ 4] = 1.75;
		calTA[ 5] = 1.95;
		calTA[ 6] = 2.2;
		calTA[ 7] = 0.75;
		calTA[ 8] = 6.;
		calTA[ 9] = 1.7;
		calTA[10] = 1.4;
		calTA[11] = -0.3;

		/* moved to constants declaration in the namespace above
		rA = 12.5;
		rB = 75.;
		deltaR = rB - rA;
		lightSpeed = 30.; // cm/ns
		dEAcuspMeV = 8.5;
		dEBcuspMeV = 24.;
		tofGlobalOffset = 2.2;
		zTargetCenter = -48.;*/

		fLastDecodedEvent = -1;
		fLastClusteredEvent = -1;
		fLastDT0SimulatedEvent = -1;
}


void
RPD::Reset()
{
	fLastDecodedEvent = -1;

	fBestMCProtonTrack = -1;
	fRPDMCTrack.clear();
	dEA_MC.clear();
	dEB_MC.clear();

	memset(fDataSentAupt, 0, sizeof(fDataSentAupt));
	memset(fDataSentAdot, 0, sizeof(fDataSentAdot));
	memset(fDataSentBupt, 0, sizeof(fDataSentBupt));
	memset(fDataSentBdot, 0, sizeof(fDataSentBdot));

	memset(fDataSentAupql, 0, sizeof(fDataSentAupql));
	memset(fDataSentAdoql, 0, sizeof(fDataSentAdoql));
	memset(fDataSentBupql, 0, sizeof(fDataSentBupql));
	memset(fDataSentBdoql, 0, sizeof(fDataSentBdoql));
	memset(fDataSentAupqh, 0, sizeof(fDataSentAupqh));
	memset(fDataSentAdoqh, 0, sizeof(fDataSentAdoqh));
	memset(fDataSentBupqh, 0, sizeof(fDataSentBupqh));
	memset(fDataSentBdoqh, 0, sizeof(fDataSentBdoqh));

	memset(fAdcAupl, 0, sizeof(fAdcAupl));
	memset(fAdcAdol, 0, sizeof(fAdcAdol));
	memset(fAdcBupl, 0, sizeof(fAdcBupl));
	memset(fAdcBdol, 0, sizeof(fAdcBdol));
	memset(fAdcAuph, 0, sizeof(fAdcAuph));
	memset(fAdcAdoh, 0, sizeof(fAdcAdoh));
	memset(fAdcBuph, 0, sizeof(fAdcBuph));
	memset(fAdcBdoh, 0, sizeof(fAdcBdoh));

	memset(fTdcAup, 0, sizeof(fTdcAup));
	memset(fTdcAdo, 0, sizeof(fTdcAdo));
	memset(fTdcBup, 0, sizeof(fTdcBup));
	memset(fTdcBdo, 0, sizeof(fTdcBdo));

	// always reset the other stuff as well
	ResetClustering();
	ResetDT0Simulation();
}


void
RPD::ResetClustering()
{
	fLastClusteredEvent = -1;

	fBestProtonTrack = -1;
	fRPDTrack.clear();
	fzT.clear();
	ftT.clear();
	dEA.clear();
	dEB.clear();
	zA_vec.clear();
	zB_vec.clear();
	fCalibratedEvent.clear();
	hits.clear();
	nTracks = 0;

	something = false;
}


void
RPD::ResetDT0Simulation()
{
	fLastDT0SimulatedEvent = -1;

	fDT0Components.clear();
}


void
RPD::Dump() const
{
	std::cout << std::endl << std::endl
	          << "######################################### RPD Event Dump #########################################" << std::endl
	          << std::endl;

	if (fLastDecodedEvent != -1) {
		std::cout << "======>  Last event decoded " << fLastDecodedEvent << std::endl
		          << std::endl;

		for (int iA = 0; iA < 12; iA++) {
			if (fDataSentAupt[iA] > 0 || fDataSentAdot[iA] > 0 ||
			    fDataSentAupql[iA] || fDataSentAdoql[iA] ||
			    fDataSentAupqh[iA] || fDataSentAdoqh[iA]) {
				std::cout << "ring A element " << iA << ":" << std::endl;
				if (fDataSentAupt[iA] > 0) {
					std::cout << "upstream TDC info:";
					for (int iHitAup = 0; iHitAup < fDataSentAupt[iA]; iHitAup++)
						std::cout << " " << fTdcAup[iA][iHitAup];
					std::cout << std::endl;
				}
				if (fDataSentAdot[iA] > 0) {
					std::cout << "downstream TDC info:";
					for (int iHitAdo = 0; iHitAdo < fDataSentAdot[iA]; iHitAdo++)
						std::cout << " " << fTdcAdo[iA][iHitAdo];
					std::cout << std::endl;
				}
				if (fDataSentAupql[iA] || fDataSentAupqh[iA]) {
					std::cout << "upstream ADC info:";
					if (fDataSentAupql[iA])
						std::cout << " low range: " << fAdcAupl[iA];
					if (fDataSentAupqh[iA])
						std::cout << " high range: " << fAdcAuph[iA];
					std::cout << std::endl;
				}
				if (fDataSentAdoql[iA] || fDataSentAdoqh[iA]) {
					std::cout << "downstream ADC info:";
					if (fDataSentAdoql[iA])
						std::cout << " low range: " << fAdcAdol[iA];
					if (fDataSentAdoqh[iA])
						std::cout << " high range: " << fAdcAdoh[iA];
					std::cout << std::endl;
				}
				std::cout << std::endl;
			}
		}

		for (int iB = 0; iB < 24; iB++) {
			if (fDataSentBupt[iB] > 0 || fDataSentBdot[iB] > 0 ||
			    fDataSentBupql[iB] || fDataSentBdoql[iB] ||
			    fDataSentBupqh[iB] || fDataSentBdoqh[iB]) {
				std::cout << "ring B element " << iB << ":" << std::endl;
				if (fDataSentBupt[iB] > 0) {
					std::cout << "upstream TDC info:";
					for (int iHitBup = 0; iHitBup < fDataSentBupt[iB]; iHitBup++)
						std::cout << " " << fTdcBup[iB][iHitBup];
					std::cout << std::endl;
				}
				if (fDataSentBdot[iB] > 0) {
					std::cout << "downstream TDC info:";
					for (int iHitBdo = 0; iHitBdo < fDataSentBdot[iB]; iHitBdo++)
						std::cout << " " << fTdcBdo[iB][iHitBdo];
					std::cout << std::endl;
				}
				if (fDataSentBupql[iB] || fDataSentBupqh[iB]) {
					std::cout << "upstream ADC info:";
					if (fDataSentBupql[iB])
						std::cout << " low range: " << fAdcBupl[iB];
					if (fDataSentBupqh[iB])
						std::cout << " high range: " << fAdcBuph[iB];
					std::cout << std::endl;
				}
				if (fDataSentBdoql[iB] || fDataSentBdoqh[iB]) {
					std::cout << "downstream ADC info:";
					if (fDataSentBdoql[iB])
						std::cout << " low range: " << fAdcBdol[iB];
					if (fDataSentBdoqh[iB])
						std::cout << " high range: " << fAdcBdoh[iB];
					std::cout << std::endl;
				}
				std::cout << std::endl;
			}
		}
	}

	if (fLastClusteredEvent != -1) {
		std::cout << "======>  Last event reconstructed " << fLastClusteredEvent << std::endl
		          << std::endl;

		std::cout << "number of reconstructed recoil particles: " << nTracks << std::endl
		          << std::endl;

		if (nTracks > 0) {
			for (int i=0; i<nTracks; ++i) {
				std::cout << "Track " << i << " Momentum : " << fRPDTrack[i].Vect().Mag();
				if (i == fBestProtonTrack)
					std::cout << "     (best proton track)";
				std::cout << std::endl
				          << "        Px       : " << fRPDTrack[i].Px()    << std::endl
				          << "        Py       : " << fRPDTrack[i].Py()    << std::endl
				          << "        Pz       : " << fRPDTrack[i].Pz()    << std::endl
				          << "        E        : " << fRPDTrack[i].E()     << std::endl
				          << "        Pt       : " << fRPDTrack[i].Pt()    << std::endl
				          << "        Theta    : " << fRPDTrack[i].Theta() << std::endl
				          << "        Phi      : " << fRPDTrack[i].Phi()   << std::endl
				          << "        Beta     : " << fRPDTrack[i].Beta()  << std::endl
				          << "        PV_Z     : " << fzT[i]               << std::endl
				          << "        PV_T     : " << ftT[i]               << std::endl
				          << "        dE(A)    : " << dEA[i]               << std::endl
				          << "        dE(B)    : " << dEB[i]               << std::endl
				          << "        Pos(A)   : " << zA_vec[i]            << std::endl
				          << "        Pos(B)   : " << zB_vec[i]            << std::endl
				          << "        Hit A in Element: " << hits[i].first << std::endl
				          << "        Hit B in Element: " << hits[i].second << std::endl
				          << std::endl;
			}
		}
	}

	if (fLastDT0SimulatedEvent != -1) {
		std::cout << "======>  Last event with DT0 simulation " << fLastDT0SimulatedEvent << std::endl
		          << std::endl;

		std::cout << "Individual components:" << std::endl;
		for (std::map<DT0Components, bool>::const_iterator component=fDT0Components.begin(); component!=fDT0Components.end(); ++component) {
			switch (component->first) {
			case DT0Component_RPD:
				std::cout << "recoil detected (RPD)          : ";
				break;
			case DT0Component_aBT:
				std::cout << "alternative beam trigger (aBT) : ";
				break;
			case DT0Component_BC:
				std::cout << "beam counter (BC)              : ";
				break;
			case DT0Component_FI01X:
				std::cout << "FI01 X dynode                  : ";
				break;
			case DT0Component_BK:
				std::cout << "combination of beam killers    : ";
				break;
			case DT0Component_BK1:
				std::cout << "beam killer 1 (BK1)            : ";
				break;
			case DT0Component_BK2:
				std::cout << "beam killer 2 (BK2)            : ";
				break;
			case DT0Component_SW:
				std::cout << "Sandwich veto (SW)             : ";
				break;
			case DT0Component_HodVeto:
				std::cout << "hodoscope vetos                : ";
				break;
			}
			std::cout << component->second << std::endl;
		}

		std::cout << std::endl
		          << "complete DT0 trigger           : " << fDT0Trigger << std::endl
		          << std::endl;
	}
}


void
RPD::Decode(const PaEvent& event)
{
	// prevent decoding the same event several times to save CPU time
	// FIXME for MC the unique event number is not really unique, e.g., analyzing
	// FIXME mDSTs with only one event will not work, as the unique event number is
	// FIXME always the same (run == 1 and spill == 1 and event number == 1; the
	// FIXME first two are always true for MC)
	if (event.UniqueEvNum() == fLastDecodedEvent)
		return;

	// clear the previous event
	Reset();

	// look for new calibration constants
	FindInDB(event.RunNum());

	// get the digits and decode them
	const std::vector<PaDigit>& digits = Get_Digits(event);
	DecodeChipDigits(digits);

	// store event number
	fLastDecodedEvent = event.UniqueEvNum();
}


void
RPD::Search(const PaEvent& e, const PaVertex& vertex, const bool correct_momentum, const bool correct_angles)
{
	// decode (and simulate) current event
	Decode(e);

	// only call Clusterize for a new event, or if one of the options changed, otherwise save the time
	// FIXME the same as in Decode for unique event number applies here
	// FIXME vertex might also change
	const bool doClusterize = (e.UniqueEvNum() != fLastClusteredEvent) || (correct_momentum != fLastCorrectMomentum) || (correct_angles != fLastCorrectAngles);

	if (doClusterize) {
		// clear the previous event or the result of the previous options
		ResetClustering();

		// do the clustering
		Clusterize(vertex, correct_momentum, correct_angles);

		// store event number and Clusterize options
		fLastClusteredEvent = e.UniqueEvNum();
		fLastCorrectMomentum = correct_momentum;
		fLastCorrectAngles = correct_angles;
	}
}


void
RPD::Search(const PaVertex& vertex, const bool correct_momentum, const bool correct_angles)
{
	// Searches the data for RPD infos and calculates the tracks
	Search(Phast::Ref().event, vertex, correct_momentum, correct_angles);
}


void
RPD::Search(const PaEvent& /*e*/, const bool /*correct_momentum*/, const double /*PV_x*/, const double /*PV_y*/, const double /*PV_z*/)
{
	std::cerr<<"Warning this method is not supported anymore" <<std::endl;
	std::cerr<<"Please use Search(const PaEvent& event, const PaVertex& vertex, const bool correct_momentum, const bool correct_angles);" <<std::endl;
	std::cerr<<"to keep the functionality set correct_angles to false"<<std::endl;
	exit(1);
}


void
RPD::Search(const bool correct_momentum, const double PV_x, const double PV_y, const double PV_z)
{
	// Searches the data for RPD infos and calculates the tracks
	Search(Phast::Ref().event, correct_momentum, PV_x, PV_y, PV_z);
}


void
RPD::Calibrate(const PaEvent& e, const int pass, const std::string& /*path*/)
{
	static bool first(true);
	static Calib_const _calib_const;
	if (first){
		// retrieve the calibration constants for the current run
		FindInDB(e.RunNum());
		// copy them
		_calib_const = rpd_constants;
	}
	if (pass == 1){// first pass of calibration (for example t0 determination)

		static int counter(0);
		counter++;
		if (counter % 10000) {

		}
	}
	if (pass == 2){ // second step in calibration (using t0 calibrate something else)

		static int counter(0);
		counter++;
		if (counter % 10000) {
			// if you gathered enough statistics
			// to perform fits, etc. you may calibrate
			// your constants
			_calib_const.tofGlobalOffset = 2.2;
			// if you want you may modify the run number as well
			//_calib_const.runnumb = e.RunNum();

			// either you
			// put the modified values back to the DB
			Set_DB_entry(_calib_const);
			// and save the complete DB into a file
			Save_DB();

			// or save only this sample in a file
			Write_constants(_calib_const);
		}
	}
	if (pass == 3){ // third step in calibration (do anything what is needed)

		static int counter(0);
		counter++;
		if (counter % 10000) {

		}
	}
}


void
RPD::Calibrate(const PaEvent& e, const std::string& path)
{
	static int last_runnumb(-1);

	static TH1F* hist_t0_calib[2][24][2]={{{NULL,},},};
	if (!hist_t0_calib[0][0][0]){
		std::cout << "\nBooking t0 histograms for calibration... " << std::flush;
		for (int iring = 0; iring < 2; iring++){
			for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
				for (int ipm = 0; ipm < 2; ipm++){
					std::stringstream title;
					title << "hist_t0_calib_" << iring << "_" << ielement << "_" << ipm;
					hist_t0_calib[iring][ielement][ipm] = new TH1F(title.str().c_str(),title.str().c_str(), 200, -50-970, 50-970);
				}
			}
		}
		std::cout << "done." << std::endl << std::flush;
	}

	std::stringstream spread, spread0;
	spread  << "t0_ch_spread_run"       << e.RunNum();
	spread0 << "t0_ch_spread_sigma_run" << e.RunNum();

	static std::map<int, TH1F*> t0_ch_spread, t0_ch_spread_sigma;
	if (t0_ch_spread.empty()) t0_ch_spread[e.RunNum()] = new TH1F (spread.str().c_str(),spread.str().c_str(), 200, -50-970, 50-970);
	if (t0_ch_spread_sigma.empty()) t0_ch_spread_sigma[e.RunNum()] = new TH1F (spread0.str().c_str(),spread0.str().c_str(), 50, -0, 10);

	// draw all histograms from time to time
	static TCanvas* t0_window(NULL);
	if (!t0_window){
		t0_window = new TCanvas("t0_window_search", "RPD helper t0 calibration window", 1200, 800);
		t0_window->Divide(12,6);
		std::cout << "t0 window" << std::flush;
	}

	/*static int counter(0);
	counter++;
	if ((counter % 10000) == 0){
		int ipad(0);
		for (int iring = 0; iring < 2; iring++){
			for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
				for (int ipm = 0; ipm < 2; ipm++){
					ipad++;
					t0_window->cd(ipad);
					gPad->Clear();
					hist_t0_calib[iring][ielement][ipm]->Draw();
					gPad->Update();
				}
			}
		}
	}*/

	// create a calibration file after one run completed
	if (last_runnumb != e.RunNum()){
		if (last_runnumb != -1){
			t0_ch_spread[last_runnumb] = new TH1F (spread.str().c_str(),spread.str().c_str(), 200, -50-970, 50-970);
			t0_ch_spread_sigma[last_runnumb] = new TH1F (spread0.str().c_str(),spread0.str().c_str(), 50, -0, 10);
			std::cout << " creating t0 calibration file for run " << last_runnumb << std::endl;
			std::stringstream filename, filename0, filename1;
			filename << path << "t0_calib_const_run_" << last_runnumb << ".txt";
			filename0 << path << "t0_sigmas_run_" << last_runnumb << ".txt";
			filename1 << path << "t0_calib_const_all.txt";
			ofstream t0_calib_file;
			ofstream t0_calib_file_all;
			ofstream t0_sigma_file;
			t0_calib_file_all.open(filename1.str().c_str(), ios::app);
			t0_calib_file.open(filename.str().c_str());
			t0_sigma_file.open(filename0.str().c_str());
			t0_calib_file << last_runnumb << ", " << std::endl;
			t0_calib_file_all << last_runnumb << ", " << std::endl;
			t0_sigma_file << last_runnumb << ", " << std::endl;
			int ipad(0);
			bool reset(true);
			for (int iring = 0; iring < 2; iring++){
				for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
					for (int ipm = 0; ipm < 2; ipm++){
						ipad++;
						t0_window->cd(ipad);
						gPad->Clear();
						if (hist_t0_calib[iring][ielement][ipm]){
							hist_t0_calib[iring][ielement][ipm]->Draw();
							std::stringstream title;
							title << "gaus_" << iring << "_" << ielement << "_" << ipm;
							//float mean = hist_t0_calib[iring][ielement][ipm]->GetMean();
							float maximum = hist_t0_calib[iring][ielement][ipm]->GetBinCenter(hist_t0_calib[iring][ielement][ipm]->GetMaximumBin());
							float rms  = hist_t0_calib[iring][ielement][ipm]->GetRMS();
							TF1* gaus = new TF1(title.str().c_str(), "gaus(0)",maximum-5, maximum+5);
							gaus->SetParameters(hist_t0_calib[iring][ielement][ipm]->GetMaximum()/rms,
									maximum,
									rms);
							hist_t0_calib[iring][ielement][ipm]->Fit(gaus, "RB");
							gaus = hist_t0_calib[iring][ielement][ipm]->GetFunction(title.str().c_str());
							t0_calib_file << - gaus->GetParameter(1) << ", ";
							t0_calib_file_all << - gaus->GetParameter(1) << ", ";
							t0_sigma_file << gaus->GetParameter(2) << ", ";
							t0_ch_spread[last_runnumb]->Fill(gaus->GetParameter(1));
							t0_ch_spread_sigma[last_runnumb]->Fill(gaus->GetParameter(2));

							gPad->Update();
							if (hist_t0_calib[iring][ielement][ipm]->GetEntries() < 1000){
								reset = false;
							}
						}
					}
				}

				t0_calib_file << std::endl;
				t0_calib_file_all << std::endl;
				t0_sigma_file << std::endl;
			}
			if (reset){
				for (int iring = 0; iring < 2; iring++){
					for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
						for (int ipm = 0; ipm < 2; ipm++){
							hist_t0_calib[iring][ielement][ipm]->Reset();
						}
					}
				}
				t0_calib_file_all << std::endl << std::endl;
				//t0_ch_spread[e.RunNum()]->Reset();
				//t0_ch_spread_sigma[e.RunNum()]->Reset();
			} else {
				t0_calib_file_all << "\\*  low number of events (< 1000) !  *\\"<< std::endl;
				t0_calib_file << "\\*  low number of events (< 1000) !  *\\"<< std::endl;
			}
			t0_calib_file.close();
			t0_calib_file_all.close();
			t0_sigma_file.close();
		}
		last_runnumb = e.RunNum();
	}

	const std::vector<PaDigit>& vdd = Get_Digits(e);
	for (int i = 0; i < int(vdd.size()); i++) { // loop over DAQ digits
		const PaDigit& d = vdd[i];

		std::string detnam = d.DecodeMapName();
		int ctr = 0;
		const std::string planeName = ChannelToPMT(detnam.data(), int(d.DigInfo(0)), ctr);

		// decoding of SADC info
		/*
		if (detnam.find("RP01Ql") == 0) {
			if (planeName.find("Au") != std::string::npos) {
				fDataSentAupq[ctr] = true;
				fAdcAup[ctr] = treatADC(d);
				if (hist_qApmup) hist_qApmup->Fill(fAdcAup[ctr]);
			}
			if (planeName.find("Ad") != std::string::npos) {
				fDataSentAdoq[ctr] = true;
				fAdcAdo[ctr] = treatADC(d);
				if (hist_qApmdo) hist_qApmdo->Fill(fAdcAdo[ctr]);
			}
			if (planeName.find("Bu") != std::string::npos) {
				fDataSentBupq[ctr] = true;
				fAdcBup[ctr] = treatADC(d);
				if (hist_qBpmup) hist_qBpmup->Fill(fAdcBup[ctr]);
			}
			if (planeName.find("Bd") != std::string::npos) {
				fDataSentBdoq[ctr] = true;
				fAdcBdo[ctr] = treatADC(d);
				if (hist_qBpmdo) hist_qBpmdo->Fill(fAdcBdo[ctr]);
			}
		}*/

		// decoding of F1 info

		if (detnam.find("RP01TA") == 0) {
			if (planeName.find("Au") != std::string::npos) { //
					hist_t0_calib[0][ctr][0]->Fill(d.DigInfo(4));
			}
			if (planeName.find("Ad") != std::string::npos) { //
					hist_t0_calib[0][ctr][1]->Fill(d.DigInfo(4));
			}
		}

		if (detnam.find("RP01TBl") == 0) {
			if (planeName.find("Bu") != std::string::npos) { //
					hist_t0_calib[1][ctr][0]->Fill(d.DigInfo(4));
			}
			if (planeName.find("Bd") != std::string::npos) { //
					hist_t0_calib[1][ctr][1]->Fill(d.DigInfo(4));
			}
		}

	}// end of loop over DAQ digits
}


const std::vector<PaDigit>&
RPD::Construct_MC_Digits(const PaEvent& e)
{
	static std::vector<PaDigit> result;
	result.clear();

	static TTree* mc_values(NULL);
	static float X(-1), Y(-1), Z(-1), dE(-1), phi(-1), T(-1), beta(-1);
	static float z0(-1); // position in the slap in respect to the middle
	static float tdctime[2] = {0.,0.}; // in units of tdc channels
	static TLorentzVector best_mc_proton;

	static int ring(-1), element(-1), pm(-1);

	// remember the position of the last PaDigit pushed back to the vector
	// [ring:<A|B>][element:<0..11|0..24>][pm:<up|down>]
	int last_tdc_entry[2][24][2];
	int last_adc_entry[2][24][2];
	// initilize to -1
	memset(&last_tdc_entry, -1, 2*24*2*sizeof(int));
	memset(&last_adc_entry, -1, 2*24*2*sizeof(int));

	fRPDMCTrack.clear();
	double energy_comp = -1;

	// find the recoiling proton from mc just for comparison
	for (int ivtx = 0; ivtx < e.NMCvertex(); ivtx++){
		const PaMCvertex& mcvertex = e.vMCvertex(ivtx);
		for (int itrk = 0; itrk < mcvertex.NMCtrack(); itrk++){
			const PaMCtrack& mctrack = e.vMCtrack(mcvertex.iMCtrack(itrk));
			if (mctrack.Pid() != 14) continue;
			if (mctrack.E() > 10) continue; // this cut has to be tuned in order to identify slow protons as recoil protons
			TLorentzVector recoil_proton = mctrack.LzVec();
			//if (recoil_proton.Theta() < 0.01) continue; // set a cut on the opening angle of the reoil proton
			fRPDMCTrack.push_back(recoil_proton);
			if (energy_comp < 0 || recoil_proton.E() < energy_comp){
				fBestMCProtonTrack = fRPDMCTrack.size()-1;
			}
		}
	}
	if (fRPDMCTrack.size() > 0 && fBestMCProtonTrack > -1){
		best_mc_proton = fRPDMCTrack[fBestMCProtonTrack];
		beta = best_mc_proton.Beta();
	}

	if (debug){
		if (!mc_values){
			std::cout << " Debuging for RPD::Construct_MC_Digits() is on: creating tree " << std::endl;
			mc_values = new TTree("mc_values", "mc_values");
			mc_values->Branch("X", &X, "X/F");
			mc_values->Branch("Y", &Y, "Y/F");
			mc_values->Branch("Z", &Z, "Z/F");
			mc_values->Branch("dE", &dE, "dE/F");
			mc_values->Branch("phi", &phi, "phi/F");
			mc_values->Branch("ring", &ring, "ring/I");
			mc_values->Branch("element", &element, "element/I");
			mc_values->Branch("pm", &pm, "pm/I");
			mc_values->Branch("T", &T, "T/F");
			mc_values->Branch("tpmup", &tdctime[0], "tpmup/F");
			mc_values->Branch("tpmdo", &tdctime[1], "tpmdo/F");
			mc_values->Branch("z0", &z0, "z0/F");
			mc_values->Branch("beta", &beta);
			//mc_values->Branch("lzVecProton", &best_mc_proton);

			if (verbose)
			// test the GetChannel output
			for (int iring = 0; iring < 2; iring++)
				for (int ielement = 0; ielement < 24; ielement++){
					if ((iring == 0) && (ielement >= 12)) continue;
					int channel_pmup, channel_pmdo, ipm;
					std::string rpdname = GetChannel(iring, ielement, true, channel_pmup, channel_pmdo);
					std::cout << " ring " << iring << " element " << ielement << " ";
					std::cout << rpdname;
					std::cout << " : tdc ch PM up " << channel_pmup << " tdc ch PM down " << channel_pmdo << std::endl;

					std::cout << " check: " << ChannelToPMT(rpdname.c_str(), channel_pmup, ipm);
					std::cout << " ipm up " << ipm << " ";
					std::cout << ChannelToPMT(rpdname.c_str(), channel_pmdo, ipm);
					std::cout << " ipm down " << ipm << std::endl;

					rpdname = GetChannel(iring, ielement, false, channel_pmup, channel_pmdo);
					std::cout << rpdname;
					std::cout << " : sadc ch PM up " << channel_pmup << " sadc ch PM down " << channel_pmdo << std::endl;
					std::cout << " check: " << ChannelToPMT(rpdname.c_str(), channel_pmup, ipm);
					std::cout << " ipm up " << ipm << " ";
					std::cout << ChannelToPMT(rpdname.c_str(), channel_pmdo, ipm);
					std::cout << " ipm down " << ipm << std::endl;
				}
		}
		//std::cout << " size of MChits " << e.MChits().size() << std::endl;
	}
	// store the most likely properties (hits next to a given proton track
	float best_dEA(-1);
	float best_dEB(-1);
	float best_phi_diff_A(-1);
	float best_phi_diff_B(-1);
	float best_TA(-1);
	float best_TB(-1);
	TVector3 best_hit_A(0.,0.,0.);
	TVector3 best_hit_B(0.,0.,0.);
	// find the RPD hits
	for (int i = 0; i < (int) e.MChits().size(); i++){
		if	(strncmp(e.MChits()[i].DetRef().Name().c_str(), "RP", 2 ) == 0 ){
			X = e.MChits()[i].X();
			Y = e.MChits()[i].Y();
			Z = e.MChits()[i].Z();
			dE= e.MChits()[i].ELos();
			T = e.MChits()[i].Time();
			phi= atan2(Y, X);
			ring = -1; //ring 0 == Ring A and ring 1 == Ring B
			// actual description
			if (strncmp(e.MChits()[i].DetRef().Name().c_str(), "RP01R1__", 8 ) == 0 ){
				ring = 0;
			}
			if (strncmp(e.MChits()[i].DetRef().Name().c_str(), "RP01R2__", 8 ) == 0 ){
				ring = 1;
			}
			// old description in ComGeant 7.04
			if (strncmp(e.MChits()[i].DetRef().Name().c_str(), "RP01RA__", 8 ) == 0 ){
				ring = 0;
			}
			if (strncmp(e.MChits()[i].DetRef().Name().c_str(), "RP01RB__", 8 ) == 0 ){
				ring = 1;
			}
			if (ring < 0){
				std::cout << "\n ******* Error in RPD::Construct_MC_Digits(): ";
				std::cout << "unknown detector name containing RP ";
				std::cout << e.MChits()[i].DetRef().Name() << "***********" << std::endl;
			}
			// apply some spatial cuts since ComGeant is not fully correct yet
			/*
			double R = sqrt(X*X+Y*Y);
			if (	(ring == 0 && !( 11 < R && R < 13 && -72 < Z && Z < -21 )) ||
					(ring == 1 && !( 72 < R && R < 78 ))
				){
				if (verbose) std::cout << " skipping due to spatial cut " << std::endl;
				continue;
			}*/ // obsolete after having fixed this issue by Pawel Schnajder
			// store some MC values
			float phi_diff = fabs(best_mc_proton.Phi()-phi);
			if (ring == 0 && (best_dEA < 0 || best_phi_diff_A < 0 || best_phi_diff_A > phi_diff)){
				best_phi_diff_A = phi_diff;
				best_dEA		= dE;
				best_TA			= T;
				best_hit_A.SetXYZ(X,Y,Z);
			}
			if (ring == 1 && (best_dEB < 0 || best_phi_diff_B < 0 || best_phi_diff_B > phi_diff)){
				best_phi_diff_B = phi_diff;
				best_dEB		= dE;
				best_TB			= T;
				best_hit_B.SetXYZ(X,Y,Z);
			}

			element = Phi_to_element(phi, ring);

			if (verbose){
				std::cout << e.MChits()[i].DetRef().Name() << std::endl;
				std::cout << " x " << X << " y " << Y << " z " << Z;
				std::cout << " dE " << dE*1e3 << " t " << T << std::endl << std::endl;
				double radius = (double) TMath::Sqrt(X*X + Y*Y);
				std::cout << e.MChits()[i].DetRef().Name() << " " << radius << " in Ring " << ring << std::endl;
			}
			// first determine the position in the element in relation to the element center
			z0 = Z - zTargetCenter - ELEMENTOFFSET[ring]; // see determination of zA in Search()
			double pmdist[2] = {0.,0.};
			for (int ipm = 0; ipm < 2; ipm++){
				// half of the element length is the minimum way of light to go to the PM
				// The upstream PM (ipm == 0) has to go more by z0 the downsteam PM less by z0
				pmdist[ipm] = (ELEMENTLENGTH[ring]/2.+pow(-1.,ipm)*z0);
			}
			// check whether a MC hit was already generated in this element
			// but it is not treated yet
			static bool doublehits_warning(false);
			if (!doublehits_warning) {
				if (	(last_adc_entry[ring][element][0]!=-1)||(last_adc_entry[ring][element][1]!=-1)||
					(last_tdc_entry[ring][element][0]!=-1)||(last_tdc_entry[ring][element][1]!=-1)){
					doublehits_warning = true;
					std::cout << "\n ***** Warning in RPD_Helper::Construct_MC_Digits(): Double hits are not treated correctly yet! ***** " << std::endl;
				}
			}

			// ************** ADC simulation **************
			int nsamples = 32;
			int sadcch[2] = {0,0}; // PM upstream and downstream
			// get the correct channel
			std::string sadcname = GetChannel(ring, element, false, sadcch[0], sadcch[1]);
			// create the 2 PM digits containing the simulated information
			for (int ipm = 0; ipm < 2; ipm++){
				// creating a vector to store the simulated RawData
				std::vector<float> sadcinfo;
				sadcinfo.resize(nsamples+2);
				// the first entry is the channel number
				sadcinfo[0] = sadcch[ipm];
				sadcinfo[1] = -1;
				// simulate an ADC signal in terms of channels, ignoring all losses
				double fdE(0);
				// take the attenuation length into account
				fdE = dE*TMath::Exp(-pmdist[ipm]/ATTENUATIONLENGTH[ring])*CORRECTION_FACTOR[ring];
				if (ring == 0)
					fdE *= ((double)dEAcusp[element])/dEAcuspMeV*1e3; // Energyloss in terms of samples in MeV
				else
					fdE *= ((double)dEBcusp[element])/dEBcuspMeV*1e3;
				// expand the sum of Charge over SADC's samples
				// first 4 entries after the info are used for the pedestal
				// calculation and therefore kept here 0
				for (int isam = 2; isam < nsamples +2; isam++){
					// note the integral is only 0.95 in this range
					sadcinfo[isam]=	TMath::Landau(isam, 11, 1., true)*1./0.95*fdE;
				}
				PaDigit sadc_digit(sadcch[ipm], sadcinfo);
				sadc_digit.StoreName(sadcname);
				if (verbose)
					std::cout << " detnam MC " << sadcname << " ch " << sadc_digit.DigInfo(0) << " time " << sadc_digit.DigInfo(4) << " pm number " << sadcch[ipm] << " Energy deposit " << fdE << std::endl;
				result.push_back(sadc_digit);
				// store the last entry
				last_adc_entry[ring][element][ipm] = result.size()-1;
			}
			// ************** TDC simulation ***************

			int tdcch[2] = {0,0}; // PM upstream and downstream
			std::string tdcname = GetChannel(ring, element, true, tdcch[0], tdcch[1]);
			tdctime[0] = 0.; tdctime[1] = 0.;
			// create a digit containing the simulated information
			// compensate the time unit in the mapping file and the time cut in Search()
			double tdctime0 = T - T0_ALL; // T/tdcresolution - 970.;
			// propagate the time of impact down to the PM's by knowing
			// the speed of light in the element and the element dimensions/position
			for (int ipm = 0; ipm < 2; ipm++){
				// the time depends on the distance of the impact to the pm
				double _tdctime = pmdist[ipm];
				if (ring == 0){
					// divided by the light speed in the element...
					_tdctime /= lightSpeedInA[element];
				}
				else {
					_tdctime /= lightSpeedInB[element];
				}
				// in terms of tdc channels
				//_tdctime /= tdcresolution; no, digits are already decoded
				// and added by the impact time in terms of channels is the pm time
				tdctime[ipm] = (_tdctime+tdctime0);
				// add the smearing by knowing the measured timeresolution
				tdctime[ipm] = gRandom->Gaus(tdctime[ipm], TIME_SMEARING[ring]);

				if (verbose){
					std::cout << " tdc time pm " << ipm << " is " << tdctime[ipm] << std::endl;
				}
				std::vector<float> tdcinfo;
				tdcinfo.resize(5);
				tdcinfo[0] = (float) tdcch[ipm];
				tdcinfo[1] = 0.;
				tdcinfo[2] = 0.;
				tdcinfo[3] = 0.;
				tdcinfo[4] = tdctime[ipm];
				PaDigit tdc_digit(tdcch[ipm], tdcinfo);
				tdc_digit.StoreName(tdcname);
				result.push_back(tdc_digit);
				// store the last entry
				last_tdc_entry[ring][element][ipm] = result.size()-1;
				if (verbose)
					std::cout << " detnam MC " << tdcname << " ch " << tdc_digit.DigInfo(0) << " time " << tdc_digit.DigInfo(4) << " pm number " << tdcch[ipm] << std::endl;
			}
			if (mc_values){
				mc_values->Fill();
			}
		}
	}
	// must be less then slab resolution pointing away from the MC correct solution
	// only the best solution is stored with hits in ring A and B
	dEA_MC.clear();
	dEB_MC.clear();
	if (best_dEA > 0 && best_dEB > 0 && best_phi_diff_A < 0.53 && best_phi_diff_B < 0.27){
		float hit_distance = (best_hit_B - best_hit_A).Mag();
		float beta = hit_distance / (best_TB-best_TA) / lightSpeed;
		if (beta > 0 && beta < 1.){
			dEA_MC.push_back(best_dEA);
			dEB_MC.push_back(best_dEB);
		}
	}

	static int counter(0);
	//std::cout << counter << std::endl;
	counter++;
	if (mc_values && draw && ((counter%1000)==0)){
		static TCanvas* debug_window(NULL);
		if (!debug_window){
			debug_window = new TCanvas("debug_window", "RPD helper MC ouput", 600, 600);
			debug_window->Divide(2,2);
		}
		debug_window->cd(1);
		gPad->Clear();
		mc_values->Draw("Y:X");
		//mc_values->Draw("z0", "ring == 0");
		gPad->Update();
		debug_window->cd(2);
		gPad->Clear();
		//mc_values->Draw("phi:Z");
		//mc_values->Draw("z0", "ring == 1");
		gPad->Update();
		debug_window->cd(3);
		gPad->Clear();
		mc_values->Draw("tpmdo:tpmup", "ring == 0");
		gPad->Update();
		debug_window->cd(4);
		gPad->Clear();
		mc_values->Draw("tpmdo:tpmup", "ring == 1");
		gPad->Update();
	}

	return result;
}


void
RPD::CheckMCresolutions() const
{
	if (!debug) return;
	// find the best matching MC recoil track
	// corresponding to the best reconstructed track
	const Phast& phast = Phast::Ref();
	// only MC data for comparison
	if (!phast.event.IsMC()) return;
	if (fBestMCProtonTrack < 0) return;
	if (fBestProtonTrack   < 0) return;
	if (dEA_MC.size() < 1) return;
	if (dEB_MC.size() < 1) return;
	TLorentzVector proton_mc = fRPDMCTrack[fBestMCProtonTrack];
	TLorentzVector proton    = fRPDTrack  [fBestProtonTrack];
	//const std::vector<PaDigit>& digits = Get_Digits(phast.event);
	int ibestprimvtx = phast.event.iBestPrimaryVertex();
	if (ibestprimvtx < 0) return;
	const PaVertex& vertex = phast.event.vVertex(ibestprimvtx);
	// assuming to have only one primary MC vertex
	int iMCprimvtx = -1;
	for (int i = 0; i < phast.event.NMCvertex(); i++){
		if (phast.event.vMCvertex(i).IsPrimary()){
			iMCprimvtx = i;
			continue;
		}
	}
	if (iMCprimvtx < 0) return;
	// take only the first MC vertex
	const PaMCvertex& mcvertex = phast.event.vMCvertex(iMCprimvtx);
	static float vtx_z_RPD;
	static float vtx_z_MC;
	static float vtx_z_Spec;
	// beta_mc was calculated from 2 best MC hits in Construct_MC_digits
	static float beta_MC;
	static float beta_RPD;
	static float dE_A;
	static float dE_B;
	static float dE_A_MC;
	static float dE_B_MC;
	static TTree* tree_rpd_resolution(NULL);
	if (!tree_rpd_resolution){
		tree_rpd_resolution = new TTree("tree_rpd_resolution", "tree_rpd_resolution");
		tree_rpd_resolution->Branch("vtx_z_RPD",&vtx_z_RPD);
		tree_rpd_resolution->Branch("vtx_z_MC",&vtx_z_MC);
		tree_rpd_resolution->Branch("vtx_z_Spec",&vtx_z_Spec);
		tree_rpd_resolution->Branch("beta_MC",&beta_MC);
		tree_rpd_resolution->Branch("beta_RPD",&beta_RPD);
		tree_rpd_resolution->Branch("dE_A",&dE_A);
		tree_rpd_resolution->Branch("dE_B",&dE_B);
		tree_rpd_resolution->Branch("dE_A_MC",&dE_A_MC);
		tree_rpd_resolution->Branch("dE_B_MC",&dE_B_MC);
	}
	vtx_z_RPD   = fzT[iBestProtonTrack()];
	vtx_z_MC    = mcvertex.Pos(2);
	vtx_z_Spec  = vertex.Z();
	beta_MC  = proton_mc.Beta();
	beta_RPD = proton.Beta();
	dE_A = dEA[iBestProtonTrack()];
	dE_B = dEB[iBestProtonTrack()];
	dE_A_MC = dEA_MC[0]; // only one hit per MC event is expected and written out
	dE_B_MC = dEB_MC[0];
	tree_rpd_resolution->Fill();
}


const std::vector<PaDigit>&
RPD::Get_Digits(const PaEvent& e)
{
	if (!e.IsMC()){
		isMC = false;
		return e.RawDigits();
	}
	else {
		isMC = true;
		return Construct_MC_Digits(e);
	}
}


int
RPD::Phi_to_element(const double phi, const int ring) const
{
	int result(-1);
	if ((ring < 0) || (ring > 1)){
		std::cout << " Error in RPD::Phi_to_element(): wrong input for RPD Ring: " << ring << std::endl;
		return result;
	}
	if ((phi < -TMath::Pi()) || (TMath::Pi() < phi)) {
		std::cout << " Error in RPD::Phi_to_element(): wrong input for phi angle: " << phi << std::endl;
		return result;
	}
	// switch to [0, 2pi[ range
	double angle = (phi < 0.) ? (TMath::TwoPi()+phi):phi;
	// changing the [0, 2pi[ range to [0,1[ range. In addition:
	// phi = 0 is the middle of a element. Therefore one has to
	// tilt the angle by 1/2 * 1/nelements
	float  element	 = (angle/TMath::TwoPi()+1./(2.*NELEMENTS[ring]))*NELEMENTS[ring];
	// quantize the element down to integer
	result = (int) element;
	// the last element in this calculation is the first one
	result = (result >= NELEMENTS[ring]) ? (result - NELEMENTS[ring]) : result;
	// test the output since I (P.Jasinski) don't trust my skills ;)
	if ((0 > result) || (result > (NELEMENTS[ring]-1))){
		std::cout << " Unexpected error in RPD:Phi_to_element(): computed element is wrong! Please report." << std::endl;
		return -1;
	}
	return result;
	/* if you would like to see this code working properly, voila:
	for (int i = 0; i < 360; i++){
		float angle = ((float)i)/360. * TMath::TwoPi();
		float x = cos(angle);
		float y = sin(angle);
		float computed_angle = atan2(y,x);
		computed_angle = (computed_angle < 0.) ? (TMath::TwoPi()+computed_angle):computed_angle;
		float element = (computed_angle/TMath::TwoPi()+1./(2.*nelements))*nelements;
		element = (element >= nelements) ? (element - nelements) : element;
		int qelement = (int) element;
		std::cout << " step \t angle \t x \t y \t atan2 \t element \t" << std::endl;
		std::cout << i << " " << angle << " " << x  << " " << y << " " << computed_angle  << " " << element << " = " << qelement << std::endl;
	}*/
}


void
RPD::DecodeChipF1Digit(const PaDigit& /*digit*/)
{
	;
}


void
RPD::DecodeChipSADCDigit(const PaDigit& /*digit*/)
{
	;
}


void
RPD::DecodeChipDigit(const PaDigit& /*digit*/)
{
	;
}


void
RPD::DecodeChipDigits(const std::vector<PaDigit>& digits)
{
	// Searches the data for RPD infos and calculates the tracks

	static TH1F* hist_tApmup(NULL);
	static TH1F* hist_tApmdo(NULL);
	static TH1F* hist_tBpmup(NULL);
	static TH1F* hist_tBpmdo(NULL);
	static TH1F* hist_qApmupl(NULL);
	static TH1F* hist_qApmdol(NULL);
	static TH1F* hist_qBpmupl(NULL);
	static TH1F* hist_qBpmdol(NULL);
	static TH1F* hist_qApmuph(NULL);
	static TH1F* hist_qApmdoh(NULL);
	static TH1F* hist_qBpmuph(NULL);
	static TH1F* hist_qBpmdoh(NULL);
	static TH1F* hist_t0_calib[2][24][2];
	if (debug && !hist_tApmup){
		std::cout << " booking histograms for debugging " << std::endl;
		hist_tApmup = new TH1F("hist_tApmup", "time for ring A pm upstream", 100, -20, 20);
		hist_tApmdo = new TH1F("hist_tApmdo", "time for ring A pm downstream", 100, -20, 20);
		hist_tBpmup = new TH1F("hist_tBpmup", "time for ring B pm upstream", 100, -20, 20);
		hist_tBpmdo = new TH1F("hist_tBpmdo", "time for ring B pm downstream", 100, -20, 20);
		hist_qApmupl = new TH1F("hist_dEApmupl", "energy deposit for ring A pm upstream", 1000, 0, 3000);
		hist_qApmdol = new TH1F("hist_dEApmdol", "energy deposit for ring A pm downstream", 1000, 0, 3000);
		hist_qBpmupl = new TH1F("hist_dEBpmupl", "energy deposit for ring B pm upstream", 1000, 0, 3000);
		hist_qBpmdol = new TH1F("hist_dEBpmdol", "energy deposit for ring B pm downstream", 1000, 0, 3000);
		hist_qApmuph = new TH1F("hist_dEApmuph", "energy deposit for ring A pm upstream large signals", 1000, 0, 3000);
		hist_qApmdoh = new TH1F("hist_dEApmdoh", "energy deposit for ring A pm downstream large signals", 1000, 0, 3000);
		hist_qBpmuph = new TH1F("hist_dEBpmuph", "energy deposit for ring B pm upstream large signals", 1000, 0, 3000);
		hist_qBpmdoh = new TH1F("hist_dEBpmdoh", "energy deposit for ring B pm downstream large signals", 1000, 0, 3000);
		for (int iring = 0; iring < 2; iring++){
			for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
				for (int ipm = 0; ipm < 2; ipm++){
					std::stringstream title;
					title << "hist_t0_calib_" << iring << "_" << ielement << "_" << ipm;
					hist_t0_calib[iring][ielement][ipm] = new TH1F(title.str().c_str(),title.str().c_str(), 200, -50-970, 50-970);
				}
			}
		}
	}

	const std::vector<PaDigit>& vdd = digits;
	for (int i = 0; i < int(vdd.size()); i++) { // loop over DAQ digits
		const PaDigit& d = vdd[i];

		std::string detnam = d.DecodeMapName();
		int ctr = 0;
		const std::string planeName = ChannelToPMT(detnam.data(), int(d.DigInfo(0)), ctr);

		if (verbose){
			std::cout << " detnam " << detnam.data() << " ch " << d.DigInfo(0) << " time " << d.DigInfo(4) << " pm number " << ctr << std::endl;
		}
		// decoding of SADC info

		// SADC for normal signals
		if (detnam.find("RP01Ql") == 0) {
			if (planeName.find("Au") != std::string::npos) {
				fDataSentAupql[ctr] = true;
				fAdcAupl[ctr] = treatADC(d);
				if (hist_qApmupl) hist_qApmupl->Fill(fAdcAupl[ctr]);
			}
			if (planeName.find("Ad") != std::string::npos) {
				fDataSentAdoql[ctr] = true;
				fAdcAdol[ctr] = treatADC(d);
				if (hist_qApmdol) hist_qApmdol->Fill(fAdcAdol[ctr]);
			}
			if (planeName.find("Bu") != std::string::npos) {
				fDataSentBupql[ctr] = true;
				fAdcBupl[ctr] = treatADC(d);
				if (hist_qBpmupl) hist_qBpmupl->Fill(fAdcBupl[ctr]);
			}
			if (planeName.find("Bd") != std::string::npos) {
				fDataSentBdoql[ctr] = true;
				fAdcBdol[ctr] = treatADC(d);
				if (hist_qBpmdol) hist_qBpmdol->Fill(fAdcBdol[ctr]);
			}
		}
		
		// SADC for large signals
		if (detnam.find("RP01Qh") == 0) {
			if (planeName.find("Au") != std::string::npos) {
				fDataSentAupqh[ctr] = true;
				fAdcAuph[ctr] = treatADC(d);
				if (hist_qApmuph) hist_qApmuph->Fill(fAdcAuph[ctr]);
			}
			if (planeName.find("Ad") != std::string::npos) {
				fDataSentAdoqh[ctr] = true;
				fAdcAdoh[ctr] = treatADC(d);
				if (hist_qApmdoh) hist_qApmdoh->Fill(fAdcAdoh[ctr]);
			}
			if (planeName.find("Bu") != std::string::npos) {
				fDataSentBupqh[ctr] = true;
				fAdcBuph[ctr] = treatADC(d);
				if (hist_qBpmuph) hist_qBpmuph->Fill(fAdcBuph[ctr]);
			}
			if (planeName.find("Bd") != std::string::npos) {
				fDataSentBdoqh[ctr] = true;
				fAdcBdoh[ctr] = treatADC(d);
				if (hist_qBpmdoh) hist_qBpmdoh->Fill(fAdcBdoh[ctr]);
			}
		}

		// decoding of F1 info

		if (detnam.find("RP01TA") == 0) {
			double f1time = d.DigInfo(4) + rpd_constants.t0[0][ctr][0];
			if (planeName.find("Au") != std::string::npos && fDataSentAupt[ctr] < 4 && TMath::Abs(f1time) < 40/*15*/) { //
				fTdcAup[ctr][fDataSentAupt[ctr]] = f1time;
				if (hist_tApmup){
					hist_tApmup->Fill(fTdcAup[ctr][fDataSentAupt[ctr]]);
					hist_t0_calib[0][ctr][0]->Fill(d.DigInfo(4));
				}
				fDataSentAupt[ctr]++;
			}
			f1time = d.DigInfo(4) + rpd_constants.t0[0][ctr][1];
			if (planeName.find("Ad") != std::string::npos && fDataSentAdot[ctr] < 4 && TMath::Abs(f1time) < 40/*15*/) { //
				fTdcAdo[ctr][fDataSentAdot[ctr]] = f1time;
				if (hist_tApmdo){
					hist_tApmdo->Fill(fTdcAdo[ctr][fDataSentAdot[ctr]]);
					hist_t0_calib[0][ctr][1]->Fill(d.DigInfo(4));
				}
				fDataSentAdot[ctr]++;
			}
		}

		if (detnam.find("RP01TBl") == 0) {
			double f1time = d.DigInfo(4) + rpd_constants.t0[1][ctr][0];
			if (planeName.find("Bu") != std::string::npos && fDataSentBupt[ctr] < 4 && TMath::Abs(f1time) < 40/*15*/) { //
				fTdcBup[ctr][fDataSentBupt[ctr]] = f1time;
				if (hist_tBpmup){
					hist_tBpmup->Fill(fTdcBup[ctr][fDataSentBupt[ctr]]);
					hist_t0_calib[1][ctr][0]->Fill(d.DigInfo(4));
				}
				fDataSentBupt[ctr]++;
			}
			f1time = d.DigInfo(4) + rpd_constants.t0[1][ctr][1];
			if (planeName.find("Bd") != std::string::npos && fDataSentBdot[ctr] < 4 && TMath::Abs(f1time) < 40/*15*/) { //
				fTdcBdo[ctr][fDataSentBdot[ctr]] = f1time;
				if (hist_tBpmdo){
					hist_tBpmdo->Fill(fTdcBdo[ctr][fDataSentBdot[ctr]]);
					hist_t0_calib[1][ctr][1]->Fill(d.DigInfo(4));
				}
				fDataSentBdot[ctr]++;
			}
		}

	}// end of loop over DAQ digits

	static int counter(0);
	counter++;
	static TCanvas* debug_window(NULL);
	//static TCanvas* t0_window(NULL);
	if ((debug) && draw && ((counter % 1000) == 0)){
		if (!debug_window){
			debug_window = new TCanvas("debug_window_search", "RPD helper MC ouput", 1000, 600);
			debug_window->Divide(4,2);
		}
		debug_window->cd(1);
		gPad->Clear();
		hist_tApmup->Draw();
		gPad->Update();
		debug_window->cd(2);
		gPad->Clear();
		hist_tApmdo->Draw();
		gPad->Update();
		debug_window->cd(3);
		gPad->Clear();
		hist_tBpmup->Draw();
		gPad->Update();
		debug_window->cd(4);
		gPad->Clear();
		hist_tBpmdo->Draw();
		gPad->Update();
		debug_window->cd(5);
		gPad->Clear();
		hist_qApmupl->Draw();
		gPad->Update();
		debug_window->cd(6);
		gPad->Clear();
		hist_qApmdol->Draw();
		gPad->Update();
		debug_window->cd(7);
		gPad->Clear();
		hist_qBpmupl->Draw();
		gPad->Update();
		debug_window->cd(8);
		gPad->Clear();
		hist_qBpmdol->Draw();
		gPad->Update();
		/*
		if (!t0_window){
			t0_window = new TCanvas("t0_window_search", "RPD helper t0 calibration", 1200, 800);
			t0_window->Divide(12,6);
		}
		int ipad(0);
		for (int iring = 0; iring < 2; iring++){
			for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
				for (int ipm = 0; ipm < 2; ipm++){
					ipad++;
					t0_window->cd(ipad);
					gPad->Clear();
					hist_t0_calib[iring][ielement][ipm]->Draw();
					gPad->Update();
				}
			}
		}*/
	}
}


void
RPD::Clusterize(const PaVertex& vertex, const bool correct_momentum, const bool correct_angles)
{
	int iB = 0;
	for (int iA = 0; iA < 12; iA++) {
		for (int ilB = -1; ilB < 2; ilB++) {

			if (iA == 0 && ilB == -1) {
				iB = 23;
			} else {
				iB = 2 * iA + ilB;
			}

			for (int iHitAup = 0; iHitAup < fDataSentAupt[iA]; iHitAup++) {
				for (int iHitAdo = 0; iHitAdo < fDataSentAdot[iA]; iHitAdo++) {
					for (int iHitBup = 0; iHitBup < fDataSentBupt[iB]; iHitBup++) {
						for (int iHitBdo = 0; iHitBdo < fDataSentBdot[iB]; iHitBdo++) {

							if (fAdcAdol[iA] < 10 || fAdcAupl[iA] < 10) {
								if (fAdcAdoh[iA] < 10 || fAdcAuph[iA] < 10 ) continue;
								}
							if (fAdcBdol[iB] < 10 || fAdcBupl[iB] < 10) {
								if (fAdcBdoh[iB] < 10 || fAdcBuph[iB] < 10) continue;
								}

							const float dtA = fTdcAup[iA][iHitAup] - fTdcAdo[iA][iHitAdo];
							float zA = dtA * lightSpeedInA[iA] / 2.;
							if (!isMC) zA += positionOffsetForA[iA];
							const float tA = (fTdcAup[iA][iHitAup] + fTdcAdo[iA][iHitAdo]) / 2.;

							const float dtB = fTdcBup[iB][iHitBup] - fTdcBdo[iB][iHitBdo];
							float zB = dtB * lightSpeedInB[iB] / 2.;
							if (!isMC) zB += positionOffsetForB[iB];
							const float tB = (fTdcBup[iB][iHitBup] + fTdcBdo[iB][iHitBdo]) / 2.;

							zA = zTargetCenter + zA - ELEMENTOFFSET_UP[0];
							zB = zTargetCenter + zB - ELEMENTOFFSET_UP[1];

							const float dof = TMath::Sqrt((zA - zB) * (zA - zB) + deltaR * deltaR);
							float beta;
							if (!isMC)
								beta = ((tB - tA + tofOffset[iA][iB] + tofGlobalOffset) > deltaTmin) ? dof / (tB - tA + tofOffset[iA][iB] + tofGlobalOffset) / lightSpeed : 1.;
							else
								beta = ((tB - tA) > deltaTmin) ? dof / (tB - tA) / lightSpeed : 1.;
							beta = (beta > 0.99995) ? 0 : beta; // introduced cutoff for not physical momenta
							beta = (beta < 0.) ? 0. : beta;
							float fdEA = TMath::Sqrt(fAdcAdol[iA] * fAdcAupl[iA]);
							float fdEB = TMath::Sqrt(fAdcBdol[iB] * fAdcBupl[iB]);
							fdEA = fdEA * dEAcuspMeV / dEAcusp[iA];
							fdEB = fdEB * dEBcuspMeV / dEBcusp[iB];


							if (verbose){
								std::cout << " measured pos z ring A " << zA << " ring B " << zB << std::endl;
								std::cout << " measured edep ring A " << fdEA << " ring B " << fdEB << std::endl;
							}

							something = true;

							if (!isMC){
								fzT.push_back(zA - rA * (zB - zA) / deltaR + zVertexOffset[iA][iB]);
								ftT.push_back((beta > 0) ? tA - dof * rA / rB / (beta * lightSpeed) - calTA[iA] : -777);
							}
							else {
								fzT.push_back(zA - rA * (zB - zA) / deltaR);
								ftT.push_back((beta > 0) ? tA - dof * rA / rB / (beta * lightSpeed): -777);
							}

							const double fCosTheta = (zB - zA) / dof;
							const double fSinTheta = deltaR / dof;
							double fTheta = TMath::ACos((zB - zA) / dof);
							double fPhi = 2. * TMath::Pi() * 15. / 360. * (iB + 1e-3);
							if (fPhi > TMath::Pi()) fPhi -= 2. * TMath::Pi();
							if (ilB) fPhi -= 2. * TMath::Pi() * 3.75 / 360. * ilB;
							if (correct_angles) {
								// Correct angles for vertex position
								fPhi = CorrectPhi(fPhi, vertex);
								//correct for beam divergence (later RPD alignment also possible to include via origin )
								CorrectAngles(vertex, fTheta, fPhi);
							}
							double fMomentum = beta * mP / TMath::Sqrt(1. - beta * beta);
							if (correct_momentum && fMomentum > 1e-6)// && fMomentum < 0.596)
								fMomentum = correctEnergyLoss(fMomentum, fSinTheta, fCosTheta, fPhi, vertex.X(), vertex.Y(), vertex.Z());// correct for energy loss in target material
							const double fMom_x = (fMomentum * fSinTheta * TMath::Cos(fPhi) > 1e-6 || fMomentum * fSinTheta * TMath::Cos(fPhi) < -1e-6) ? fMomentum * fSinTheta * TMath::Cos(fPhi) : 0;
							const double fMom_y = (fMomentum * fSinTheta * TMath::Sin(fPhi) > 1e-6 || fMomentum * fSinTheta * TMath::Sin(fPhi) < -1e-6) ? fMomentum * fSinTheta * TMath::Sin(fPhi) : 0;
							const double fMom_z = (fMomentum * fCosTheta > 1e-6 || fMomentum * fCosTheta < -1e-6) ? fMomentum * fCosTheta : 0;
							hits.push_back(std::make_pair(iA, iB));
							nTracks++;
							TLorentzVector fRPD4Vector; // assigning proton mass to the 4-vector
							fRPD4Vector.SetXYZM(fMom_x, fMom_y, fMom_z, mP); // assigning proton mass to the 4-vector
							fRPDTrack.push_back(fRPD4Vector);
							dEA.push_back(fdEA);
							dEB.push_back(fdEB);
							zA_vec.push_back(zA);
							zB_vec.push_back(zB);
							fCalibratedEvent.push_back(true);
						}
					}
				}
			}
		}
	}

	if (!something) {
		fzT.push_back(-777);
		ftT.push_back(-777);
		hits.push_back(std::make_pair(-1, -1));
		fCalibratedEvent.push_back(false);
		fRPDTrack.push_back(TLorentzVector(0., 0., 0., 0.));
		dEA.push_back(0);
		dEB.push_back(0);
		zA_vec.push_back(-777);
		zB_vec.push_back(-777);
	}

	// guess best proton track, still to be modified! only 2008/2009 hydrogen data!
	fBestProtonTrack = -1;
	for (unsigned int i = 0; i < fRPDTrack.size(); i++)
		if (	(fRPDTrack[i].Vect().Mag() > 0.) &&
				(fabs(fzT[i] + 48.) < 40.) && // rough cut around the target
				fCalibratedEvent[i])
			fBestProtonTrack = i;
	CheckMCresolutions();
}


Calib_const
RPD::Load_DB_entry(const float db_entries[]) const
{
	Calib_const result;

	int pos(0);
	result.runnumb = (int) db_entries[pos];
	pos++;

	for (int iring = 0; iring < 2; iring++){
		for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
			for (int ipm = 0; ipm < 2; ipm++){
				result.t0[iring][ielement][ipm] = db_entries[pos];
				pos++;
			}
		}
	}

	for (int iring = 0; iring < 2; iring++){
		for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
				result.positionOffset[iring][ielement] = db_entries[pos];
				pos++;
		}
	}

	for (int ielementA = 0; ielementA < NELEMENTS[0]; ielementA++){
		for (int ielementB = 0; ielementB < NELEMENTS[1]; ielementB++){
			result.zVertexOffset[ielementA][ielementB] = db_entries[pos];
			pos++;
		}
	}

	for (int iring = 0; iring < 2; iring++){
		for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
			result.dEcusp[iring][ielement] = db_entries[pos];
			pos++;
		}
	}

	for (int ielement = 0; ielement < NELEMENTS[0]; ielement++){
		result.calTA[ielement] = db_entries[pos];
		pos++;
	}

	result.tofGlobalOffset = db_entries[pos];
	pos++;

	return result;
}


bool
RPD::Load_DB()
{
	bool result(false);

	std::cout << "\nLoading RPD database";
	calib_constants_db.clear();

	/********************** old part ***********************/

	int nentries = sizeof(const_t0)/sizeof(float);
	// we expect 1 runnumber, ring A and B
	// with 12 and 24 elements times
	// 2 pms(pm up and pm down)
	int ndata = 1+NELEMENTS[0]*2+NELEMENTS[1]*2;
	// has the db the correct size?
	assert(nentries%ndata == 0);
	// number of runs that are stored
	int nruns = (int) (nentries/ndata);
	// fill the mapped db
	int pos_const_t0(-1);
	for (int ientry = 0; ientry < nruns; ientry++){
		pos_const_t0++;
		// the first entry is the run number
		int runnumb =  (int) const_t0[pos_const_t0];
		// looks strange since the runnumber is also stored inside
		// the struct as well as it is the key for the map used here
		calib_constants_db[runnumb].runnumb = runnumb;
		// Initialize the t0 entries if needed
		memset(&calib_constants_db[runnumb].t0, 0, (ndata-1)*sizeof(float));
		// the rest of the entries are the pm t0's
		for (int iring = 0; iring < 2; iring++){
		std::cout << "." << std::flush;
			for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
				for (int ipm = 0; ipm < 2; ipm++){
					//if ((iring == 0) && (ielement >= 12)) continue;
					pos_const_t0++;
					calib_constants_db[runnumb].t0[iring][ielement][ipm] = const_t0[pos_const_t0];
					// test the filling procedure
					if (verbose){
						std::cout << " run number " << calib_constants_db[runnumb].runnumb;
						std::cout << " ring " << iring << " element " << ielement;
						std::cout << " pm " << ipm;
						std::cout << " t0 const " << calib_constants_db[runnumb].t0[iring][ielement][ipm];
						std::cout << std::endl;
					}
				}
			}
		}
	}

	std::cout << " done.\n" << std::endl;

	//Save_DB(); // only needed once if the structure of the DB was changed

	result = true;

	last_DB_runnumb = -1;

	/********************** end of old part *************************/

	/********************** new part *******************************/
	/*
	int nentries = sizeof(const_db)/sizeof(float);
	std::cout << nentries << " / " << NENTRIES_CALIB_DB << std::endl;
	// has the db the correct size?
	assert(nentries%NENTRIES_CALIB_DB == 0);
	// number of runs that are stored
	int nruns = (int) (nentries/NENTRIES_CALIB_DB);
	std::cout << nruns << std::endl;
	// fill the mapped db
	for (int ientry = 0; ientry < nruns; ientry++){
		// give the reference to the n'th entry to load
		Calib_const calib_const = Load_DB_entry(&const_db[ientry*NENTRIES_CALIB_DB]);
		// save the constants to the DB map
		calib_constants_db[calib_const.runnumb]=calib_const;
	}

	std::cout << "done " << std::endl;
	result = true;

	last_DB_runnumb = -1;

	//Save_DB();
	*/
	/********************************* end of new part ************/

/*
	std::map<int, Calib_const>::const_iterator im;

		std::cout << "testing db first time" << std::endl;
		if (Calib_constants_db.begin()==Calib_constants_db.end())
			std::cout << " hm. something wrong! " << std::endl;
		for (im = Calib_constants_db.begin(); im != Calib_constants_db.end(); im++){
			std::cout << (*im).first << " has " << (*im).second.runnumb << std::endl;
		}
*/
	return result;
}


void
RPD::FindInDB(const int runnumb)
{
/*
	std::map<int, Calib_const>::const_iterator im;

		std::cout << "testing db second time" << std::endl;
		if (Calib_constants_db.begin()==Calib_constants_db.end())
			std::cout << " hm. something wrong! " << std::endl;
		for (im = Calib_constants_db.begin(); im != Calib_constants_db.end(); im++){
			std::cout << (*im).first << " has " << (*im).second.runnumb << std::endl;
		}
		std::cout << " size is " << Calib_constants_db.size() << std::endl;
*/
	// check for the last checked run number to
	// enhance the performance a bit
	if (last_DB_runnumb == runnumb){
		return;
	}
	last_DB_runnumb = runnumb;
	// in case of MC initialize to t0 calibration stored as
	// run number 1
	if (isMC){
		rpd_constants = calib_constants_db[1];
		//int ndata = NELEMENTS[0]*2+NELEMENTS[1]*2;
		if (rpd_constants.runnumb != 1)
			std::cout << " Error in RPD::FindInDB(): loading of MC t0 calibration constants failed! " << std::endl;
		/*
		for (int iring = 0; iring < 2; iring++){
			for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
				for (int ipm = 0; ipm < 2; ipm++){
					rpd_t0_constants.t0[iring][ielement][ipm] = T0_ALL;
				}
			}
		}*/
		//memset(&rpd_t0_constants.t0, T0_ALL, ndata*sizeof(float)); wrong!
		if (verbose) std::cout << " New runnumber " << runnumb << ". MC t0 calibration values set " << std::endl;
		return;
	}
	// in all other cases search in DB
	std::map<int, Calib_const>::const_iterator im;

	// find the element grater than the run number
	{
		im = calib_constants_db.upper_bound(runnumb);
		if(im != calib_constants_db.end()) { // found
			// point to the element before
			if (im != calib_constants_db.begin()) {
				im--;
				//out = (*im).second;
				rpd_constants = (*im).second;
			}
		}
		else{ // the last element is taken
			rpd_constants = (*im).second;
		}
	}
	if (verbose){
		std::cout << " New runnumber " << runnumb;
		std::cout << ": using t0 calibration constants from run ";
		std::cout << rpd_constants.runnumb << "." << std::endl;
	}
}


bool
RPD::Save_DB(const std::string& filename) const
{
	bool result(false);

	std::map<int, Calib_const>::const_iterator db_it;
	for (db_it=calib_constants_db.begin() ; db_it!= calib_constants_db.end(); db_it++ ){
		Write_constants(db_it->second, filename);
	}
	return result;
}


bool
RPD::Write_constants(const Calib_const& calib_const, const std::string& filename) const
{
	bool result(false);

	// if filename not given, create one
	std::stringstream _filename;
	if (filename == ""){
		_filename << "calib_const_run_" << calib_const.runnumb << ".txt";
	} else {
		_filename << filename;
	}
	std::cout << " writing calibration entry for run " << calib_const.runnumb << std::endl;
	std::cout << " in " << _filename.str() << std::endl;

	ofstream calib_file;
	calib_file.open(_filename.str().c_str(), ios::app);

	calib_file << std::endl;
	calib_file << " /***********************************/" << std::endl;
	calib_file << calib_const.runnumb << ", " << std::endl;
	calib_file << " // t0 " << std::endl;
	for (int iring = 0; iring < 2; iring++){
		for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
			for (int ipm = 0; ipm < 2; ipm++){
				calib_file << calib_const.t0[iring][ielement][ipm] << ", ";
			}
		}
		calib_file << std::endl;
	}
	calib_file << std::endl;

	calib_file << " // positionOffset" << std::endl;
	for (int iring = 0; iring < 2; iring++){
		for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
				calib_file << calib_const.positionOffset[iring][ielement] << ", ";
		}
		calib_file << std::endl;
	}
	calib_file << std::endl;

	calib_file << " // zVertexOffset" << std::endl;
	for (int ielementA = 0; ielementA < NELEMENTS[0]; ielementA++){
		for (int ielementB = 0; ielementB < NELEMENTS[1]; ielementB++){
				calib_file << calib_const.zVertexOffset[ielementA][ielementB] << ", ";
		}
		calib_file << std::endl;
	}

	calib_file << " // dEcusp" << std::endl;
	for (int iring = 0; iring < 2; iring++){
		for (int ielement = 0; ielement < NELEMENTS[iring]; ielement++){
				calib_file << calib_const.dEcusp[iring][ielement] << ", ";
		}
		calib_file << std::endl;
	}
	calib_file << std::endl;

	calib_file << " // calTA" << std::endl;
	for (int ielement = 0; ielement < NELEMENTS[0]; ielement++){
		calib_file << calib_const.calTA[ielement] << ", ";
	}
	calib_file << std::endl;

	calib_file << " // tofGlobalOffset" << std::endl;
	calib_file << calib_const.tofGlobalOffset << ", ";
	calib_file << std::endl;

	calib_file.close();

	result = true;
	return result;
}


bool
RPD::Set_DB_entry(const Calib_const& calib_const)
{
	bool result(false);
	calib_constants_db[calib_const.runnumb] = calib_const;
	return result;
}


int
RPD::MapPhiElement(const double phi, const std::string& ring) const
{
	const size_t foundA = ring.find("A");
	const size_t foundB = ring.find("B");
	if (foundA != std::string::npos && foundB == std::string::npos) return Phi_to_element(phi, 0);
	if (foundA == std::string::npos && foundB != std::string::npos) return Phi_to_element(phi, 1);
	return -1;
}


void
RPD::CorrectAngles(const PaVertex& vertex, double& fTheta, double& fPhi) const
{
	// rotation of reference system in direction of the beam
	const Phast& phast = Phast::Ref();
	TVector3 beam = phast.event.vParticle(vertex.InParticle()).ParInVtx(vertex.MyIndex()).Mom3();
	TVector3 vec_direction = beam.Unit();
	TVector3 vec_origin(0.,0.,1.); // this to be changed if somebody does a RPD alignment
	// get the angle of the vector
	double angle = vec_origin.Angle(vec_direction);
	// get the rotation axis perpendicular to the plane between these both
	TVector3 vec_rotation = vec_origin.Cross(vec_direction);
	vec_rotation = vec_rotation.Unit();
	// rotate around this axis by the given angle
	TVector3 proton_direction(0.,0.,0.);
	proton_direction.SetPtThetaPhi(1.,fTheta,fPhi);
	proton_direction.Rotate(-angle, vec_rotation);
	fTheta=proton_direction.Theta();
	fPhi=proton_direction.Phi();
}


bool
RPD::IsDT0(const PaEvent& event)
{
	// check whether current PaEvent is triggered by the DT0 trigger

	// for real data use the trigger mask
	if (!event.IsMC())
		return ((event.TrigMask()&0x1) == 0x1); 

	// only Monte Carlo data arrives here

	// only do the simulation for a new event
	// FIXME the same as in Decode for unique event number applies here
	if (event.UniqueEvNum() == fLastDT0SimulatedEvent)
		return fDT0Trigger;

	// FIXME Hodoscope Veto System (Vtot) to be treated! missing!
	const bool hit_HodVeto(false);

	bool hit_BC(false);
	bool hit_FI01X(false);
	bool hit_BK1(false);
	bool hit_BK2(false);
	bool hit_SW(false);

	//loop over all MC hits of this event
	for (std::vector<PaMChit>::const_iterator hit=event.MChits().begin(); hit!=event.MChits().end(); ++hit) {
		// search for alternative beam trigger (aBT) (beam counter (HK03X1__) and FI01 X dynode (FI01X1__))
		if (hit->DetRef().Name() == "HK03X1__")
			hit_BC = true;
		if (hit->DetRef().Name() == "FI01X1__")
			hit_FI01X = true;

		// search for hits in either beam killer (BK1 (HK01X1__) or BK2 (HK02X1__))
		if (hit->DetRef().Name() == "HK01X1__")
			hit_BK1 = true;
		if (hit->DetRef().Name() == "HK02X1__")
			hit_BK2 = true;

		// search for hits in the Sandwich detector (HH01X1__ to HH12X1__)
		if (hit->DetRef().Name().substr(0, 2) == "HH" && hit->DetRef().Name().substr(4) == "X1__")
			hit_SW = true;
	}

	// decode (and simulate) current event
	Decode(event);

	// search for RPD trigger
	bool hit_RPD(false);
	for (int iA = 0; iA < 12; iA++) {
		for (int ilB = -1; ilB < 2; ilB++) {
			int iB = 2 * iA + ilB;
			if (iA == 0 && ilB == -1)
				iB = 23;

			for (int iHitAup = 0; iHitAup < fDataSentAupt[iA]; iHitAup++) {
				for (int iHitAdo = 0; iHitAdo < fDataSentAdot[iA]; iHitAdo++) {
					for (int iHitBup = 0; iHitBup < fDataSentBupt[iB]; iHitBup++) {
						for (int iHitBdo = 0; iHitBdo < fDataSentBdot[iB]; iHitBdo++) {
							// FIXME thresholds to be tuned
							if (fAdcAdol[iA] < 10 || fAdcAupl[iA] < 10)
								if (fAdcAdoh[iA] < 10 || fAdcAuph[iA] < 10)
									continue;
							if (fAdcBdol[iB] < 10 || fAdcBupl[iB] < 10)
								if (fAdcBdoh[iB] < 10 || fAdcBuph[iB] < 10)
									continue;
							hit_RPD = true;
						}
					}
				}
			}
		}
	}

	// combination for alternative beam trigger
	const bool hit_aBT(hit_BC && hit_FI01X);

	// combination of beam killers (BK1 && BK2)
	const bool hit_BK(hit_BK1 && hit_BK2);

	// final combination for DT0 (aBT && RPD && !(BK || SW || Veto))
	// FIXME take time coincidence into account, for the moment simulations
	// FIXME with pile-up will not trigger DT0 properly
	fDT0Trigger = hit_aBT && hit_RPD && !(hit_BK || hit_SW || hit_HodVeto);

	// store event number
	fLastDT0SimulatedEvent = event.UniqueEvNum();

	// store decision for individual components
	fDT0Components[DT0Component_RPD]     = hit_RPD;
	fDT0Components[DT0Component_aBT]     = hit_aBT;
	fDT0Components[DT0Component_BC]      = hit_BC;
	fDT0Components[DT0Component_FI01X]   = hit_FI01X;
	fDT0Components[DT0Component_BK]      = hit_BK;
	fDT0Components[DT0Component_BK1]     = hit_BK1;
	fDT0Components[DT0Component_BK2]     = hit_BK2;
	fDT0Components[DT0Component_SW]      = hit_SW;
	fDT0Components[DT0Component_HodVeto] = hit_HodVeto;

	return fDT0Trigger;
}


bool
RPD::IsDT0()
{
	return IsDT0(Phast::Ref().event);
}


double
RPD::CorrectPhi(const double Phi, const PaVertex& /*vertex*/) const
{
//   needs to be corrected / deactivated for the moment

//   double r = 75.;	//Radius = 750 mm
//   double x = vertex.X(), y = vertex.Y();
//   double u = TMath::Cos(Phi), v = TMath::Sin(Phi);
// 
//   double t = (r*r-y*y)*u*u+(r*r-x*x)*v*v+2*u*v*x*y;
//   t = TMath::Sqrt(t);
//   t = -t - u*x - v*y;
//   t = t/(u*u+v*v);
//   
//   TVector3 vec(x+t*u, y+t*v, 0);
//
//   return vec.Phi();
   return Phi;
}


double
RPD::AngleWeight(double Phi) const
{
   TH1D *Angles = new TH1D("Angles","Angle",108,-3.141593,3.141593);
   Angles->SetBinContent(1,0.9629771);
   Angles->SetBinContent(2,0.9657779);
   Angles->SetBinContent(3,0.9246591);
   Angles->SetBinContent(4,0.8815469);
   Angles->SetBinContent(5,0.9038115);
   Angles->SetBinContent(6,0.9521089);
   Angles->SetBinContent(7,0.9682258);
   Angles->SetBinContent(8,0.9728625);
   Angles->SetBinContent(9,0.9705696);
   Angles->SetBinContent(10,0.9676951);
   Angles->SetBinContent(11,0.9618725);
   Angles->SetBinContent(12,0.9032147);
   Angles->SetBinContent(13,0.8395902);
   Angles->SetBinContent(14,0.8594779);
   Angles->SetBinContent(15,0.9276961);
   Angles->SetBinContent(16,0.9750057);
   Angles->SetBinContent(17,0.9817806);
   Angles->SetBinContent(18,0.9836317);
   Angles->SetBinContent(19,0.9862065);
   Angles->SetBinContent(20,0.9949239);
   Angles->SetBinContent(21,0.985653);
   Angles->SetBinContent(22,0.9531195);
   Angles->SetBinContent(23,0.924263);
   Angles->SetBinContent(24,0.9483279);
   Angles->SetBinContent(25,0.9692644);
   Angles->SetBinContent(26,0.9798989);
   Angles->SetBinContent(27,0.9794876);
   Angles->SetBinContent(28,0.9844519);
   Angles->SetBinContent(29,0.9992052);
   Angles->SetBinContent(30,0.9940022);
   Angles->SetBinContent(31,0.9749219);
   Angles->SetBinContent(32,0.9275844);
   Angles->SetBinContent(33,0.9443691);
   Angles->SetBinContent(34,0.9859272);
   Angles->SetBinContent(35,0.9996039);
   Angles->SetBinContent(36,1);
   Angles->SetBinContent(37,0.9892842);
   Angles->SetBinContent(38,0.9826287);
   Angles->SetBinContent(39,0.9759224);
   Angles->SetBinContent(40,0.9532084);
   Angles->SetBinContent(41,0.8910261);
   Angles->SetBinContent(42,0.8669815);
   Angles->SetBinContent(43,0.9172901);
   Angles->SetBinContent(44,0.9547904);
   Angles->SetBinContent(45,0.9562657);
   Angles->SetBinContent(46,0.9532871);
   Angles->SetBinContent(47,0.9557096);
   Angles->SetBinContent(48,0.9540514);
   Angles->SetBinContent(49,0.9468983);
   Angles->SetBinContent(50,0.9105025);
   Angles->SetBinContent(51,0.8879155);
   Angles->SetBinContent(52,0.9091719);
   Angles->SetBinContent(53,0.9396917);
   Angles->SetBinContent(54,0.9391559);
   Angles->SetBinContent(55,0.9355146);
   Angles->SetBinContent(56,0.9442752);
   Angles->SetBinContent(57,0.9436759);
   Angles->SetBinContent(58,0.9373734);
   Angles->SetBinContent(59,0.9086793);
   Angles->SetBinContent(60,0.8873289);
   Angles->SetBinContent(61,0.8977248);
   Angles->SetBinContent(62,0.9399736);
   Angles->SetBinContent(63,0.9541835);
   Angles->SetBinContent(64,0.9491278);
   Angles->SetBinContent(65,0.9463752);
   Angles->SetBinContent(66,0.9384678);
   Angles->SetBinContent(67,0.9096468);
   Angles->SetBinContent(68,0.883751);
   Angles->SetBinContent(69,0.8828547);
   Angles->SetBinContent(70,0.9070059);
   Angles->SetBinContent(71,0.9218684);
   Angles->SetBinContent(72,0.9212234);
   Angles->SetBinContent(73,0.9198598);
   Angles->SetBinContent(74,0.9241106);
   Angles->SetBinContent(75,0.909949);
   Angles->SetBinContent(76,0.8540083);
   Angles->SetBinContent(77,0.8001981);
   Angles->SetBinContent(78,0.8303319);
   Angles->SetBinContent(79,0.885135);
   Angles->SetBinContent(80,0.8861685);
   Angles->SetBinContent(81,0.8778523);
   Angles->SetBinContent(82,0.8810518);
   Angles->SetBinContent(83,0.8901831);
   Angles->SetBinContent(84,0.8826414);
   Angles->SetBinContent(85,0.8608339);
   Angles->SetBinContent(86,0.8330946);
   Angles->SetBinContent(87,0.8742693);
   Angles->SetBinContent(88,0.8956578);
   Angles->SetBinContent(89,0.8979102);
   Angles->SetBinContent(90,0.8976334);
   Angles->SetBinContent(91,0.9013763);
   Angles->SetBinContent(92,0.9081257);
   Angles->SetBinContent(93,0.8969198);
   Angles->SetBinContent(94,0.8598106);
   Angles->SetBinContent(95,0.8537366);
   Angles->SetBinContent(96,0.9115462);
   Angles->SetBinContent(97,0.9320307);
   Angles->SetBinContent(98,0.9442675);
   Angles->SetBinContent(99,0.9472436);
   Angles->SetBinContent(100,0.9551941);
   Angles->SetBinContent(101,0.9529265);
   Angles->SetBinContent(102,0.9310556);
   Angles->SetBinContent(103,0.8905183);
   Angles->SetBinContent(104,0.8991417);
   Angles->SetBinContent(105,0.9373784);
   Angles->SetBinContent(106,0.9545542);
   Angles->SetBinContent(107,0.9563013);
   Angles->SetBinContent(108,0.9557502);
   Angles->SetBinContent(109,4.062873e-05);
   Angles->SetMinimum(0);
   Angles->SetMaximum(1.1);

   while(Phi>2* TMath::Pi()) Phi = Phi - 2*TMath::Pi();
   while(Phi<0) Phi = Phi + 2*TMath::Pi();

   Int_t bin = Phi/(2*TMath::Pi())*108;

   return Angles->GetBinContent(bin);
}


double
RPD::AngleWeight(const TLorentzVector& vec) const
{
  return AngleWeight(vec.Phi());
}
