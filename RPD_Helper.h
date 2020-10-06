#ifndef Phast_RPD_Helper_h
#define Phast_RPD_Helper_h
/*
 RPD Helper for Phast

 comments & questions to
 Johannes Bernhard <Johannes.Bernhard@cern.ch>

 Enjoy!

 changes:

 09-02-12 18:50 jbernhar : Changed defintion of phi according to root definition (i.e. -pi < phi <= pi)
 09-02-18 00:18 jbernhar : Implemented Etienne's corrections (energy loss in the target, calibration constants, vertex correction)
 09-02-21 18:30 jbernhar : Introduced multiple RPD Tracks
 09-02-23 18:50 jbernhar : guarded several variables against division by 0
 09-02-26 23:53 jbernhar : changed class structure to be according to Phast
 09-04-01 16:29 jbernhar : fixed bug in energy corrections, Tobi removed RPD:: in front of member functions
 09-04-15 13:19 jasinski : formatted appearance, added default constructor, created new method: RPD::Initialize() to be called by constructors, Reset() introduced and called in Search() for the case that RPD is static
 09-04-21 09:12 jasinski : minor changes to some methods
 09-04-21 21:59 jbernhar : some changes to Reset(), splitting into RPD_Helper.h and RPD_Helper.cc, energy correction now true on default, doxygen doc added
 09-04-22 14:45 jbernhar : fixed possible memset bug in Initialize()
 09-04-22 16:00 schluter : Fixed all memsets, cleanup header use.
 09-04-22 18:45 jbernhar : added vertex time + calibration constants
 09-04-24 15:18 schluter : Fix boundary cases in energy correction.
 09-05-05 12:01 schluter : Make return types const& where it makes sense,
                           move one-line functions to header.
 09-05-08 21:25 schluter : Do away with requirement of only one static instance.
 09-11-11 22:20 jbernhar : Changed constructor and iBestProtonTrack like proposed by Sergei Gerassimov
 09-12-06 14:00 jasinski : Including simple MC simulation for MChits
 10-01-11 16:00 jasinski : effective light speed in the elements is taken into account
 10-01-14 10:00 jasinski : bug in getMCdigits found. MC digits of last events were kept
 10-02-02 12:00 jasinski : merging code pieces in both Search() methods
 10-02-03 12:00 jasinski : tdc smearing included
 10-02-08 12:00 jasinski : db for t0 calibration constants included
 10-02-12 12:00 jasinski : replacing "slat" by "slab", several changes due to new ComGeant RPD description
 10-02-25 12:00 jasinski : introducing singleton structure. Multiple intances will be kept allowed.
 14-07-01 16:00 hubers   : introduce polynomial description of energyloss
 */

/*!
  \class RPD
  \brief RPD Helper Class

  Class to treat COMPASS RPD information.

  \author Johannes.Bernhard@cern.ch
*/

#include <map>
#include <string>
#include <vector>

#include <TLorentzVector.h>

class PaDigit;
class PaEvent;
class PaVertex;

struct Calib_const {
	// [ring:<A|B>][element:<0..11|0..24>][pm:<up|down>]
	int	  runnumb;
	float t0[2][24][2];
	// [ring:<A|B>][element:<0..11|0..24>]
	float positionOffset[2][24];
	// [element:<0..11>][element:<0..24>]
	float zVertexOffset[12][24];
	// [ring:<A|B>][element:<0..11|0..24>]
	float dEcusp[2][24];
	// [element:<0..11>]
	float calTA[12];
	float tofGlobalOffset;

	/******** changing the db structure needs also changes in **********/
	/*
	 *
	 * Load_DB_entry();
	 * Write_constants();
	 * RPD_DB/calib.db
	 * NENTRIES_CALIB_DB
	 *
	 * Recommended way to change the layout of a DB
	 *
	 * 1. Add the variables you need to this struct
	 * 2. Change Write_constants() so it writes the DB according to the new structure
	 * 	- add the new variables to the output stream
	 *  - remove old ones if not needed anymore
	 * 3. Run the program by loading the old DB and call Save_DB() once
	 *  - "database.db" is written containing the DB with the new structure
	 * 4. Replace the old "RPD_DB/calib.db" by the entries of "database.db"
	 * 5. Change Load_DB_entry() according to the new layout of your DB
	 * 	- have a look into Write_constants() to see the order in the array of floats
	 * 6. Remove the values in the struct that are not needed anymore
	 *  - The same values you removed from the output stream in (2.)
	 * 7. Change NENTRIES_CALIB_DB by counting all values in this struct
	 * 8. Try your new DB out!
	 *
	 * Questions? -> jasinski@kph.uni-mainz.de (or Promme@web.de if not valid anymore)
	 */
};

// this number must fit the number of entries
// in struct Calib_const in order
// to perform some db checks before loading it
// this is the number of entries that will be read
// per calibration entry
// nota bene that sometimes the ring A(0) has a float with 24 entries
// but only 12 calibration constants will be read
const int NENTRIES_CALIB_DB = 1 + 2*12 + 2*24 + 12 + 24 + 12*24 + 12 + 24 + 12 + 1;

class RPD {

public:

	// individual DT0 components
	enum DT0Components {
		DT0Component_RPD,
		DT0Component_aBT, DT0Component_BC, DT0Component_FI01X,
		DT0Component_BK, DT0Component_BK1, DT0Component_BK2,
		DT0Component_SW,
		DT0Component_HodVeto
	};

	//! Get a reference to the RPD singleton
	// usage: RPD& rpd = RPD::Instance();
	static RPD& Instance();

private:

	// pointer to the class in case of singleton use
	static RPD* instance;

	//! Constructor
	// private so that a RPD object cannot be created directly, but can only be used as a singleton
	RPD();

	//! Reset all stored data
	void Reset();

	//! Decode an event
	// In case of real data, get the DAQ digits, and decode those, in case
	// of Monte Carlo, simulate the detector response and then decode the
	// result.
	void Decode(const PaEvent& event);

	// event number of event decoded last, only decode each event once
	long long fLastDecodedEvent;

	// decoded TDC data
	// up to five times can be measured per slab
	int fDataSentAupt[12];
	float fTdcAup[12][5];
	int fDataSentAdot[12];
	float fTdcAdo[12][5];
	int fDataSentBupt[24];
	float fTdcBup[24][5];
	int fDataSentBdot[24];
	float fTdcBdo[24][5];

	// decoded ADC data
	bool fDataSentAupql[12];
	float fAdcAupl[12];
	bool fDataSentAdoql[12];
	float fAdcAdol[12];
	bool fDataSentBupql[24];
	float fAdcBupl[24];
	bool fDataSentBdoql[24];
	float fAdcBdol[24];
	bool fDataSentAupqh[12];
	float fAdcAuph[12];
	bool fDataSentAdoqh[12];
	float fAdcAdoh[12];
	bool fDataSentBupqh[24];
	float fAdcBuph[24];
	bool fDataSentBdoqh[24];
	float fAdcBdoh[24];

	//! Reset the reconstructed data
	void ResetClustering();

	// event number of event clustered last, only cluster each event once with the same options
	long long fLastClusteredEvent;

	// whether the momenta were corrected during the last call to Clusterize
	bool fLastCorrectMomentum;

	// whether the angles were corrected during the last call to Clusterize
	bool fLastCorrectAngles;

	//! Reset the DT0 trigger simulation
	void ResetDT0Simulation();

	// event number of last event the DT0 trigger was simulated for
	long long fLastDT0SimulatedEvent;

	// DT0 trigger decision
	bool fDT0Trigger;

	// trigger decision for individual components
	std::map<DT0Components, bool> fDT0Components;

	int nTracks;

	bool isMC; // is the last reconstructed event an MC simulated one?

	// database for RPD calibration values
	std::map<int, Calib_const> calib_constants_db;
	// usually it should be in RPD_Helper.cc as a global variable
	// but it does not work! I've been fighting 5 hours and it looks
	// like a problem of the gcc. The DB is getting lost!
	// Ask P.Jasinski(me) for details

	// loading all DB constants given by an array of values
	bool Load_DB();

	// load one DB constants entry starting from one
	// entry in a array of values
	Calib_const Load_DB_entry(const float db_entries[]) const;

	// save DB to file
	bool Save_DB(const std::string& filename = "database.db") const;

	// append calibration constants to a file
	// if no filename is given a new
	// file containing the runnumber is created
	bool Write_constants(const Calib_const& calib_const, const std::string& filename = "") const;

	// Write calibration constants to DB
	// if calib_const.runnumb exists in the DB the entry is overwritten
	// if calib_const.runnumb does not exist a new entry is created
	bool Set_DB_entry(const Calib_const& calib_const);

	// search in DB for constants
	void FindInDB(const int runnumb);

	std::string ChannelToPMT(const std::string& RpdPlanes, const int chId, int &thisPMT) const;

	// Simulation of Digits via MC requires to know the channel number
	// of SADC and TDC equipment. This is given by the properties
	// in:
	// ring: 0 or 1 for A or B
	// element: 0..11 or 0..23 for ring 0 or 1
	// tdc : true if false : sadc
	//
	// out:
	// corresponding channel_pmup (upstream) and channel_pmdo (downstream)
	//
	// the detector name is returned
	std::string GetChannel(const int ring, const int element, const bool tdc, int &channel_pmup, int &channel_pmdo) const;

	// FIXME make pRpd const
	double correctEnergyLoss(double pRpd, const double sinTheta, const double cosTheta, const double phi, const double xVertex, const double yVertex, const double zVertex, const bool debug=false) const;

	double correction_energy(const double p, const int mat, const double thickness, const bool debug) const;
	double correction_energyPol(const double p, const int mat, const double thickness, const bool debug) const;

	double terpol3(const double xx1, const double xx2, const double xx3, const double yy1, const double yy2, const double yy3, const double x) const;

	float treatADC(const PaDigit& d) const;

	// create MC Digits out of MCHits in the RPD
	// TDC as well as ADC digits are created
	// this part will be moved to Coral sooner or later
	const std::vector<PaDigit>& Construct_MC_Digits(const PaEvent& e);

	// this is needed to distinguish between MC and real data
	// in case of MC one needs to construct the signals
	const std::vector<PaDigit>& Get_Digits(const PaEvent& e);

	// Since MC does not return the element that was hit yet,
	// the corresponding element is calculated by knowing the
	// RPD ring and the phi angle
	// Ring A = 0; Ring B = 1
	// if angle is determined by x and y position, one
	// should call atan2(y,x) to determine the angle phi
	// thus phi's range is expected to be [0,pi] ]-pi,0[
	// wrong input returns -1
	int Phi_to_element(const double phi, const int ring) const;

	// function to correct for rotated RPD system with respect to the beam direction
	void CorrectAngles(const PaVertex& vertex, double& fTheta, double& fPhi) const;

	int fBestProtonTrack;
	int fBestMCProtonTrack;

	std::vector<TLorentzVector> fRPDTrack; // RPD Tracks are stored here
	std::vector<TLorentzVector> fRPDMCTrack; // MC RPD Tracks are stored here

	std::vector<double> fzT; // z-component of the primary vertex seen by the RPD for every track

	std::vector<double> ftT; // RPD vertex time

	std::vector<double> dEA; // energy loss in A for every track

	std::vector<double> dEB;// energy loss in B for every track

	std::vector<double> dEA_MC; // MC energy loss in A for the best matching hit

	std::vector<double> dEB_MC; // MC energy loss in B for the best matching hit

	std::vector<double> zA_vec; //z Position of Ring A hit

	std::vector<double> zB_vec; //z Position of Ring B hit

	std::vector<bool> fCalibratedEvent; // Helper to see whether beta is between 0 and 1 (for every track)

	std::vector<std::pair<int, int> > hits; // pair of A and B hits for every track

	bool something; // RPD saw something

/********************** coral part ************************
 *  starting to introduce some coral like structure
 *  to keep the code of the helper class similar to it
 *  that way one day we may move easier to coral
 */


	/* decoding a F1 type digit
	 * In		: copy of ChipF1 digit
	 * Out		:
	 * returns	:
	 * comments	:
	 *
	 */
	void DecodeChipF1Digit(const PaDigit& digit);//(const CS::ChipF1::Digit &digit);

	/* decoding a SADC type digit
	 * In		: copy of ChipSADC digit
	 * Out		:
	 * returns	:
	 * comments	:
	 *
	 */
	void DecodeChipSADCDigit(const PaDigit& digit);//(const CS::ChipSADC::Digit &digit);

	/* decoding one digit
	 * In		: copied digit
	 * Out		:
	 * returns	:
	 * comments	:
	 * - calling the decoding of F1 TDC digits or SADC digits
	 *   depending on the type of this Digit
	 * - implementation as public since it is given to be public
	 *   as virtual method in base class. Coral will probably call
	 *   directly DecodeChipDigits?!
	 */
	void DecodeChipDigit(const PaDigit& digit);//const CS::Chip::Digit &digit);

	/* decoding a bunch of Digits
	 * In		: copied vector of Digits
	 * Out		:
	 * returns	:
	 * comments	:
	 */
	void DecodeChipDigits(const std::vector<PaDigit>& digits);//const CS::Chip::Digits &digits);

	/* clusterization of decoded data per event
	 * In		:
	 * Out		:
	 * returns	:
	 * comments	:
	 * Most probably it won't be easy to move this part to
	 * coral since the corresponding primary vertex is not given
	 * there.
	 * In coral Clusterize() has no arguments!
	 * So we will have to separate the vertex part and keep it
	 * in the helper class.
	 */

	void Clusterize(const PaVertex& vertex, const bool correct_momentum, const bool correct_angles);

	/* in case of debugging mode and MC data
	 * best reconstructed proton and vertex
	 * is compared to MC primary vertex and MC proton track
	 */
	void CheckMCresolutions() const;

/************************** end of coral part *****************/

public:

	//! Vector of all RPD tracks
	const std::vector<TLorentzVector>& vTracks() const { return fRPDTrack; }

	//! std::vector of all MC RPD tracks
	const std::vector<TLorentzVector>& vMCTracks() const { return fRPDMCTrack; }

	//! Number of RPD tracks
	int nTrack() const { return nTracks; }

	//! Z component of the primary vertex seen by the RPD
	const std::vector<double>& PV_Z() const { return fzT; }

	//! RPD Vertex time
	const std::vector<double>& PV_T() const { return ftT; }

	//! Vector of pairs of Ring A and B hits
	const std::vector<std::pair<int, int> >& Hits() const { return hits; } // returns pairs of A and B hits

	//! Returns whether the RPD saw something
	bool HasTracks() const { return something; }

	//! Number of the best RPD track
	int iBestProtonTrack() const { return fBestProtonTrack; }

	//! Number of the best RPD MC track (not the track in PaEvent::vMCtrack())
	int iBestProtonMCTrack() const { return fBestMCProtonTrack; }

	//! Energy loss in Ring A for every track
	const std::vector<double>& dE_A() const { return dEA; }

	//! Energy loss in Ring B for every track
	const std::vector<double>& dE_B() const { return dEB; }

	//! Position of hits in Ring A
	const std::vector<double>& ZA() const { return zA_vec; }

	//! Position of hits in Ring B
	const std::vector<double>& ZB() const { return zB_vec; }

	//! Search for RPD tracks. PV_x PV_y PV_z are the x and y coordinates for a given primary vertex if one wishes to use the energy correction function (z coordinate necessary for nuclear targets)
	void Search(const PaEvent& e, const PaVertex& vertex, const bool correct_momentum = true, const bool correct_angles = true);
	void Search(const PaVertex& vertex, const bool correct_momentum = true, const bool correct_angles = true);
	void Search(const PaEvent& e, const bool correct_momentum = true, const double PV_x = 0., const double PV_y = 0., const double PV_z = -777.);
	void Search(const bool correct_momentum = true, const double PV_x = 0., const double PV_y = 0., const double PV_z = -777.);

	//! produce a calibration file on a run by run basis
	// the produced file has to be merged with the existing one in the RPD_DB folder
	// last run is never calibrated (for 1 run no calibration)
	void Calibrate(const PaEvent& e, const std::string& path = "");

	//! example of how to write calibration constants out into a DB
	// by using previously calibrated constants
	// use phast -f option to run several times over the same events
	void Calibrate(const PaEvent& e, const int pass = 1, const std::string& path = "");

	//! Dumps Track information on stdout
	void Dump() const;

	//! Check if tracks have a meaningful calibration, obsolete
	const std::vector<bool>& CheckCalibration() const { return fCalibratedEvent; }

	//! Provides Element Number for a given phi angle and Ring ("A","B")
	int MapPhiElement(const double phi, const std::string& ring) const;

	//! Check for DT0 trigger in the current event (both MC and RD)
	bool IsDT0(const PaEvent& event);
	bool IsDT0();

	//! Get the decision for the individual DT0 comonents
	const std::map<DT0Components, bool>& GetDT0Components() const {return fDT0Components; }

	//! Retrieve Phi angle corrected for vertex position
	double CorrectPhi(const double Phi, const PaVertex& vertex) const;

	//! Return scaling value for a certain phi angle (using a lookup table from a histogram)
	// FIXME make Phi const
	double AngleWeight(double Phi) const;
	double AngleWeight(const TLorentzVector& vec) const;

};

#endif
