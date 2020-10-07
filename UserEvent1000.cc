# include <cassert>

#include "TCanvas.h"
#include "TTree.h"

#include "ECALDefinition.h"
#include "BeamDefinition.h"
#include "PhotonDefinition.h"
#include "RecoilDefinition.h"
#include "RPDDefinition.h"
#include "ScatteredDefinition.h"
#include "SelectionHelper.h"
#include "VertexDefinition.h"


void
UserEvent1000(PaEvent& event)
{
	// Tree
	static TTree* tree(NULL);

	// Scattered particles
	static int C_NUMBER_PARTICLES_CHARGED = 3;
	static int C_NUMBER_PARTICLES_NEUTRAL = 2;

	// Target region
	static double C_TARGET_UP   = -260;  // upstream boundary [cm]
	static double C_TARGET_DOWN =  160;  // downstream boundary [cm]

	// Tools
	static SelectionHelper*     selectionHelper;
	static VertexDefinition*    vertexDefinition;
	static ECALDefinition*      ecalDefinition;
	static RPDDefinition*       rpdDefinition;
	static BeamDefinition*      beamDefinition;
	static ScatteredDefinition* scatteredDefinition;
	static PhotonDefinition*    photonDefinition;
	static RecoilDefinition*    recoilDefinition;

	static bool first(true);

	// Initial configuration
	if (first) {
		first = false;
		tree  = new TTree("UserEvent1000", "Events");
		assert(tree);

		selectionHelper     = new SelectionHelper    (*tree);
		vertexDefinition    = new VertexDefinition   (*tree);
		ecalDefinition      = new ECALDefinition     (*tree);
		rpdDefinition       = new RPDDefinition      (*tree);
		beamDefinition      = new BeamDefinition     (*tree);
		scatteredDefinition = new ScatteredDefinition(*tree, C_NUMBER_PARTICLES_CHARGED, C_NUMBER_PARTICLES_NEUTRAL);
		photonDefinition    = new PhotonDefinition   (*tree);
		recoilDefinition    = new RecoilDefinition   (*tree);

		Phast::Ref().HistFileDir("UserEvent");
	}
	selectionHelper->reset();

	// --- MC truth ---
	if (event.IsMC()) {
		for (int iMCvertex = 0; iMCvertex < event.NMCvertex(); ++iMCvertex) {
			const PaMCvertex& vertexMC = event.vMCvertex(iMCvertex);
			if (vertexMC.IsPrimary() and not vertexMC.IsPileup()) {
				vector<PaMCtrack> tracksMC    = scatteredDefinition->getScatteredParticleMCtracks(event, vertexMC);
				PaMCtrack         beamTrackMC = beamDefinition->getBeamTrackMC                   (event, vertexMC);

				vertexDefinition->fillMC   (event   , vertexMC);
				beamDefinition->fillMC     (beamTrackMC);
				scatteredDefinition->fillMC(tracksMC, vertexMC);
				photonDefinition->fillMC   (event   , vertexMC);
				recoilDefinition->fillMC   (event   , vertexMC);
			}
		}
	}

	// Statistics
	selectionHelper->fillStatistics("Total");

	// DT0 trigger
	if (not selectionHelper->isDT0()) {
		if (event.IsMC()) {
			selectionHelper->handleMCEvent(*tree, false, event);
		}
		return;
	}
	selectionHelper->fillStatistics("DT0");

	// At least one primary vertex in event
	const int numberPrimaryVertex = vertexDefinition->getNumberPrimaryVertices(event);
	if (numberPrimaryVertex > 1) {
		if (event.IsMC()) {
			selectionHelper->handleMCEvent(*tree, false, event);
		}
		return;
	}
	selectionHelper->fillStatistics("One vertex");

	// Number of neutral clusters in event
	//TODO don't we need minimum 4 photons?
	if (ecalDefinition->getNumberNeutralClusters(event) < 2) {
		if (event.IsMC()) {
			selectionHelper->handleMCEvent(*tree, false, event);
		}
		return;
	}
	selectionHelper->fillStatistics("ECAL Cluster");

	// Loop over vertices
	// fixed multiple vertices
	bool accepted = false;
	const PaVertex& vertex = event.vVertex(event.iBestPrimaryVertex());
	
	// Vertex cut around target
	if (not (C_TARGET_UP <= vertex.Z() and vertex.Z() <= C_TARGET_DOWN)) {
		return;
	}
	selectionHelper->fillStatistics("Target");


	// Number of outgoing particles in primary vertex
	if (vertex.NOutParticles() != C_NUMBER_PARTICLES_CHARGED) {
		return;
	}
	selectionHelper->fillStatistics("Outgoing");

	// Tracks in RPD
	rpdDefinition->initRPD(event);
	if (not rpdDefinition->cutHasTracks()) {
		return;
	}
	selectionHelper->fillStatistics("RPD tracks");

	// Best proton in RPD
	if (not rpdDefinition->cutHasBestProton()) {
		return;
	}
	selectionHelper->fillStatistics("RPD proton");

	const PaParticle& beamParticle = beamDefinition->getBeamParticle(event, vertex);

	// Scattered particle
	const vector<PaParticle>& scatteredParticles = scatteredDefinition->getScatteredParticles(event, vertex);

	// Check charge conservation
	if (not selectionHelper->isCorrectCharge(beamParticle, scatteredParticles)) {
		return;
	}
	selectionHelper->fillStatistics("Charge");


	// --- FILL ACCEPTED EVENTS ---

	// Accepted
	accepted = true;

	selectionHelper->fill    (event);
	vertexDefinition->fill   (event, vertex);
	ecalDefinition->fill     (event);
	rpdDefinition->fill      (event, vertex);
	scatteredDefinition->fill(event       , scatteredParticles, vertex);
	beamDefinition->fill     (event);
	

	// Always fill tree if MC - use isAcceptedMC branch to distinguish
	if (event.IsMC()) {
		selectionHelper->handleMCEvent(*tree, accepted, event);
	} else if (accepted) {
		// Fill real data tree if event is accepted
		tree->Fill();
	}

	Phast::Ref().h_file = tree->GetCurrentFile();
}
