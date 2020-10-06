#include "VertexDefinition.h"


VertexDefinition::VertexDefinition(TTree& tree)
	: _phast(Phast::Ref()),
	  _vertexVector    (NULL),
	  _vertexVectorMC  (NULL),
	  _numberVertices  (0),
	  _numberOutgoing  (0),
	  _vertexX         (0),
	  _vertexY         (0),
	  _vertexZ         (0),
	  _numberMCVertices(0),
	  _numberMCOutgoing(0),
	  _vertexMCX       (0),
	  _vertexMCY       (0),
	  _vertexMCZ       (0)
{
	_phast.h_file->cd();

	tree.Branch("Vertex_Number",      &_numberVertices,   "numberVertices/I");
	tree.Branch("Vertex_Outgoing",    &_numberOutgoing,   "numberOutgoing/I");
	tree.Branch("Vertex_X",           &_vertexX,          "vertexX/D");
	tree.Branch("Vertex_Y",           &_vertexY,          "vertexY/D");
	tree.Branch("Vertex_Z",           &_vertexZ,          "vertexZ/D");
	tree.Branch("Vertex_MC_Number",   &_numberMCVertices, "numberMCVertices/I");
	tree.Branch("Vertex_MC_Outgoing", &_numberMCOutgoing, "numberMCOutgoing/I");
	tree.Branch("Vertex_MC_X",        &_vertexMCX,        "vertexMCX/D");
	tree.Branch("Vertex_MC_Y",        &_vertexMCY,        "vertexMCY/D");
	tree.Branch("Vertex_MC_Z",        &_vertexMCZ,        "vertexMCZ/D");
}


VertexDefinition::~VertexDefinition()
{
	if (_vertexVector) {
		delete _vertexVector;
	}
	if (_vertexVectorMC) {
		delete _vertexVectorMC;
	}
}


void
VertexDefinition::fill(const PaEvent&  event,
                       const PaVertex& vertex)
{
	_numberVertices = getNumberPrimaryVertices(event);
	_numberOutgoing = getNumberOutgoing       (vertex);
	_vertexX        = vertex.X();
	_vertexY        = vertex.Y();
	_vertexZ        = vertex.Z();
	_vertexVector   = new TVector3(_vertexX, _vertexY, _vertexZ);
}


void
VertexDefinition::fillMC(const PaEvent&    eventMC,
                         const PaMCvertex& vertexMC)
{
	if (eventMC.IsMC()) {
		_numberMCVertices = getNumberPrimaryMCVertices(eventMC);
		_numberMCOutgoing = getNumberMCOutgoing       (vertexMC);
		_vertexMCX        = vertexMC.Pos(0);
		_vertexMCY        = vertexMC.Pos(1);
		_vertexMCZ        = vertexMC.Pos(2);
		_vertexVectorMC   = new TVector3(_vertexMCX, _vertexMCY, _vertexMCZ);
	}
}


int
VertexDefinition::getNumberPrimaryVertices(const PaEvent& event) const
{
	int numberPrimaryVertices = 0;
	for (int i = 0; i < event.NVertex(); ++i) {
		const PaVertex& paVertex = event.vVertex(i);
		if (paVertex.IsPrimary()) {
			numberPrimaryVertices++;
		}
	}
	return numberPrimaryVertices;
}


int
VertexDefinition::getNumberPrimaryMCVertices(const PaEvent& event) const
{
	int numberPrimaryVertices = 0;
	const vector<PaMCvertex>& mcVertices = event.vMCvertex();
	for (int i = 0; i < event.NMCvertex(); ++i) {
		const PaMCvertex& paMCvertex = mcVertices[i];
		if (paMCvertex.IsPrimary()) {
			numberPrimaryVertices++;
		}
	}
	return numberPrimaryVertices;
}
