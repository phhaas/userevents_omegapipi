#ifndef VERTEXDEFINITION_H
#define VERTEXDEFINITION_H

#include "TH1D.h"
#include "TH2D.h"

#include "Phast.h"

class VertexDefinition
{

public:

	VertexDefinition(TTree& tree);
	~VertexDefinition();

	void fill  (const PaEvent& event,   const PaVertex&   vertex);
	void fillMC(const PaEvent& eventMC, const PaMCvertex& vertexMC);

	int    getNumberPrimaryVertices(const PaEvent&  event)  const;


private:

	int    getNumberOutgoing       (const PaVertex& vertex) const { return vertex.NOutParticles(); }
	double getVertexX() const { return _vertexX; }
	double getVertexY() const { return _vertexY; }
	double getVertexZ() const { return _vertexZ; }

	int    getNumberPrimaryMCVertices(const PaEvent&    event)    const;
	int    getNumberMCOutgoing       (const PaMCvertex& vertexMC) const { return vertexMC.NMCtrack(); }
	double getVertexMCX() const { return _vertexMCX; }
	double getVertexMCY() const { return _vertexMCY; }
	double getVertexMCZ() const { return _vertexMCZ; }

	const TVector3& getVertexVector()   const { return *_vertexVector;   }
	const TVector3& getVertexVectorMC() const { return *_vertexVectorMC; }
	Phast& _phast;

	TVector3* _vertexVector;
	TVector3* _vertexVectorMC;

	// output tree variables
	int    _numberVertices;
	int    _numberOutgoing;
	double _vertexX;
	double _vertexY;
	double _vertexZ;
	int    _numberMCVertices;
	int    _numberMCOutgoing;
	double _vertexMCX;
	double _vertexMCY;
	double _vertexMCZ;

};

#endif  // VERTEXDEFINITION_H
