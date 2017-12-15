


#include "diff.h"
#include <GU/GU_Detail.h>
#include <GEO/GEO_Hedge.h>
#include <GEO/GEO_HedgeInterface.h>
#include <SYS/SYS_Math.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <OP/OP_Node.h>
#include <UT/UT_DSOVersion.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstdio>
#include <PRM/PRM_Template.h>
#include <OP/OP_AutoLockInputs.h>
#include <GEO/GEO_Point.h>
#include <GEO/GEO_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <iterator>
#include <vector>
#include <iostream>
#include <algorithm>	
#include <mesh.h>



static PRM_Name scalarfield("Field", "Field");

//____________________________________________________________________________________________________________________________ boringHdkstuff
void
newSopOperator(OP_OperatorTable *table)
{
	OP_Operator *op;
	op = new OP_Operator("ddg_TrigClosed", //internal name
		"ddg_TrigClosed.v1", 				 //UI Name
		DIFF_SOP::myConstructor,		 //How to build the SOP	
		DIFF_SOP::MyTemplateList,		 //list of parameters		
		1,								 //min number of sources
		1,								 //max number of inputs
		0,								 //local variables
		OP_FLAG_GENERATOR);				 //flag as a generator type		

	table->addOperator(op);
}


OP_Node * DIFF_SOP::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {

	return new DIFF_SOP(net, name, op);
}



static PRM_Name field("fieldop", "FieldOperation");
static PRM_Name fielparam[] =
{
	PRM_Name("grad",    "Gradient"),
	PRM_Name("div",     "Divergence"),
	PRM_Name("curl",    "Laplacian"),
	PRM_Name("hod",     "VectorDecomp"),
	PRM_Name("vecsource", "VectorSource"),
	PRM_Name(0)
};
enum
{
	DIFF_SOP_grad, DIFF_SOP_div, DIFF_SOP_lap, DIFF_SOP_hod, DIFF_SOP_vecsource
};

static PRM_ChoiceList fieldmenu(PRM_CHOICELIST_SINGLE, fielparam);

PRM_Template DIFF_SOP::MyTemplateList[] = {
	PRM_Template(PRM_STRING ,1, &scalarfield),
	PRM_Template(PRM_ORD, 1, &field,0, &fieldmenu),
	PRM_Template()
};

DIFF_SOP::DIFF_SOP(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op) {}

DIFF_SOP::~DIFF_SOP() {}


int
DIFF_SOP::fielchoice()
{
	return evalInt(field.getToken(), 0, 0.0f);
}


//_________________________________________________________________________________________________________________________________ VectorDecomp
void  DIFF_SOP::Laplacian(GA_ROHandleV3 vectorfield, OP_Context &context) {
	//http://courses.cms.caltech.edu/cs177/hmw/Hmw3.pdf

	GA_Attribute *vec1 = gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "a", 3);
	GA_RWHandleV3 hedgvec1(vec1);

	GA_Attribute *vec2 = gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "b", 3);
	GA_RWHandleV3 hedgvec2(vec2);
	GA_Attribute *vec3 = gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "c", 3);
	GA_RWHandleV3 hedgvec3(vec3);
	//_______________________________________________________________________ FaceIndicesEdges 
	int ptCount = gdp->getNumPoints();
	int faceCount = gdp->getNumPrimitives();

	//______________Vertices Matrix 
	Eigen::MatrixXd  V(ptCount, 3);
	int r = 0;
	for (GA_Iterator it(gdp->getPointRange(NULL)); !it.atEnd(); ++it)
	{
		UT_Vector3 p = gdp->getPos3(*it);
		V.row(r) = Eigen::Vector3d(p.x(), p.y(), p.z());
		r++;
	}
	//______________Indices Matrix and Edge List 
	Eigen::MatrixXi F(faceCount, 3);
	UT_Array< const GA_Primitive * >prims;
	gdp->getPrimitivesOfType(GA_PRIMPOLY, prims);
	r = 0;

	// init sort/stack arrays

	edgeList.setSize(0);
	UT_Array <UT_Vector2>stack;
	stack.setSize(0);

	// Loop trough faces, if twinhalfedge and findable in stack -> pass, 
	//if not. mark as visited by adding in stack array,
	//swap for corresponding twinhalfedge
	// 
	for (int i = 0; i<prims.size(); ++i) {
		const GA_Primitive *prim = prims(i);
		F.row(r) = Eigen::Vector3i(prim->getPointIndex(0), prim->getPointIndex(1), prim->getPointIndex(2));
		GEO_Hedge curr = polinterface->polyHedge(i);
		GEO_Hedge next = polinterface->nextPrimitiveHedge(curr);
		GEO_Hedge prev = polinterface->prevPrimitiveHedge(curr);
		UT_Array<int>::iterator it;
		if (polinterface->isPrimary(prev) != 1) {
			int dst = polinterface->dstPoint(prev);
			int src = polinterface->srcPoint(prev);
			bool it = stack.find(UT_Vector2(dst, src));
			if (it != 1) {
				stack.append(UT_Vector2(src, dst));
				edgeList.append(UT_Vector2(dst, src));
			}
		}
		else {
			edgeList.append(UT_Vector2(polinterface->dstPoint(prev), polinterface->srcPoint(prev)));

			it = find(primhedgelist.begin(), primhedgelist.end(), polinterface->hedgePoly(prev));
			if (it == primhedgelist.end())
				primhedgelist.append(polinterface->hedgePoly(prev));
		}
		if (polinterface->isPrimary(next) != 1) {
			int dst = polinterface->dstPoint(next);
			int src = polinterface->srcPoint(next);
			bool it = stack.find(UT_Vector2(dst, src));
			if (it != 1) {
				stack.append(UT_Vector2(src, dst));
				edgeList.append(UT_Vector2(dst, src));
			}
		}
		else {
			edgeList.append(UT_Vector2(polinterface->dstPoint(next), polinterface->srcPoint(next)));
			it = find(primhedgelist.begin(), primhedgelist.end(), polinterface->hedgePoly(next));
			if (it == primhedgelist.end())
				primhedgelist.append(polinterface->hedgePoly(next));
		}
		if (polinterface->isPrimary(curr) != 1) {
			int dst = polinterface->dstPoint(curr);
			int src = polinterface->srcPoint(curr);
			bool it = stack.find(UT_Vector2(dst, src));
			if (it != 1) {
				stack.append(UT_Vector2(src, dst));
				edgeList.append(UT_Vector2(dst, src));
			}
		}
		else {
			edgeList.append(UT_Vector2(polinterface->dstPoint(curr), polinterface->srcPoint(curr)));
			it = find(primhedgelist.begin(), primhedgelist.end(), polinterface->hedgePoly(curr));
			if (it == primhedgelist.end())
				primhedgelist.append(polinterface->hedgePoly(curr));
		}
		r++;
	}



	SparseMat	d0(edgeList.size(), gdp->getNumPoints());
	d0 = DIFF_SOP::deriative0(d0);

	SparseMat	d1(gdp->getNumPrimitives(), edgeList.size());
	d1 = DIFF_SOP::deriative1(d1);


	SparseMat	hstar2(edgeList.size(), edgeList.size());
	hstar2 = DIFF_SOP::star2(hstar2);


	SparseMat	hstar2invert(edgeList.size(), edgeList.size());
	hstar2invert = DIFF_SOP::star2inv(hstar2invert);


	//___________________________________________________________________________ 1-form
	//
	UT_Array<fpreal>oneform;
	covector(oneform, vectorfield);

	// __________________________________________________________________________ δdα = δω

	Eigen::Map<EigenMat>w(oneform.data(), oneform.size(), 1);

	SparseMat Laplace;
	EigenMat da(gdp->getNumPoints(), 1);
	EigenMat oneformscalar(edgeList.size(), 1);


	Eigen::SimplicialLDLT<SparseMat> solver;

	// δ =  d0dual *1
	// ∆ =  δ * doprimal
	Laplace = (d0.transpose() * hstar2 * d0);

	// 0-form L-matrix v*v
	solver.compute(Laplace);

	// ∆α =δω
	da = solver.solve(d0.transpose() * hstar2 *w);

	//δdα
	oneformscalar = d0 * da;

	// _________________________________________________________________________ dδβ = dw
	//
	SparseMat Laplace2;
	EigenMat twoformscalar(gdp->getNumPrimitives(), 1);
	EigenMat beta(edgeList.size(), 1);
	// 1/diagonalElement


	Eigen::SparseLU<SparseMat> lu;


	//δ = *1(-1) d1dual
	//∆ = δd1primal
	Laplace2 = (d1 * hstar2invert * d1.transpose());

	// 2-form L-matrix f*f
	lu.compute(Laplace2);

	//∆β = dw
	twoformscalar = lu.solve(d1*w);

	//δβ
	beta = hstar2invert * d1.transpose() * twoformscalar;



	// ____________________________________________________________________ y

	//EigenMat harmonic;
	EigenMat harmonic;
	harmonic = w - da - beta;


	//___________________________________________________________________________ Vectorfieldinterpolation

	UT_Array<fpreal>curflfree;
	curflfree.setSize(oneformscalar.size());
	UT_Array<fpreal>da2;
	da2.setSize(beta.size());

	UT_Array<fpreal>da3;
	da3.setSize(harmonic.size());
	Eigen::Map<EigenMat>(curflfree.data(), oneformscalar.rows(), oneformscalar.cols()) = oneformscalar;
	Eigen::Map<EigenMat>(da2.data(), beta.rows(), beta.cols()) = beta;
	Eigen::Map<EigenMat>(da3.data(), harmonic.rows(), harmonic.cols()) = harmonic;


	for (int i = 0; i < gdp->getNumPrimitives(); i++) {
		hedgvec1.set(i, oneformtofield(curflfree, i));
		hedgvec2.set(i, oneformtofield(da2, i));
		hedgvec3.set(i, oneformtofield(da3, i));

	}

}
//_______________________________________________________________________________ DecOperators

//_______________________deriative0form 
//edge-node incidence matrix, providing a gradient operator 

//E*P 
DIFF_SOP::SparseMat DIFF_SOP::deriative0(SparseMat d0) {

	std::vector<coeffs> tripletList;
	tripletList.reserve(edgeList.size()*gdp->getNumPoints());

	for (int edge = 0; edge < edgeList.size(); ++edge) {
		GA_Offset minus = gdp->pointOffset(edgeList.data()[edge].x());
		GA_Offset plus = gdp->pointOffset(edgeList.data()[edge].y());
		tripletList.push_back(coeffs(edge, minus, -1));
		tripletList.push_back(coeffs(edge, plus, 1));
	}
	d0.setFromTriplets(tripletList.begin(), tripletList.end());
	return d0;
}
//_______________________deriative1form 
// important condition wich needs to be verified is, for an edge shared by two faces, 
// for each face it should have an opposite orientation

//F*E
DIFF_SOP::SparseMat DIFF_SOP::deriative1(SparseMat d1) {

	fpreal val;
	int edgeindex;
	std::vector<coeffs> tripletList;

	tripletList.reserve(gdp->getNumPrimitives()* edgeList.size());

	for (int prim = 0; prim < gdp->getNumPrimitives(); ++prim) {
		GEO_Hedge hedge = polinterface->polyHedge(prim);
		do {
			hedge = polinterface->nextPrimitiveHedge(hedge);
			int ix = DIFF_SOP::getedgeix(hedge);
			UT_Vector2 a = edgeList.data()[ix];
			GEO_Hedge ahedge = polinterface->findHedgeWithEndpoints(a.y(), a.x());
			UT_Vector2 tmpa = UT_Vector2(polinterface->dstPoint(hedge), polinterface->srcPoint(hedge));
			if (a == tmpa)
				val = 1.0;
			else val = -1.0;
			tripletList.push_back(coeffs(prim, ix, val));

		} while (hedge != polinterface->polyHedge(prim));

	}
	d1.setFromTriplets(tripletList.begin(), tripletList.end());
	return d1;
}
//P*P
DIFF_SOP::SparseMat DIFF_SOP::star1(SparseMat star1) {
	int ptCount = gdp->getNumPoints();

	std::vector<coeffs> tripletList;
	tripletList.reserve(gdp->getNumPoints() * gdp->getNumPoints());

	UT_Array<int> adjacentfaces;
	fpreal	dualarea = 0;

	UT_ValArray< GA_OffsetArray> neighbourArray;
	gdp->buildRingZeroPoints(neighbourArray, NULL);
	for (int pt = 0; pt < ptCount; ++pt) {
		GA_OffsetArray neighbours = GA_OffsetArray(neighbourArray[pt]);

		for (int j = 0; j < neighbours.size(); j++) {

			GEO_Hedge hedge = polinterface->findHedgeWithEndpoints(pt, neighbours[j]);
			GEO_Hedge nt = polinterface->nextEquivalentHedge(hedge);

			UT_Array<int>::iterator it;
			it = find(adjacentfaces.begin(), adjacentfaces.end(), polinterface->hedgePoly(hedge));
			if (it == adjacentfaces.end()) {
				adjacentfaces.append(polinterface->hedgePoly(hedge));
			}
			it = find(adjacentfaces.begin(), adjacentfaces.end(), polinterface->hedgePoly(nt));
			if (it == adjacentfaces.end()) {
				adjacentfaces.append(polinterface->hedgePoly(nt));
			}
		}
		for (int i = 0; i < adjacentfaces.size(); i++) {
			GEO_Primitive *prim = gdp->getGEOPrimitive(adjacentfaces.data()[i]);
			dualarea += prim->calcArea();
		}
		dualarea = (dualarea / 3.0);
		tripletList.push_back(coeffs(pt, pt, dualarea));
	}
	star1.setFromTriplets(tripletList.begin(), tripletList.end());
	return star1;

}
//E*E
DIFF_SOP::SparseMat DIFF_SOP::star2(SparseMat star2) {

	std::vector<coeffs> tripletList;
	tripletList.reserve(edgeList.size()* edgeList.size());

	for (int edge = 0; edge < edgeList.size(); ++edge) {
		GA_Offset pt1 = gdp->pointOffset(edgeList.data()[edge].x());
		GA_Offset pt2 = gdp->pointOffset(edgeList.data()[edge].y());
		GEO_Hedge cotAlphaHe = polinterface->findHedgeWithEndpoints(pt1, pt2);
		GEO_Hedge cotBetaHe = polinterface->nextEquivalentHedge(cotAlphaHe);
		fpreal cotW = (DIFF_SOP::cotangent(cotAlphaHe) + DIFF_SOP::cotangent(cotBetaHe)) / 2.0;
		tripletList.push_back(coeffs(edge, edge, cotW));

	}
	star2.setFromTriplets(tripletList.begin(), tripletList.end());
	return star2;
}

DIFF_SOP::SparseMat DIFF_SOP::star2inv(SparseMat star2) {

	std::vector<coeffs> tripletList;
	tripletList.reserve(edgeList.size()* edgeList.size());

	for (int edge = 0; edge < edgeList.size(); ++edge) {
		GA_Offset pt1 = gdp->pointOffset(edgeList.data()[edge].x());
		GA_Offset pt2 = gdp->pointOffset(edgeList.data()[edge].y());
		GEO_Hedge cotAlphaHe = polinterface->findHedgeWithEndpoints(pt1, pt2);
		GEO_Hedge cotBetaHe = polinterface->nextEquivalentHedge(cotAlphaHe);
		fpreal cotW = (DIFF_SOP::cotangent(cotAlphaHe) + DIFF_SOP::cotangent(cotBetaHe)) / 2.0;
		if (cotW != 0) cotW = 1 / cotW;
		tripletList.push_back(coeffs(edge, edge, cotW));

	}
	star2.setFromTriplets(tripletList.begin(), tripletList.end());
	return star2;
}

//F*F
DIFF_SOP::SparseMat DIFF_SOP::star3(SparseMat star3) {
	std::vector<coeffs> tripletList;
	tripletList.reserve(gdp->getNumPrimitives()* gdp->getNumPrimitives());
	;
	for (int prim = 0; prim < gdp->getNumPrimitives(); ++prim) {
		GEO_Primitive *Geoprim = gdp->getGEOPrimitive(prim);
		fpreal area = 1 / Geoprim->calcArea();
		tripletList.push_back(coeffs(prim, prim, area));
	}
	star3.setFromTriplets(tripletList.begin(), tripletList.end());
	return star3;
}
//_____________________________________________________________________________________________ Helpers


void DIFF_SOP::covector(UT_Array <fpreal>&oneformlist, GA_ROHandleV3 &vectorfield) {

	UT_Vector3 vfield;
	UT_Vector3 vfieldtwin;
	UT_Vector3 Plus;
	fpreal oneform;

	for (int i = 0; i < edgeList.size(); i++) {

		GA_Offset pt1 = gdp->pointOffset(edgeList.data()[i].x());
		GA_Offset pt2 = gdp->pointOffset(edgeList.data()[i].y());
		GEO_Hedge temp = polinterface->findHedgeWithEndpoints(pt1, pt2);
		GEO_Hedge next = polinterface->nextEquivalentHedge(temp);
		UT_Vector3 df = gdp->getPos3(pt2) - gdp->getPos3(pt1);

		if (polinterface->isBoundaryHedge(temp) != 1) {
			vfield = vectorfield.get(polinterface->hedgePoly(temp));
		}
		else vfield = 0;

		if (polinterface->isBoundaryHedge(next) != 1) {
			vfieldtwin = vectorfield.get(polinterface->hedgePoly(next));
		}
		else vfieldtwin = 0;

		Plus = (vfield + vfieldtwin)* .5;
		oneform = dot(Plus, df);
		oneformlist.append(oneform);

	}
}

UT_Vector3 DIFF_SOP::oneformtofield(UT_Array<fpreal> &oneformlist, GA_Offset primoffset) {

	GEO_Primitive *prim = gdp->getGEOPrimitive(primoffset);
	GEO_Hedge Hedge1 = polinterface->polyHedge(primoffset);
	GEO_Hedge hnext = polinterface->nextPrimitiveHedge(Hedge1);
	GEO_Hedge hprev = polinterface->prevPrimitiveHedge(Hedge1);

	UT_Vector3 prev = gdp->getPos3(polinterface->srcPoint(hprev)) - gdp->getPos3(polinterface->dstPoint(hprev));
	UT_Vector3 next = gdp->getPos3(polinterface->srcPoint(hnext)) - gdp->getPos3(polinterface->dstPoint(hnext));
	UT_Vector3 cur = gdp->getPos3(polinterface->srcPoint(Hedge1)) - gdp->getPos3(polinterface->dstPoint(Hedge1));

	UT_Vector3 np = cur - next;
	UT_Vector3 pc = prev - cur;
	UT_Vector3 nc = next - prev;

	UT_Vector3 N = prim->computeNormal();
	N.normalize();
	fpreal area = prim->calcArea();

	int ixhprev = getedgeix(hprev);
	int ixhnext = getedgeix(hnext);
	int hedge1 = getedgeix(Hedge1);

	fpreal edge1 = oneformlist.data()[ixhprev];
	fpreal edge2 = oneformlist.data()[ixhnext];
	fpreal edge3 = oneformlist.data()[hedge1];

	UT_Vector2 a = edgeList.data()[ixhprev];
	UT_Vector2 b = edgeList.data()[ixhnext];
	UT_Vector2 c = edgeList.data()[hedge1];

	UT_Vector2 tmpa = UT_Vector2(polinterface->dstPoint(hprev), polinterface->srcPoint(hprev));
	UT_Vector2 tmpb = UT_Vector2(polinterface->dstPoint(hnext), polinterface->srcPoint(hnext));
	UT_Vector2 tmpc = UT_Vector2(polinterface->dstPoint(Hedge1), polinterface->srcPoint(Hedge1));

	if (a != tmpa) edge1 *= -1;
	if (b != tmpb) edge2 *= -1;
	if (c != tmpc) edge3 *= -1;

	return cross(N, np * edge1 + pc *edge2 + nc * edge3) / ((2.0 * area));

}

//pretty weak TODO 

/*
int DIFF_SOP::getedgeix(GEO_Hedge h, UT_Array<UT_Vector2> &edgeList) {


int primnum;
if (polinterface->isPrimary(h) != 1) {
GEO_Hedge twin = polinterface->sym(h);
primnum = polinterface->hedgePoly(twin);
}
else {

for (auto it = primhedgelist.begin(); it != primhedgelist.end(); it++) {
if (polinterface->hedgePoly(h) == *it);
//std::cout << "y" << (it)-primhedgelist.begin() << std::endl;
return (it)-primhedgelist.begin();
}
};
bool cond = 0;
UT_Array<UT_Vector2>::iterator start;

start = edgeList.begin() + (primnum * 3);
UT_Array<UT_Vector2>::iterator end = start + 3;

UT_Vector2 tmp = UT_Vector2(polinterface->dstPoint(h), polinterface->srcPoint(h));
UT_Vector2 swp = UT_Vector2(polinterface->srcPoint(h), polinterface->dstPoint(h));
//std::cout << tmp << std::endl;
int t = 0;

do {
int primtmp = primnum - t;
start = edgeList.begin() + (primtmp * 3);
for (auto it = start; it != end; it++) {
//	std::cout << *it << "test" << std::endl;
if (tmp == *it) {
cond = 1;
//	std::cout << ((it)-start) + (primtmp * 3) << "edge" << std::endl;
return ((it)-start) + (primtmp * 3);
}
if (swp == *it) {
cond = 1;
//	std::cout << ((it)-start) + (primtmp * 3) << "edge2" << std::endl;
return (it - start) + (primtmp * 3);
}
}
t++;
} while (cond = 1);
}

*/
int DIFF_SOP::getedgeix(GEO_Hedge h) {

	UT_Vector2 tmp = UT_Vector2(polinterface->dstPoint(h), polinterface->srcPoint(h));
	UT_Vector2 swp = UT_Vector2(polinterface->srcPoint(h), polinterface->dstPoint(h));

	if (edgeList.find(tmp) == -1)
		return	edgeList.find(swp);
	else
		return edgeList.find(tmp);


}


fpreal DIFF_SOP::cotangent(GEO_Hedge halfedge) {


	GEO_Hedge hnext = polinterface->nextPrimitiveHedge(halfedge);
	GEO_Hedge hprev = polinterface->prevPrimitiveHedge(halfedge);
	UT_Vector3F nexta = polinterface->hedgeVector(hnext);
	nexta.negate();
	UT_Vector3F prev = polinterface->hedgeVector(hprev);
	return dot(prev, nexta) / cross(nexta, prev).length();


}

//_____________________________________________________________________________________________________________________________Gradient

void DIFF_SOP::gradient(GA_ROHandleF scalarfield, GA_Offset primOffset) {

	GA_Attribute *vec1 = gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "Gradient", 3);
	GA_RWHandleV3 hedgvec1(vec1);
	GEO_Hedge Hedge1 = polinterface->polyHedge(primOffset);
	GA_Offset vtx0 = polinterface->srcPoint(Hedge1);
	GEO_Hedge  hnext = polinterface->nextPrimitiveHedge(Hedge1);
	GA_Offset vtx1 = polinterface->srcPoint(hnext);
	GEO_Hedge  hprev = polinterface->prevPrimitiveHedge(Hedge1);
	GA_Offset vtx2 = polinterface->srcPoint(hprev);
	UT_Vector3F prev = polinterface->hedgeVector(hprev);
	UT_Vector3F next = polinterface->hedgeVector(hnext);
	UT_Vector3F curr = polinterface->hedgeVector(Hedge1);

	UT_Vector3F	   n = cross(next, prev); n.normalize();
	fpreal area = gdp->getGEOPrimitive(primOffset)->calcArea();

	UT_Vector3F  grad =
		cross(n, next) *  scalarfield.get(vtx0, 0) +
		cross(n, prev) *  scalarfield.get(vtx1, 0) +
		cross(n, curr) *  scalarfield.get(vtx2, 0);

	grad = grad / (2.0 * area);
	hedgvec1.set(primOffset, grad);
}
//_____________________________________________________________________________________________________________________________Divergence

void DIFF_SOP::divergence(GA_ROHandleV3 vectorfield, GA_Offset primOffset) {

	GA_Attribute *div = gdp->addFloatTuple(GA_ATTRIB_POINT, "divergence", 1);
	GA_RWHandleF divhandle(div);
	fpreal t = 0;
	GEO_Primitive *prim = gdp->getGEOPrimitiveByIndex(primOffset);
	std::vector<int> primpoints;
	UT_Vector3F grad = vectorfield(primOffset);
	for (int i = 0; i < 3; i++) {
		GA_Offset ptnum = prim->getPointOffset(i);
		primpoints.push_back(ptnum);
	}
	std::vector<int>  a = primpoints;
	std::vector<int>  b = primpoints;
	std::rotate(a.begin(), a.begin() + 1, a.end());
	primpoints.insert(primpoints.end(), a.begin(), a.end());
	std::rotate(b.rbegin(), b.rbegin() + 1, b.rend());
	primpoints.insert(primpoints.end(), b.begin(), b.end());

	for (auto it = primpoints.begin(); it != primpoints.end(); advance(it, 3)) {
		UT_Vector3F curr = gdp->getPos3(*it) - gdp->getPos3(*(it + 1));
		UT_Vector3F next = gdp->getPos3(*(it + 2)) - gdp->getPos3(*(it + 1));
		UT_Vector3F prev = gdp->getPos3(*it) - gdp->getPos3(*(it + 2));
		fpreal cot1 = dot(prev, -next) / cross(prev, -next).length();
		fpreal cot2 = dot(-curr, next) / cross(-curr, next).length();
		t = cot1 * (dot(curr, grad)) + cot2 * (dot(-prev, grad));
		divhandle.add(*it, t);
	}
}

//__________________________________________________________________________________________________________________________COOK
OP_ERROR
DIFF_SOP::cookMySop(OP_Context &context) {
	OP_Node::flags().timeDep = 1;
	int fields = fielchoice();
	OP_Node::flags().timeDep = 1;

	fpreal now = context.getTime();
	OP_AutoLockInputs inputs(this);
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();
	duplicateSource(0, context);
	polinterface = new GEO_PolyInterface(gdp);
	UT_String scalarname;
	SCALAR(scalarname, now);

	if (!scalarname.isstring()) { addMessage(SOP_ATTRIBUTE_INVALID, "field parameter is not defined"); }
	else {
		const GA_Attribute *ah = gdp->findAttribute(GA_ATTRIB_POINT, scalarname);
		const GA_Attribute *dh = gdp->findAttribute(GA_ATTRIB_PRIMITIVE, scalarname);
		if (!ah) {
			addWarning(SOP_ATTRIBUTE_INVALID, "no field values found");
		}
		else {
			if (gdp->findPointAttribute(scalarname)->getTupleSize() < 3) {
				sclar = gdp->findFloatTuple(GA_ATTRIB_POINT, scalarname, 1);
			}
			if (gdp->findPointAttribute(scalarname)->getTupleSize() >= 3) {
				vectorfield = gdp->findFloatTuple(GA_ATTRIB_POINT, scalarname, 3);
			}
		}
		if (dh) {
			vectorfield = gdp->findFloatTuple(GA_ATTRIB_PRIMITIVE, scalarname, 3);
			DIFF_SOP::clearErrors();
		}
	}
	if (sclar.isValid()) {

		if (fields == DIFF_SOP_grad) {
			for (int j = 0; j < gdp->getNumPrimitives(); j++) { gradient(sclar, j); }
		}
		else { gdp->destroyAttribute(GA_ATTRIB_PRIMITIVE, "Gradient", 0); }

	}
	if (vectorfield.isValid()) {

		if (fields == DIFF_SOP_div) {
			for (int j = 0; j < gdp->getPrimitiveRange().getEntries(); j++) { divergence(vectorfield, j); }
		}
		else { gdp->destroyAttribute(GA_ATTRIB_POINT, "Divergece", 0); }

		if (fields == DIFF_SOP_hod) {
			Laplacian(vectorfield, context);
		}
		else { gdp->destroyAttribute(GA_ATTRIB_POINT, "HodgeDecomp", 0); }
	}

	return error();
}







