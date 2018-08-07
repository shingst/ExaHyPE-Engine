// Implementing all of the slicers in a single file to decrease
// compilation time, for the moment.

#include "Slicer.h"
#include "RegionSlicer.h"
#include "CartesianSlicer.h"

#include "tarch/logging/Log.h"
#include "exahype/parser/Parser.h"

#include <sstream>

using namespace exahype::plotters;
using namespace std;

// storage for the RegionSlicer limits
const double RegionSlicer::defaultLeftBottomFront = -std::numeric_limits<double>::max();
const double RegionSlicer::defaultRightTopBack    = +std::numeric_limits<double>::max();

static tarch::logging::Log _log("exahype::plotters::Slicer");

CartesianSlicer::CartesianSlicer(const dvec& _req, const ivec& _active, int _baseDim) : 
		targetDim(_baseDim - tarch::la::sum(_active)),
		baseDim(_baseDim),
		req(_req),
		active(_active),
		activeAxes(-1),
		runningAxes(-1) {
		
	for(int i=0; i<DIMENSIONS; i++) {
		activeAxes(i) = disabled;
		runningAxes(i) = disabled;
		
		// This algorithm is crazy. Needed a lot of debugging with standalone
		// examples, but now its tested for DIM<=3.
		
		for(int j=i; j<DIMENSIONS; j++) { // forward check for actives
			if(active(j)) {
				activeAxes(i)=j;
				for(int k=0; k<i; k++) { // backward check if not already included
					if(activeAxes(k)==j)
						activeAxes(i)=disabled;
				}
				if(activeAxes(i)!=disabled)
					break;
			}
		}
		
		for(int j=i; j<DIMENSIONS; j++) { // forward check for actives
			if(!active(j)) {
				runningAxes(i)=j;
				for(int k=0; k<i; k++) { // backward check if not already included
					if(runningAxes(k)==j)
						runningAxes(i)=disabled;
				}
				if(runningAxes(i)!=disabled)
					break;
			}
		}
	}
}
	
CartesianSlicer* CartesianSlicer::fromSelectionQuery(const exahype::parser::ParserView select) {
	dvec r; ivec v;
	v(0) = select.isValueValidDouble("x");
	r(0) = select.getValueAsDouble("x");
	v(1) = select.isValueValidDouble("y");
	r(1) = select.getValueAsDouble("y");
	#if DIMENSIONS==3
	v(2) = select.isValueValidDouble("z");
	r(2) = select.getValueAsDouble("z");
	#endif

	// v(i) == true == 1 means that value was provided, v(i) == false == 0 means that not.
	return new CartesianSlicer(r, v);
}

/**
 * A variant of tarch::la::Vector::toString which replcaes infVal with repl.
 **/
std::string valueReplPrinter(const tarch::la::Vector<DIMENSIONS, double>& vec, const double infVal, const std::string& repl) {
	stringstream s;
	s << "[";
	for(int i=0; i < DIMENSIONS; i++) {
		if(vec[i] == infVal) s << repl;
		else s << vec[i];
		if(i + 1 < DIMENSIONS) s << ",";
	}
	s << "]";
	return s.str();
}

std::string RegionSlicer::toString() const {
	// a beautiful infinity character in utf8. Could also just use "inf".
	// Unicode works in most Linux and Mac terminals but not Windows. (https://stackoverflow.com/a/12020179)
	std::string inf = "\u221E", plus = "+", minus = "-";
	
	stringstream s;
	s << "RegionSlicer("
	  << valueReplPrinter(_regionOfInterestLeftBottomFront, defaultLeftBottomFront, minus + inf)
	  << ","
	  << valueReplPrinter(_regionOfInterestRightTopBack, defaultRightTopBack, plus + inf)
	  << ")";
	
	return s.str();
}

std::string CartesianSlicer::toString() const {
	stringstream s;
	s << "CartesianSlicer(Dim["<<baseDim<<" -> "<<targetDim<<"], req="<<req<<"="<<planeLabel()<<")";
	return s.str();
}

// debugging stuff, should not be operator<< but be named like "debugString" or so.
//std::ostream& operator<<(std::ostream &s,const CartesianSlicer& c) {
std::string CartesianSlicer::debugVerbose() {
	stringstream s;
	s << "CartesianSlicer, Reducing Dimension " << baseDim << " to " << targetDim << ":\n";
	s << "   req = " << req << "\n";
	s << "   active = " << active << "\n";
	s << "   activeAxes = " << activeAxes << "\n";
	s << "   runningAxes = " << runningAxes << "\n";
	return s.str();
}

std::string CartesianSlicer::planeLabel() const {
	// 1D cutting:
	bool active_2 = DIMENSIONS == 3 ? active(2) : true;
	if(!active(0) &&  active(1) &&  active_2) return "x";
	if( active(0) && !active(1) &&  active_2) return "y";
	if( active(0) &&  active(1) && !active_2) return "z";
	
	// 2D cutting, returns "xy", "xz" or "yz"
	if(!active(0) && !active(1) &&  active_2) return "xy";
	if( active(0) && !active(1) && !active_2) return "yz";
	if(!active(0) &&  active(1) &&  active_2) return "xz";
	return "unknowns";
}

RegionSlicer* RegionSlicer::fromSelectionQuery(const exahype::parser::ParserView select) {
	dvec regionOfInterestLeftBottomFront, regionOfInterestRightTopBack;
	double x;
	
	x = select.getValueAsDouble("left");
	regionOfInterestLeftBottomFront(0) = select.getValueAsDouble("left") ? x : defaultLeftBottomFront; // "-", min
	x = select.getValueAsDouble( "bottom" );
	regionOfInterestLeftBottomFront(1) = select.getValueAsDouble("bottom") ? x : defaultLeftBottomFront; // "-", min
	#if DIMENSIONS==3
	x = select.getValueAsDouble( "front" );
	regionOfInterestLeftBottomFront(2) = select.getValueAsDouble("front") ? x : defaultLeftBottomFront; // "-", min
	#endif
	
	x = select.getValueAsDouble( "right" );
	regionOfInterestRightTopBack(0) = select.getValueAsDouble("right") ? x : defaultRightTopBack;
	x = select.getValueAsDouble( "top" );
	regionOfInterestRightTopBack(1) = select.getValueAsDouble("top") ? x : defaultRightTopBack;
	#if DIMENSIONS==3
	x = select.getValueAsDouble( "back" );
	regionOfInterestRightTopBack(2) = select.getValueAsDouble("back") ? x : defaultRightTopBack;
	#endif
	
	return new RegionSlicer(regionOfInterestLeftBottomFront, regionOfInterestRightTopBack);
}

Slicer* Slicer::bestFromSelectionQuery(const exahype::parser::ParserView select) {
	logInfo("bestFromSelectionQuery", "Scanning plotting selection query '"<<select.dump()<<"'");
	
	// Build up the registry:
	Slicer *a = CartesianSlicer::fromSelectionQuery(select);
	Slicer *b = RegionSlicer::fromSelectionQuery(select);

	if(a->clips() && b->clips()) {
		logError("bestFromSelectionQuery", "Warning: Several slicing strategies apply to the given arguments '"<<select.dump()<<"'. I choose " << a->getIdentifier());
	}

	if(a->clips()) { delete b; return a;}
	if(b->clips()) { delete a; return b; }
	
	// nothing clips
	delete a; delete b;
	return new NonSlicer;
}

	
