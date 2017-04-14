/*
 * FluidNetwork.hpp
 *
 *  Created on: Apr 23, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef FLUIDNETWORK_HPP_
#define FLUIDNETWORK_HPP_
#include <iterator>
#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "ConnectLat.hpp"
#include "Vertices.hpp"
#include "ProProcess.hpp"
#include "MiscTools.hpp"
#include <eigen3/Eigen/Sparse>

typedef Eigen::SparseMatrix<double> 	SpMat;
typedef std::set<CustomPair<unsigned,unsigned> > LatfNodeSet;

/*
enum {
	ActiveFNode,
	InactiveFNode,
};
*/
class FluidNetwork {
	friend class Vertcies;
public:
	FluidNetwork();
	FluidNetwork	(	const double		_minApecture);

	~FluidNetwork();

	void updateFluidNetwork					(	Vertices*						p_verti,
												ConnectLat*						p_conLat,
												const double                    time);

	bool putPriInjectNode					(	const unsigned					LatticeID,
												const double					injectRate,
												const double					pressure,
												Vertices*						p_verti,
												ConnectLat*						p_conLat);

	bool putInjectNode						(	const unsigned					LatticeID,
												const double					injectRate,
												Vertices*						p_verti,
												ConnectLat*						p_conLat);

	bool removeMultiFNode					(	std::vector<unsigned>			latList);

	bool reactivateMultiFNode				(	std::vector<unsigned>			latList);

	void setInitialList                     (   ConnectLat*                     p_conLat);

	void printConnectivity                  (   const unsigned                  step);

	std::vector<Point>							getAllFNode	();

	std::vector<IDpair> 						getAllPipe	();

	std::vector<double> 						getAllPipeWidth	();

	std::vector<unsigned> getPerimeter      (   const Point                     center,
	                                            const double                    lavg);

	std::vector<unsigned>   getLatIDList    ();

protected:
	double						maxCvDist;
	unsigned					maxCoordNo;
	double						minApecture;
	LatfNodeSet					latFNodeSet;		//Store LatticeID and corresponding fNode to facilitate searching using LatticeID
	std::set<unsigned>			activeFNodeList;
	std::set<unsigned>			inactiveLatID;		//First is fNodeID, second is LatticeID
    std::vector<unsigned>       initCrackList;
	std::set<unsigned>			conFNode;	//Connected Active FNode
	std::vector<unsigned>		dualFNode;	//Temporary store DOF which is both injection Node and neighbor of Primary injection node
	struct FNode {
		FNode (	unsigned 	_LatticeID,
				Point 		_centroid,
				double      _Cl,
				double      _t0)
		{
			LatticeID = _LatticeID;
			centroid = _centroid;
			old_p	= 0.0;
			t_exp = 0.0;
			leakVol = 0.0;
			Cl = _Cl;
			path	= -1;
			t0 = _t0;
			_isActive = true;
		}
		struct Nb {
			Nb (	unsigned	_fNodeID,
					unsigned	_pipeID) {
				fNodeID = _fNodeID;
				pipeID = _pipeID;
			}
			unsigned fNodeID;
			unsigned pipeID;
			bool operator< 	(	const Nb& rhs) const { return fNodeID < rhs.fNodeID; }
			bool operator== (	const Nb& rhs) const { return fNodeID == rhs.fNodeID; }
			bool operator== (	const unsigned rhs) const { return fNodeID == rhs; }
		};
		unsigned                        LatticeID;
		Point							centroid;
		double							old_p;
		double                          current_p;
		double                          new_p;
		double                          Cl;
		double                          t_exp;
		double                          leakVol;
		double                          t0;
		std::vector<Nb>	 				nbActive;
		std::vector<Nb>					nbInactive;
		int								path;

		void putNb 		(	const unsigned 		nbFNodeID,
							const unsigned 		nbPipeID,
							const bool 			__isActive) {
			Nb nb(nbFNodeID,nbPipeID);
			if (__isActive) {
				nbActive.push_back(nb);
			} else {
				nbInactive.push_back(nb);
			}
		}
		void setIsActive	( bool __isActive) {
			_isActive = __isActive;
		}
		void setInactiveNb (	unsigned fNodeID) {
			std::vector<Nb>::iterator it = std::find(nbActive.begin(),nbActive.end(),fNodeID);
			if (it!=nbActive.end()) {
				nbInactive.push_back(*it);
				remove(nbActive,*it);
			} else {
				tripOut << "WARNING!!! Input fNodeID = " << fNodeID << " cannot be found!!" << std::endl;
			}
		}
		void setActiveNb (	unsigned fNodeID) {
			std::vector<Nb>::iterator it = std::find(nbInactive.begin(),nbInactive.end(),fNodeID);
			if (it!=nbInactive.end()) {
				nbActive.push_back(*it);
				remove(nbInactive,*it);
			} else {
				tripOut << "WARNING!!! Input fNodeID = " << fNodeID << " cannot switch to active!!" << std::endl;
			}
		}
		bool isActive() {
			return _isActive;
		}
		bool operator< 	(	const FNode& rhs) const { return LatticeID < rhs.LatticeID; }
		bool operator== (	const FNode& rhs) const { return LatticeID == rhs.LatticeID; }
		bool operator== (	const unsigned rhs) const { return LatticeID == rhs; }
	private:
		bool							_isActive;
	};

	struct INode {
		INode (	unsigned 	_fNodeID,
				double		injectRate) {
			fNodeID = _fNodeID;
			qin		= injectRate;
		}
		unsigned 	fNodeID;
		double 		qin;
		bool operator< 	(	const INode& rhs) const { return fNodeID < rhs.fNodeID; }
		bool operator== (	const INode& rhs) const { return fNodeID == rhs.fNodeID; }
		bool operator== (	const unsigned rhs) const { return fNodeID == rhs; }
	};

	struct PriINode {
		PriINode () {};
		PriINode (	unsigned 	_fNodeID,
					double		injectRate,
					double		pressure) {
			fNodeID = _fNodeID;
			qin		= injectRate;
			p		= pressure;
		}
		unsigned 	fNodeID;
		double 		qin;
		double		p;
		bool operator< 	(	const PriINode& rhs) const { return fNodeID < rhs.fNodeID; }
		bool operator== (	const PriINode& rhs) const { return fNodeID == rhs.fNodeID; }
		bool operator== (	const unsigned rhs) const { return fNodeID == rhs; }
	};

	struct Pipe {
		Pipe (	unsigned FNodeID1, unsigned FNodeID2, double length1, double length2, double _width)
		{
			nbFNodeID[0] 	= FNodeID1;
			nbFNodeID[1] 	= FNodeID2;
			halfLength[0] 	= length1;
			halfLength[1] 	= length2;
			width			= _width;
			isActive		= true;
		}
		std::array<unsigned,2> 		nbFNodeID;
		std::array<double,2>		halfLength;
		double	 					width;
		bool						isActive;

		void setIsActive	( bool _isActive) {
			isActive = _isActive;
		}
	};

	std::vector<FNode>			fNodeList;
	std::vector<Pipe>			pipeList;
	std::vector<INode>			iNodeList;
	std::vector<PriINode>		priINodeList;
	PriINode					priINode0;
	std::vector<unsigned> getPriINodeNbList	();
	void searchNbFNode						(	Vertices*						p_verti,
												const Surface					&v1);

	void putFNode							(	const unsigned 					LatticeID,
												Vertices*						p_verti,
												const ConnectLat*				p_conLat,
												const double                    time);

	void putPipe							(	const unsigned					fNodeID1,
												const unsigned					fNodeID2,
												Vertices*						p_verti);

	bool setInactiveFNode					(	unsigned						LatticeID);

	bool setActiveFNode						(	unsigned						LatticeID);

	void setPermActiveFNode                 ();

	std::vector<unsigned>   getAllLatIDList ();

	void checkInactiveNbFNode				();

	void connectedFNode						();

	void connectedFNodeSet					();

	std::vector<unsigned>	getActiveLatIDList	();
	void updateMaxCoordNo					(	unsigned						fNodeID) {
		if (fNodeList[fNodeID].nbActive.size()>maxCoordNo) {
			maxCoordNo = fNodeList[fNodeID].nbActive.size();
			dualOut << "Fluid node max coordination no has increased to " << maxCoordNo << std::endl;
		}
	}
	Point  getMainCrackCenter               ();

	bool checkInActiveFNode					();

	bool checkNbActive 						();

	void printAllID							(	const unsigned 				step);
};

#endif /* FLUIDNETWORK_HPP_ */
