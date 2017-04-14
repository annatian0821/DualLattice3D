/*
 * ProProcess.hpp
 *
 *  Created on: Oct 25, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef PROPROCESS_HPP_
#define PROPROCESS_HPP_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "MiscTools.hpp"
//#include "PCG_Eigen.hpp"
#include <eigen3/Eigen/Dense>

class GLatForce
{
	friend class PCG_Eigen;
public:
    GLatForce() {};

    ~GLatForce() {};
    double 	calculate			();

    double  getLatElongation    (   const unsigned          LatticeID);

    double  getLatExtension     (   const unsigned          LatticeID);

    double 	getMaxLatNormForce (	unsigned 				maxBreak);

    double 	calAllLatForce		(	const bool				isCal_inLat,
                                    const bool				isCal_bLat,
                                    const bool				isCal_vLat);

    void 	calNormEnergy		();

    double	getCrackOpening		( 	const unsigned			LatticeID);

    double  getLatTensileStress (   const unsigned          LatticeID);

    double 	getEnergyRelease	( 	const unsigned			LatticeID)
    {
    	return 0.5*LatTab[LatticeID].k[0]*LatTab[LatticeID].et[0]*LatTab[LatticeID].et[0];
    }

    std::array<double,DofPerNode> 	getResidual			(	const unsigned			NodeID);

    double 						getMaxLatForce			(	const unsigned			NodeID);

    double 						getMaxLatMoment			(	const unsigned			NodeID);
private:
    void 	calSingleLatForce	( 	const unsigned			LatticeID,
                                    double*					max);
    Eigen::VectorXd getMemForce	(	const unsigned			NodeID,
    								const unsigned			k);
};





class RmLattice
{
public:

private:

};



#endif /* PROPROCESS_HPP_ */
