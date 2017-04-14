/*
 * ProProcess.cpp
 *
 *  Created on: Oct 25, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "ProProcess.hpp"

using namespace std;

double GLatForce::calculate	()
{
    double max;

    Clock clock_cal;
    clock_cal.start("LatCal");

    max = calAllLatForce(true,true,true);

    clock_cal.get();
    return max;
}

double GLatForce::calAllLatForce	(	const bool				isCal_inLat,
                                        const bool				isCal_bLat,
                                        const bool				isCal_vLat)
{
    double max_e = -Huge;

    if (isCal_inLat) {
        for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
            max_e = std::max(getLatElongation(LatticeID),max_e);
        }
    }

    if (isCal_bLat) {
        for (unsigned LatticeID = tLatticeNum; LatticeID < tLatticeNum+tbLatticeNum; LatticeID++) {
            max_e = std::max(getLatElongation(LatticeID),max_e);
        }
    }

    if (isCal_vLat) {
        for (unsigned LatticeID = tLatticeNum+tbLatticeNum; LatticeID < tLatticeNum+tbLatticeNum+tvLatticeNum; LatticeID++) {
            max_e = std::max(getLatElongation(LatticeID),max_e);
        }
    }

    return max_e;
}

double  GLatForce::getLatElongation    (   const unsigned          LatticeID)
{
    double sum = 0.0;
    unsigned NodeID1 = LatTab[LatticeID].nb[0];
    unsigned NodeID2 = LatTab[LatticeID].nb[1];
    double diff;

    for (unsigned i = 0; i < Dim; i++) {
        diff = GNodeTab[NodeID1].coord[i] - GNodeTab[NodeID2].coord[i]
                 + GNodeTab[NodeID1].d[i] - GNodeTab[NodeID2].d[i];

        sum += diff*diff;
    }
    return sqrt(sum)- LatTab[LatticeID].length;
}

double  GLatForce::getLatTensileStress    (   const unsigned          LatticeID)
{
    double sum = 0.0;
    unsigned NodeID1 = LatTab[LatticeID].nb[0];
    unsigned NodeID2 = LatTab[LatticeID].nb[1];
    double diff;
    for (unsigned i = 0; i < Dim; i++) {
        diff = GNodeTab[NodeID1].coord[i] - GNodeTab[NodeID2].coord[i]
                 + GNodeTab[NodeID1].d[i] - GNodeTab[NodeID2].d[i];

        sum += diff*diff;
    }
    double e =  sqrt(sum)- LatTab[LatticeID].length;
    return LatTab[LatticeID].k[0]*e/LatTab[LatticeID].area;
}


double  GLatForce::getLatExtension    (   const unsigned          LatticeID)
{
    unsigned NodeID1 = LatTab[LatticeID].nb[0];
    unsigned NodeID2 = LatTab[LatticeID].nb[1];

    std::array<double,Dim> delta;
    for (unsigned d = 0; d < Dim; d++) {
        delta[d] = GNodeTab[NodeID1].d[d] - GNodeTab[NodeID2].d[d]
                    + GNodeTab[NodeID1].d0[d] - GNodeTab[NodeID2].d0[d];
    }
    Geometry geo;
    return geo.dot(delta, LatTab[LatticeID].axes[0]);
}

double GLatForce::getCrackOpening	( 	const unsigned			LatticeID)
{
	unsigned NodeID1 = LatTab[LatticeID].nb[0];
	unsigned NodeID2 = LatTab[LatticeID].nb[1];
	double sum = 0.0;
	for (unsigned i = 0; i < Dim; i++) {
		double diff = GNodeTab[NodeID1].coord[i] - GNodeTab[NodeID2].coord[i]
		     + GNodeTab[NodeID1].d[i] - GNodeTab[NodeID2].d[i];
		sum += diff*diff;
	}
	double aperture = sqrt(sum) - LatTab[LatticeID].length - LatTab[LatticeID].initOpening;
	aperture = std::max(aperture,MinApecture);
}

std::array<double,DofPerNode> GLatForce::getResidual			(	const unsigned			NodeID)
{
	std::array<double,DofPerNode> residual;
	Geometry geo;
	for (unsigned d=0; d<DofPerNode; d++) {
		residual[d] = 0.0;
	}
	for (unsigned k=0; k<GNodeTab[NodeID].nbNodeID.size(); k++) {
		unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[k];
		unsigned LatticeID = GNodeTab[NodeID].nbLatticeID[k];
		std::array<double,Dim> delta;
		for (unsigned d=0; d<Dim; d++) {
			delta[d] = GNodeTab[nbNodeID].d[d] - GNodeTab[NodeID].d[d] ;
		}
		double f0 = geo.dot(delta,LatTab[LatticeID].axes[0]);
		for (unsigned d=0; d<Dim; d++) {
			residual[d] += f0*LatTab[LatticeID].axes[0][d]*LatTab[LatticeID].k[0];
		}
		if (KNum>1) {
			double f1 = geo.dot(delta,LatTab[LatticeID].axes[1]);
			double f2 = geo.dot(delta,LatTab[LatticeID].axes[2]);
			for (unsigned d=0; d<Dim; d++) {
				residual[d] += f1*LatTab[LatticeID].axes[1][d]*LatTab[LatticeID].k[1];
				residual[d] += f2*LatTab[LatticeID].axes[2][d]*LatTab[LatticeID].k[1];
			}

		}
		if (KNum >2) {
			std::array<double,Dim> rDelta,rf;
			for (unsigned d=0; d<Dim; d++) {
				rDelta[d] = GNodeTab[nbNodeID].d[d+Dim] - GNodeTab[NodeID].d[d+Dim];
			}
			for (unsigned d=0; d<Dim; d++) {
				rf[d] = geo.dot(rDelta,LatTab[LatticeID].axes[d]);
			}
			for (unsigned d=0; d<Dim; d++) {
				residual[d+Dim] += rf[d]*LatTab[LatticeID].axes[0][d]*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]);
				residual[d+Dim] += rf[d]*LatTab[LatticeID].axes[1][d]*LatTab[LatticeID].k[2];
				residual[d+Dim] += rf[d]*LatTab[LatticeID].axes[2][d]*LatTab[LatticeID].k[3];
			}
		}
	}
	for (unsigned d=0; d<DofPerNode; d++) {
		residual[d] += GNodeTab[NodeID].extF[d];
	}
	return residual;
}

double GLatForce::getMaxLatForce			(	const unsigned			NodeID)
{
	double maxForce = 0.0;
	Geometry geo;
	for (unsigned k=0; k<GNodeTab[NodeID].nbLatticeID.size(); k++) {
		unsigned LatticeID = GNodeTab[NodeID].nbLatticeID[k];
		unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[k];
		std::array<double,Dim> delta;
		for (unsigned d=0; d<Dim; d++) {
			delta[d] = GNodeTab[nbNodeID].d[d] - GNodeTab[NodeID].d[d];
		}
		double resultantForce = 0.0;
		if (KNum==1) {
			resultantForce = geo.dot(delta,LatTab[LatticeID].axes[0])*LatTab[LatticeID].k[0];
		} else {
			double f0 = geo.dot(delta,LatTab[LatticeID].axes[0])*LatTab[LatticeID].k[0];
			double f1 = geo.dot(delta,LatTab[LatticeID].axes[1])*LatTab[LatticeID].k[1];
			double f2 = geo.dot(delta,LatTab[LatticeID].axes[2])*LatTab[LatticeID].k[1];
			resultantForce = std::sqrt (f0*f0+f1*f1+f2*f2);
		}
		if (maxForce < resultantForce) {
			maxForce = resultantForce;
		}
	}
	if (maxForce<Tiny) {
		maxForce = 1.0;
	}
	return maxForce;
}

double GLatForce::getMaxLatMoment			(	const unsigned			NodeID)
{
	double maxMoment = 0.0;
	Geometry geo;
	for (unsigned k=0; k<GNodeTab[NodeID].nbLatticeID.size(); k++) {
		unsigned LatticeID = GNodeTab[NodeID].nbLatticeID[k];
		unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[k];
		std::array<double,Dim> delta;
		for (unsigned d=0; d<Dim; d++) {
			delta[d] = GNodeTab[nbNodeID].d[d+Dim] - GNodeTab[NodeID].d[d+Dim];
		}
		double resultantMoment = 0.0;
		double f0 = geo.dot(delta,LatTab[LatticeID].axes[0])*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]);
		double f1 = geo.dot(delta,LatTab[LatticeID].axes[1])*LatTab[LatticeID].k[2];
		double f2 = geo.dot(delta,LatTab[LatticeID].axes[2])*LatTab[LatticeID].k[3];
		resultantMoment = std::sqrt (f0*f0+f1*f1+f2*f2);
		if (maxMoment < resultantMoment) {
			maxMoment = resultantMoment;
		}
	}
	if (maxMoment<Tiny) {
		maxMoment = 1.0;
	}
	return maxMoment;
}



