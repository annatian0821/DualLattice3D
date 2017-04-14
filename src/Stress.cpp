/*
 * Stress.cpp
 *
 *  Created on: 31 Aug 2015
 *      Author: john
 */

#include "Stress.hpp"

void Stress::calBasic			()
{
//    unsigned 	nbNodeID,LatticeID;
//    double		idx,idy,idz,jdx,jdy,jdz,F,lj;
//	double		stress[6];
    Geometry	geo;

    nNum = tNodeNum + tbNodeNum + tvNodeNum;

    stress.clear();
    stress.resize(nNum);
//    GLatForce l;
    PCG_Eigen	pcg;

#pragma omp parallel for if (Use_OpenMP) shared (GNodeTab)
    for (unsigned NodeID = 0; NodeID < nNum; NodeID++) {
        for (unsigned k=0; k<GNodeTab[NodeID].n; k++) {
            unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[k];
            double idx = GNodeTab[nbNodeID].centroid[0] - GNodeTab[NodeID].centroid[0];
            double idy = GNodeTab[nbNodeID].centroid[1] - GNodeTab[NodeID].centroid[1];
            double idz = GNodeTab[nbNodeID].centroid[2] - GNodeTab[NodeID].centroid[2];
            Tensor memForce = pcg.getMemForce(NodeID,k);

            stress[NodeID][0] += idx * memForce[0];
            stress[NodeID][1] += idy * memForce[1];
            stress[NodeID][2] += idz * memForce[2];
            stress[NodeID][3] += idx * memForce[1];
            stress[NodeID][4] += idx * memForce[2];
            stress[NodeID][5] += idy * memForce[2];
        }
        for (unsigned d=0; d<6; d++) {
            stress[NodeID][d] /= 2.0*GNodeTab[NodeID].v;
        }
    }
}

void Stress::calPriStress			()
{
    calBasic();
    pStress.clear();
    pStress.resize(nNum);
    pDir.resize(nNum);

#pragma omp parallel for if (Use_OpenMP)
    for (unsigned i = 0; i < stress.size() ; i++) {
        //	Stress invarient
        double i1 = stress[i][0] + stress[i][1] + stress[i][2];
        double i2 = stress[i][0]*stress[i][1] + stress[i][1]*stress[i][2] + stress[i][0]*stress[i][2]
                - 	stress[i][3]*stress[i][3] - stress[i][4]*stress[i][4] - stress[i][5]*stress[i][5];
        double i3 = stress[i][0]*stress[i][1]*stress[i][2] + 2.*stress[i][3]*stress[i][4]*stress[i][5]
                -	stress[i][0]*stress[i][5]*stress[i][5] - stress[i][1]*stress[i][4]*stress[i][4] - stress[i][2]*stress[i][3]*stress[i][3];

        double A = 2./3.*sqrt(i1*i1-3.*i2);
        double phi = 1./3.*acos((2.*i1*i1*i1-9.*i1*i2+27*i3)/2/pow(i1*i1-3*i2+tiny,1.5));
        std::vector<iDouble>        tmpStress;
        for (unsigned d=0; d<3; d++) {
            iDouble newitem(i1/3. + A*cos(phi + d*2.*Pi/3.),d);
            tmpStress.push_back(newitem);
        }

        std::sort (tmpStress.begin(),tmpStress.end());

        for (unsigned d=0; d<3; d++) {
            pStress[i][d] = tmpStress[d].data;
        }

        calPDir (i,&tmpStress);
        tmpStress.clear();
    }
}

void Stress::calPDir				(	const unsigned				i,
                                        std::vector<iDouble>*		p_tmpStress)
{
    for (unsigned d=0; d<3; d++) {
        double  a[3];
        a[0] = (stress[i][1]- pStress[i][d])*(stress[i][2]- pStress[i][d])
               - stress[i][5]*stress[i][5];
        a[1] = -(stress[i][3]*(stress[i][2]-pStress[i][d]) - stress[i][4]*stress[i][5]);
        a[2] = stress[i][3]*stress[i][5] - (stress[i][1]-pStress[i][d])*stress[i][4];

        double ABC = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])+tiny;

        for (unsigned k = 0; k<3; k++) {
            pDir[i][3*(*p_tmpStress)[d].index+k] = a[k]/ABC;
        }
    }
}

double Stress::getStress			(	const unsigned NodeID,
                                        const unsigned k)
{
    return stress[NodeID][k];
}

double Stress::getPriStress			(	const unsigned NodeID,
                                        const unsigned k)
{
    return pStress[NodeID][k];
}

double Stress::getStress_p          (   const unsigned NodeID)
{
    return (pStress[NodeID][0]+pStress[NodeID][1]+pStress[NodeID][2])/3.0;
}

double Stress::getStress_q         (    const unsigned NodeID)
{
    double s12 = (pStress[NodeID][0]-pStress[NodeID][1]);
    double s13 = (pStress[NodeID][0]-pStress[NodeID][2]);
    double s23 = (pStress[NodeID][1]-pStress[NodeID][2]);
    return std::sqrt((s12*s12+s13*s13+s23*s23)/2.0);
}

double Stress::getPriDir			(	const unsigned NodeID,
                                        const unsigned d,
                                        const unsigned k)
{
    return pDir[NodeID][3*d+k];
}

Vec Stress::getPriDir               (   const unsigned NodeID,
                                        const unsigned d)
{
    Vec vec;
    vec[0] = pDir[NodeID][3*d];
    vec[1] = pDir[NodeID][3*d+1];
    vec[2] = pDir[NodeID][3*d+2];
    return vec;
}
