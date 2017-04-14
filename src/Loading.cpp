/*
 * Loading.cpp
 *
 *  Created on: Oct 25, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "Loading.hpp"

using namespace std;

void Loading::addLoad				(	const double			maxNorm,
            							const unsigned 			face,
            							const unsigned 			dir)
{
	double t_applyLoad = incrRatio/maxNorm - 1.0;
	dualOut << "maxNorm = " << maxNorm << std::endl;
	deltaF = sumF * t_applyLoad;
	sumF += deltaF;
	applyGUniPressureOnBoundary(deltaF,face,dir);
	dualOut << "Load increase, deltaF = " << deltaF << " sumF = " << sumF << std::endl;

}

void Loading::addLoad               (   const double            maxNorm)
{
    double applyDisp = incrRatio/maxNorm;
    dualOut << "maxNorm = " << maxNorm << std::endl;
    for (unsigned i=0; i<nodeLists.constrain.size(); i++) {
        nodeLists.constrain[i].conDisp[0] *= applyDisp;
    }
    dualOut << "Load increase, Disp = " << nodeLists.constrain[0].conDisp[0] << std::endl;
}

void Loading::addLoad				(	const double			maxNorm,
            							ConnectLat*				p_conLat) {
	cout << "maxNorm = " << maxNorm << " oldMaxNorm = " << oldMaxNorm << " deltaF = " << deltaF << endl;
	double newDeltaF = (incrRatio/maxNorm - 1.0)*sumF;
	deltaF = std::min(std::max(newDeltaF , (incrRatio-1.0)*sumF/2.0),(incrRatio-1.0)*sumF*20);
	oldMaxNorm = maxNorm;
	sumF += deltaF;
	std::vector<unsigned>	connectLat = p_conLat->getGroup(0);
	applyConstFluidPressure (connectLat,sumF);
	dualOut << "Load increase, deltaF = " << deltaF << " sumF = " << sumF << std::endl;

}

void Loading::reduceLoad			(	const double            maxNorm,
                                        const unsigned 			face,
            							const unsigned 			dir)
{
	deltaF = sumF *(std::min(decrRatio/maxNorm-1.0,decrRatio-1.0));
	sumF += deltaF;
	applyGUniPressureOnBoundary(deltaF,face,dir);
	dualOut << "Load decrease, deltaF = " << deltaF << " sumF = " << sumF << endl;
}

void Loading::reduceLoad               (   const double            maxNorm)
{
    double applyDisp = std::min(decrRatio/maxNorm,decrRatio);
    dualOut << "maxNorm = " << maxNorm << std::endl;
    for (unsigned i=0; i<nodeLists.constrain.size(); i++) {
        nodeLists.constrain[i].conDisp[0] *= applyDisp;
    }
    dualOut << "Load increase, Disp = " << nodeLists.constrain[0].conDisp[0] << std::endl;
}

void Loading::reduceLoad			(	const double			maxNorm,
            							ConnectLat*				p_conLat)
{
	deltaF = std::min(sumF *(decrRatio/maxNorm-1.0),0.0);
	sumF += deltaF;
	oldMaxNorm = maxNorm;
	vector<unsigned>	connectLat = p_conLat->getGroup(0);
	applyConstFluidPressure (connectLat,sumF);
	dualOut << "Fluid pressure decrease, deltaF = " << deltaF << " sumF = " << sumF << endl;
}

// Initial load at boundary
void Loading::initLoad				(	const double			_deltaF,
            							const unsigned 			face,
            							const unsigned 			dir)
{
	deltaF = _deltaF;
	applyGUniPressureOnBoundary(deltaF,face,dir);
	sumF += deltaF;
	dualOut << "deltaF = " << deltaF << " sumF = " << sumF << endl;
}

void Loading::initialDisp           (   const unsigned           face,
                                        const unsigned           dir,
                                        const double             disp)
{
    for (unsigned i=0; i<nodeLists.boundary[face].size(); i++) {
        Constrain con(nodeLists.boundary[face][i],dir,disp);
        nodeLists.constrain.push_back(con);
    }
}

// Apply initial fluid pressure
void Loading::initLoad				(	const double			_deltaF,
            							ConnectLat*				p_conLat)
{
	deltaF = _deltaF;
	vector<unsigned>	connectLat = p_conLat->getGroup(0);
	applyConstFluidPressure (connectLat,deltaF);
//	p_fCal-> updateUniformP(deltaF);
	dualOut << "Initial fluid pressure, p = " << deltaF  << std::endl;
}

void Loading::setOldMaxNorm			(	const unsigned			maxBreak)
{
	GLatRemoval	glatRemoval;
	oldMaxNorm = glatRemoval.getMaxLatNormForce (maxBreak);
	dualOut << "Initial p_max is set to be " << oldMaxNorm << std::endl;
}

void Loading::assignPresDisp		(	const unsigned 		NodeID,
                                        const double 		dx,
                                        const double 		dy,
                                        const double 		dz,
                                        const bool* 		flag)
{
    if (flag[0]) {
        GNodeTab[NodeID].d[0] = dx;
    }
    if (flag[1]) {
        GNodeTab[NodeID].d[1] = dy;
    }
    if (flag[2]) {
        GNodeTab[NodeID].d[2] = dz;
    }
}


void Loading::assignNodeExtForce	(	const unsigned 			NodeID,
                                        const double 			fx,
                                        const double 			fy,
                                        const double 			fz)
{
    GNodeTab[NodeID].extdF[0] = fx;
    GNodeTab[NodeID].extdF[1] = fy;
    GNodeTab[NodeID].extdF[2] = fz;

    GNodeTab[NodeID].extF[0] += fx;
    GNodeTab[NodeID].extF[1] += fy;
    GNodeTab[NodeID].extF[2] += fz;
}
void Loading::applyGUniPressureOnBoundary ( const double                pressure,
                                            const unsigned              face,
                                            const unsigned              dir)
{
    dualOut << "Apply uniform pressure on boundaries..."<< std::endl;
    dualOut << " Pressure = " << pressure << " , face = " << face << " , dir = " << dir << std::endl;

    double force;
    double totalForce = 0.0;
    switch (dir) {
    case 0:
        for (unsigned i=0; i<nodeLists.boundary[face].size(); i++) {
            force = pressure*nodeLists.boundaryArea[face][i];
            totalForce += force;
            assignNodeExtForce(nodeLists.boundary[face][i],force, 0.0, 0.0);
        }
        break;
    case 1:
        for (unsigned i=0; i<nodeLists.boundary[face].size(); i++) {
            force = pressure*nodeLists.boundaryArea[face][i];
            totalForce += force;
            assignNodeExtForce(nodeLists.boundary[face][i],0.0, force, 0.0);
        }
        break;
    case 2:
        for (unsigned i=0; i<nodeLists.boundary[face].size(); i++) {
            force = pressure*nodeLists.boundaryArea[face][i];
            totalForce += force;
            assignNodeExtForce(nodeLists.boundary[face][i],0.0, 0.0, force);
        }
        break;
    default:
        dualOut << "Loading::applyUniPressureOnBoundary -- wrong dir input!!" << endl;
    }
    std::cout << "Total force = " << totalForce << std::endl;
}

void Loading::outputForceNodes              (   const unsigned              face)
{
    dualOut << "Output force nodes ..."<< std::endl;
    char fileName[255]; //filename
    sprintf(fileName,"%s/forceNodes.vtk",OutputFolder);
    ofstream file(fileName, ios::out);
    file << "NodeID " << '\t' << "area" << '\n';
    double totalForce = 0.0;
    for (unsigned i=0; i<nodeLists.boundary[face].size(); i++) {
        file << nodeLists.boundary[face][i] << '\t' << nodeLists.boundaryArea[face][i] << '\n';
    }
    std::cout << "Output force nodes finished!" << std::endl;
}

void Loading::setForceZero ()
{
    unsigned k;
    unsigned NodeID = 0;
    unsigned size = GNodeTab.size();

    dualOut << "Set all external force on node to be zero..."<< endl;

    for (NodeID = 0; NodeID < size; NodeID++)
    {
        for (k=0; k<Dim; k++)
        {
            GNodeTab[NodeID].extdF[k] = 0.0;
            GNodeTab[NodeID].extF[k] = 0.0;
        }
    }
}

void Loading::applyConstFluidPressure 		(	std::vector<unsigned> 	latList,
												const double			sumP)
{
	for (unsigned i=0; i<latList.size(); i++) {
		applySingleConstFluidPressure (latList[i],sumP);
	}
	dualOut << "Constant fluid pressure " << sumP << " is applied on fracture surface" << endl;
}

void Loading::applyFluidPressure			(	FluidCal*				p_fCal,
                                                ConnectLat*             p_conLat)
{
	dualOut << "Updating fluid pressure..." << std::endl;
	std::vector<unsigned>   group0 = p_conLat ->getGroup(0);
	for (unsigned i=0; i<group0.size(); i++) {
	    applySingleConstFluidPressure (group0[i],0.0);
	}
	std::vector<IDAndInfo> pairList = p_fCal -> getIDAndPressure();
	std::vector<double>    s  = p_fCal -> getStorageRate   ();

	for (unsigned i=0; i<pairList.size(); i++) {
		if (i<10) {
			dualOut << pairList[i].first << ' ' << pairList[i].second << ' ' << s[i] << std::endl;
		}
		applySingleConstFluidPressure (pairList[i].first,pairList[i].second);
	}
}

void Loading::applySingleConstFluidPressure (	const unsigned		LatticeID,
												const double		deltaP,
												const double		sumP)
{
    unsigned NodeID[2] = { 	LatTab[LatticeID].nb[0],
    						LatTab[LatticeID].nb[1] };
    std::array<double,Dim>	delta;
    double					dist = 0.0;
    for (unsigned i=0; i<Dim; i++) {
    	delta[i] = GNodeTab[NodeID[0]].coord[i] -GNodeTab[NodeID[1]].coord[i];
    	dist += delta[i]*delta[i];
    }
    dist = std::sqrt(dist);
    double deltaF 	= deltaP*LatTab[LatticeID].area;
    double sumF 	= sumP*LatTab[LatticeID].area;

    Geometry geo;
    for (unsigned n=0; n<2; n++) {
    	double sign = 1.0-2.0*(n%2);
    	if (GNodeTab[NodeID[n]].type!=Type_Unstable) {
    	    for (unsigned i=0; i<Dim; i++) {
    	        GNodeTab[NodeID[n]].extdF[i] = sign*deltaF*delta[i]/dist;
    	        GNodeTab[NodeID[n]].extF[i] += sign*deltaF*delta[i]/dist;
    	    }
		}
    }
}

void Loading::applySingleConstFluidPressure (	const unsigned		LatticeID,
												const double		sumP)
{
    unsigned NodeID[2] = { 	LatTab[LatticeID].nb[0],
    						LatTab[LatticeID].nb[1] };
    std::array<double,Dim>	delta;
    double				dist = 0.0;
    for (unsigned i=0; i<Dim; i++) {
    	delta[i] = GNodeTab[NodeID[0]].coord[i] -GNodeTab[NodeID[1]].coord[i] ;
    	dist += delta[i]*delta[i];
    }
    dist = sqrt(dist);
    double old_sumF;
    double sumF 	= sumP*LatTab[LatticeID].area;

    for (unsigned n=0; n<2; n++) {
    	double sign = 1.0-2.0*(n%2);
    	if (GNodeTab[NodeID[n]].type!=Type_Unstable) {
    		for (unsigned i=0; i<Dim; i++) {
    			old_sumF = GNodeTab[NodeID[n]].extF[i];
    			GNodeTab[NodeID[n]].extF[i] = sign*sumF*delta[i]/dist;
    			GNodeTab[NodeID[n]].extdF[i] = GNodeTab[NodeID[n]].extF[i] - old_sumF;
			}
		}
    }
}

void Loading::logSumF	(	const unsigned 		step )
{
	static bool isFirst = true;
	char	fileName[255];
	sprintf(fileName,"%s/%sappliedStressLog.txt.txt",OutputSubFolder,OutputFilePrefix);
	ofstream file(fileName, ios::out | ios::app);
	if (isFirst) {
		file << "# Total Step" << '\t' << "applied stress (kPa)" << '\n';
	}
	file << step << ' ' << sumF << '\n';
	isFirst = false;
}

double Loading::getFacePressure   (     const unsigned      face,
                                        const unsigned      dir)
{
    double sumF = 0.0;
    double sumArea = 0.0;
    for (unsigned i=0; i<nodeLists.boundary[face].size(); i++) {
        sumF += GNodeTab[nodeLists.boundary[face][i]].extF[dir];
        sumArea += nodeLists.boundaryArea[face][i];
    }
    return sumF/sumArea;
}

void Loading::setConForceZero   ()
{
    for (unsigned i=0; i<nodeLists.constrain.size(); i++) {
        unsigned NodeID = nodeLists.constrain[i].NodeID;
        for (unsigned k=0; k<nodeLists.constrain[i].conDir.size(); k++) {
            GNodeTab[NodeID].extF[nodeLists.constrain[i].conDir[k]] = 0.0;
        }
    }
}

void Loading::clearAllLoading	()
{
	for (unsigned i = 0; i<GNodeTab.size(); i++) {
		for (unsigned j=0; j<Dim; j++) {
		    GNodeTab[i].extF[j] = GNodeTab[i].extdF[j] = 0.0;
		}
	}
}
