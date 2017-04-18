/*
 * FluidNetwork.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "FluidNetwork.hpp"

FluidNetwork::FluidNetwork() {
	maxCvDist = 0.0;
	maxCoordNo = 0;
}


FluidNetwork::FluidNetwork	(	const double		_minApecture) {
	maxCvDist = 0.0;
	maxCoordNo = 0;
	minApecture = _minApecture;
}


FluidNetwork::~FluidNetwork() {
}



void FluidNetwork::updateFluidNetwork					(	Vertices*				p_verti,
															ConnectLat*				p_conLat,
															const double            time)
{
	dualOut << "updating FluidNetwork..." << std::endl;
	//To obtain the LatticeID of existing Fluid Node and Group 0 in ConnectLat and sort them. Extract the difference between two list
	//which is the lattice to be added in the fluid network
/*	std::vector<unsigned>	latGroup0 = p_conLat->getGroupActive(0,minApecture);
	std::cout <<  "latGroup0.size() = " << latGroup0.size() << std::endl;
	std::vector<unsigned>	oldfList = getAllLatIDList();
	std::vector<unsigned>	putLatList;
	setDiffVector(latGroup0,oldfList,putLatList);
	unsigned cnt = 0, fNodeID;
	for (unsigned i=0; i<putLatList.size(); i++) {
		if (std::binary_search(latGroup0.begin(),latGroup0.end(),putLatList[i])) {
		    putFNode(putLatList[i],p_verti,p_conLat,time);
		    cnt++;
		}
	}*/
	std::vector<unsigned>   latGroup0 = p_conLat->getGroupActive(0,minApecture);
	latGroup0.insert(latGroup0.end(),initCrackList.begin(),initCrackList.end());
	uniqueVector(&latGroup0);
	std::cout <<  "latGroup0.size() = " << latGroup0.size() << std::endl;
	std::vector<unsigned>   oldfList = getActiveLatIDList();
	std::vector<unsigned>   putLatList;
	setDiffVector(latGroup0,oldfList,putLatList);
	unsigned cnt = 0, fNodeID;
	for (unsigned i=0; i<putLatList.size(); i++) {
	    if (std::binary_search(latGroup0.begin(),latGroup0.end(),putLatList[i])) {
	        if (std::find(inactiveLatID.begin(),inactiveLatID.end(),putLatList[i])==inactiveLatID.end()) {
	            putFNode(putLatList[i],p_verti,p_conLat,time);
	            cnt++;
	        }
	        else {
	            if (setActiveFNode  (putLatList[i])) {
	                fNodeID = std::find(latFNodeSet.begin(),latFNodeSet.end(),putLatList[i])->i2;
	            }
	        }
	    } else {
	        if (setInactiveFNode(putLatList[i])) {
	            fNodeID = std::find(latFNodeSet.begin(),latFNodeSet.end(),putLatList[i])->i2;
	        }
	    }
	}
	connectedFNodeSet();
	dualOut << "No of lattice becomes fluid nodes = " << cnt << " inActive fluid node in this step = " << putLatList.size()-cnt << std::endl;
	dualOut << "Total fluid node = " << fNodeList.size() << " Total pipe no = " << pipeList.size() << std::endl;
}

void FluidNetwork::connectedFNode						()
{
	dualOut << "connectedFNode starts...." << std::endl;
	std::vector<unsigned>	current;
	conFNode.clear();
	unsigned fNodeID;
	for (unsigned i=0; i<priINodeList.size(); i++) {
	    current.push_back(fNodeID = priINodeList[i].fNodeID);
	}
	unsigned nNum;
	std::vector<unsigned>	front;
	std::vector<unsigned>	back;
	std::vector<unsigned>	tmp;
	while (!current.empty()) {
		std::copy(current.begin(),current.end(),std::inserter (conFNode,conFNode.end()));
		for (unsigned i=0; i<current.size(); i++) {
			nNum = fNodeList[current[i]].nbActive.size();
			for (unsigned j=0; j<nNum; j++) {
				front.push_back(fNodeList[current[i]].nbActive[j].fNodeID);
			}
		}
		uniqueVector(&front);
		tmp = current;
		back.insert(back.end(),current.begin(),current.end());
		uniqueVector(&back);
		current = setDiffVector(front,back);
		back = tmp;
		front.clear();
	}
	std::for_each(conFNode.begin(),conFNode.end(),printData<unsigned>);
	dualOut << "connectedFNode ended...." << std::endl;
}

void FluidNetwork::connectedFNodeSet						()
{
	std::cout << "connectedFNodeSet starts...." << std::endl;
	std::set<unsigned>	current;
	conFNode.clear();
	unsigned fNodeID;
	for (unsigned i=0; i<priINodeList.size(); i++) {
		fNodeID = priINodeList[i].fNodeID;
		current.insert(fNodeID);
		if (fNodeList[fNodeID].isActive()) {
			for (unsigned j=0; j<fNodeList[fNodeID].nbActive.size(); j++)
			current.insert(fNodeList[fNodeID].nbActive[j].fNodeID);
		}
	}
	unsigned nNum;
	std::set<unsigned>	front;
	std::set<unsigned>	back;
	while (!current.empty()) {
		for (std::set<unsigned>::iterator it=current.begin(); it!=current.end(); ++it) {
			nNum = fNodeList[*it].nbActive.size();
			for (unsigned j=0; j<nNum; j++) {
				if (fNodeList[*it].isActive()) {
					front.insert(fNodeList[*it].nbActive[j].fNodeID);
				}
			}
		}
		conFNode.insert(current.begin(),current.end());
		current.clear();
		std::set_difference(front.begin(),front.end(),conFNode.begin(),conFNode.end(),std::inserter(current,current.end()));
		front.clear();
	}
	dualOut << "conFNode.size = " << conFNode.size() << std::endl;
	std::cout << "connectedFNodeSet ended...." << std::endl;
}

std::vector<unsigned> FluidNetwork::getPriINodeNbList	()
{
	using Eigen::VectorXd;
	std::vector<unsigned>	nbPriINodeList;
	unsigned fNodeID;
	for (unsigned i=0; i<priINodeList.size(); i++) {
		fNodeID = priINodeList[i].fNodeID;
		dualOut << "fNodeList[fNodeID].nbActive[j].fNodeID...";
		for (unsigned j=0; j< fNodeList[fNodeID].nbActive.size(); j++) {
//			if (std::find(priINodeList.begin(),priINodeList.end(),fNodeList[fNodeID].nbActive[j].fNodeID)==priINodeList.end()) {
		    dualOut << fNodeList[fNodeID].nbActive[j].fNodeID << std::endl;
				nbPriINodeList.push_back(fNodeList[fNodeID].nbActive[j].fNodeID);
//			}
		}
	}
	uniqueVector(&nbPriINodeList);
	return nbPriINodeList;
}

void FluidNetwork::putFNode								(	const unsigned 					LatticeID,
															Vertices*						p_verti,
															const ConnectLat*				p_conLat,
															const double                    t0)
{
	Geometry geo;
	activeFNodeList.insert(fNodeList.size());
	CustomPair<unsigned,unsigned> cp(LatticeID,fNodeList.size());
	latFNodeSet.insert(cp);
	Surface v = p_verti-> getVertiCoordArray (LatticeID);
	FNode newNode (LatticeID,LatTab[LatticeID].centroid,LeakOffCoeff,t0);
	fNodeList.push_back(newNode);
	searchNbFNode (p_verti,v);
}

// It should be follow after fNodeList.push_back
void FluidNetwork::searchNbFNode					(	Vertices*			p_verti,
														const Surface		&v1)
{
	std::vector<unsigned>	nbFNodeIDList;
	unsigned fNodeID = fNodeList.size()-1;
	Geometry geo;
	double cvDist1 = geo.getMaxCentroidVertiDist(fNodeList[fNodeID].centroid,v1);
	double diff,cvDist2;
	Surface	v2;
	//The last FNode is itself
	bool isNbActive;
	for (unsigned nbFNodeID=0; nbFNodeID<fNodeList.size()-1; nbFNodeID++) {
		diff = geo.dist(fNodeList[nbFNodeID].centroid,fNodeList[fNodeID].centroid)-cvDist1;
		if (diff<maxCvDist) {
			v2 = p_verti->getVertiCoordArray	(fNodeList[nbFNodeID].LatticeID);
			cvDist2 = geo.getMaxCentroidVertiDist(fNodeList[nbFNodeID].centroid,v2);
			if (diff<cvDist2) {
				if (p_verti->isConnectedExact(fNodeList[fNodeID].LatticeID,fNodeList[nbFNodeID].LatticeID)) {
					isNbActive = fNodeList[nbFNodeID].isActive();
					putPipe(nbFNodeID,fNodeID,p_verti);
					fNodeList[fNodeID].putNb(nbFNodeID,pipeList.size()-1,isNbActive);
					fNodeList[nbFNodeID].putNb(fNodeID,pipeList.size()-1,true);
					updateMaxCoordNo(nbFNodeID);
				}
			}
		}
	}
	if (cvDist1>maxCvDist) {
		maxCvDist = cvDist1;
	}
	updateMaxCoordNo(fNodeID);
}

void FluidNetwork::checkInactiveNbFNode					()
{
	if (activeFNodeList.size()<=1) {
		return;
	}
	unsigned fNodeID;
	for (std::set<unsigned>::iterator it=activeFNodeList.begin(); it!=activeFNodeList.end(); ++it) {
		fNodeID = *it;
		if (fNodeList[fNodeID].nbActive.empty()) {
			setInactiveFNode(fNodeList[fNodeID].LatticeID);
		}
	}
}

void FluidNetwork::putPipe						(	const unsigned		fNodeID1,
													const unsigned		fNodeID2,
													Vertices*			p_verti)
{
	Geometry geo;
	Line commonEdge = p_verti->getCommonEdge(fNodeList[fNodeID1].LatticeID,fNodeList[fNodeID2].LatticeID);
	double halfLen1 = geo.linePointDist(fNodeList[fNodeID1].centroid,commonEdge);
	double halfLen2 = geo.linePointDist(fNodeList[fNodeID2].centroid,commonEdge);
	double width = geo.dist(commonEdge[0],commonEdge[1]);
	Pipe newPipe(fNodeID1,fNodeID2,halfLen1,halfLen2,width);
	pipeList.push_back(newPipe);
}

bool FluidNetwork::putPriInjectNode		(	const unsigned				LatticeID,
											const double				injectRate,
											const double				pressure,
											Vertices*					p_verti,
											ConnectLat*					p_conLat)
{
	unsigned fNodeID;
	std::vector<FNode>::iterator it = std::find(fNodeList.begin(),fNodeList.end(),LatticeID);
	if (it==fNodeList.end()) {
		fNodeID = fNodeList.size();
		putFNode (LatticeID,p_verti,p_conLat,0.0);
	} else {
		fNodeID = it - fNodeList.begin();
	}
	PriINode  priInjectNode (fNodeID,injectRate,pressure);
	priINodeList.push_back(priInjectNode);
	return true;
}

bool FluidNetwork::putInjectNode		(	const unsigned				LatticeID,
											const double				injectRate,
											Vertices*					p_verti,
											ConnectLat*					p_conLat)
{
	unsigned fNodeID;
	LatfNodeSet::iterator it = std::find(latFNodeSet.begin(),latFNodeSet.end(),LatticeID);
	if (it==latFNodeSet.end()) {
		fNodeID = fNodeList.size();
		putFNode (LatticeID,p_verti,p_conLat,0.0);
	} else {
		fNodeID = it->i2;
	}
	INode newInjectNode (fNodeID,injectRate);
	iNodeList.push_back(newInjectNode);
	return true;
}

void FluidNetwork::setInitialList       (       ConnectLat*             p_conLat)
{
    std::vector<unsigned>    connectLat = p_conLat->getGroup(0);
    dualOut << "FluidCal::setInitialList: set initial list... " << std::endl;
    for (unsigned i=0; i<connectLat.size(); i++) {
        initCrackList.push_back(connectLat[i]);
    }
    std::sort(initCrackList.begin(),initCrackList.end());
}

bool FluidNetwork::removeMultiFNode		(	std::vector<unsigned>			latList)
{
	bool allSuccess = true;
	for (unsigned i=0; i<latList.size(); i++) {
		allSuccess = (allSuccess && setInactiveFNode(latList[i]));
	}
	return allSuccess;
}
bool FluidNetwork::setInactiveFNode		(	unsigned						LatticeID)
{
	if (activeFNodeList.size()<=1) {
		return false;
	}
	std::vector<unsigned>	rmFNode;
	std::vector<FNode>::iterator it = std::find (fNodeList.begin(),fNodeList.end(),LatticeID);
	unsigned fNodeID = it-fNodeList.begin();
	if (it!=fNodeList.end()) {
		if (it -> isActive()) {
			it -> setIsActive(false);
			if (activeFNodeList.erase(fNodeID)==0) {
				tripOut << "fNodeID= " << fNodeID << " cannot be removed from active list!" << std::endl;
				return false;
			}
			inactiveLatID.insert(LatticeID);
			for (unsigned i=0; i<(it->nbActive).size(); i++) {
				fNodeList[it ->nbActive[i].fNodeID].setInactiveNb(fNodeID);
			}
			for (unsigned i=0; i<(it->nbInactive).size(); i++) {
				fNodeList[it ->nbInactive[i].fNodeID].setInactiveNb(fNodeID);
//				}
			}
			fNodeList[fNodeID].old_p = 0.0;
//			}
			return true;
		} else {
			tripOut << "[FluidNetwork::setInactiveFNode] fNodeID = " << fNodeID << " is inactive already so nothing is done" << std::endl;
			return false;
		}
	} else {
		tripOut << "[FluidNetwork::setInactiveFNode] fNodeID = " << fNodeID << " is not a fluid node so no fluid node is removed" << std::endl;
		return false;
	}
}

bool FluidNetwork::reactivateMultiFNode			(	std::vector<unsigned>			latList)
{
	bool allSuccess = true;
	for (unsigned i=0; i<latList.size(); i++) {
		allSuccess = (allSuccess && setActiveFNode(latList[i]));
	}
	return allSuccess;
}

bool FluidNetwork::setActiveFNode		(	unsigned						LatticeID)
{
	std::vector<FNode>::iterator it = std::find (fNodeList.begin(),fNodeList.end(),LatticeID);
	unsigned fNodeID = it-fNodeList.begin();
	if (it!=fNodeList.end()) {
		if (!(it->isActive())) {
			if ((it->nbActive).empty()) {
				return false;
			}
			activeFNodeList.insert(fNodeID);
			if (inactiveLatID.erase(LatticeID)==0) {
				tripOut << "fNodeID = " << fNodeID << " cannot be removed from inactive list!" << std::endl;
			}
			for (unsigned i=0; i<(it->nbActive).size(); i++) {
				fNodeList[it->nbActive[i].fNodeID].setActiveNb(fNodeID);
			}
			for (unsigned i=0; i<(it->nbInactive).size(); i++) {
				fNodeList[it->nbInactive[i].fNodeID].setActiveNb(fNodeID);
			}
			it-> setIsActive(true);
			return true;
		} else {
			tripOut << "[FluidNetwork::setActiveFNode] fNodeID = " << fNodeID
					<< " is active already so nothing is done" << std::endl;
			return false;
		}
	} else {
		tripOut << "[FluidNetwork::setActiveFNode] fNodeID = " << fNodeID
				<< " is not a fluid node so no fluid node is reactivate" << std::endl;
		return false;
	}
}

bool FluidNetwork::checkInActiveFNode		()
{
	GLatForce lforce;
	bool isFound = false;
	for (std::set<unsigned>::iterator it=activeFNodeList.begin(); it!=activeFNodeList.end(); ++it) {
		if (lforce.getCrackOpening(fNodeList[*it].LatticeID)<Tiny) {
			setInactiveFNode(fNodeList[*it].LatticeID);
			isFound = true;
			tripOut << "WARNING!!! fNodeID = " << *it << " is in activeFNodeList but it has -ve crack opening, pls check!!" << std::endl;
		}
	}
	return isFound;
}

std::vector<unsigned>	FluidNetwork::getLatIDList	() {
	std::vector<unsigned> latList;
	for (unsigned i=0; i<fNodeList.size(); i++) {
		latList.push_back(fNodeList[i].LatticeID);
	}
	return latList;
}

std::vector<unsigned>	FluidNetwork::getActiveLatIDList	() {
	std::vector<unsigned> latList;
	for (std::set<unsigned>::iterator it=activeFNodeList.begin(); it!=activeFNodeList.end(); ++it) {
		latList.push_back(fNodeList[*it].LatticeID);
	}
	return latList;
}

std::vector<unsigned>   FluidNetwork::getAllLatIDList    () {
    std::vector<unsigned> latList;
    for (unsigned fNodeID = 0; fNodeID<fNodeList.size(); fNodeID++) {
        latList.push_back(fNodeList[fNodeID].LatticeID);
    }
    return latList;
}

std::vector<Point>	FluidNetwork::getAllFNode			() {
	std::vector<Point> allFNode;
	for (unsigned i=0; i<fNodeList.size(); i++) {
		allFNode.push_back(fNodeList[i].centroid);
	}
	return allFNode;
}



std::vector<IDpair> FluidNetwork::getAllPipe	()
{
	std::vector<IDpair> allPipe;
	IDpair pipeEndID;
	for (unsigned i=0; i<pipeList.size(); i++) {
		pipeEndID.first = pipeList[i].nbFNodeID[0];
		pipeEndID.second = pipeList[i].nbFNodeID[1];
		allPipe.push_back(pipeEndID);
	}
	return allPipe;
}

std::vector<double> FluidNetwork::getAllPipeWidth	()
{
	std::vector<double> 	allWidth;
	for (unsigned i=0; i<pipeList.size(); i++) {
		allWidth.push_back(pipeList[i].width);
	}
	return allWidth;
}

Point  FluidNetwork::getMainCrackCenter      ()
{
    std::vector<Point>  pList;

    for (auto it = fNodeList.begin();it != fNodeList.end(); ++it) {
        pList.push_back(it->centroid);
    }
    Point center;
    center[0] = center[1] = center[2] = 0.0;
    for (unsigned i=0; i<pList.size(); i++) {
        center[0] += pList[i][0];
        center[1] += pList[i][1];
        center[2] += pList[i][2];
    }
    center[0] /=pList.size();
    center[1] /=pList.size();
    center[2] /=pList.size();
    return center;
}

std::vector<unsigned>   FluidNetwork::getPerimeter (    const Point         center,
                                                        const double        lavg)
{
    Geometry geo;
    std::vector<unsigned>   perimeterList;
    std::vector<std::array<double,2> >     distAngList;
    for (auto it = fNodeList.begin();it != fNodeList.end(); ++it) {
        std::array<double,2> distAng;
        distAng[0]=geo.dist(it->centroid,center);
        distAng[1] = std::atan2(it->centroid[1] - center[1],it->centroid[0] - center[0]);
        if (distAng[1] < 0.0) {
            distAng[1] += 2*Pi;
        }
        distAngList.push_back(distAng);
    }
    double N = fNodeList.size();
    double adj = 1.1;
    unsigned n = 3;
    unsigned res = adj*4.0*(std::sqrt(N)-1)/n;
    if (res==0) {
        tripOut << "res = 0, quit!" << std::endl;
        return perimeterList;
    }
    double dl = 2*Pi/res;
    std::vector<std::vector<std::pair<double,unsigned> > > dirPartition;
    dirPartition.resize(res);

    unsigned i=0;
    for (auto it = fNodeList.begin();it != fNodeList.end(); ++it) {
        unsigned aPart = (distAngList[i][1]-Tiny)/dl;
        if (aPart>=dirPartition.size()) {
            tripOut << "aPart>=dirPartition.size(), aPart = " << aPart << std::endl;
        }
        std::pair<double,unsigned>  rNodeID;
        rNodeID.first = distAngList[i][0];
        rNodeID.second = it->LatticeID;
        dirPartition[aPart].push_back(rNodeID);
        i++;
    }
    double R = (sqrt(N)-1.0)/2.0*lavg;
    for (unsigned i=0; i<res; i++) {
        std::sort(dirPartition[i].begin(), dirPartition[i].end());
        double rmax = (dirPartition[i][0].first+dirPartition[i][1].first+dirPartition[i][2].first)/3.0;
        unsigned nMax = std::min((rmax/R)*n,3.0);
        for (unsigned j=0; j<nMax; j++) {
            perimeterList.push_back(dirPartition[i][j].second);
        }
    }
    return perimeterList;
}

bool  FluidNetwork::checkNbActive ()
{
	if (activeFNodeList.size()<=1) {
		return true;
	}
	GLatForce lforce;
	bool isFound = false;
	for (std::set<unsigned>::iterator it=activeFNodeList.begin(); it!=activeFNodeList.end(); ++it) {
		for (unsigned j=0; j<fNodeList[*it].nbActive.size(); j++) {
			if (!fNodeList[fNodeList[*it].nbActive[j].fNodeID].isActive()) {
				tripOut << "WARNING!!! fNodeID = " << *it << " contain inactive nbFNode" << std::endl;
				isFound = true;
			} else if (lforce.getCrackOpening(fNodeList[fNodeList[*it].nbActive[j].fNodeID].LatticeID)<Tiny){
				tripOut << "WARNING!!! fNodeID = " << *it << " contain nbFNode with -ve crack opening" << std::endl;
				isFound = true;
			}
		}
	}
	return isFound;
}

void FluidNetwork::printAllID	(	const unsigned 		step)
{
	char fileName[255];
	sprintf(fileName,"%s/%sFluidID%04d.vtk",OutputSubFolder,OutputFilePrefix,step);
	std::ofstream file(fileName, std::ios::out | std::ios::app);
	file << "activeFNodeList : ";
	for (std::set<unsigned>::iterator it= activeFNodeList.begin(); it!=activeFNodeList.end(); ++it) {
		file << *it << ' ';
	}
	file << '\n';
	file << "inactiveLatID : ";
	for (std::set<unsigned>::iterator it= inactiveLatID.begin(); it!=inactiveLatID.end(); ++it) {
		file << std::find(latFNodeSet.begin(),latFNodeSet.end(),*it)->i2 << ' ';
	}
	file << '\n';
	file << "conFNode : ";
	for (std::set<unsigned>::iterator it= conFNode.begin(); it!=conFNode.end(); ++it) {
		file << *it << ' ';
	}
	file << '\n';
	for (unsigned fNodeID = 0; fNodeID < fNodeList.size(); fNodeID++) {
		file << "fNodeID = " << fNodeID << " active? " << fNodeList[fNodeID].isActive() << " activeNbFNodeID : ";
		for (unsigned i=0; i < fNodeList[fNodeID].nbActive.size(); i++) {
			file << "("<< fNodeList[fNodeID].nbActive[i].fNodeID << "," << fNodeList[fNodeID].nbActive[i].pipeID << ") ,";
		}
		file << " inactiveNbFNodeID : ";
		for (unsigned i=0; i < fNodeList[fNodeID].nbInactive.size(); i++) {
			file << "("<< fNodeList[fNodeID].nbInactive[i].fNodeID << "," << fNodeList[fNodeID].nbInactive[i].pipeID << ") ,";
		}
		file << '\n';
	}
	for (unsigned pipeID = 0; pipeID < pipeList.size(); pipeID++) {
		file << "pipeID = " << pipeID << " nb0 = " << pipeList[pipeID].nbFNodeID[0] << " nb1 = "  << pipeList[pipeID].nbFNodeID[1] << '\n';
	}
}

void FluidNetwork::printConnectivity   (   const unsigned      step)
{
    char fileName[255];
    sprintf(fileName,"%s/%sFluidConnectivity%04d.vtk",OutputSubFolder,OutputFilePrefix,step);
    std::ofstream file(fileName, std::ios::out);
    Geometry geo;
    std::vector<double> distList;
    std::vector<double> offsetList;
    std::vector<double> activeNodeRatioList;
    std::vector<unsigned> coordList;
    std::vector<double> timeList;
    std::vector<double> apectureList;
    Point center;
    center[0] = Nx*UnitLength/2.0;
    center[1] = Ny*UnitLength/2.0;
    center[2] = Nz*UnitLength/2.0;
    const double radius = CrackRadius*Nx*UnitLength;
    Vec     normVec;
    normVec[0] = std::sin(Theta*Pi/180.0);        normVec[1] = 0.0;       normVec[2] = std::cos(Theta*Pi/180.0);
    for (unsigned i=0; i<fNodeList.size(); i++) {
        if (true) {
            std::array<double,2>    distOffset = geo.distOffsetProj(fNodeList[i].centroid,center,normVec);
            distList.push_back(distOffset[0]-radius);
            offsetList.push_back(distOffset[1]);
            coordList.push_back(fNodeList[i].nbActive.size()+fNodeList[i].nbInactive.size());
            activeNodeRatioList.push_back(((double) fNodeList[i].nbActive.size())/(fNodeList[i].nbActive.size()+fNodeList[i].nbInactive.size()));
            timeList.push_back(fNodeList[i].t0);
        }
    }
    file << "# time " << '\t' << "dist" << '\t' << "offset" << '\t' << "CoordNo" << '\t' << "activeNodeRatio" << std::endl;
    for (unsigned i=0; i<timeList.size(); i++) {
        file << timeList[i] << '\t' << distList[i] << '\t' << offsetList[i] << '\t'
             << coordList[i] << '\t' << activeNodeRatioList[i] << std::endl;
    }
}

