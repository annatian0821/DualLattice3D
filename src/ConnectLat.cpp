/*
 * ConnectLat.cpp
 *
 *  Created on: Jan 30, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "ConnectLat.hpp"

//using namespace std;

ConnectLat::ConnectLat() {
    // TODO Auto-generated constructor stub
}

ConnectLat::~ConnectLat() {
    // TODO Auto-generated destructor stub
}

void	ConnectLat::putLattice		(	const unsigned 			LatticeID,
										const unsigned			step,
										const double			time,
										const unsigned          type,
                                        Vertices*				p_verti,
                                        ConnectLat*				p_self)
{
    bool isPut = false;
    static std::vector<unsigned> printList;
    std::vector<LatNetwork::iterator> itListToMerge;
    LatNetwork::iterator it_Group;
    LatGroup::iterator it_Element;
    std::list<FLatData> newList;

    FLatData	fLatData (LatticeID,step,time,type,false);

    for (it_Group = latNetwork.begin(); it_Group != latNetwork.end(); ++it_Group) {
        for (it_Element=it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
            if (p_verti->isConnectedExact(LatticeID,it_Element->LatticeID)) {
                if (!isPut) {
                    it_Group->push_back(fLatData);
                    isPut = true;
                }
                itListToMerge.push_back(it_Group);
                break;
            }
        }
    }
    if (!isPut)
    {
        newList.clear();
        newList.push_back(fLatData);
        latNetwork.push_back(newList);
    }
    else
    {
        mergeLatList(&itListToMerge,p_verti);
    }
    rmInnerCrack ();
    if (isLogOffset) {
    	putOffset(LatticeID,p_verti);
    }
    disConList.insert(LatticeID);
}

void	ConnectLat::putMultiConLat			(	std::vector<unsigned>* 		p_latList,
												const unsigned				step,
												const double				time)
{
	LatGroup newList;
	for (	unsigned i = 0; i < p_latList->size(); i++) {
		FLatData newLat ((*p_latList)[i],step,time,false);
		newList.push_back(newLat);
	}
	latNetwork.push_back(newList);
	rmInnerCrack ();
}

void	ConnectLat::rmInnerCrack			()
{
	unsigned rmCnt = 0;
	for (LatNetwork::iterator it_Group=latNetwork.begin(); it_Group!=latNetwork.end(); it_Group++ ) {
		rmCnt += it_Group->size();
		it_Group->remove_if (isBothEndUnstable());
		rmCnt -= it_Group->size();
	}
	latNetwork.remove_if(isEmptyList());

	if (rmCnt!=0) {
		dualOut << rmCnt << " inner crack found and removed" << std::endl;
	}
}

unsigned ConnectLat::reConnectLat	()
{
	static unsigned cntAcc = 0;
	unsigned		cnt = 0;
    GLatForce		glforce;
    LatNetwork::iterator it_Group=latNetwork.begin();
    for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
        if ((glforce.getCrackOpening(it_Element->LatticeID)<minApecture)
        		&&(!isOneEndUnstable(*it_Element))
        		&&(it_Element->isReconnect ==false)) {
       		it_Element->isReconnect = true;
       		reConnectLatSingle(it_Element->LatticeID);
        	reConList.insert(it_Element->LatticeID);
        	disConList.erase(it_Element->LatticeID);
        	cnt++;
        }
    }
    ++it_Group;
    for (; it_Group!=latNetwork.end(); ++it_Group ) {
        for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
           	if ((glforce.getCrackOpening(it_Element->LatticeID)<0.0)
           			&&(!isOneEndUnstable(*it_Element))
           			&&(it_Element->isReconnect ==false)) {
           		it_Element->isReconnect = true;
           		reConnectLatSingle(it_Element->LatticeID);
           		reConList.insert(it_Element->LatticeID);
           		disConList.erase(it_Element->LatticeID);
            	cnt++;
        	}
        }
    }
    cntAcc +=cnt;
    if (cnt!=0) {
    	dualOut << "No of broken lattice re-connected in this step = " << cnt << " Total re-connected lattice = " << reConList.size() << std::endl;
    } else {
    	dualOut << "NO broken lattice re-connected in this step" << std::endl;
    }
    return cnt;
}

void ConnectLat::reConnectLatSingle	(	const unsigned		LatticeID)
{
	unsigned nbNodeID,nbNum;
	LatTab[LatticeID].k[0] = NormLatStiffness*LatTab[LatticeID].area/LatTab[LatticeID].length/10.0;
	for (unsigned i=0; i<2; i++) {
		nbNodeID = LatTab[LatticeID].nb[i];
		if (GNodeTab[nbNodeID].n-GNodeTab[nbNodeID].ne>=1) {
			GNodeTab[nbNodeID].ne++;
		} else {
			dualOut << "[ConnectLat::reConnectLatSingle] - LatticeID = " << LatticeID << " nbNodeID = "
					<< nbNodeID << " n - ne < 1, please check!!" << std::endl;
		}
	}
}
bool ConnectLat::setIsReconnect			(	const unsigned			LatticeID,
    										const bool				inBool)
{
	bool isFound = false;
	LatGroup::iterator it;
	FLatData target(LatticeID,0,0.0,true);
	for (LatNetwork::iterator it_Group=latNetwork.begin(); it_Group!=latNetwork.end(); it_Group++ ) {
		it = find(it_Group->begin(),it_Group->end(),target);
		if (it!=it_Group->end()) {
			if (it->isReconnect!=inBool) {
				it->isReconnect = inBool;
				isFound = true;
				break;
			} else {
				dualOut << "ConnectLat::setIsReconnect - The lattice is found but input of isReconnect is same as "
							"current one which should not be the case, continue to search..." << std::endl;
			}
		}
	}
	return isFound;
}
double		ConnectLat::getTotalFluidVol			()
{
	double tVol = 0.0;
	LatNetwork::iterator it_Group=latNetwork.begin();
	unsigned LatticeID;
	double crackOpening;
	GLatForce lforce;
	for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
		if (isBothEndStable	(*it_Element)) {
			LatticeID = it_Element->LatticeID;
			crackOpening = lforce.getCrackOpening(LatticeID);
			if (crackOpening>0.0) {
				tVol += LatTab[LatticeID].area * crackOpening;
			}
		}
	}
    return tVol;
}

unsigned	ConnectLat::getListSize			()
{
    return latNetwork.size();
}

unsigned	ConnectLat::getLatNum			()
{
    unsigned cnt=0;
    LatNetwork::iterator it_Group;
    for (it_Group=latNetwork.begin(); it_Group!=latNetwork.end(); it_Group++ )
    {
        cnt+=it_Group->size();
    }
    return cnt;
}


std::vector<unsigned>  			ConnectLat::getGroup			(	const unsigned 		groupID)
{
    std::vector<unsigned>				outVec;

    if (groupID<latNetwork.size()) {
    	LatNetwork::iterator it_Group = std::next(latNetwork.begin(),groupID);
    	outVec.reserve(it_Group->size());
        for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
        	outVec.push_back(it_Element->LatticeID);
        }
    } else {
    	dualOut << "The group input does not exist, an empty vector is returned" << std::endl;
    }
    return outVec;
}

std::vector<unsigned>  			ConnectLat::getGroupActive			(	const unsigned 			groupID,
																		const double			minApecture)
{
    std::vector<unsigned>				outVec;
    GLatForce lforce;
    if (groupID<latNetwork.size()) {
    	LatNetwork::iterator it_Group = std::next(latNetwork.begin(),groupID);
    	outVec.reserve(it_Group->size());
        for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
        	if ((lforce.getCrackOpening(it_Element->LatticeID)>minApecture)){
        		outVec.push_back(it_Element->LatticeID);
        	}
        }
    } else {
    	dualOut << "The group input does not exist, an empty vector is returned" << std::endl;
    }
    return outVec;
}

void ConnectLat::sortLatGroup ()
{
    latNetwork.sort(compareLatGroupReverse());
}

std::vector<unsigned>   ConnectLat::getClusterList  ()
{
    sortLatGroup ();
    std::vector<unsigned>   clusterList;
    for (auto it = latNetwork.begin(); it!= latNetwork.end(); ++it) {
        clusterList.push_back(it->size());
    }
    return clusterList;
}

Point  ConnectLat::getMainCrackCenter      ()
{
    std::vector<Point>  pList;

    for (auto it = latNetwork.begin()->begin();it != latNetwork.begin()->end(); ++it) {
        pList.push_back(LatTab[it->LatticeID].centroid);
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

std::vector<std::vector<unsigned> >  ConnectLat::getPerimeter    (    const Point                   center,
                                                                      const double                  lavg,
                                                                      const std::array<unsigned,3>  xyz)
{
    dualOut << "ConnectLat::getPerimeter..." << std::endl;
    Geometry geo;
    std::vector<std::vector<unsigned> >  perimeterList;
    std::vector<std::array<double,2> >     distAngList;
    for (auto it = latNetwork.begin()->begin();it != latNetwork.begin()->end(); ++it) {
        Point centroid = LatTab[it->LatticeID].centroid;
        std::array<double,2> distAng;
        distAng[0]=geo.dist(centroid,center);
        distAng[1] = std::atan2(centroid[xyz[1]] - center[xyz[1]],centroid[xyz[0]] - center[xyz[0]]);
        if (distAng[1] < 0.0) {
            distAng[1] += 2*Pi;
        }
        distAngList.push_back(distAng);
    }
    double N = distAngList.size();
    double adj = 1.5;
    unsigned n = 3;
    unsigned res = adj*4.0*(std::sqrt(N)-1)/n;
    if (res<3) {
        tripOut << "res < 0, quit!" << std::endl;
        return perimeterList;
    }
    std::cout << "N = " << N << " res = " << res << std::endl;
    double dl = 2*Pi/res;
    std::vector<std::vector<std::pair<double,unsigned> > > dirPartition;
    dirPartition.resize(res);
    perimeterList.resize(res);

    unsigned i=0;
    for (auto it = latNetwork.begin()->begin();it != latNetwork.begin()->end(); ++it) {
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
    double R = 0.0;
    unsigned cnt = 0;
    for (unsigned i=0; i<res; i++) {
        std::sort(dirPartition[i].begin(), dirPartition[i].end());
        if (dirPartition[i].size()>=n) {
            for (unsigned j=0; j<n; j++) {
                R += dirPartition[i][dirPartition[i].size()-1-j].first;
                cnt++;
            }
        }
    }
    R/=cnt;
    for (unsigned i=0; i<res; i++) {
        std::sort(dirPartition[i].begin(), dirPartition[i].end());
        if (dirPartition[i].size()>=n) {
            double rmax = 0.0;
            if (dirPartition[i].size()>=n) {
                for (unsigned j=0; j<n; j++) {
                    rmax += dirPartition[i][dirPartition[i].size()-1-j].first;
                }
            }
            rmax /=n;
            unsigned nMax = (rmax/R)*n;
            if (nMax<n) {
                nMax =n;
            }
            for (unsigned j=0; j<nMax; j++) {
                perimeterList[i].push_back(dirPartition[i][dirPartition[i].size()-1-j].second);
            }
        }
    }
    return perimeterList;
}

std::vector<std::vector<unsigned> > 	ConnectLat::getAllGroup		()
{
    std::vector<unsigned>			tmpVec;
    std::vector<std::vector<unsigned> >	outAllVec;
    LatNetwork::iterator it_Group;
    sortLatGroup ();
    for (it_Group=latNetwork.begin(); it_Group!=latNetwork.end(); ++it_Group)
    {
    	tmpVec.clear();
    	for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
    		tmpVec.push_back(it_Element->LatticeID);
    	}
        outAllVec.push_back(tmpVec);
    }

    return outAllVec;
}

std::vector<unsigned> 					ConnectLat::getAllGroupVec		()
{
    std::vector<unsigned>			outVec;
    LatNetwork::iterator it_Group;
    for (it_Group=latNetwork.begin(); it_Group!=latNetwork.end(); ++it_Group) {
    	for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
    		outVec.push_back(it_Element->LatticeID);
    	}
    }
    return outVec;
}

std::vector<unsigned>                   ConnectLat::getAllLat      (    const unsigned      typeID)
{
    std::vector<unsigned>           outVec;
    LatNetwork::iterator it_Group;
    for (it_Group=latNetwork.begin(); it_Group!=latNetwork.end(); ++it_Group) {
        for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
            if (it_Element->type==typeID) {
                outVec.push_back(it_Element->LatticeID);
            }
        }
    }
    return outVec;
}

std::vector<unsigned>                   ConnectLat::getReConLat()
{
    std::vector<unsigned> outLat { std::begin(reConList),std::end(reConList)};
    return outLat;
}

std::vector<unsigned>                   ConnectLat::getDisConLat()
{
    std::vector<unsigned> outLat { std::begin(disConList),std::end(disConList)};
    return outLat;
}


std::vector<std::vector<unsigned> >		ConnectLat::getAllStep			()
{
	std::vector<unsigned>					tmpVec;
	std::vector<std::vector<unsigned> >		outAllVec;
	LatNetwork::iterator it_Group;
	for (it_Group=latNetwork.begin(); it_Group!=latNetwork.end(); ++it_Group)
	{
	    tmpVec.clear();
	    for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
	    	tmpVec.push_back(it_Element->step);
	    }
	    outAllVec.push_back(tmpVec);
	}
	return outAllVec;
}

std::vector<std::vector<double> >		ConnectLat::getAllTime			()
{
	std::vector<double>					tmpVec;
	std::vector<std::vector<double> >		outAllVec;
	LatNetwork::iterator it_Group;
	for (it_Group=latNetwork.begin(); it_Group!=latNetwork.end(); ++it_Group)
	{
	    tmpVec.clear();
	    for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
	    	tmpVec.push_back(it_Element->time);
	    }
	    outAllVec.push_back(tmpVec);
	}
	return outAllVec;
}

std::vector<std::vector<unsigned> >       ConnectLat::getAllType          ()
{
    std::vector<unsigned>                     tmpVec;
    std::vector<std::vector<unsigned> >       outAllVec;
    LatNetwork::iterator it_Group;
    for (it_Group=latNetwork.begin(); it_Group!=latNetwork.end(); ++it_Group) {
        tmpVec.clear();
        for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
            tmpVec.push_back(it_Element->type);
        }
        outAllVec.push_back(tmpVec);
    }
    return outAllVec;
}

void	ConnectLat::mergeLatList	(	std::vector<LatNetwork::iterator>* 	p_itListToMerge,
                                        const Vertices*						p_verti)
{
    if (p_itListToMerge->size() >= 2) {
        for (unsigned i=p_itListToMerge->size()-1; i>0; i--) {
            (*p_itListToMerge)[0]->merge(*(*p_itListToMerge)[i]);
        }
    }
    latNetwork.remove_if(isEmptyList());
}

void	ConnectLat::printAll		(	const unsigned			step)
{
    char fileName[255];

    sprintf(fileName,"%s/%snetworkAll%04d.vtk",OutputFolder,OutputFilePrefix,step);
    std::ofstream cLfile(fileName, std::ios::out);

    unsigned i = 0;

    for (LatNetwork::iterator it_Group = latNetwork.begin(); it_Group != latNetwork.end(); ++it_Group) {
        cLfile << i << " - ";
        for (LatGroup::iterator it_Element=it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
            cLfile << it_Element->LatticeID << " ";
        }
        cLfile << '\n';
        i++;
    }
}

void	ConnectLat::putOffset			(	const unsigned		LatticeID,
											Vertices*			p_verti)
{
	Surface v = p_verti-> getVertiCoordArray (LatticeID);
	Geometry geo;
	Point centroid = geo.getCentroid(LatticeID,v);
	double dz = centroid[2] - centre[2];
	cntStep++;cntTotal++;
	sumDzStep+=dz,sumDzTotal+=dz;
	sumSqDzStep+=dz*dz,sumSqDzTotal+=dz*dz;
	return;
}

void	ConnectLat::logOffset			(	const double		time)
{
	static bool isFirstLog = true;
	char fileName[255];
	std::sprintf(fileName,"%s/%soffsetLog.txt",OutputSubFolder,OutputFilePrefix);
	std::ofstream offsetLog(fileName, std::ios::out | std::ios::app);
	if (isFirstLog) {
		offsetLog 	<< "# time" << '\t' << "offset - Step" << '\t' << "offset - Total " << '\t'
					<< "offset RMS - Step " << '\t' << "offset RMS - Total " << '\n';
		isFirstLog = false;
	}
	if (cntStep>0) {
		offsetLog 	<< time << '\t' << sumDzStep/cntStep << '\t' << sumDzTotal/cntTotal << '\t'
					<< std::sqrt((sumSqDzStep-sumDzStep*sumDzStep/cntStep)/cntStep) << '\t'
					<<  std::sqrt((sumSqDzTotal-sumDzTotal*sumDzTotal/cntTotal)/cntTotal) << '\n';
	}
	sumDzStep = 0.0; sumSqDzStep = 0.0;
	cntStep = 0;
}

void	ConnectLat::logOffset			(	const unsigned		step)
{
	static bool isFirstLog = true;
	char fileName[255];
	std::sprintf(fileName,"%s/%soffsetLog.txt",OutputSubFolder,OutputFilePrefix);
	std::ofstream offsetLog(fileName, std::ios::out | std::ios::app);
	if (isFirstLog) {
		offsetLog 	<< "# step" << '\t' << "offset - Step" << '\t' << "offset - Total " << '\t'
					<< "offset RMS - Step " << '\t' << "offset RMS - Total " << '\n';
		isFirstLog = false;
	}
	if (cntStep>0) {
		offsetLog 	<< step << '\t' << sumDzStep/cntStep << '\t' << sumDzTotal/cntTotal << '\t'
					<< std::sqrt((sumSqDzStep-sumDzStep*sumDzStep/cntStep)/cntStep) << '\t'
					<<  std::sqrt((sumSqDzTotal-sumDzTotal*sumDzTotal/cntTotal)/cntTotal) << '\n';
	}
	sumDzStep = 0.0; sumSqDzStep = 0.0;
	cntStep = 0;
}

void	ConnectLat::logConOffset	(	const double		time,
										Vertices*			p_verti)
{
	static bool isFirstLog = true;
	static double	oldTime = 0.0;
	static unsigned stepCnt = 0;
	static unsigned totalCnt = 0;
	static double	stepSumDz = 0.0;
	static double	stepSumSqDz = 0.0;
	static double	totalSumDz = 0.0;
	static double	totalSumSqDz = 0.0;
	char fileName[255];
	std::sprintf(fileName,"%s/%soffsetConLog.txt",OutputSubFolder,OutputFilePrefix);
	std::ofstream offsetConLog(fileName, std::ios::out | std::ios::app);
	if (isFirstLog) {
		offsetConLog 	<< "# time" << '\t' << "offset - Step" << '\t' << "offset - Total " << '\t'
					<< "offset RMS-Step " << '\t' << "offset RMS-Total " << '\n';
		isFirstLog = false;
	}
	double latFailTime,dz;
	unsigned LatticeID;
	LatNetwork::iterator it_Group = latNetwork.begin();
	for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
		if (isBothEndStable	(*it_Element)) {
			latFailTime = it_Element->time;
			if (latFailTime > oldTime + Tiny) {
				LatticeID = it_Element->LatticeID;
				Surface v = p_verti-> getVertiCoordArray (LatticeID);
				Geometry geo;
				Point centroid = geo.getCentroid(LatticeID,v);
				dz = centroid[2] - centre[2];

				stepCnt++;totalCnt++;
				stepSumDz +=dz; totalSumDz +=dz;
				stepSumSqDz +=dz*dz; totalSumSqDz +=dz*dz;
			}
		}
	}
	if (stepCnt>0) {
		offsetConLog 	<< time << '\t' << stepSumDz/stepCnt << '\t' << totalSumDz/totalCnt << '\t'
						<< std::sqrt((stepSumSqDz-stepSumDz/stepCnt)/stepCnt) << '\t'
						<< std::sqrt((totalSumSqDz-totalSumDz/totalCnt)/totalCnt) << '\n';
	}
	stepSumDz = 0.0; stepSumSqDz = 0.0;
	stepCnt = 0;
	oldTime = time;
}

void	ConnectLat::logConOffset	(	const unsigned		step,
										Vertices*			p_verti)
{
	static bool isFirstLog = true;
	static unsigned	oldStep = 0.0;
	static unsigned stepCnt = 0;
	static unsigned totalCnt = 0;
	static double	stepSumDz = 0.0;
	static double	stepSumSqDz = 0.0;
	static double	totalSumDz = 0.0;
	static double	totalSumSqDz = 0.0;
	char fileName[255];
	std::sprintf(fileName,"%s/%soffsetConLog.txt",OutputSubFolder,OutputFilePrefix);
	std::ofstream offsetConLog(fileName, std::ios::out | std::ios::app);
	if (isFirstLog) {
		offsetConLog 	<< "# Step" << '\t' << "offset - Step" << '\t' << "offset - Total " << '\t'
					<< "offset RMS-Step " << '\t' << "offset RMS-Total " << '\n';
		isFirstLog = false;
	}
	double latFailStep,dz;
	unsigned LatticeID;
	LatNetwork::iterator it_Group = latNetwork.begin();
	for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
		if (isBothEndStable	(*it_Element)) {
			latFailStep = it_Element->step;
			if (latFailStep > oldStep) {
				LatticeID = it_Element->LatticeID;
				Surface v = p_verti-> getVertiCoordArray (LatticeID);
				Geometry geo;
				Point centroid = geo.getCentroid(LatticeID,v);
				dz = centroid[2] - centre[2];

				stepCnt++;totalCnt++;
				stepSumDz +=dz; totalSumDz +=dz;
				stepSumSqDz +=dz*dz; totalSumSqDz +=dz*dz;
			}
		}
	}
	if (stepCnt>0) {
		offsetConLog 	<< step << '\t' << stepSumDz/stepCnt << '\t' << totalSumDz/totalCnt << '\t'
						<< std::sqrt((stepSumSqDz-stepSumDz/stepCnt)/stepCnt) << '\t'
						<< std::sqrt((totalSumSqDz-totalSumDz/totalCnt)/totalCnt) << '\n';
	}
	stepSumDz = 0.0; stepSumSqDz = 0.0;
	stepCnt = 0;
	oldStep = step;
}

void	ConnectLat::logFLat				(	const double		newTime,
											const unsigned 		rmLatAcc)
{
	static unsigned oldGroupSize = getGroup(0).size();
	static double	oldTime = InitTime;
	static bool isFirstLog = true;
	unsigned nowGroupSize = getGroup(0).size();
	char fileName[255];
	std::sprintf(fileName,"%s/%sfLatLog.vtk",OutputSubFolder,OutputFilePrefix);
	std::ofstream fLatLog(fileName, std::ios::out | std::ios::app);
	if (isFirstLog) {
		fLatLog << "# time" << '\t' << "fLatNum" << '\t' << "Diff in conLat "<< '\n';
		isFirstLog = false;
	}
	fLatLog << oldTime << ' ' << rmLatAcc << ' ' << nowGroupSize - oldGroupSize << '\n';
	oldTime = newTime;
	oldGroupSize = nowGroupSize;
	return;
}

void ConnectLat::logFractureVol     (    const double            sumF,
                                         const unsigned          totalStep)
{
    std::cout << "Writing fracture volume..." << '\n';
    static bool isFirstLog = true;
    char fileName[255];
    std::sprintf(fileName,"%s/%smainFracVol.txt",OutputSubFolder,OutputFilePrefix);
    std::ofstream file(fileName, std::ios::out | std::ios::app);
    if (isFirstLog) {
        file   << "# totalStep" << '\t' << "sumF" << '\t' << "Total Vol" << '\n';
        isFirstLog = false;
    }
    std::vector<unsigned>   latList = getGroupActive(0,0.0);
    double totVol = 0.0;
    GLatForce lforce;
    LatNetwork::iterator it_Group = std::next(latNetwork.begin(),0);
    for (LatGroup::iterator it_Element = it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
        unsigned LatticeID = it_Element->LatticeID;
        totVol += std::max(lforce.getCrackOpening(LatticeID),0.0)*LatTab[LatticeID].area;
    }
    file << totalStep << '\t' << sumF << '\t' << totVol << std::endl;

}

unsigned    ConnectLat::setType             (   const unsigned          _LatticeID,
                                            const unsigned          _type)
{
    for (LatNetwork::iterator it_Group = latNetwork.begin(); it_Group != latNetwork.end(); ++it_Group) {
        for (LatGroup::iterator it_Element=it_Group->begin(); it_Element!=it_Group->end(); ++it_Element) {
            if (it_Element->LatticeID == _LatticeID) {
                unsigned oldType = it_Element->type;
                it_Element->type = _type;
                return oldType;
            }
        }
    }
    return 999;
}

void	ConnectLat::swapGroupFrontBack	()
{
	std::swap(latNetwork.front(),latNetwork.back());
}
void	ConnectLat::reset				()
{
	latNetwork.clear();
	reConList.clear();
}
