/*
 * Statistics.cpp
 *
 *  Created on: Nov 12, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "Statistics.hpp"

using namespace std;


Statistics::Statistics() {
    // TODO Auto-generated constructor stub

}

Statistics::~Statistics() {
    // TODO Auto-generated destructor stub
}

void    EnergyCal::fillForceNodeList       (   const unsigned  faceID,
                                               const unsigned  dir)
{
    for (unsigned i=0; i<nodeLists.boundary[faceID].size(); i++) {
        forceNodeIDList.push_back(nodeLists.boundary[faceID][i]);
        oldDisp.push_back(0.0);
    }
    loadDir = dir;
}

void    EnergyCal::record      (    FluidCal*                   p_fCal)
{
    oldApecture.clear();
    fLatID = p_fCal->getLatID();
    oldPressure = p_fCal->getPressure(0.0);
    GLatForce   glf;
    Geometry geo;
    for (unsigned i=0; i<fLatID.size(); i++) {
        oldApecture.push_back(glf.getCrackOpening(fLatID[i]));
    }
    oldTotalStrainEnergy = 0.0;
    for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
    	oldTotalStrainEnergy += getTotalStrainEnergy(LatticeID);
    }
}

void    EnergyCal::record      ()
{
    for (unsigned i=0; i<forceNodeIDList.size(); i++) {
        oldDisp[i] = GNodeTab[forceNodeIDList[i]].d[loadDir];
    }
    oldTotalStrainEnergy = 0.0;
    for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
        oldTotalStrainEnergy += getTotalStrainEnergy(LatticeID);
    }
}

void    EnergyCal::compute     (    double                      time)
{
    static bool isFirst = true;
    char    fileName[255];
    std::sprintf(fileName,"%s/energyLog.txt",OutputSubFolder);
    std::ofstream file(fileName, std::ios::out | std::ios::app);
    if (isFirst) {
        file << "#time" << '\t' << "totalFracArea" << '\t' <<  "ERR" << '\t' << "ERR2" << '\t' << "accFracArea" << '\t'
                << "accERR" << std::endl;
        isFirst = false;
    }
    if (totalFracArea<Tiny) {
        dualOut << "[EnergyCal::compute] : There is no fracturing, no statistics for ERR is done!" << std::endl;
        return;
    }
    double ERR = ER / totalFracArea;
    double ERR2 = totalBrokenLatEnergy/totalFracArea;
    file << '\t' << totalFracArea << '\t' << ERR << '\t' << ERR2 << '\t' << ERR/ERR2 << std::endl;
    file << time << '\t' << totalFracArea << '\t' << ERR << '\t' << ERR2 << '\t'  << accFracArea << '\t' << accER/accFracArea << std::endl;
    accER += ER;
    ER = 0.0;
    accFracArea += totalFracArea;
    totalFracArea = 0.0;
    totalBrokenLatEnergy = 0.0;
}

void    EnergyCal::accumulate     ()
{
    double  work=0.0;
    GLatForce   glf;
    if (ProblemID==1) {
        for (unsigned i=0; i<fLatID.size(); i++) {
            work += oldPressure[i]*LatTab[fLatID[i]].area*(glf.getCrackOpening(fLatID[i])-oldApecture[i]);
        }
    } else {
        for (unsigned i=0; i<forceNodeIDList.size(); i++) {
            work += GNodeTab[forceNodeIDList[i]].extF[loadDir]*(GNodeTab[forceNodeIDList[i]].d[loadDir]-oldDisp[i]);
        }
    }

    double newTotalStrainEnergy = 0.0;
    for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
        newTotalStrainEnergy += getTotalStrainEnergy(LatticeID);
    }
    double diffSE = oldTotalStrainEnergy - newTotalStrainEnergy;
    ER += work + diffSE;
}

void    EnergyCal::compute     (    unsigned                        step,
                                    double                          sumF)
{
    static bool isFirst = true;
    char    fileName[255];
    std::sprintf(fileName,"%s/energyLog.txt",OutputSubFolder);
    std::ofstream file(fileName, std::ios::out | std::ios::app);
    if (isFirst) {
        file << "#Step" << '\t' << "sumF" << '\t' << "totalFracArea" << '\t' <<  "ERR" << '\t' << "ERR2" << '\t' << "accFracArea" << '\t'
                << "accERR" << std::endl;
        isFirst = false;
    }
    if (totalFracArea<Tiny) {
        dualOut << "[EnergyCal::compute] : There is no fracturing, no statistics for ERR is done!" << std::endl;
        return;
    }
    double ERR = ER / totalFracArea;
    double ERR2 = totalBrokenLatEnergy/totalFracArea;
    file << step << '\t' << sumF << totalFracArea << '\t' << ERR << '\t' << ERR2 << '\t' << accFracArea << '\t' << accER/accFracArea << std::endl;
    accER += ER;
    ER = 0.0;
    accFracArea += totalFracArea;
    totalFracArea = 0.0;
    totalBrokenLatEnergy = 0.0;
}

double EnergyCal::getTotalStrainEnergy	(	const unsigned		LatticeID)
{
	Eigen::MatrixXd     lSM(2*DofPerNode, 2*DofPerNode);
	PCG_Eigen	pcg;
	if (KNum==1) {
		Eigen::MatrixXd     lSM2 = pcg.getSpringLocalSM(LatticeID);
		lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
		lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
		lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -lSM2;
		lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = lSM2;
	} else if (KNum==2) {
		Eigen::MatrixXd     lSM2 = pcg.getShearLocalSM(LatticeID);
		lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
		lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
		lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -lSM2;
		lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = lSM2;
	} else {
		lSM = pcg.getFullLocalSM2(LatticeID);
	}
	Eigen::VectorXd		vd(2*DofPerNode);
	unsigned NodeID = LatTab[LatticeID].nb[0];
	unsigned nbNodeID = LatTab[LatticeID].nb[1];
	for (unsigned d=0; d<DofPerNode; d++) {
		vd[d] = GNodeTab[NodeID].d[d];
		vd[d+DofPerNode] = GNodeTab[nbNodeID].d[d];
	}
	return 0.5*vd.transpose()*lSM*vd;
}


void EnergyCal::setAreaZero ()
{
    totalFracArea =0;
}
void EnergyCal::accArea     (   double      fracArea)
{
    totalFracArea +=fracArea;
}
void EnergyCal::accBrokenLatEnergy  (   double              energy)
{
    totalBrokenLatEnergy += energy;
}
void EnergyCal::setBrokenLatEnergyZero  ()
{
    totalBrokenLatEnergy = 0.0;
}

void Statistics::Fluid_pressure     (   FluidCal*               p_fCal,
                                        unsigned                step)
{
    dualOut << "Writing fluid pressure distribution..." << '\n';

    char    fileName[255];
    sprintf(fileName,"%s/%sstat-fluidP%04d.txt",path,prefix,step);
    std::ofstream file(fileName, std::ios::out);

    std::vector<Point>      fNode = p_fCal->getAllFNode();
    std::vector<double>     fPressure = p_fCal->getPressure(0.0);
    std::vector<double>     fStorageRate = p_fCal-> getStorageRate   ();
    std::vector<unsigned>     fLatticeID = p_fCal-> getLatIDList();
    Geometry geo;
    unsigned tDispfNode = fNode.size();
    Point center = fNode[0];
    GLatForce lforce;
    double total_storage = 0.0;
    double total_storageRate = 0.0;
    double time = p_fCal->getTime();
    double dt = p_fCal-> getTimeStep();
    file << "R" << '\t' << "p" << '\t' << "crack opening" << '\t' << "q" << '\t' << "Area" << std::endl;
    for (unsigned i=0; i < tDispfNode; i++)
    {
        file << geo.dist(center,fNode[i]) << '\t' << std::max(fPressure[i],0.0) << '\t'
                << lforce.getCrackOpening(fLatticeID[i]) << '\t'
                << fStorageRate[i] << '\t'
                << LatTab[fLatticeID[i]].area << '\n';
        if (fPressure[i]>0.0) {
            double tStorage = lforce.getCrackOpening(fLatticeID[i])*LatTab[fLatticeID[i]].area;
            file << tStorage<< std::endl;
            total_storage += tStorage;
        } else {
            file << 0.0 << std::endl;
        }
        total_storageRate += fStorageRate[i];
    }
    file << "Time = " << '\t' << time << std::endl;
    file << "total_storage/time = " << '\t' << total_storage/(time) << '\n';
}

void Statistics::compute	(	const unsigned 			resolution)
{
    Clock	stClock;
    unitAng = Pi/((double)resolution);
    unitAngInDeg = unitAng*180/Pi;

    stClock.start("Stat");
    dualOut << "Computing statistical data..." << endl;
    fill_eLatIDListSimple();

    glatDir("xy",resolution);
    glatDir("xz",resolution);
    glatDir("yz",resolution);

    stClock.get();
}

void Statistics::fillRefineLatList ()
{
    dualOut << "Starting Statistics::fillRefineLatList..."<< std::endl;
    if (refineNode[0]>=refineNode[1]) {
        dualOut << "refineNode[0]>=refineNode[1], nothing is done!" << std::endl;
        return;
    }
    for (unsigned NodeID = refineNode[0]; NodeID < refineNode[1]; NodeID++) {
        for (unsigned k=0; k<GNodeTab[NodeID].nbLatticeID.size(); k++) {
            unsigned LatticeID = GNodeTab[NodeID].nbLatticeID[k];
            unsigned nb0 = LatTab[LatticeID].nb[0];
            unsigned nb1 = LatTab[LatticeID].nb[1];
            if ((nb0 >=refineNode[0])&&(nb1 <=refineNode[1])) {
                refineLatList.push_back(LatticeID);
            }
        }
    }
    std::sort(refineLatList.begin(), refineLatList.end());
    uniqueVector(&refineLatList);
    dualOut << "Total refineLat No = " << refineLatList.size() << std::endl;
}

std::vector<Statistics::RLat> Statistics::getRLatList ( std::vector<std::vector<unsigned> >      &periList,
                                                        const Point                              center)
{
    dualOut << "Starting Statistics::getRLatList..."<< std::endl;
    std::vector<RLat>       rLatList = getRLatList (periList,&innerLatList,center);
    std::vector<unsigned>   latList;
    for (unsigned i=0; i<rLatList.size(); i++) {
        latList.push_back(rLatList[i].LatticeID);
    }
    std::sort(latList.begin(),latList.end());
    std::vector<unsigned>   removeList;
    std::set_symmetric_difference (latList.begin(),latList.end(),innerLatList.begin(),innerLatList.end(),std::back_inserter(removeList));
    for (unsigned i=0; i<removeList.size(); i++) {
        eraseValue(innerLatList,removeList[i]);
    }
    std::sort (rLatList.begin(),rLatList.end());
    dualOut << "Total rLat = " << rLatList.size() << std::endl;
    return rLatList;
}

std::vector<Statistics::RLat> Statistics::getRNodeList ( std::vector<std::vector<unsigned> >      &periList,
                                                         std::vector<unsigned>                    &nodeList,
                                                         const Point                              center)
{
    dualOut << "Starting Statistics::getRNodeList..."<< std::endl;

    std::vector<RLat>       rNodeList;
    if (nodeList.size()<500) {
        tripOut << "The nodeList given is too small, nodeList.size() = " << nodeList.size() << std::endl;
        return rNodeList;
    }
    unsigned res = periList.size();
    double dl = 2*Pi/res;
    Geometry geo;
    for (unsigned i=0; i<nodeList.size(); i++) {
        unsigned NodeID = nodeList[i];
        Point centroid = GNodeTab[NodeID].coord;
        double angle = std::atan2(centroid[xyz[1]]-center[xyz[1]],centroid[xyz[0]]-center[xyz[0]]);
        if (angle<0.0) {
            angle += 2*Pi;
        }
        unsigned aPart = (angle-Tiny)/dl;
        if (aPart>=res) {
            aPart = res -1;
        }
        double R = 0.0;
        for (unsigned p=0; p<periList[aPart].size(); p++) {
            R += geo.dist(LatTab[periList[aPart][p]].centroid,center);
        }
        R /= periList[aPart].size();
        std::array<double,2> diffC;
        diffC[0] = centroid[xyz[0]] - center[xyz[0]];
        diffC[1] = centroid[xyz[1]] - center[xyz[1]];
        double mag = std::sqrt(diffC[0]*diffC[0]+diffC[1]*diffC[1]);
        if (R > mag+100*Tiny) {
            RLat    rLat(-1,-1,0.0,NodeID);
            rNodeList.push_back(rLat);
            continue;
        }
        double minR = Huge;
        double dz;
        double a;
        bool isFound = false;

        for (unsigned p=0; p<periList[aPart].size(); p++) {
            std::array<double,2> diff;
            diff[0] = (centroid[xyz[0]] - LatTab[periList[aPart][p]].centroid[xyz[0]]);
            diff[1] = (centroid[xyz[1]] - LatTab[periList[aPart][p]].centroid[xyz[1]]);
            double r = (diffC[0]*diff[0]+diffC[1]*diff[1])/mag;
            double d = geo.dist(LatTab[periList[aPart][p]].centroid,centroid);
            double lavg =  std::sqrt(LatTab[periList[aPart][p]].area);
            if (1.15*lavg>d) {
                r /= 2.0;
            } else {
                r -= lavg/2.0;
            }
            if ((r<minR)&&(r>Tiny)) {
                minR = r;
                a = mag-r;
                dz = centroid[xyz[2]] - LatTab[periList[aPart][p]].centroid[xyz[2]];
                isFound = true;
            }
        }
        if (isFound) {
            RLat    rLat(minR,a,dz,NodeID);
            rNodeList.push_back(rLat);
        } else {
            RLat    rLat(-1,-1,0,NodeID);
            rNodeList.push_back(rLat);
        }
    }
    std::sort (rNodeList.begin(),rNodeList.end());
    dualOut << "Total rNode = " << rNodeList.size() << std::endl;
    return rNodeList;
}

std::vector<Statistics::RLat> Statistics::getRLatList ( std::vector<std::vector<unsigned> >      &periList,
                                                        std::vector<unsigned>*                   p_latList,
                                                        const Point                              center)
{
    dualOut << "Starting Statistics::getRLatList..."<< std::endl;

    std::vector<RLat>       rLatList;
    if (p_latList->empty()) {
        return rLatList;
    }
    unsigned res = periList.size();
    double dl = 2*Pi/res;
    Geometry geo;
    for (unsigned i=0; i<p_latList->size(); i++) {
        unsigned LatticeID = (*p_latList)[i];
        Point centroid = LatTab[LatticeID].centroid;
        double angle = std::atan2(centroid[xyz[1]]-center[xyz[1]],centroid[xyz[0]]-center[xyz[0]]);
        if (angle<0.0) {
            angle += 2*Pi;
        }
        unsigned aPart = (angle-Tiny)/dl;
        if (aPart>=res) {
            aPart = res -1;
        }
        double R = 0.0;
        for (unsigned p=0; p<periList[aPart].size(); p++) {
            R += geo.dist(LatTab[periList[aPart][p]].centroid,center);
        }
        R /= periList[aPart].size();
        std::array<double,2> diffC;
        diffC[0] = centroid[xyz[0]] - center[xyz[0]];
        diffC[1] = centroid[xyz[1]] - center[xyz[1]];
        double mag = std::sqrt(diffC[0]*diffC[0]+diffC[1]*diffC[1]);
        if (R > mag+100*Tiny) {
            RLat    rLat(-1,-1,0.0,LatticeID);
            rLatList.push_back(rLat);
            continue;
        }
        double minR = Huge;
        double dz;
        double a;
        bool isFound = false;
        for (unsigned p=0; p<periList[aPart].size(); p++) {
            std::array<double,2> diff;
            diff[0] = (centroid[xyz[0]] - LatTab[periList[aPart][p]].centroid[xyz[0]]);
            diff[1] = (centroid[xyz[1]] - LatTab[periList[aPart][p]].centroid[xyz[1]]);
            double r = (diffC[0]*diff[0]+diffC[1]*diff[1])/mag;
            double d = geo.dist(LatTab[periList[aPart][p]].centroid,centroid);
            double lavg =  std::sqrt(LatTab[periList[aPart][p]].area);
            if (1.15*lavg>d) {
                r /= 2.0;
            } else {
                r -= lavg/2.0;
            }
            if ((r<minR)&&(r>Tiny)) {
                minR = r;
                a = mag-r;
                dz = centroid[xyz[2]] - LatTab[periList[aPart][p]].centroid[xyz[2]];
                isFound = true;
            }
        }
        if (isFound) {
            RLat    rLat(minR,a,dz,LatticeID);
            rLatList.push_back(rLat);
        } else {
            RLat    rLat(-1,-1,0,LatticeID);
            rLatList.push_back(rLat);
        }
    }
    std::sort (rLatList.begin(),rLatList.end());
    dualOut << "Total rLat = " << rLatList.size() << std::endl;
    return rLatList;
}

void Statistics::writeRLatStats     (    ConnectLat*        p_conLat,
                                         PCG_Eigen*         p_EigenPCG,
                                         const double       sumF,
                                         const unsigned     step)
{
    dualOut << "Writing rLatStats..." << '\n';
    Clock sClock;
    sClock.start("rLatStats");

    char    fileName[255];
    sprintf(fileName,"%s/%sstat-rLat%04d.txt",path,prefix,step);
    std::ofstream file(fileName, std::ios::out);
    file << "# sumF = " << sumF << std::endl;
    file << "# r" << '\t' << "a" << '\t' << "r/a" << '\t' << "dz/avgD" << '\t' << "Sigma_n" << '\t' << "Sigma_s" << std::endl;
    Point center = p_conLat->getMainCrackCenter();

    if (CrackRadius>Tiny) {
        fillInnerLatList(0.02,getLimits(center),p_conLat);
    } else {
        fillInnerLatList(0.02,p_conLat);
    }
    std::vector<std::vector<unsigned> > periList = p_conLat->getPerimeter(center,avgD,xyz);
    std::vector<RLat> rLatList = getRLatList(periList,center);
    std::vector<std::pair<double,std::pair<double,double> > > s_nList;
    std::vector<std::pair<double,std::pair<double,double> > > s_sList;
    for (unsigned i=0; i<rLatList.size(); i++) {
        if ((rLatList[i].a>0.0)&&(rLatList[i].r>0.0)) {
            unsigned LatticeID = rLatList[i].LatticeID;
            Tensor memForce = p_EigenPCG->getMemForce(LatticeID);
            double s_n = memForce[0]/LatTab[LatticeID].area;
            double s_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
            std::pair<double,double> pos2D(rLatList[i].r/rLatList[i].a,rLatList[i].dz/avgD);
            std::pair<double,std::pair<double,double> > tri(s_n/sumF,pos2D);
            s_nList.push_back(tri);
            tri.first = s_s/sumF;
            s_sList.push_back(tri);
            file << rLatList[i].r << '\t' << rLatList[i].a << '\t' << rLatList[i].r/rLatList[i].a << '\t'
                 << rLatList[i].dz/avgD << '\t' << s_n/sumF << '\t' << s_s/sumF << std::endl;
        }
    }
    printData2D(s_nList,"s_n",sumF,step,100,10);
    printData2D(s_sList,"s_s",sumF,step,100,10);
    sClock.get();
}

void Statistics::fill_eLatIDList()
{
    unsigned LatticeID;

    for (LatticeID=0; LatticeID < tLatticeNum; LatticeID++)
    {
        if ((LatTab[LatticeID].k[0] > 2*Tiny)&&(LatTab[LatticeID].k[0] < Huge/2.0))
        {
            eLatIDList.push_back(LatticeID);
        }
    }
}

void Statistics::fill_eLatIDListSimple()
{
    unsigned LatticeID;

    eLatIDList.clear();

    for (LatticeID=0; LatticeID < tLatticeNum; LatticeID++)
    {
        eLatIDList.push_back(LatticeID);
    }
}

void Statistics::getPlaneInfo 	(	const string 		plane,
                                    unsigned* 			xyz0,
                                    unsigned* 			xyz1,
                                    unsigned* 			ID)
{
    if (plane=="xy") {
        *xyz0 = 0;
        *xyz1 = 1;
        *ID = 0;
    }
    else if (plane=="xz")
    {
        *xyz0 = 0;
        *xyz1 = 2;
        *ID = 1;
    }
    else if (plane=="yz")
    {
        *xyz0 = 1;
        *xyz1 = 2;
        *ID = 2;
    }
    else
    {
        dualOut << "Statistics::latDir - Wrong input of plane!" << endl;
    }
}

void Statistics::getInnerBoundaryNodeList ()
{
	unsigned nNum;
	std::vector<unsigned>	inNodeList;
	std::cout << "Starting Statistics::getInnerBoundaryNodeList... \n";
	for (unsigned NodeID=0; NodeID < tInNodeNum; NodeID++) {
		nNum = GNodeTab[NodeID].nbNodeID.size();
		for (unsigned k=0; k<nNum; k++) {
			if (GNodeTab[NodeID].nbNodeID[k]>=tInNodeNum) {
				inNodeList.push_back(NodeID);
			}
		}
	}
	nNum = inNodeList.size();
	double avgDist = pow(((double) Nx*Ny*Nz*UnitLength*UnitLength*UnitLength)/tInNodeNum,1/Dim);
	double xUplimit = Nx*UnitLength-avgDist;
	double yUplimit = Ny*UnitLength-avgDist;
	double zUplimit = Nz*UnitLength-avgDist;
	unsigned NodeID;
	for (unsigned i=0; i<nNum; i++) {
		NodeID = inNodeList[i];
		if (GNodeTab[NodeID].coord[0]<avgDist) 		{ inBoundary[Face5].push_back(NodeID); }
		if (GNodeTab[NodeID].coord[0]>xUplimit) 	{ inBoundary[Face3].push_back(NodeID); }
		if (GNodeTab[NodeID].coord[1]<avgDist) 		{ inBoundary[Face2].push_back(NodeID); }
		if (GNodeTab[NodeID].coord[1]>yUplimit) 	{ inBoundary[Face4].push_back(NodeID); }
		if (GNodeTab[NodeID].coord[2]<avgDist) 		{ inBoundary[Face1].push_back(NodeID); }
		if (GNodeTab[NodeID].coord[2]>zUplimit) 	{ inBoundary[Face6].push_back(NodeID); }
	}
}

void Statistics::updateStat 	(	const double 		var,
                                    const unsigned 		resolution,
                                    const double 		min,
                                    const double 		max,
                                    vector<unsigned>* 	p_Stat)
{
    unsigned i;
    double range = max - min;
    double step = range/resolution;

    for (i=0; i<resolution; i++)
    {
        if ((var > min+i*step)&&(var < min+(i+1)*step))
        {
            (*p_Stat)[i]++;
            break;
        }
    }
}

unsigned Statistics::fillAngList (  const unsigned      resolution,
                                    ConnectLat*         p_conLat)
{
    const unsigned minRes = 10;
    const unsigned minAvgSize = 12;
    unsigned res = resolution;
    if (res < minRes) {
        res = minRes;
    }
    if (innerLatList.size()/res <minAvgSize) {
        res = innerLatList.size()/minAvgSize;
    }
    unsigned N = innerLatList.size();
    for (unsigned d=0; d<3; d++) {
        latAngList[d].clear();
        latAngList[d].resize(res);
    }
    std::array<unsigned,3> xyz1,xyz2;
    xyz1[0] = 1; xyz2[0] = 0;
    xyz1[1] = 2; xyz2[1] = 0;
    xyz1[2] = 2; xyz2[2] = 1;
    double dl = 2.0*Pi/res;
    unsigned half = res/2;
    fillInnerLatList(0.02,p_conLat);
    for (unsigned i=0; i<innerLatList.size(); i++) {
        unsigned LatticeID = innerLatList[i];
        for (unsigned d=0; d<3; d++) {
            double angle = std::atan2(LatTab[LatticeID].axes[0][xyz1[d]],LatTab[LatticeID].axes[0][xyz2[d]]);
            if ((angle > Pi/2.0)&&(angle <= Pi+Tiny)) {
                angle = Pi - angle;
            } else if ((angle < 0.0)&&( angle >= -Pi/2.0)) {
                angle = -angle;
            } else if (angle < -Pi/2.0) {
                angle = Pi+angle;
            }
            unsigned aPart = (angle-Tiny)/dl;
            if (aPart>=latAngList[d].size()) {
                tripOut << "aPart>=latAngList[d].size(), aPart = " << aPart << std::endl;
                aPart = latAngList[d].size()-1;
            }
            latAngList[d][aPart].push_back(LatticeID);
            latAngList[d][half-1-aPart].push_back(LatticeID);
            latAngList[d][half+aPart].push_back(LatticeID);
            latAngList[d][res-1-aPart].push_back(LatticeID);
        }
    }

    char    fileName[255];
    std::sprintf(fileName,"%s/%sstat-LatAng.txt",path,prefix);
    std::ofstream file(fileName);
    std::time_t  now     = time(0);
    char*   dt      = ctime(&now);
    const double scaling = ((double) res)/minAvgSize/N/4;

    file  << "# File - Statistics of lattice direction distribution " << std::endl;
    file  << "# Saved in " << fileName << '\t' <<  "created on " << dt << std::endl;
    file  << "# " << "Angle" << '\t' << "Frequency_xy" <<  '\t' << "Frequency_xz" << '\t' << "Frequency_yz" << std::endl;
    for (unsigned i=0; i<res; i++) {
        file << (i)*dl*180.0/Pi <<  '\t'
             << ((double) latAngList[0][i].size())*scaling <<  '\t'
             << ((double) latAngList[1][i].size())*scaling <<  '\t'
             << ((double) latAngList[2][i].size())*scaling <<  std::endl
             << (i+1)*dl*180.0/Pi <<  '\t'
             << ((double) latAngList[0][i].size())*scaling <<  '\t'
             << ((double) latAngList[1][i].size())*scaling <<  '\t'
             << ((double) latAngList[2][i].size())*scaling <<  std::endl;
    }
    file << "N = " << N << std::endl;
    return res;
}

void Statistics::writeLatCluster    (   ConnectLat*         p_conLat,
                                        const unsigned      resolution,
                                        const unsigned      step)
{
    std::vector<unsigned>   clusterList = p_conLat->getClusterList();
    if (clusterList.empty()) {
        return;
    }
    printCluster (&clusterList,"Cluster",step,resolution);
}

void Statistics::writeLatClusterAgg    (   ConnectLat*         p_conLat,
                                           const unsigned      step)
{
    static bool isFirst = true;
    std::vector<unsigned>   clusterList = p_conLat->getClusterList();
    if (clusterList.empty()) {
        return;
    }
    char    fileName[255];
    std::sprintf(fileName,"%s/%sstat_ClusterAgg.txt",path,prefix);
    std::ofstream inFile(fileName, std::ios::out | std::ios::app);


    if (isFirst) {
        inFile  << "# File - Statistics of ClusterAgg " << std::endl;
        inFile  << "# step" << '\t' << "Total FailLatNum" << '\t'
                << "1st Cluster LatNum" << '\t' << "2nd Cluster LatNum"  << '\t' << "3rd Cluster LatNum" << '\t'
                << "4th Cluster LatNum"<< '\t' << "5th Cluster LatNum"<< '\t' << "6th Cluster LatNum" << '\t'
                << "7th Cluster LatNum"<< '\t' << "8th Cluster LatNum"<< '\t' << "9th Cluster LatNum" << '\t'
                << "10th Cluster LatNum" << '\t'
                << "Isolated LatNum 1" << '\t' << "Isolated LatNum 2" << '\t' << "Isolated LatNum 3" << '\t'
                << "Isolated LatNum 4" << '\t' << "Isolated LatNum 5" << '\t' << "Isolated LatNum 6" << '\t'
                << "Isolated LatNum 7" << '\t' << "Isolated LatNum 8" << '\t' << "Isolated LatNum 9" << '\t'
                << "Isolated LatNum 10" << std::endl;
        isFirst = false;
    }
    unsigned isoLat = 10;
    std::vector<unsigned> isoLatList(isoLat,0);
    std::vector<unsigned> cluster;
    std::vector<std::pair<unsigned,unsigned> >   Cnt;
    int i;
    unsigned currentLatNum = 1;
    for (i=clusterList.size()-1; i>=0; i--) {
        if (clusterList[i]>isoLat) {
            break;
        } else if (clusterList[i]==currentLatNum) {
            isoLatList[currentLatNum-1]+=currentLatNum;
        } else {
            currentLatNum = clusterList[i];
            isoLatList[currentLatNum-1]+=currentLatNum;
        }
    }
    inFile << step << '\t' << p_conLat->getLatNum() << '\t';
//    unsigned outClusterNum = clusterList.size();
    for (unsigned d=0; d<10; d++) {
        if (d<clusterList.size()) {
            inFile << clusterList[d] << '\t';
        } else {
            inFile << 0 << '\t';
        }
    }
    for (unsigned d=0; d<isoLat; d++) {
        inFile << isoLatList[d] << '\t';
    }
    inFile << std::endl;
}
LatPriDirList Statistics::fillLatPriDirList (  const unsigned      resolution,
                                               const unsigned      step)
{
    dualOut << "Starting Statistics::fillLatPriDirList... \n";
    Clock   sClock;
    sClock.start("stat_fillLatPriDirList");
    LatPriDirList       latPriDirList;
    const unsigned minRes = 10;
    const unsigned minAvgSize = 12;
    unsigned res = resolution;
    if (res < minRes) {
        res = minRes;
    }
    std::vector<unsigned> nList { std::begin(nodeLists.free),std::end(nodeLists.free) };
    for (unsigned d=0; d<6; d++) {
        for (unsigned i=0; i<nodeLists.boundary[d].size(); i++) {
            eraseValue(nList,nodeLists.boundary[d][i]);
        }
    }
    unsigned N = nList.size();
    if (nList.size()/res < minAvgSize) {
        res = nList.size() / minAvgSize;
    }
    for (unsigned dd=0; dd<3; dd++) {
        for (unsigned d=0; d<3; d++) {
            latPriDirList[dd][d].clear();
            latPriDirList[dd][d].resize(resolution);
        }
    }
    std::array<unsigned,3> xyz1,xyz2;
    xyz1[0] = 1; xyz2[0] = 0;
    xyz1[1] = 2; xyz2[1] = 0;
    xyz1[2] = 2; xyz2[2] = 1;
    double dl = Pi/2.0/res;
    Stress stress;
    stress.calBasic();
    stress.calPriStress();

    for (unsigned i=0; i<nList.size(); i++) {
        unsigned NodeID = nList[i];
        std::array<Vec,3> v;
        v[0] = stress.getPriDir(NodeID,0);
        v[1] = stress.getPriDir(NodeID,1);
        v[2] = stress.getPriDir(NodeID,2);
        for (unsigned d=0; d<3; d++) {
            for (unsigned dd=0; dd<3; dd++) {
                double angle = std::atan2(v[dd][xyz1[d]],v[dd][xyz2[d]]);
                if ((angle > Pi/2.0)&&(angle <= Pi+Tiny)) {
                    angle = Pi - angle;
                } else if ((angle < 0.0)&&( angle >= -Pi/2.0)) {
                    angle = -angle;
                } else if (angle < -Pi/2.0) {
                    angle = Pi+angle;
                }
                unsigned aPart = (angle-Tiny)/dl;
                if (aPart>=latPriDirList[dd][d].size()) {
                    tripOut << "aPart>=latPriDirList[dd][d].size(), aPart = " << aPart << std::endl;
                    aPart = latPriDirList[dd][d].size() -1;
                }
                latPriDirList[dd][d][aPart].push_back(NodeID);
            }
        }
    }

    char    fileName[255];
    std::sprintf(fileName,"%s/%sstat-LatPriDirAng_%04d.txt",path,prefix,step);
    std::ofstream file(fileName);
    std::time_t  now     = time(0);
    char*   dt      = ctime(&now);
    const double scaling = ((double) res)/minAvgSize/N;
    file  << "# File - Statistics of node principal directions distribution " << std::endl;
    file  << "# Saved in " << fileName << '\t' <<  "created on " << dt << std::endl;
    file  << "# " << "Angle" << '\t' << "PD1 - Frequency_xy" <<  '\t' << "PD1 - Frequency_xz" << '\t' << "PD1 - Frequency_yz" << '\t'
                                     << "PD2 - Frequency_xy" <<  '\t' << "PD2 - Frequency_xz" << '\t' << "PD2 - Frequency_yz" << '\t'
                                     << "PD3 - Frequency_xy" <<  '\t' << "PD3 - Frequency_xz" << '\t' << "PD3 - Frequency_yz" << std::endl;

    //1st quarter
    for (unsigned i=0; i<res; i++) {
        file << (i)*dl*180.0/Pi <<  '\t'
             << ((double) latPriDirList[0][0][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][2][i].size())*scaling << '\t'
             << ((double) latPriDirList[1][0][i].size())*scaling << '\t'
             << ((double) latPriDirList[1][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[1][2][i].size())*scaling << '\t'
             << ((double) latPriDirList[2][0][i].size())*scaling << '\t'
             << ((double) latPriDirList[2][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[2][2][i].size())*scaling << std::endl
             << (i+1)*dl*180.0/Pi <<  '\t'
             << ((double) latPriDirList[0][0][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][2][i].size())*scaling << '\t'
             << ((double) latPriDirList[1][0][i].size())*scaling << '\t'
             << ((double) latPriDirList[1][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[1][2][i].size())*scaling << '\t'
             << ((double) latPriDirList[2][0][i].size())*scaling << '\t'
             << ((double) latPriDirList[2][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[2][2][i].size())*scaling << std::endl;
    }
    //2nd quarter
    for (unsigned i=0; i<res; i++) {
        file << (i)*dl*180.0/Pi+90.0 <<  '\t'
             << ((double) latPriDirList[0][0][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][2][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[1][0][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[1][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[1][2][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[2][0][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[2][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[2][2][res-1-i].size())*scaling << std::endl
             << (i+1)*dl*180.0/Pi+90.0 <<  '\t'
             << ((double) latPriDirList[0][0][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][2][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[1][0][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[1][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[1][2][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[2][0][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[2][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[2][2][res-1-i].size())*scaling << std::endl;
    }
    //3rd quarter
    for (unsigned i=0; i<res; i++) {
        file << (i)*dl*180.0/Pi+180.0 <<  '\t'
             << ((double) latPriDirList[0][0][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][2][i].size())*scaling << '\t'
             << ((double) latPriDirList[1][0][i].size())*scaling << '\t'
             << ((double) latPriDirList[1][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[1][2][i].size())*scaling << '\t'
             << ((double) latPriDirList[2][0][i].size())*scaling << '\t'
             << ((double) latPriDirList[2][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[2][2][i].size())*scaling << std::endl
             << (i+1)*dl*180.0/Pi+180.0 <<  '\t'
             << ((double) latPriDirList[0][0][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][2][i].size())*scaling << '\t'
             << ((double) latPriDirList[1][0][i].size())*scaling << '\t'
             << ((double) latPriDirList[1][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[1][2][i].size())*scaling << '\t'
             << ((double) latPriDirList[2][0][i].size())*scaling << '\t'
             << ((double) latPriDirList[2][1][i].size())*scaling <<  '\t'
             << ((double) latPriDirList[2][2][i].size())*scaling << std::endl;
    }
    //4th quarter
    for (unsigned i=0; i<res; i++) {
        file << (i)*dl*180.0/Pi+270.0 <<  '\t'
             << ((double) latPriDirList[0][0][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][2][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[1][0][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[1][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[1][2][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[2][0][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[2][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[2][2][res-1-i].size())*scaling << std::endl
             << (i+1)*dl*180.0/Pi+270.0 <<  '\t'
             << ((double) latPriDirList[0][0][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[0][2][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[1][0][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[1][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[1][2][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[2][0][res-1-i].size())*scaling << '\t'
             << ((double) latPriDirList[2][1][res-1-i].size())*scaling <<  '\t'
             << ((double) latPriDirList[2][2][res-1-i].size())*scaling << std::endl;
    }
    file << "N = " << N << std::endl;
    sClock.get();
    return latPriDirList;
}

bool Statistics::writeFailLatPolar (    const unsigned             resolution,
                                        std::vector<unsigned>      latList,
                                        const unsigned             step,
                                        const double               sumF,
                                        const char*                vName,
                                        PCG_Eigen*                 p_EigenPCG)
{
    dualOut << "Writing stat_writeFailLatPolar-" << vName << "..." << '\n';
    Clock   sClock;
    sClock.start("stat_writeFailLatPolar");
    std::array<std::vector<std::vector<unsigned> >,3>    latAng;
    const unsigned minRes = 8;
    const unsigned minSize = 12;
    unsigned N = latList.size();
    if (N <minRes*minSize) {
        dualOut << "The number of failure lattice is too small N = " << N
                << ". No stat is outputted!" << std::endl;
        return false;
    }
    unsigned res = resolution*4/4;
    if (latList.size()/res<minSize) {
        res = N/minSize*4/4;
    }
    for (unsigned d=0; d<3; d++) {
        latAng[d].resize(res);
    }
    std::array<unsigned,3> xyz1,xyz2;
    xyz1[0] = 1; xyz2[0] = 0;
    xyz1[1] = 2; xyz2[1] = 0;
    xyz1[2] = 2; xyz2[2] = 1;
    double dl = 2.0*Pi/res;
    unsigned half = res/2;
    for (unsigned i=0; i<N; i++) {
        unsigned LatticeID = latList[i];
        for (unsigned d=0; d<3; d++) {
            double angle = std::atan2(LatTab[LatticeID].axes[0][xyz1[d]],LatTab[LatticeID].axes[0][xyz2[d]]);
            if ((angle > Pi/2.0)&&(angle <= Pi+Tiny)) {
                angle = Pi - angle;
            } else if ((angle < 0.0)&&( angle >= -Pi/2.0)) {
                angle = -angle;
            } else if (angle < -Pi/2.0) {
                angle = Pi+angle;
            }
/*            if (angle < 0.0) {
                angle += 2*Pi;
            }*/
            unsigned aPart = (angle-Tiny)/dl;
            if (aPart>=latAng[d].size()) {
                tripOut << "aPart>=latAng[d].size(), aPart = " << aPart << std::endl;
                aPart = latAng[d].size()-1;
            }
            latAng[d][aPart].push_back(LatticeID);
            latAng[d][half-1-aPart].push_back(LatticeID);
            latAng[d][half+aPart].push_back(LatticeID);
            latAng[d][res-1-aPart].push_back(LatticeID);
        }
    }

    char    fileName[255];
    std::sprintf(fileName,"%s/%sstat-failLatAng_%s-%04d.txt",path,prefix,vName,step);
    std::ofstream file(fileName);
    std::time_t  now     = time(0);
    char*   dt      = ctime(&now);

    file  << "# File - Statistics of failure lattice direction distribution , step = " << step << std::endl;
    file  << "# Saved in " << fileName << '\t' <<  "created on " << dt << std::endl;
    file  << "# " << "Angle" << '\t' << "Frequency_xy" <<  '\t' << "Frequency_xz" << '\t' << "Frequency_yz"
          << '\t' << "Frequency_xy/N" <<  '\t' << "Frequency_xz/N" << '\t' << "Frequency_yz/N"<< std::endl;
    N *=4;
    for (unsigned i=0; i<res; i++) {
        file << (i)*dl*180.0/Pi <<  '\t' << latAng[0][i].size()*res/minSize <<  '\t' << latAng[1][i].size()*res/minSize <<  '\t'
             << latAng[2][i].size()*res/minSize <<  '\t' << ((double) latAng[0][i].size())/((double) N) << '\t'
             << ((double) latAng[1][i].size())/((double) N) <<  '\t' << ((double) latAng[2][i].size())/ ((double) N) <<  std::endl;
        file << (i+1)*dl*180.0/Pi <<  '\t' << latAng[0][i].size()*res/minSize <<  '\t' << latAng[1][i].size()*res/minSize <<  '\t'
             << latAng[2][i].size()*res/minSize <<  '\t' << ((double) latAng[0][i].size())/((double) N) << '\t'
             << ((double) latAng[1][i].size())/((double) N) <<  '\t' << ((double) latAng[2][i].size())/ ((double) N) <<  std::endl;
    }
    sClock.get();
    return true;
}

void Statistics::writelatPolar (  const unsigned      resolution,
                                  const unsigned      step,
                                  const double        sumF,
                                  ConnectLat*         p_conLat,
                                  PCG_Eigen*          p_EigenPCG)
{
    dualOut << "Writing stat_writelatPolar..." << '\n';
    Clock   sClock;
    sClock.start("stat_writelatPolar");
    unsigned res = fillAngList(resolution,p_conLat);

    if (res==0) {
        tripOut << "Too few lattice, no stats is outputted!" << std::endl;
        return;
    }

    std::vector<std::vector<double> > snList;
    std::vector<std::vector<double> > ssList;
    std::vector<std::vector<double> > snRatioList;
    std::vector<std::vector<double> > ssRatioList;
    std::vector<std::vector<double> > ss_snList;
    snList.resize(3);
    ssList.resize(3);
    ss_snList.resize(3);
    snRatioList.resize(3);
    ssRatioList.resize(3);

    for (unsigned d=0; d<3; d++) {
        snList[d].resize(res);
        ssList[d].resize(res);
        ss_snList[d].resize(res);
        snRatioList[d].resize(res);
        ssRatioList[d].resize(res);
    }

    for (unsigned d=0; d<3; d++) {
        for (unsigned i=0; i<res; i++) {
            double sum_sn = 0.0;
            double sum_ss = 0.0;
            double sum_snR = 0.0;
            double sum_ssR = 0.0;
            for (unsigned k=0; k<latAngList[d][i].size(); k++) {
                unsigned LatticeID = latAngList[d][i][k];
                Tensor memForce = p_EigenPCG -> getMemForce(LatticeID);
                double s_n = memForce[0]/LatTab[LatticeID].area;
                double s_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
                {
                    sum_sn += s_n;
                    sum_ss += s_s;

                    sum_snR += std::max(s_n/LatTab[LatticeID].et[0],s_n/LatTab[LatticeID].et[1]);
                    sum_ssR += s_s/LatTab[LatticeID].shearStr;
                }
            }
            if (latAngList[d][i].size()!=0) {
                sum_sn /= latAngList[d][i].size();
                sum_ss /= latAngList[d][i].size();

                sum_snR /= latAngList[d][i].size();
                sum_ssR /= latAngList[d][i].size();

                ss_snList[d][i] = sum_ss/sum_sn;
            } else {
                ss_snList[d][i] = 0.0;
            }
            snList[d][i] = sum_sn/sumF;
            ssList[d][i] = sum_ss/sumF;

            snRatioList[d][i] = sum_snR;
            ssRatioList[d][i] = sum_ssR;
        }
    }
    char    fileName[255];
    std::sprintf(fileName,"%s/%sstat_LatPolar-%04d.txt",path,prefix,step);
    std::ofstream file(fileName); //, ios::out | ios::app);

    time_t  now     = time(0);
    char*   dt      = ctime(&now);

    file  << "# File - Statistics of Lattice data in Polar form " << std::endl;
    file  << "# Saved in " << fileName << '\t' << "Step = " << step << '\t' <<  "created on " << dt << std::endl;
    file  << "# " << "Angle" << '\t' << "Sn_xy" <<  '\t' << "Sn_xz" << '\t' << "Sn_yz" << '\t'
            << "Ss_xy" << '\t' << "Ss_xz" << '\t' << "Ss_yz" << '\t'
            << "SnR_xy" <<  '\t' << "SnR_xz" << '\t' << "SnR_yz" << '\t'
            << "SsR_xy" << '\t'  << "SsR_xz" << '\t' << "SsR_yz" << '\t'
            << "Ss/sn_xy" << '\t' << "Ss/sn_xz" << '\t' << "Ss/sn_yz"<< std::endl;
    double dl = 2*Pi/res;
    // 1st quarter
    for (unsigned i=0; i<res; i++) {
        file << (i)*dl*180.0/Pi <<  '\t'
             << snList[0][i] <<  '\t' << snList[1][i] <<  '\t' << snList[2][i] <<  '\t'
             << ssList[0][i] <<  '\t' << ssList[1][i] <<  '\t' << ssList[2][i] <<  '\t'
             << snRatioList[0][i] <<  '\t' << snRatioList[1][i] <<  '\t' << snRatioList[2][i] <<  '\t'
             << ssRatioList[0][i] <<  '\t' << ssRatioList[1][i] <<  '\t' << ssRatioList[2][i] <<  '\t'
             << ss_snList[0][i] <<  '\t' << ss_snList[1][i] <<  '\t' << ss_snList[2][i] <<  std::endl;
        file << (i+1)*dl*180.0/Pi <<  '\t'
             << snList[0][i] <<  '\t' << snList[1][i] <<  '\t' << snList[2][i] <<  '\t'
             << ssList[0][i] <<  '\t' << ssList[1][i] <<  '\t' << ssList[2][i] <<  '\t'
             << snRatioList[0][i] <<  '\t' << snRatioList[1][i] <<  '\t' << snRatioList[2][i] <<  '\t'
             << ssRatioList[0][i] <<  '\t' << ssRatioList[1][i] <<  '\t' << ssRatioList[2][i] <<  '\t'
             << ss_snList[0][i] <<  '\t' << ss_snList[1][i] <<  '\t' << ss_snList[2][i] <<  std::endl;
    }
    sClock.get();
}

void Statistics::writelatPolar (  const unsigned      resolution,
                                  const unsigned      nc,
                                  const unsigned      step,
                                  const double        sumF,
                                  ConnectLat*         p_conLat,
                                  PCG_Eigen*          p_EigenPCG)
{
    dualOut << "Writing stat_writelatPolarFull..." << '\n';
    Clock   sClock;
    sClock.start("stat_writelatPolar");
    unsigned res = fillAngList(resolution,p_conLat);
    unsigned res0 = res;
    if (res==0) {
        tripOut << "[Statistics::writelatPolar] Too few lattice, no stats is outputted!" << std::endl;
        return;
    }
    if (res%4!=0) {
        tripOut << "[Statistics::writelatPolar] Resolution entered is not divisible by 4, process the data with care" << std::endl;
    }
    res /= 4;
    if (res%nc!=0) {
        tripOut << "[Statistics::writelatPolar] Resolution entered is not divisible by coarseness = " << nc << " , process the data with care" << std::endl;
    }
    res /= nc;
    std::array<std::vector<std::vector<unsigned> >,3> latList;
    std::array<std::vector<std::vector<double> >,3> snList;
    std::array<std::vector<std::vector<double> >,3> ssList;
    std::array<std::vector<std::vector<double> >,3> snRatioList;
    std::array<std::vector<std::vector<double> >,3> ssRatioList;

    for (unsigned d=0; d<3; d++) {
        latList[d].resize(res);
        snList[d].resize(res);
        ssList[d].resize(res);
        snRatioList[d].resize(res);
        ssRatioList[d].resize(res);
    }
    unsigned quad = res0/4;
    for (unsigned d=0; d<3; d++) {
        for (unsigned i=0; i<res; i++) {
            for (unsigned n=0; n<nc; n++) {
                latList[d][i].insert(latList[d][i].end(),latAngList[d][i*nc+n].begin(),latAngList[d][i*nc+n].end());
                latList[d][i].insert(latList[d][i].end(),latAngList[d][2*quad-1-(i*nc+n)].begin(),latAngList[d][2*quad-1-(i*nc+n)].end());
                latList[d][i].insert(latList[d][i].end(),latAngList[d][2*quad+i*nc+n].begin(),latAngList[d][2*quad+i*nc+n].end());
                latList[d][i].insert(latList[d][i].end(),latAngList[d][4*quad-1-(i*nc+n)].begin(),latAngList[d][4*quad-1-(i*nc+n)].end());
            }
        }
    }
    for (unsigned d=0; d<3; d++) {
        for (unsigned i=0; i<res; i++) {
            for (unsigned k=0; k<latList[d][i].size(); k++) {
                unsigned LatticeID = latList[d][i][k];
                Tensor memForce = p_EigenPCG -> getMemForce(LatticeID);
                double s_n = memForce[0]/LatTab[LatticeID].area;
                snList[d][i].push_back(s_n/sumF);
                snRatioList[d][i].push_back(std::max(s_n/LatTab[LatticeID].et[0],s_n/LatTab[LatticeID].et[1]));
                if (ShearToAxialStiffnessRatio>Tiny) {
                    double s_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
                    ssList[d][i].push_back(s_s/sumF);
                    ssRatioList[d][i].push_back(s_s/LatTab[LatticeID].shearStr);
                }
            }
        }
    }

    printPDF (snList[0],"Polar_sn_xy",step,res0,tailLen);
    printPDF (snList[1],"Polar_sn_xz",step,res0,tailLen);
    printPDF (snList[2],"Polar_sn_yz",step,res0,tailLen);

    printPDF (snRatioList[0],"Polar_snRatio_xy",step,res0,0.0);
    printPDF (snRatioList[1],"Polar_snRatio_xz",step,res0,0.0);
    printPDF (snRatioList[2],"Polar_snRatio_yz",step,res0,0.0);
    if (ShearToAxialStiffnessRatio>Tiny) {
        printPDF (ssList[0],"Polar_ss_xy",step,res0,tailLen);
        printPDF (ssList[1],"Polar_ss_xz",step,res0,tailLen);
        printPDF (ssList[2],"Polar_ss_yz",step,res0,tailLen);

        printPDF (ssRatioList[0],"Polar_ssRatio_xy",step,res0,0.0);
        printPDF (ssRatioList[1],"Polar_ssRatio_xz",step,res0,0.0);
        printPDF (ssRatioList[2],"Polar_ssRatio_yz",step,res0,0.0);
    }
    sClock.get();

}

void Statistics::latDir		(	const char* 		plane,
                                const unsigned 		resolution)
{
    unsigned i,LatticeID,cnt;
    unsigned nbNodeID[2];
    unsigned pID = 0;

    unsigned xyz[2];

    double angle,cntP;

    getPlaneInfo(plane,&xyz[0],&xyz[1],&pID);

    statLatAng[pID].resize(resolution,0);
    mean_dir[pID] = 0.0;
    sd_dir[pID] = 0.0;

    for (unsigned i = 0; i < eLatIDList.size(); i++) {
        LatticeID = eLatIDList[i];
        nbNodeID[0] = LatTab[LatticeID].nb[0];
        nbNodeID[1] = LatTab[LatticeID].nb[1];

        angle = atan2(GNodeTab[nbNodeID[1]].coord[xyz[1]]-GNodeTab[nbNodeID[0]].coord[xyz[1]],
                      GNodeTab[nbNodeID[1]].coord[xyz[0]]-GNodeTab[nbNodeID[0]].coord[xyz[0]]);
        if (angle < 0.0) {
            angle += Pi;
        }

        eLatAngList[pID].push_back(angle);
        mean_dir[pID] +=angle;

        updateStat(angle,resolution,0.0,Pi,&statLatAng[pID]);
    }

    mean_dir[pID] /= eLatAngList[pID].size();

    //Calculation the sample standard deviation
    for (i = 0; i < eLatIDList.size(); i++)
    {
        sd_dir[pID] += pow(eLatAngList[pID][i] - mean_dir[pID],2);
    }

    sd_dir[pID] = sqrt(sd_dir[pID]/eLatAngList[pID].size());

    cnt 	= 0;
    cntP	= 0.0;

    print (&statLatAng[pID], resolution,0.0,180.0,pID, "dir", plane);

    for (i=0; i<resolution; i++)
    {
        cnt += statLatAng[pID][i];
        cntP += statLatAng[pID][i]*100.0/((double) eLatAngList[pID].size());
    }

}

void Statistics::glatDir		(	const char* 		plane,
                                    const unsigned 		resolution)
{
    unsigned i,LatticeID;
    unsigned nbNodeID[2];
    unsigned pID = 0;

    unsigned xyz[2];

    double angle;

    getPlaneInfo(plane,&xyz[0],&xyz[1],&pID);

    cnt_dir[pID].clear();

    if (lastRecordStep==0)
    {
        statLatAng[pID].resize(resolution,0);
    }
    mean_dir[pID] = 0.0;
    sd_dir[pID] = 0.0;

    for (i = lastLatNum; i < eLatIDList.size(); i++)
    {
        LatticeID = eLatIDList[i];
        nbNodeID[0] = LatTab[LatticeID].nb[0];
        nbNodeID[1] = LatTab[LatticeID].nb[1];

        angle = atan2(GNodeTab[nbNodeID[1]].coord[xyz[1]]-GNodeTab[nbNodeID[0]].coord[xyz[1]],
                      GNodeTab[nbNodeID[1]].coord[xyz[0]]-GNodeTab[nbNodeID[0]].coord[xyz[0]]);
        if (angle < 0.0)
        {
            angle += Pi;
        }

        eLatAngList[pID].push_back(angle);
        mean_dir[pID] +=angle;

        updateStat(angle,resolution,0.0,Pi,&statLatAng[pID]);
    }

    mean_dir[pID] /= eLatAngList[pID].size();

    //Calculation the sample standard deviation
    for (i = 0; i < eLatIDList.size(); i++)
    {
        sd_dir[pID] += pow(eLatAngList[pID][i] - mean_dir[pID],2);
    }

    sd_dir[pID] = sqrt(sd_dir[pID]/eLatAngList[pID].size());

    print (&statLatAng[pID], resolution,0.0,Pi,pID, "dir", plane);

}

void Statistics::glatDir		(	const char* 		plane,
                                    const unsigned 		resolution,
                                    const unsigned 		step)
{
    unsigned i,LatticeID;
    unsigned nbNodeID[2];
    unsigned pID = 0;

    unsigned xyz[2];

    double angle;

    char	varName[20];

    getPlaneInfo(plane,&xyz[0],&xyz[1],&pID);
    cnt_dir[pID].clear();

    if (lastRecordStep==0)
    {
        statLatAng[pID].resize(resolution,0);
    }
    mean_dir[pID] = 0.0;
    sd_dir[pID] = 0.0;

    for (i = lastLatNum; i < eLatIDList.size(); i++)
    {
        LatticeID = eLatIDList[i];
        nbNodeID[0] = LatTab[LatticeID].nb[0];
        nbNodeID[1] = LatTab[LatticeID].nb[1];

        angle = atan2(GNodeTab[nbNodeID[1]].coord[xyz[1]]-GNodeTab[nbNodeID[0]].coord[xyz[1]],
                      GNodeTab[nbNodeID[1]].coord[xyz[0]]-GNodeTab[nbNodeID[0]].coord[xyz[0]]);
        if (angle < 0.0)
        {
            angle += Pi;
        }

        eLatAngList[pID].push_back(angle);
        mean_dir[pID] +=angle;

        updateStat(angle,resolution,0.0,Pi,&statLatAng[pID]);
    }

    mean_dir[pID] /= eLatAngList[pID].size();

    //Calculation the sample standard deviation
    for (i = 0; i < eLatIDList.size(); i++)
    {
        sd_dir[pID] += pow(eLatAngList[pID][i] - mean_dir[pID],2);
    }

    sd_dir[pID] = sqrt(sd_dir[pID]/eLatAngList[pID].size());

    sprintf(varName,"%04d-dir",step);

    print (&statLatAng[pID], resolution,0.0,Pi,pID, varName, plane);
}


void Statistics::print 		(	const vector<unsigned>* 	p_stat,
                                const unsigned 				resolution,
                                const double 				min,
                                const double 				max,
                                const unsigned 				pID,
                                const char* 				varName1,
                                const char* 				varName2)
{
    unsigned i;
    unsigned cnt = 0;
    double cntP = 0.0;
    double step = fabs(max-min)/resolution;
    double num;

    time_t 	now 	= time(0);
    char* 	dt 		= ctime(&now);
    char	fileName[255];

    sprintf(fileName,"%s/%sstat-%s-%s.txt",path,prefix,varName1,varName2);

    ofstream statFile(fileName); //, ios::out | ios::app);

    statFile << "# File - " << fileName << '\t' << varName1 << '\t'<< varName2 << '\t' <<  "created on " << dt << endl;
    statFile << "# Angle" << '\t' << "stat_no" <<  '\t' << "stat_percent" << endl;

    cnt =0;
    cntP=0.0;

    if (eLatIDList.size() ==0)
    {
        num = 1.0;
    }
    else
    {
        num = ((double) eLatIDList.size());
    }

    for (i= 0; i < resolution; i++)
    {
        statFile << min + (i+0.5)*step << '\t' << (*p_stat)[i] << '\t' << (*p_stat)[i]*100.0/((double) eLatIDList.size()) << endl;
        cnt += (*p_stat)[i];
        cntP += (*p_stat)[i]*100.0/num;
    }
    statFile << "# Total = " << cnt << " Total Percent = " << cntP << endl;
    statFile << "# Mean = " << mean_dir[pID]*180/Pi << " SD = " << sd_dir[pID]*180/Pi << endl;
}

void Statistics::latStiffness	(	const char* 		plane,
                                    const unsigned 		resolution)
{
    unsigned 		i, LatticeID;
    unsigned 		pID = 0;
    unsigned 		xyz[2];
    vector<double>	stiList;

    double ke;
    double min = Huge;
    double max = Tiny;

    getPlaneInfo(plane,&xyz[0],&xyz[1],&pID);
    statLatStiffness[pID].resize(resolution,0.0);

    for (i = 0; i < eLatIDList.size(); i++)
    {
        LatticeID = eLatIDList[i];
        ke = LatTab[LatticeID].k[0]*fabs(cos(eLatAngList[pID][i]));
        if (ke < min)
        {
            min = ke;
        }
        if (ke > max)
        {
            max = ke;
        }
        stiList.push_back(ke);
    }

    for (i = 0; i < eLatIDList.size(); i++)
    {
        updateStat(stiList[i],resolution,min-Tiny,max+Tiny,&statLatStiffness[pID]);
    }
    print (&statLatStiffness[pID], resolution,min,max, pID,"stiffness",plane);
}

void FailStats::write			(	const unsigned 			step,
                                    const unsigned			resolution)
{
    unsigned	i,j;

    unitAng = Pi/((double)resolution);
    unitAngInDeg = unitAng*180/Pi;

    dualOut << "Generating failure lattice statistical data..." << " step = "<< step << endl;

//	eLatIDList.clear();
    for (i=lastRecordStep; i<fLatPerStep.size(); i++)
    {
        for (j=0; j<fLatPerStep[i].size(); j++)
        {
            eLatIDList.push_back(fLatPerStep[i][j]);
        }
    }

    glatDir("xy",resolution,step);
    glatDir("xz",resolution,step);
    glatDir("yz",resolution,step);

    if (cnt%3==0) {
        fillInnerNodeList ();
        writeCoordNo(step,innerNodeList);
        writeLatLen(step,2*resolution);
        writeLatArea(step,10*resolution,false);
        writeLatBreakDisp(step,10*resolution,false);
    }
    lastRecordStep = step;
    lastLatNum = eLatIDList.size();
    cnt++;
}

void Statistics::writeStress	(	const std::vector<unsigned>*    p_weakPlaneLatID,
									const double                    sumF,
									const unsigned				    timeStep,
									const double				    time)
{
	static bool		isFirstTime = true;
	dualOut << "Writing stat_Stress..." << '\n';
	Clock   sClock;
	sClock.start("stat_Stress");
	char	fileName[255];
	sprintf(fileName,"%s/%sstat_Stress-%0d.txt",path,prefix,timeStep);
	ofstream fileS(fileName, ios::out | ios::app);
	sprintf(fileName,"%s/%sstat_PriStress-%0d.txt",path,prefix,timeStep);
	ofstream filePS(fileName, ios::out | ios::app);
	sprintf(fileName,"%s/%sstat_LatSurStress-%0d.txt",path,prefix,timeStep);
	ofstream fileLSS(fileName, ios::out | ios::app);
	sprintf(fileName,"%s/%sstat_stressConCell-%0d.txt",path,prefix,timeStep);
	ofstream fileSC2(fileName, ios::out | ios::app);
	sprintf(fileName,"%s/%sstat_stressAgg.txt",path,prefix);
	ofstream fileAgg(fileName, ios::out | ios::app);
	Stress stress;
	stress.calPriStress();
	sClock.get();
	double sqs_s11 = 0.0, sqs_s22 = 0.0, 	sqs_s33 = 0.0;
	double sqs_s1 = 0.0, sqs_s2 = 0.0, 		sqs_s3 = 0.0;
	double sum_s1 = 0.0, sum_s2 = 0.0,		sum_s3 = 0.0;
	double sum_s11 = 0.0, sum_s22 = 0.0, 	sum_s33 = 0.0;
	unsigned cntInNode = 0;
	std::vector<double>   s1List;
	std::vector<double>   s2List;
	std::vector<double>   s3List;
	std::vector<double>   s11List;
	std::vector<double>   s22List;
	std::vector<double>   s33List;
	std::vector<double>   s12List;
	std::vector<double>   s13List;
	std::vector<double>   s23List;

	std::vector<double>   pList;
	std::vector<double>   qList;

	Point centre;
	centre[0] = UnitLength*Nx*0.5;
	centre[1] = UnitLength*Ny*0.5;
	centre[2] = UnitLength*Nz*0.5;
	Geometry geo;
	const double radius = CrackRadius*Nx*UnitLength;
	double dist;
	fileS << "# sumF = " << '\t' <<  sumF << std::endl;
	filePS << "# sumF = " << '\t' << sumF << std::endl;
	fillInnerNodeList ();
	std::vector<unsigned> nList = innerNodeList;

	for (std::vector<unsigned>::iterator it = nList.begin(); it != nList.end(); ++it) {
		if (*it>=tInNodeNum) {
			break;
		}
		double s1 = stress.getPriStress(*it,0);	double s2 = stress.getPriStress(*it,1);	double s3 = stress.getPriStress(*it,2);
		double s11 = stress.getStress(*it,0);	double s22 = stress.getStress(*it,1);	double s33 = stress.getStress(*it,2);
		double s12 = stress.getStress(*it,3);   double s13 = stress.getStress(*it,4);   double s23 = stress.getStress(*it,5);
		s11List.push_back(-s11/sumF);       s22List.push_back(-s22/sumF);       s33List.push_back(-s33/sumF);
		s12List.push_back(-s12/sumF);       s13List.push_back(-s13/sumF);       s23List.push_back(-s23/sumF);
		s1List.push_back(-s1/sumF);        s2List.push_back(-s2/sumF);         s3List.push_back(-s3/sumF);

		pList.push_back(-(s1+s2+s3)/3.0/sumF);
		qList.push_back(std::sqrt(((s1-s2)*(s1-s2)+(s2-s3)*(s2-s3)+(s3-s1)*(s3-s1))/2.0)/std::fabs(sumF));

		if (!p_weakPlaneLatID->empty()) {
            if ((*it>=crackNodeNum)&&(*it<crackNodeNum+outerCrackNodeNum)) {
                dist = geo.dist(GNodeTab[*it].coord,centre);
                fileSC2 << dist - radius << '\t' << s33 << '\n';
            }
		}
		fileS 	<< s11 << ' ' << s22 << ' ' << s33 << '\n';
		filePS 	<< s1  << ' ' << s2  << ' ' << s3  << '\n';
		sum_s11 += s11; 		sum_s22 += s22; 		sum_s33 += s33;
		sqs_s11 += s11*s11; 	sqs_s22 += s22*s22; 	sqs_s33 += s33*s33;
		sum_s1  += s1;			sum_s2  += s2;			sum_s3  += s3;
		sqs_s1  += s1*s1;		sqs_s2  += s2*s2;		sqs_s3  += s3*s3;
		cntInNode++;
	}
	printPDF(&s11List,"s11",timeStep,100,tailLen);
	printPDF(&s22List,"s22",timeStep,100,tailLen);
	printPDF(&s33List,"s33",timeStep,100,tailLen);
	printPDF(&s12List,"s12",timeStep,100,tailLen);
	printPDF(&s13List,"s13",timeStep,100,tailLen);
	printPDF(&s23List,"s23",timeStep,100,tailLen);
	printPDF(&s1List,"s1",timeStep,100,tailLen);
	printPDF(&s2List,"s2",timeStep,100,tailLen);
	printPDF(&s3List,"s3",timeStep,100,tailLen);
	printPDF(&pList,"p",timeStep,100,tailLen);
	printPDF(&qList,"q",timeStep,100,tailLen);

	GLatForce l;
	double lss = 0.0;
	double lssSum = 0.0;
	double lssSqSum = 0.0;
	unsigned cntLat = 0;
	for (unsigned LatticeID = 0; LatticeID <tLatticeNum; LatticeID++) {
		if (!LatTab[LatticeID].isBreak) {
			if ((LatTab[LatticeID].nb[0] < tInNodeNum) && (LatTab[LatticeID].nb[1] < tInNodeNum)) {
				lss = l.getLatElongation(LatticeID)*LatTab[LatticeID].k[0]/LatTab[LatticeID].area;
				fileLSS << lss << '\n';
				lssSum += lss; lssSqSum += lss*lss;
				cntLat++;
			}
		}
	}
	if (isFirstTime) {
		fileAgg << "Time" << '\t' << "Sigma_xx (mean,SD)" << '\t' << "Sigma_yy (mean,SD)" << '\t' << "Sigma_zz (mean,SD)" << '\t'
				<< "Sigma_1 (mean,SD)" << '\t' << "Sigma_2 (mean,SD)" << '\t' << "Sigma_3 (mean,SD)" << '\t'
				<< "Sigma_suf (mean,SD)" << '\n';
		isFirstTime = false;
	}
	fileAgg << time << ' '
			<< sum_s11/cntInNode << ' ' << std::sqrt((sqs_s11 - sum_s11*sum_s11/cntInNode)/cntInNode) << ' '
			<< sum_s22/cntInNode << ' ' << std::sqrt((sqs_s22 - sum_s22*sum_s22/cntInNode)/cntInNode) << ' '
			<< sum_s33/cntInNode << ' ' << std::sqrt((sqs_s33 - sum_s33*sum_s33/cntInNode)/cntInNode) << ' '
			<< sum_s1/cntInNode << ' ' << std::sqrt((sqs_s1 - sum_s1*sum_s1/cntInNode)/cntInNode) << ' '
			<< sum_s2/cntInNode << ' ' << std::sqrt((sqs_s2 - sum_s2*sum_s2/cntInNode)/cntInNode) << ' '
			<< sum_s3/cntInNode << ' ' << std::sqrt((sqs_s3 - sum_s3*sum_s3/cntInNode)/cntInNode) << ' '
			<< lssSum/cntLat 	<< ' ' << std::sqrt((lssSqSum - lssSum*lssSum/cntLat)/cntLat) << '\n';
	sClock.get();
	dualOut << "Writing stat_Stress finished" << '\n';
}

void Statistics::writeStress    (   const std::vector<unsigned>*    p_weakPlaneLatID,
                                    ConnectLat*                     p_conLat,
                                    const double                    sumF,
                                    const unsigned                  timeStep,
                                    const double                    time)
{
    static bool     isFirstTime = true;
    dualOut << "Writing stat_Stress..." << '\n';
    Clock   sClock;
    sClock.start("stat_Stress");
    Stress stress;
    stress.calPriStress();
    sClock.get();
    unsigned cntInNode = 0;
    std::vector<double>   s1List;
    std::vector<double>   s2List;
    std::vector<double>   s3List;
    std::vector<double>   s11List;
    std::vector<double>   s22List;
    std::vector<double>   s33List;
    std::vector<double>   s12List;
    std::vector<double>   s13List;
    std::vector<double>   s23List;

    std::vector<double>   pList;
    std::vector<double>   qList;




    Point center = p_conLat->getMainCrackCenter();
    fillInnerNodeList (getLimits(center));
    std::vector<unsigned> nList = innerNodeList;

    std::vector<std::vector<unsigned> > periList = p_conLat->getPerimeter(center,avgD,xyz);
    std::vector<std::pair<double,std::pair<double,double> > >
        s1_2DList,s2_2DList,s3_2DList,s11_2DList,s22_2DList,s33_2DList,s12_2DList,s13_2DList,s23_2DList,p_2DList,q_2DList;

    std::vector<RLat> rNodeList = getRNodeList (periList,nList,center);
    if (rNodeList.size()==0) {
        return;
    }

    char    fileName[255];
    std::sprintf(fileName,"%s/%sstat_stressConCell-%0d.txt",path,prefix,timeStep);
    std::ofstream fileSC2(fileName, std::ios::out | std::ios::app);
    std::sprintf(fileName,"%s/%sstat_stressRaw-%0d.txt",path,prefix,timeStep);
    std::ofstream fileRaw(fileName, std::ios::out | std::ios::app);
    fileRaw << "# r/a" << '\t' << "dz/avgD" << '\t' << "s1" << '\t' << "s2" << '\t' << "s3"
            << '\t' << "s11" << '\t' << "s22" << '\t' << "s33" << '\t' << "s12" << '\t' << "s13" << '\t' << "s23"
            << '\t' << "p" << '\t' << "q" << std::endl;
    unsigned i = 0;
    Geometry geo;
    for (std::vector<unsigned>::iterator it = nList.begin(); it != nList.end(); ++it) {
        if (*it>=tInNodeNum) {
            break;
        }
        double s1 = stress.getPriStress(*it,0); double s2 = stress.getPriStress(*it,1); double s3 = stress.getPriStress(*it,2);
        double s11 = stress.getStress(*it,0);   double s22 = stress.getStress(*it,1);   double s33 = stress.getStress(*it,2);
        double s12 = stress.getStress(*it,3);   double s13 = stress.getStress(*it,4);   double s23 = stress.getStress(*it,5);

        s1 /=-sumF;                         s2 /=-sumF;                         s3 /=-sumF;
        s11 /=-sumF;                        s22 /=-sumF;                        s33 /=-sumF;
        s12 /=-sumF;                        s13 /=-sumF;                        s23 /=-sumF;

        double p = (s1+s2+s3)/3.0;
        double q = std::sqrt(((s1-s2)*(s1-s2)+(s2-s3)*(s2-s3)+(s3-s1)*(s3-s1))/2.0);

        s11List.push_back(s11);       s22List.push_back(s22);       s33List.push_back(s33);
        s12List.push_back(s12);       s13List.push_back(s13);       s23List.push_back(s23);
        s1List.push_back(s1);         s2List.push_back(s2);          s3List.push_back(s3);

        pList.push_back(p);
        qList.push_back(q);

        if ((rNodeList[i].a>0.0)&&(rNodeList[i].r>0.0)) {
            std::pair<double,double> pos2D(rNodeList[i].r/rNodeList[i].a,rNodeList[i].dz/avgD);
            std::pair<double,std::pair<double,double> > triplet(s1,pos2D);
            s1_2DList.push_back(triplet);
            triplet.first = s2; s2_2DList.push_back(triplet);
            triplet.first = s3; s3_2DList.push_back(triplet);
            triplet.first = s11; s11_2DList.push_back(triplet);
            triplet.first = s22; s22_2DList.push_back(triplet);
            triplet.first = s33; s33_2DList.push_back(triplet);
            triplet.first = s12; s12_2DList.push_back(triplet);
            triplet.first = s13; s13_2DList.push_back(triplet);
            triplet.first = s23; s23_2DList.push_back(triplet);
            triplet.first = p; p_2DList.push_back(triplet);
            triplet.first = q; q_2DList.push_back(triplet);
            fileRaw << rNodeList[i].r/rNodeList[i].a << '\t' << rNodeList[i].dz/avgD << '\t'
                    << s1 << '\t' << s2 << '\t' << s3 << '\t'
                    << s11 << '\t' << s22 << '\t'<< s33 << '\t'
                    << s12 << '\t' << s13 << '\t'<< s23 << '\t'
                    << p << '\t' << q << std::endl;
        }

        if (!p_weakPlaneLatID->empty()) {
            if ((*it>=crackNodeNum)&&(*it<crackNodeNum+outerCrackNodeNum)) {
                fileSC2 << rNodeList[i].r << '\t' << s33 << '\n';
            }
        }
        cntInNode++;
        i++;
    }
    unsigned res = 100;
    printPDF(&s11List,"s11",timeStep,res,tailLen);
    printPDF(&s22List,"s22",timeStep,res,tailLen);
    printPDF(&s33List,"s33",timeStep,res,tailLen);
    printPDF(&s12List,"s12",timeStep,res,tailLen);
    printPDF(&s13List,"s13",timeStep,res,tailLen);
    printPDF(&s23List,"s23",timeStep,res,tailLen);
    printPDF(&s1List,"s1",timeStep,res,tailLen);
    printPDF(&s2List,"s2",timeStep,res,tailLen);
    printPDF(&s3List,"s3",timeStep,res,tailLen);
    printPDF(&pList,"p",timeStep,res,tailLen);
    printPDF(&qList,"q",timeStep,res,tailLen);

    unsigned res1 = 200;
    unsigned res2 = 20;

    printData2D(s1_2DList,"s1",sumF,timeStep,res1,res2);
    printData2D(s2_2DList,"s2",sumF,timeStep,res1,res2);
    printData2D(s3_2DList,"s3",sumF,timeStep,res1,res2);
    printData2D(s11_2DList,"s11",sumF,timeStep,res1,res2);
    printData2D(s22_2DList,"s22",sumF,timeStep,res1,res2);
    printData2D(s33_2DList,"s33",sumF,timeStep,res1,res2);
    printData2D(s12_2DList,"s12",sumF,timeStep,res1,res2);
    printData2D(s13_2DList,"s13",sumF,timeStep,res1,res2);
    printData2D(s23_2DList,"s23",sumF,timeStep,res1,res2);
    printData2D(p_2DList,"p",sumF,timeStep,res1,res2);
    printData2D(q_2DList,"q",sumF,timeStep,res1,res2);

    sClock.get();
    dualOut << "Writing stat_Stress finished" << '\n';
}

void Statistics::writeWeakPlaneStress    (   ConnectLat*                     p_conLat,
                                             PCG_Eigen*                      p_EigenPCG,
                                             std::vector<unsigned>*          p_weakPlaneLatID,
                                             const double                    sumF,
                                             const unsigned                  step)
{
    GLatForce l;
    if (!p_weakPlaneLatID->empty()) {
        char    fileName[255];
        sprintf(fileName,"%s/%sstat_stressConLat-%0d.txt",path,prefix,step);
        ofstream fileSC(fileName, ios::out | ios::app);
        fileSC << "# sumF = " << '\t' << sumF << std::endl;
        fileSC << "# r" << '\t' << "a" << '\t' << "r/a" << '\t' << "Sigma_n/SumF" << "Sigma_s/SumF" << std::endl;
        Point center = p_conLat->getMainCrackCenter();
        std::vector<std::vector<unsigned> > periList = p_conLat->getPerimeter(center,avgD,xyz);
        std::vector<RLat> rLatList = getRLatList(periList,p_weakPlaneLatID,center);
        Point midPoint;
        for (unsigned i=0; i< rLatList.size(); i++) {
            unsigned LatticeID = rLatList[i].LatticeID;
            if (!LatTab[LatticeID].isBreak) {
                Tensor memForce = p_EigenPCG->getMemForce(LatticeID);
                double s_n = memForce[0]/LatTab[LatticeID].area;
                double s_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
                fileSC << rLatList[i].r << '\t' << rLatList[i].a << '\t' << rLatList[i].r/rLatList[i].a << '\t'
                       << s_n/sumF << '\t' << s_s/sumF << std::endl;
            }
        }
    }
}

void Statistics::writeLatStress    (   const double                    sumF,
                                       const unsigned                  timeStep,
                                       ConnectLat*                     p_conLat,
                                       PCG_Eigen*                      p_EigenPCG)
{
    dualOut << "Writing stat_LatStress..." << '\n';
    Clock   sClock;
    sClock.start("stat_LatStress");
    std::vector<double>   axialList,axialForceList;
    std::vector<double>   shearList,shearForceList;
    std::vector<double>   axialRatioList;
    std::vector<double>   shearRatioList;
    std::vector<double>   ratioList;
    std::vector<double>   absRatioList;

    fillInnerLatList(0.02,p_conLat);
    sClock.get();
    axialList.resize(innerLatList.size());
    axialForceList.resize(innerLatList.size());
    shearList.resize(innerLatList.size());
    shearForceList.resize(innerLatList.size());
    ratioList.resize(innerLatList.size());
    axialRatioList.resize(innerLatList.size());
    shearRatioList.resize(innerLatList.size());
    absRatioList.resize(innerLatList.size());
    if (ShearToAxialStiffnessRatio > Tiny) {
        #pragma omp parallel for if (Use_OpenMP) shared (axialList,shearList,ratioList,absRatioList)
        for (unsigned  i = 0 ; i < innerLatList.size() ; i++) {
            unsigned LatticeID = innerLatList[i];
            Tensor memForce = p_EigenPCG -> getMemForce(LatticeID);
            double s_n = memForce[0]/LatTab[LatticeID].area;
            double s_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
            {
                axialList[i] = s_n/sumF;
                axialForceList[i] = memForce[0]/sumF;
                shearList[i] = s_s/sumF;
                double tensileStr = LatTab[LatticeID].et[0]*LatTab[LatticeID].k[0]/LatTab[LatticeID].area;
                double compStr = LatTab[LatticeID].et[1]*LatTab[LatticeID].k[0]/LatTab[LatticeID].area;
                axialRatioList[i] = std::max(s_n/tensileStr,s_n/compStr);
                shearForceList[i] = s_s*LatTab[LatticeID].area/sumF;
                axialForceList[i] = memForce[0]/sumF;
                shearRatioList[i] = s_s/LatTab[LatticeID].shearStr;
                ratioList[i] = s_s/s_n;
            }
        }
        sClock.get();
        printPDF(&axialList,"lat-sn",timeStep,100,tailLen);
        printPDF(&axialForceList,"lat-Fn",timeStep,100,tailLen);
        printPDF(&shearList,"lat-ss",timeStep,100,tailLen);
        printPDF(&shearForceList,"lat-Fs",timeStep,100,tailLen);
        printPDF(&axialRatioList,"lat-snRatio",timeStep,100,tailLen);
        printPDF(&shearRatioList,"lat-ssRatio",timeStep,100,tailLen);
        printPDF(&ratioList,"lat-ss_sn",timeStep,100,tailLen);
    } else {
        #pragma omp parallel for if (Use_OpenMP) shared (axialList)
        for (unsigned  i = 0 ; i < innerLatList.size() ; i++) {
            unsigned LatticeID = innerLatList[i];
            Tensor memForce = p_EigenPCG -> getMemForce(LatticeID);
            double s_n = memForce[0]/LatTab[LatticeID].area;
            axialList[i] = s_n/sumF;
        }
        printPDF(&axialList,"lat-sn",timeStep,100,tailLen);
    }
    sClock.get();
    dualOut << "Writing stat_LatStress finished" << '\n';
}

void Statistics::setCrackMonthLatList   (   const double                width,
                                            std::vector<unsigned>*      p_latList)
{
    dualOut << "Set crack mouth lattice list..." << std::endl;
    for (unsigned i=0; i<p_latList->size(); i++) {
        unsigned LatticeID = (*p_latList)[i];
        unsigned NodeID = LatTab[LatticeID].nb[0];
        if (GNodeTab[NodeID].coord[0]<width) {
            crackMonthLatList[0].push_back(LatticeID);
        } else if (GNodeTab[NodeID].coord[0]>Nx*UnitLength - width) {
            crackMonthLatList[1].push_back(LatticeID);
        }
    }
    dualOut << "crackMonthLatList[0].size() = " << crackMonthLatList[0].size()
            << " crackMonthLatList[1].size() = " << crackMonthLatList[1].size() << std::endl;
}

void Statistics::writeStressStrain		(	const unsigned				faceDir,
											const unsigned				dir,
											const double				sumF,
											const unsigned              totalStep)
{
	static bool isFirstTime = true;
	char	fileName[255];
	cout << "Writing dispRecord..." << '\n';

	sprintf(fileName,"%s/%sdispRec_Face%0d%0d.txt",path,prefix,faceDir,dir);
	ofstream dispFile(fileName, ios::out | ios::app);

	if (isFirstTime) {
		dispFile << "#sumF \t avgDisp/avgDist \t sumF/(avgDisp/avgDist) "
				 <<	"\t avgInDisp/avgInDist \t sumF/(avgInDisp/avgInDist)" << endl;
	}

	double avgDisp = 0.0	, avgDist = 0.0;
	double avgInDisp = 0.0	, avgInDist = 0.0;
	switch (faceDir)
	{
	case 0:
		avgDisp 	= avgFaceDisp(Face3,dir,false);
		avgDist		= avgFaceCoord(Face3,faceDir,false);
		avgInDisp 	= avgFaceDisp(Face3,dir,true);
		avgInDist	= avgFaceCoord(Face3,faceDir,true);

		avgDisp 	-= avgFaceDisp(Face5,dir,false);
		avgDist		-= avgFaceCoord(Face5,faceDir,false);
		avgInDisp 	-= avgFaceDisp(Face5,dir,true);
		avgInDist	-= avgFaceCoord(Face5,faceDir,true);
		break;
	case 1:
		avgDisp 	= avgFaceDisp(Face4,dir,false);
		avgDist		= avgFaceCoord(Face4,faceDir,false);
		avgInDisp 	= avgFaceDisp(Face4,dir,true);
		avgInDist	= avgFaceCoord(Face4,faceDir,true);

		avgDisp 	-= avgFaceDisp(Face2,dir,false);
		avgDist		-= avgFaceCoord(Face2,faceDir,false);
		avgInDisp 	-= avgFaceDisp(Face2,dir,true);
		avgInDist	-= avgFaceCoord(Face2,faceDir,true);
		break;
	case 2:
		avgDisp 	= avgFaceDisp(Face6,dir,false);
		avgDist		= avgFaceCoord(Face6,faceDir,false);
		avgInDisp 	= avgFaceDisp(Face6,dir,true);
		avgInDist	= avgFaceCoord(Face6,faceDir,true);

		avgDisp 	-= avgFaceDisp(Face1,dir,false);
		avgDist		-= avgFaceCoord(Face1,faceDir,false);
		avgInDisp 	-= avgFaceDisp(Face1,dir,true);
		avgInDist	-= avgFaceCoord(Face1,faceDir,true);
		break;
	default:
		dualOut << "Wrong input of faceDir, no data is recorded" << endl;
		return;
	}

	std::cout 	<< sumF      << '\t'                 << avgDisp/avgDist 		<< '\t' << sumF/(avgDisp/avgDist)
						     << '\t' << avgInDisp/avgInDist 	<< '\t' << sumF/(avgInDisp/avgInDist)<< endl;
	dispFile	<< totalStep << '\t' << sumF << '\t' << avgDisp/avgDist 		<< '\t' << sumF/(avgDisp/avgDist)
						     << '\t' << avgInDisp/avgInDist 	<< '\t' << sumF/(avgInDisp/avgInDist)<< endl;
	isFirstTime = false;
}

void Statistics::writeStressStrain      (   const double                sumF,
                                            const unsigned              totalStep)
{
    static bool isFirstTime = true;
    char    fileName[255];
    cout << "Writing dispRecord..." << '\n';

    sprintf(fileName,"%s/%sdispRec_Face.txt",path,prefix);
    ofstream dispFile(fileName, ios::out | ios::app);

    if (isFirstTime) {
        dispFile << "#sumF" << '\t' << "avgDispZ/avgDistZ" << '\t' << "sumF/(avgDispZ/avgDistZ)" << '\t'
                << "avgDispX/avgDistX" << '\t' << "sumF/(avgDispX/avgDistX)" << '\t'
                << "avgDispY/avgDistY" << '\t' << "sumF/(avgDispY/avgDistY)" << '\t'
                << "Poisson ratio" <<  std::endl;
    }

    std::array<double,3> avgDisp,avgDist;
    avgDisp[0]     = avgFaceDisp(Face3,0,false);
    avgDist[0]     = avgFaceCoord(Face3,0,false);

    avgDisp[0]     -= avgFaceDisp(Face5,0,false);
    avgDist[0]     -= avgFaceCoord(Face5,0,false);

    avgDisp[1]     = avgFaceDisp(Face4,1,false);
    avgDist[1]     = avgFaceCoord(Face4,1,false);

    avgDisp[1]     -= avgFaceDisp(Face2,1,false);
    avgDist[1]     -= avgFaceCoord(Face2,1,false);

    avgDisp[2]     = avgFaceDisp(Face6,2,false);
    avgDist[2]     = avgFaceCoord(Face6,2,false);

    avgDisp[2]     -= avgFaceDisp(Face1,2,false);
    avgDist[2]     -= avgFaceCoord(Face1,2,false);


    std::cout   << totalStep << '\t' << sumF << '\t'
                << avgDisp[2]/avgDist[2] << '\t' << sumF/(avgDisp[2]/avgDist[2]) << '\t'
                << avgDisp[0]/avgDist[0] << '\t' << sumF/(avgDisp[0]/avgDist[0]) << '\t'
                << avgDisp[1]/avgDist[1] << '\t' << sumF/(avgDisp[1]/avgDist[1]) << '\t'
                << -(avgDisp[0]/avgDist[0]+avgDisp[1]/avgDist[1])/(2.0*avgDisp[2]/avgDist[2]) << std::endl;
    dispFile    << totalStep << '\t' << sumF << '\t'
                << avgDisp[2]/avgDist[2] << '\t' << sumF/(avgDisp[2]/avgDist[2]) << '\t'
                << avgDisp[0]/avgDist[0] << '\t' << sumF/(avgDisp[0]/avgDist[0]) << '\t'
                << avgDisp[1]/avgDist[1] << '\t' << sumF/(avgDisp[1]/avgDist[1]) << '\t'
                << -(avgDisp[0]/avgDist[0]+avgDisp[1]/avgDist[1])/(2.0*avgDisp[2]/avgDist[2]) << std::endl;
    isFirstTime = false;
}

void Statistics::writeCrackMouthDispHist      (     const double                sumF)
{
    static bool isFirstTime = true;
    char    fileName[255];
    cout << "Writing CrackMouthDispHist..." << '\n';

    sprintf(fileName,"%s/%scrackMouthDispHist.txt",path,prefix);
    ofstream dispFile(fileName, ios::out | ios::app);

    if (isFirstTime) {
        dispFile << "#sumF \t crackMouthDispMax \t crackMonthDispAvg" << endl;
    }
    double avgDisp1 = 0.0;
    GLatForce   l;
    if (crackMonthLatList[0].empty()) {
        tripOut << "[Statistics::writeCrackMouthDispHist] crackMonthLatList[0] is empty!!" << std::endl;
    } else {
        for (unsigned i=0; i<crackMonthLatList[0].size(); i++) {
            avgDisp1 += l.getCrackOpening(crackMonthLatList[0][i]);
        }
        avgDisp1 /= crackMonthLatList[0].size();
    }
    double avgDisp2 = 0.0;
    if (crackMonthLatList[1].empty()) {
       tripOut << "[Statistics::writeCrackMouthDispHist] crackMonthLatList[1] is empty!!" << std::endl;
    } else {
        for (unsigned i=0; i<crackMonthLatList[1].size(); i++) {
            avgDisp2 += l.getCrackOpening(crackMonthLatList[1][i]);
        }
        avgDisp2 /= crackMonthLatList[1].size();
    }
    std::cout   << sumF << '\t' << std::max(avgDisp1,avgDisp2) << '\t'<< (avgDisp1+avgDisp2)/2.0 << std::endl;
    dispFile    << sumF << '\t' << std::max(avgDisp1,avgDisp2) << '\t'<< (avgDisp1+avgDisp2)/2.0 << std::endl;
    isFirstTime = false;
}
void Statistics::writeMaxNorm	(	const double				maxNorm,
									const double				sumF,
									const unsigned				step)
{
	char	fileName[255];
	cout << "Writing MaxNormRecord..." << '\n';
	ofstream* p_dispFile;
	sprintf(fileName,"%s/%sMaxNormRec.txt",path,prefix);
	if (step==0) {
		p_dispFile = new ofstream(fileName, ios::out);
		*p_dispFile << "#step" << '\t' << "maxNorm" << '\t' << "sumF" << '\t' << endl;
	} else {
		p_dispFile = new ofstream (fileName, ios::out | ios::app);
	}
	cout 		<< step << '\t' << maxNorm << '\t' << sumF << endl;
	*p_dispFile << step << '\t' << maxNorm << '\t' << sumF << endl;
	p_dispFile->close();
	delete p_dispFile;
}

void Statistics::writeTotalFluidVol	(	ConnectLat*					p_conLat,
										const unsigned				step)
{
	char	fileName[255];
	cout << "Writing TotalFluidVolRecord..." << '\n';
	ofstream* p_file;
	sprintf(fileName,"%s/%stFluidVolRec.txt",path,prefix);
	if (step==0) {
		p_file = new ofstream(fileName, ios::out);
		*p_file << "#step" << '\t' << "total Fluid Volume injected" << endl;
	} else {
		p_file = new ofstream (fileName, ios::out | ios::app);
	}

	double vol = p_conLat ->getTotalFluidVol();
	cout 		<< step << '\t' << vol  << endl;
	*p_file 	<< step << '\t' << vol << endl;
	p_file->close();
	delete p_file;
}

double Statistics::avgFaceDisp		(	const unsigned			face,
										const unsigned			dir,
										const bool				isInner)
{
	unsigned NodeID;
	double sum = 0.0;
	unsigned size;
	if (isInner) {
		size = inBoundary[face].size();
		for (unsigned i = 0; i<size; i++) {
			NodeID = inBoundary[face][i];
			sum += GNodeTab[NodeID].d[dir];
		}
	} else {
		size = nodeLists.inBoundary[face].size();
		for (unsigned i = 0; i<size; i++) {
			NodeID = nodeLists.inBoundary[face][i];
			sum += GNodeTab[NodeID].d[dir];
		}
	}
	return sum/size;
}

double Statistics::avgFaceCoord		(	const unsigned		face,
										const unsigned		dir,
										const bool			isInner)
{
	unsigned NodeID;
	unsigned size = nodeLists.inBoundary[face].size();
	double sum = 0.0;

	if (isInner) {
		size = inBoundary[face].size();
		for (unsigned i = 0; i<size; i++) {
			NodeID = inBoundary[face][i];
			sum += GNodeTab[NodeID].coord[dir];
		}
	} else {
		size = nodeLists.inBoundary[face].size();
		for (unsigned i = 0; i<size; i++) {
			NodeID = nodeLists.inBoundary[face][i];
			sum += GNodeTab[NodeID].coord[dir];
		}
	}
	return sum/size;
}

void Statistics::writeCoordNo 	(	const unsigned		                step,
                                    const std::vector<unsigned>&        nodeList)
{
    vector <unsigned> nCnt;
    unsigned nMax = 0;

    char	fileName[255];

    time_t 	now 	= time(0);
    char* 	dt 		= ctime(&now);

    cout << "Writing Node Coordination No stats..." << '\n';

    sprintf(fileName,"%s/%sstat_ne-%04d.txt",path,prefix,step);

    ofstream CNoFile(fileName);

    cout << "Statistics :: maxCoordNo = " << maxCoordNo << endl;
    nCnt.resize(maxCoordNo+1);

    unsigned cnt = 0;
    for (unsigned i = 0; i < nodeList.size(); i++)
    {
        nCnt[GNodeTab[nodeList[i]].ne]++;
    }
    cout << "end" << endl;

    for (unsigned i = 0; i<maxCoordNo; i++) {
        if (nCnt[i]!=0) {
            nMax = i;
        }
    }

    CNoFile << "# File - Statistics of Coordination Number distribution" << fileName << '\t' << "Step = " << step << '\t' <<  "created on " << dt << endl;
    CNoFile << "# Coordination NO" << '\t' << "No of Node" <<  '\t' << "percent" << endl;
    for (unsigned i=0; i<nMax; i++)
    {
        CNoFile << i << '\t' << nCnt[i] << '\t' << ((double) nCnt[i])/nodeList.size()*100 << '\n';
    }
}

void Statistics::writeNodeVol        (   const unsigned                      step,
                                         const unsigned                      resolution,
                                         const std::vector<unsigned>&        nodeList)
{
    vector <unsigned> nCnt;
    unsigned nMax = 0;

    char    fileName[255];

    cout << "Writing Node Volume stats..." << '\n';

    double sumLen = 0.0;
    unsigned cntLat = 0;
    std::vector<double>   volList;
    std::vector<double>   aspectRatioList;
    std::vector<double>   aspectRatioList2;
    std::vector<double>   spericityList;
    for (unsigned i = 0; i < nodeList.size(); i++) {
        unsigned NodeID = nodeList[i];
        double vol = GNodeTab[NodeID].v;
        double area = 0.0;
        double len = 0.0;
        for (unsigned n=0; n<GNodeTab[NodeID].nbLatticeID.size(); n++) {
            area += LatTab[GNodeTab[NodeID].nbLatticeID[n]].area;
            sumLen += LatTab[GNodeTab[NodeID].nbLatticeID[n]].length;
            len += LatTab[GNodeTab[NodeID].nbLatticeID[n]].length;
            cntLat++;
        }
        len /= GNodeTab[NodeID].nbLatticeID.size();
        volList.push_back(vol);
        aspectRatioList.push_back(vol/area);
        aspectRatioList2.push_back(vol/(area*len));
        spericityList.push_back(std::pow(6*vol,2.0/3.0)*std::pow(Pi,1.0/3.0)/area);
    }
    sumLen /= cntLat;
    for (unsigned i=0; i<aspectRatioList.size(); i++) {
        aspectRatioList[i] /= sumLen;
    }
    printPDF (&volList,"InVol",step,resolution,tailLen);
    printPDF (&aspectRatioList,"InAspectRatio",step,resolution,tailLen);
    printPDF (&aspectRatioList2,"InAspectRatio2",step,resolution,tailLen);
    printPDF (&spericityList,"InSpericity",step,resolution,tailLen);
}

std::array<double,6> Statistics::getLimits (    const Point     centre)
{
    std::array<double,6> limits;
    double rRatio = 4*avgD;
    if (GeometryID==2) {
        limits[0] = -Huge;
        limits[1] = Huge;
        limits[2] = centre[1]-rRatio/2.0;
        limits[3] = centre[1]+rRatio/2.0;
        limits[4] = -Huge;
        limits[5] = Huge;
    } else {
        limits[0] = -Huge;
        limits[1] = Huge;
        limits[2] = -Huge;
        limits[3] = Huge;
        limits[4] = centre[2]-rRatio/2.0;
        limits[5] = centre[2]+rRatio/2.0;
    }
    return limits;
}

bool Statistics::fillInnerNodeList ()
{
    if (!innerNodeList.empty()) {
        return false;
    }

    for (unsigned NodeID = pxNodeNum; NodeID < tInNodeNum; NodeID++) {
        innerNodeList.push_back(NodeID);
    }
    for (unsigned d=0; d<6; d++) {
        for (unsigned i=0; i<nodeLists.boundary[d].size(); i++) {
            eraseValue(innerNodeList,nodeLists.boundary[d][i]);
        }
    }
    return true;
}

bool Statistics::fillInnerNodeList ( const std::array<double,6> limits)
{
    innerNodeList.clear();
    for (unsigned NodeID = pxNodeNum; NodeID < refineNode[1]; NodeID++) {
        if (isInside(GNodeTab[NodeID].coord,limits)) {
            innerNodeList.push_back(NodeID);
        }
    }
    for (unsigned d=0; d<6; d++) {
        for (unsigned i=0; i<nodeLists.boundary[d].size(); i++) {
            eraseValue(innerNodeList,nodeLists.boundary[d][i]);
        }
    }
    dualOut << "innerNodeList.size() = " << innerNodeList.size() << std::endl;
    return true;
}

bool Statistics::fillInnerLatList (     const double                  minAreaFactor,
                                        const std::array<double,6>    limits,
                                        ConnectLat*                   p_conLat)
{
    innerLatList.clear();
        double minArea = minAreaFactor*avgD*avgD;
        if (Use_UnevenMesh) {
            for (unsigned i=0; i<refineLatList.size(); i++) {
                unsigned LatticeID = refineLatList[i];
                if (LatTab[LatticeID].area>minArea) {
                    if (((LatTab[LatticeID].nb[0] < tInNodeNum)&&(LatTab[LatticeID].nb[0] > pxNodeNum))&&
                            ((LatTab[LatticeID].nb[1] < tInNodeNum)&&(LatTab[LatticeID].nb[1] > pxNodeNum))) {
                        if (isInside(LatticeID,limits)) {
                            innerLatList.push_back(LatticeID);
                        }
                    }
                }
            }
        } else {
            for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
                if (LatTab[LatticeID].area>minArea) {
                    if (((LatTab[LatticeID].nb[0] < tInNodeNum)&&(LatTab[LatticeID].nb[0] > pxNodeNum))&&
                        ((LatTab[LatticeID].nb[1] < tInNodeNum)&&(LatTab[LatticeID].nb[1] > pxNodeNum))) {
                        if (isInside(LatticeID,limits)) {
                            innerLatList.push_back(LatticeID);
                        }
                    }
                }
            }
        }
        rmBoundaryLat(innerLatList);

    std::vector<unsigned> failLat = p_conLat -> getAllGroupVec      ();
    for (unsigned i=0; i< failLat.size(); i++) {
        eraseValue(innerLatList,failLat[i]);
    }
    return true;
}

void Statistics::rmBoundaryLat (    std::vector<unsigned>   &latList)
{
    std::array<bool,6> face = getFaceFixity ();
    for (unsigned d=0; d<6; d++) {
        if (face[d]) {
            for (unsigned i=0; i<nodeLists.boundary[d].size(); i++) {
                unsigned NodeID = nodeLists.boundary[d][i];
                for (unsigned k=0; k<GNodeTab[NodeID].nbLatticeID.size(); k++) {
                    eraseValue(latList,GNodeTab[NodeID].nbLatticeID[k]);
                }
            }
        }
    }
}

bool Statistics::fillInnerLatList (     const double        minAreaFactor,
                                        ConnectLat*         p_conLat)
{
    double minArea = minAreaFactor*avgD*avgD;
    if (innerLatList.empty()) {
        if (Use_UnevenMesh) {
            innerLatList = refineLatList;
        } else {
            for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
                if (LatTab[LatticeID].area>minArea) {
                    if (((LatTab[LatticeID].nb[0] < tInNodeNum)&&(LatTab[LatticeID].nb[0] > pxNodeNum))&&
                        ((LatTab[LatticeID].nb[1] < tInNodeNum)&&(LatTab[LatticeID].nb[1] > pxNodeNum))) {
                        innerLatList.push_back(LatticeID);
                    }
                }
            }
            rmBoundaryLat(innerLatList);
        }
    }

    std::vector<unsigned> failLat = p_conLat -> getAllGroupVec      ();
    for (unsigned i=0; i< failLat.size(); i++) {
        eraseValue(innerLatList,failLat[i]);
    }
    return true;
}

void Statistics::limitInnerLatList (     const double                  minAreaFactor,
                                         const std::array<double,6>    limits)
{
    std::vector<unsigned>   latList = innerLatList;
    for (unsigned i=0; i<latList.size(); i++) {
        if (!isInside(latList[i],limits)){
            eraseValue(innerLatList,latList[i]);
        }
    }
    return;
}

bool Statistics::isInside       (   const unsigned                  LatticeID,
                                    const std::array<double,6>      limits)
{
    Point centroid = LatTab[LatticeID].centroid;
    if ((centroid[0]>limits[0])&&(centroid[0]<limits[1])&&
        (centroid[1]>limits[2])&&(centroid[0]<limits[3])&&
        (centroid[2]>limits[4])&&(centroid[0]<limits[5])) {
        return true;
    }
    return false;
}

bool Statistics::isInside       (   const Point                     centroid,
                                    const std::array<double,6>      limits)
{
    if ((centroid[0]>limits[0])&&(centroid[0]<limits[1])&&
        (centroid[1]>limits[2])&&(centroid[0]<limits[3])&&
        (centroid[2]>limits[4])&&(centroid[0]<limits[5])) {
        return true;
    }
    return false;
}

void Statistics::writeLatLen 	(	const unsigned		step,
                                    const unsigned		resolution)
{
	vector<double>	tmpList;	//temporary store all the values for statistics

	cout << "Writing Lattice Length stats..." << '\n';
	for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
		if (LatTab[LatticeID].isBreak == false) {
			tmpList.push_back(LatTab[LatticeID].length);
		}
	}
	printPDF (&tmpList,"LatLen",step,resolution,tailLen);
}

void Statistics::writeInLatAll 	(	const unsigned				pxNodeNum,
									const unsigned				step,
                                    const unsigned				resolution,
                                    ConnectLat*                 p_conLat)
{
	std::vector<double>	lenList;		//temporary store all length of lattice
	std::vector<double>	areaList;		//temporary store all area denoted by lattice
	std::vector<double>	tensionList;	//temporary store all tension strength of lattice

//	const double ge = 2*NormLatStiffness*ReConRatio;
	fillInnerLatList(0.02,p_conLat);
	std::cout << "Writing Inner Lattice Length stats..." << '\n';
	for (unsigned i=0; i<innerLatList.size(); i++) {
	    unsigned LatticeID = innerLatList[i];
	    lenList.push_back(LatTab[LatticeID].length);
	    areaList.push_back(LatTab[LatticeID].area);
	}
	fillInnerNodeList ();
	writeCoordNo(step,innerNodeList);
	printPDF (&lenList,"InLatLen",step,resolution,tailLen);
	writeNodeVol (step,resolution,innerNodeList);
	printPDF (&areaList,"InLatArea",step,resolution,tailLen);
}

void Statistics::writeLatArea 	(	const unsigned			step,
                                    const unsigned			resolution,
                                    const bool				isLogScale)
{
	std::vector<double>	tmpList;	//temporary store all the values for statistics

	cout << "Writing Lattice Area stats..." << '\n';
	for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
		if (LatTab[LatticeID].isBreak == false)
		{
	    	tmpList.push_back(LatTab[LatticeID].area);
	    }
	}

	printPDF (&tmpList,"Area",step,resolution,tailLen);
}

void Statistics::writeInLatArea 	(	const unsigned			step,
										const unsigned			resolution,
										const bool				isLogScale)
{
	std::vector<double>	tmpList;	//temporary store all the values for statistics

	cout << "Writing Lattice Area stats..." << '\n';
	for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
		if (LatTab[LatticeID].isBreak == false) {
			if ((LatTab[LatticeID].nb[0] < tInNodeNum)&&
				((LatTab[LatticeID].nb[1] < tInNodeNum))) {
				tmpList.push_back(LatTab[LatticeID].area);
			}
	    }
	}
	printPDF (&tmpList,"Area",step,resolution,tailLen);
}

void Statistics::writeLatBreakDisp 	(	const unsigned			step,
                                    	const unsigned			resolution,
                                    	const bool				isLogScale)
{
    vector<double>	tmpList;	//temporary store all the values for statistics

    std::cout << "Writing lattice breaking disp stats..." << '\n';
    for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
    	if (LatTab[LatticeID].isBreak == false)
    	{
    		tmpList.push_back(LatTab[LatticeID].et[0]);
    	}
    }
    printPDF (&tmpList,"Ft0",step,resolution,tailLen);
}

void Statistics::print2 				(	const vector<double>*	p_vList,
											const char*				vName,
											const unsigned			step,
                                    		const unsigned			resolution,
                                    		const bool				isLogScale)
{
	char	fileName[255];
	double	Max = -Huge;
	double	Min = Huge;

	vector<unsigned>	Cnt (resolution,0);

	sprintf(fileName,"%s/%sstat_%s-%04d.txt",path,prefix,vName,step);
	ofstream inFile(fileName); //, ios::out | ios::app);

	time_t 	now 	= time(0);
	char* 	dt 		= ctime(&now);

	inFile 	<< "# File - Statistics of " << vName << " distribution " << endl;
	inFile 	<< "# Saved in " << fileName << '\t' << "Step = " << step << '\t' <<  "created on " << dt << endl;
	inFile	<< "# " << vName << '\t' << "Frequency" <<  '\t' << "percent" << endl;

	for (unsigned i = 0; i < p_vList->size(); i++) {
		if ((*p_vList)[i]>Max) {
			Max = (*p_vList)[i];
		}
		if ((*p_vList)[i]<Min) {
			Min = (*p_vList)[i];
		}
	}

	if (isLogScale) {
	    if (Min <= 0) {
	    	dualOut << "Min value is smaller or equal to zero, a small positive value Tiny is assigned as Min value for logarithm output" << endl;
	        Min = Tiny;
	    }

	    if (Max <=0) {
	        dualOut << "All the values are non-positive, Stat output -Statistics::writeLatArea- abort!!" << endl;
	        return;
	    }
	    double dl = log(pow((Max/Min),1./((double) resolution)));

	    for (unsigned i = 0; i<p_vList->size(); i++)	{
	    	Cnt[log((*p_vList)[i]/Min)/dl]++;
	    }

	    dl = exp(dl);
	    double dlAcc = 1.0;

	    for (unsigned i=0; i<resolution; i++) {
	    	inFile << (Min*dlAcc) << '\t' << Cnt[i] << '\t' << ((double) Cnt[i])/p_vList->size()*100 << '\n';
	        dlAcc *=dl;
	    }
	} else {
	    double dl = (Max - Min)/resolution;
	    for (unsigned i = 0; i < p_vList->size(); i++) {
	    	Cnt[((*p_vList)[i]-Min)/dl]++;
	    }
	    for (unsigned i=0; i<resolution; i++) {
	        inFile << (Min+i*dl)<< '\t' << Cnt[i] << '\t' << ((double) Cnt[i])/p_vList->size()*100 << '\n';
	    }
	}

	std::array<double,4> statPara = getStatInfo(*p_vList);
	inFile  << "# " << "Mean = " << '\t' << statPara[0] << '\t' << " SD = " <<  '\t' << statPara[1] << '\t'
	        << " Skewness = " <<  '\t' << statPara[1] << '\t' << " Kurtosis = " <<  '\t' << statPara[3] << std::endl;
}

void Statistics::printPDF 				(	std::vector<double>*	    p_vList,
											const char*					vName,
											const unsigned				step,
                                    		const unsigned				resolution,
                                    		const double                tailSize)
{
    if (p_vList->size()<100) {
        tripOut << "The population size of " << vName << " is too small, N = " << p_vList->size()
                << " the statistics is not representative, no stats file is generated!" << std::endl;
        return;
    }

    unsigned res = resolution;
    if (p_vList->size()/res<20) {
        res = p_vList->size()/20;
    }
    char	fileName[255];
	double	Max = -Huge;
	double	Min = Huge;

	std::vector<unsigned>	Cnt (res,0);


	std::sprintf(fileName,"%s/%sstat_%s-%04d.txt",path,prefix,vName,step);
	std::ofstream inFile(fileName); //, ios::out | ios::app);

	time_t 	now 	= time(0);
	char* 	dt 		= ctime(&now);

	inFile 	<< "# File - Statistics of " << vName << " distribution " << endl;
	inFile 	<< "# Saved in " << fileName << '\t' << "Step = " << step << '\t' <<  "created on " << dt << endl;
	inFile	<< "# " << vName << '\t' << "Frequency" <<  '\t' << "percent" << '\t' << "PDF" << '\t' << "CDF"  << endl;

	std::sort(p_vList->begin(),p_vList->end());
	unsigned N = p_vList->size();
	unsigned tail = std::max(N / res * tailSize,1.01);

	unsigned maxPos = N - tail;
	if (tail>=N) {
	    maxPos = 0;
	}

	unsigned minPos = tail;
	Max = (*p_vList)[maxPos];
	Min = (*p_vList)[minPos];
	double dl = (Max - Min)/res;
	for (unsigned i = 0; i < p_vList->size(); i++) {
	    unsigned pos = ((*p_vList)[i]-Min)/dl;
	    if (pos >=res) {
	        pos = res - 1;
	    }
	    Cnt[pos]++;
	}
	double integral = dl*p_vList->size();
	double cdf = 0.0;
	double pdf;
	for (unsigned i=0; i<res; i++) {
		pdf = ((double) Cnt[i])/integral;
		cdf += pdf*dl;
	    inFile << (Min+(i+0.5)*dl)<< '\t' << Cnt[i] << '\t' << ((double) Cnt[i])/p_vList->size()*100 << '\t' << pdf << '\t' << cdf << '\n';
	}

	std::array<double,4> statPara = getStatInfo(*p_vList);
	inFile  << "# " << "Mean = " << '\t' << statPara[0] << '\t' << " SD = " <<  '\t' << statPara[1] << std::endl;
	inFile  << "# Skewness = " <<  '\t' << statPara[2] << '\t' << " Kurtosis = " <<  '\t' << statPara[3] << std::endl;
}

void Statistics::printData2D             (  std::vector<std::pair<double,std::pair<double,double> > >   &pList,
                                            const char*                                                 vName,
                                            const double                                                sumF,
                                            const unsigned                                              step,
                                            const unsigned                                              resolution1,
                                            const unsigned                                              resolution2)
{
    if (pList.size()<500) {
        tripOut << "The population size of " << vName << " is too small, N = " << pList.size()
                << " the statistics is not representative, no stats file is generated!" << std::endl;
        return;
    }

    unsigned res1 = resolution1;
    unsigned res2 = resolution2;
    if (pList.size()/res1/res2<5) {
        double r = std::max(std::sqrt(pList.size()/res1/res2),2.0);
        res1 /= r;
        res2 /= r;
        dualOut << "res1 becomes = " << res1 << " , res2 becomes = " << res2 << std::endl;
    }
    char    fileName[255];

    std::vector<std::vector<unsigned> >   Cnt (res1);
    std::vector<std::vector<double> >   Data (res1);
    for (unsigned i=0; i<res1; i++) {
        Cnt[i].resize(res2,0);
        Data[i].resize(res2,0.0);
    }


    std::sprintf(fileName,"%s/%sstat2D_%s-%04d.txt",path,prefix,vName,step);
    std::ofstream inFile(fileName); //, ios::out | ios::app);

    time_t  now     = time(0);
    char*   dt      = ctime(&now);

    inFile  << "# File - Statistics of " << vName << " distribution " << endl;
    inFile  << "# Saved in " << fileName << '\t' << "Step = " << step << '\t' <<  "created on " << dt << endl;
    inFile  << "# SumF = " << '\t' << sumF << std::endl;
    inFile  << "# " << "r/a" <<  '\t' << "dz/AvgD" << '\t' << vName << '\t'  << endl;

    std::vector<double>     p1List,p2List;
    for (unsigned i=0; i<pList.size(); i++) {
 //       std::cout << pList[i].first << ' ' << pList[i].second.first << ' ' <<  pList[i].second.second << std::endl;
        p1List.push_back(pList[i].second.first);
        p2List.push_back(pList[i].second.second);
    }
    std::sort(p1List.begin(),p1List.end());
    std::sort(p2List.begin(),p2List.end());
    unsigned N = pList.size();

    unsigned minPos = 0, maxPos = N-1;
    double Max1 = p1List[maxPos];
    double Min1 = p1List[minPos];
    double Max2 = p2List[maxPos];
    double Min2 = p2List[minPos];
    double dl1 = (Max1 - Min1)/res1;
    double dl2 = (Max2 - Min2)/res2;
    for (unsigned i = 0; i < pList.size(); i++) {
        unsigned pos1 = (pList[i].second.first-Min1)/dl1;
        unsigned pos2 = (pList[i].second.second-Min2)/dl2;
        if (pos1 >=res1) {
            pos1 = res1 - 1;
        }
        if (pos2 >=res2) {
            pos2 = res2 - 1;
        }
        Cnt[pos1][pos2]++;
        Data[pos1][pos2]+=pList[i].first;
    }
    inFile << "0" << '\t';
    for (unsigned j=0; j<res2; j++) {
        inFile << (Min2+(j+0.5)*dl2) << '\t';
    }
    inFile << '\n';
    for (unsigned i=0; i<res1; i++) {
        inFile << Min1 + (i+0.5)*dl1 << '\t';
        for (unsigned j=0; j<res2; j++) {
            if (Cnt[i][j]>0) {
                inFile << Data[i][j]/((double) Cnt[i][j]) << '\t';
            } else {
                double data = 0.0;
                unsigned cnt = 0;
                unsigned n=0;
                do {
                    n++;
                    unsigned iStart,jStart,iEnd,jEnd;
                    if (i<=n) {
                        iStart = 0;
                    } else {
                        iStart = i - n;
                    }
                    if (i+n>res1) {
                        iEnd = res1;
                    } else {
                        iEnd = i + n;
                    }
                    if (j<=n) {
                        jStart = 0;
                    } else {
                        jStart = j - n;
                    }
                    if (j+n>res2) {
                        jEnd = res2;
                    } else {
                        jEnd = j + n;
                    }
                    for (unsigned ii=iStart; ii<iEnd; ii++) {
                        for (unsigned jj=jStart; jj<jEnd; jj++) {
                            if (Cnt[ii][jj]>0) {
                                data += Data[ii][jj];
                                cnt += Cnt[ii][jj];
                            }
                        }
                    }
                } while (cnt ==0);
                inFile << data/((double) cnt) << '\t';
            }
        }
        inFile << '\n';
    }
}

void Statistics::printCluster           (   std::vector<unsigned>*      p_vList,
                                            const char*                 vName,
                                            const unsigned              step,
                                            const unsigned              resolution)
{
    unsigned res = resolution;
    if (p_vList->size()<resolution/4) {
        tripOut << "The population size of " << vName << " is too small, N = " << p_vList->size()
                << " and resolution = " << resolution << " .No output is done" << std::endl;
        return;
    }
    if (p_vList->size()<resolution) {
        tripOut << "The population size of " << vName << " is too small, N = " << p_vList->size()
                << " and resolution = " << resolution << " resolution is set to N" << std::endl;
        res = p_vList->size();
    }

    char    fileName[255];
    std::vector<std::pair<unsigned,unsigned> >   Cnt;
    std::sprintf(fileName,"%s/%sstat_%s-%04d.txt",path,prefix,vName,step);
    std::ofstream inFile(fileName); //, ios::out | ios::app);

    time_t  now     = time(0);
    char*   dt      = ctime(&now);

    inFile  << "# File - Statistics of " << vName << " distribution " << endl;
    inFile  << "# Saved in " << fileName << '\t' << "Step = " << step << '\t' <<  "created on " << dt << endl;
    inFile  << "# " << vName << '\t' << "Frequency"  << '\t' << "PDF" << '\t' << "CDF"  << endl;

    std::sort(p_vList->begin(),p_vList->end());
    std::vector<unsigned>   vList2;
    for (int i=p_vList->size()-1; i>=0; i--) {
        vList2.push_back((*p_vList)[i]);
    }
    unsigned N = vList2.size();
    std::pair<unsigned,unsigned> p;
    p.first = vList2[0];
    p.second = 1;
    Cnt.push_back(p);
    unsigned l=0;
    for (unsigned i=1; i<N;i++) {
        if (Cnt[l].first==vList2[i]) {
            Cnt[l].second++;
        } else {
            l++;
            if (l>=res) {
                dualOut << "The container is full, some items are discarded!" << std::endl;
                break;
            } else {
                p.first=vList2[i];
                p.second = 1;
                Cnt.push_back(p);
            }
        }
    }
    if (l<res) {
        dualOut << "All the items are taken into account but the resoulion is smaller than specificed,"
                << "l= " << l << " res = " << res << std::endl;
    }
    for (unsigned i=0; i<Cnt.size(); i++) {
        inFile << Cnt[i].first << '\t' << Cnt[i].second << '\t' << ((double) Cnt[i].second)/((double) N) <<  '\n';
    }
    inFile << "Total failLat Num = " << '\t' << N << std::endl;
}



void Statistics::printPDF               (   std::vector<std::vector<double> >                     &vList,
                                            const char*                                           vName,
                                            const unsigned                                        step,
                                            const unsigned                                        resolution,
                                            const double                                          tailSize)
{
    if (vList[0].size()<100) {
        tripOut << "The population size of " << vName << " is too small, N = " << vList[0].size()
                << " the statistics is not representative, no stats file is generated!" << std::endl;
        return;
    }

    unsigned chunkNum = vList.size();
    unsigned res = resolution;
    if (vList[0].size()/res<20) {
        res = vList[0].size()/20;
    }
    char    fileName[255];

    std::sprintf(fileName,"%s/%sstat_%s-%04d.txt",path,prefix,vName,step);
    std::ofstream inFile(fileName); //, ios::out | ios::app);

    time_t  now     = time(0);
    char*   dt      = ctime(&now);

    inFile  << "# File - Statistics of " << vName << " distribution " << endl;
    inFile  << "# Saved in " << fileName << '\t' << "Step = " << step << '\t' << " Num of chunk = " << vList.size() << '\t' << "created on " << dt << endl;
    inFile  << "# ";
    for (unsigned t=0; t<chunkNum; t++) {
        inFile  << vName << '\t' << "Frequency" << '\t' << "PDF" << '\t' << "CDF"  << '\t';
    }
    inFile << std::endl;

    std::vector<unsigned> maxPos;
    std::vector<unsigned> minPos;
    std::vector<double>  Max;
    std::vector<double>  Min;
    std::vector<double>  dl;
    std::vector<std::vector<unsigned> >  Cnt;
    std::vector<unsigned> N;

    maxPos.resize(chunkNum);
    minPos.resize(chunkNum);
    Max.resize(chunkNum);
    Min.resize(chunkNum);
    dl.resize(chunkNum);
    Cnt.resize(chunkNum);
    N.resize(chunkNum);
    unsigned totN = 0;
    std::vector<double> integral;

    for (unsigned t = 0; t < chunkNum; t++) {
        std::sort(vList[t].begin(),vList[t].end());
        N[t] = vList[t].size();
        totN += N[t];
        unsigned tail = N[t] / res * tailSize;
        maxPos[t] = N[t] - tail;
        minPos[t] = tail;
        Max[t] = vList[t][maxPos[t]];
        Min[t] = vList[t][minPos[t]];
        dl[t] = (Max[t] - Min[t])/res;
        Cnt[t].resize(res,0);
        integral.push_back(dl[t]*N[t]);
    }


    for (unsigned t = 0; t < chunkNum; t++) {
        for (unsigned i = 0; i < vList[t].size(); i++) {
            unsigned pos = (vList[t][i]-Min[t])/dl[t];
            if (pos >=res) {
                pos = res - 1;
            }
            Cnt[t][pos]++;
        }
    }

    std::vector<double> cdf (chunkNum,0.0);
    double pdf;
    for (unsigned i=0; i<res; i++) {
        for (unsigned t=0; t<chunkNum; t++) {
            pdf = ((double) Cnt[t][i])/integral[t];
            cdf[t] += pdf*dl[t];
            inFile << (Min[t]+(i+0.5)*dl[t])<< '\t' << Cnt[t][i] << '\t' << pdf*N[t]/totN << '\t' << cdf[t] << '\t';
        }
        inFile << '\n';
    }
    for (unsigned t=0; t<chunkNum; t++) {
        std::array<double,4> statPara = getStatInfo(vList[t]);
        inFile  << "# " << "Chuck = " << t << '\t' << "N[t]/totN = " << N[t]/totN << " Mean = " << '\t' << statPara[0]
                << '\t' << " SD = " <<  '\t' << statPara[1] << '\t'
                << "Skewness = " << '\t' << statPara[2] << '\t' << " Kurtosis = " <<  '\t' << statPara[3] << std::endl;
    }
    inFile << "totN = " << totN << std::endl;
}

//calculate mean and standard variation for a given list of double p_vList
void Statistics::getStatInfo 			(	const std::vector<double>*	p_vList,
											double*						p_mean,
											double*						p_SD)
{
	double sum = 0.0, sqSum = 0.0;
	for (unsigned i=0; i<p_vList->size(); i++) {
		sum += (*p_vList)[i];
		sqSum += (*p_vList)[i]*(*p_vList)[i];
	}
	double size = p_vList->size();
	*p_mean = sum / size;

	*p_SD = std::sqrt((sqSum-sum*sum/size)/size);
}

std::array<double,4> Statistics::getStatInfo            (   const std::vector<double>  &vList)
{
    std::array<double,4> statsPara;
    double sum = 0.0, sqSum = 0.0;
    for (unsigned i=0; i<vList.size(); i++) {
        sum += vList[i];
        sqSum += vList[i]*vList[i];
    }
    double size = vList.size();
    statsPara[0] = sum / size;
    statsPara[1] = std::sqrt((sqSum-sum*sum/size)/size);
    double cubicSum = 0.0, p4Sum = 0.0;
    for (unsigned i=0; i<vList.size(); i++) {
        double diff = vList[i] - statsPara[0];
        cubicSum += diff*diff*diff;
        p4Sum += diff*diff*diff*diff;
    }
    statsPara[2] = cubicSum/(size*statsPara[1]*statsPara[1]*statsPara[1]);
    statsPara[3] = p4Sum/(size*statsPara[1]*statsPara[1]*statsPara[1]*statsPara[1])-3.0;
    return statsPara;
}


