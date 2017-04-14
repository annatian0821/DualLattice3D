/*
 * InitCond.cpp
 *
 *  Created on: Oct 12, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "InitCond.hpp"

using namespace std;

unsigned				initRmNode;									//Number of node removed to form initial crack
unsigned				initRmLat;									//Number of lattice removed to form initial crack
unsigned                PriInjLattice;

InitCond::InitCond			(	GLatRemoval*			_p_glatRemoval)
{
	p_glatRemoval = _p_glatRemoval;
}

// For UCS
void 	InitCond::setInitCond		(	ConnectLat*						p_conLat,
                                        Vertices*						p_verti,
                                        PCG_Eigen*						p_EigenPCG,
                                        std::vector<unsigned>*			p_pxCrackLatList,
                                        const std::array<unsigned,2>	pxCrackPerimeterNum,
                                        const unsigned					maxCoordNo)
{
	cout << "p_pxCrackLatList->size() = " << p_pxCrackLatList ->size() << std::endl;
	initNodeAndLatState();
    applyinitialCondition(InitStress,p_verti,p_conLat,p_EigenPCG,maxCoordNo);
    applyBoundaryRestrain(BoundaryMat);
    for (std::vector<unsigned>::iterator it=p_pxCrackLatList->begin(); it!=p_pxCrackLatList->end(); ++it) {
    	unsigned LatticeID = *it;
    	if (!LatTab[LatticeID].isBreak) {
    		p_glatRemoval->removeLattice (LatticeID,p_conLat,p_verti,p_EigenPCG,0,0.0);
    		LatTab[LatticeID].k[0] = Tiny;
    		if (DofPerNode==6) {
    		    LatTab[LatticeID].k[1] = LatTab[LatticeID].k[2] = LatTab[LatticeID].k[3] = Tiny;
    		}
    	}
    	LatTab[LatticeID].et[0] = Tiny;
    }
    setZero();
}

// For hydraulic fracturing
void 	InitCond::setInitCond				(	ConnectLat*						p_conLat,
                                                GLatRemoval*                    p_glatRemoval,
        										Vertices*						p_verti,
        										FluidCal*						p_fCal,
        										PCG_Eigen*						p_EigenPCG,
        										std::vector<unsigned>*			p_pxCrackLatList,
        										const std::array<unsigned,2>	pxCrackPerimeterNum,
        										const double					initP,
        										const unsigned					maxCoordNo)
{
	initNodeAndLatState();
	if (false) {
        genPriInjectWell(p_pxCrackLatList,InjectRate,initP,p_conLat,p_glatRemoval,p_verti,p_fCal);
        rmPxCrack (p_conLat,p_verti,p_EigenPCG,p_pxCrackLatList);
        GLatRemoval     glr(1e-4);
        glr.disConnectLat(p_verti,p_conLat,p_EigenPCG,0,0.0,p_pxCrackLatList);
    }
	applyinitialCondition(InitStress,p_verti,p_conLat,p_EigenPCG,maxCoordNo);
	applyBoundaryRestrain(BoundaryMat);
	if (true) {
	    genPriInjectWell(p_pxCrackLatList,InjectRate,initP,p_conLat,p_glatRemoval,p_verti,p_fCal);
	    rmPxCrack (p_conLat,p_verti,p_EigenPCG,p_pxCrackLatList);
	    GLatRemoval     glr(1e-4);
	    glr.disConnectLat(p_verti,p_conLat,p_EigenPCG,0,0.0,p_pxCrackLatList);
	}
	setZero();
	p_EigenPCG->updateDOF (true);
}

void    InitCond::setInitCrackAperture (  std::vector<unsigned>*           p_pxCrackLatList)
{
    for (unsigned LatticeID = 0; LatticeID < LatTab.size(); LatticeID++) {
        LatTab[LatticeID].initOpening = 0.0;
    }
    GLatForce lforce;
    for (std::vector<unsigned>::iterator it=p_pxCrackLatList->begin(); it!=p_pxCrackLatList->end(); ++it) {
        LatTab[*it].initOpening=lforce.getCrackOpening(*it);
    }
}

void    InitCond::rmPxCrack (   ConnectLat*                    p_conLat,
                                Vertices*                      p_verti,
                                PCG_Eigen*                     p_EigenPCG,
                                std::vector<unsigned>*         p_pxCrackLatList)
{
    std::vector<unsigned>   rmLatID = *p_pxCrackLatList;
    unsigned LatticeID;
    for (unsigned i=0; i<rmLatID.size(); i++) {
        LatticeID = rmLatID[i];
        if (!LatTab[LatticeID].isBreak) {
            p_glatRemoval->removeLattice (LatticeID,p_conLat,p_verti,p_EigenPCG,0,0.0);
        }
        LatTab[LatticeID].et[0] = Tiny;
    }
}

void 	InitCond::initNodeAndLatState	()
{
	bool		isFind;
	//To sort neighbour node according to NodeID
	for (unsigned NodeID = 0; NodeID < GNodeTab.size(); NodeID++) {
		sortPair(&GNodeTab[NodeID].nbNodeID,&GNodeTab[NodeID].nbLatticeID);
	    isFind = false;
	    for (unsigned i=0; i<GNodeTab[NodeID].nbNodeID.size(); i++) {
	    	if (GNodeTab[NodeID].nbNodeID[i]>NodeID) {
	        	GNodeTab[NodeID].start = i;
	            isFind = true;
	            break;
	        }
	    }
	    if (!isFind) {
	    	GNodeTab[NodeID].start = GNodeTab[NodeID].nbNodeID.size();
	    }
	}
}

void InitCond::applyinitialCondition	(	const Tensor6				initStress,
											Vertices*					p_verti,
											ConnectLat*					p_conLat,
											PCG_Eigen*					p_EigenPCG,
											const unsigned				maxCoordNo) {
	bool isInSituStress = (std::fabs(initStress[0])>Tiny)||(std::fabs(initStress[1])>Tiny)||(std::fabs(initStress[2])>Tiny)
							||(std::fabs(initStress[3])>Tiny)||(std::fabs(initStress[4])>Tiny)||(std::fabs(initStress[5])>Tiny);
	if (!isInSituStress) {
	    dualOut << "[InitCond::applyinitialCondition]: No initial stress to be applied \n";
		return;
	}
	if (true) {
	    applyBoundaryStress (initStress);
	    RestrainMat restrainMat = getRestrainMat(initStress);
	    applyBoundaryRestrain (restrainMat);
	} else {
	    applyBoundaryStressFree (initStress);
	    applyMinRestrain ();
	}
	setBoreHoleConstrain ();
	stablizeInitModel (p_verti,p_conLat,p_EigenPCG,maxCoordNo,5*MaxPerStep_BrokenLattice,50);
}


void InitCond::applyBoundaryStress		( 	const Tensor6				initStress)
{
	Loading load;
	//For 3 axial stresses, sigma_xx,sigma_yy,sigma_zz
	unsigned face[3] = { Face3, Face4, Face6 };
	for (unsigned i=0; i<3; i++) {
		if (fabs(initStress[i])>Tiny) {
//			load.applyGUniPressureOnBoundary(-initStress[i],face[i],p_Di,i);
			load.applyGUniPressureOnBoundary(-initStress[i],face[i],i);
		}
	}

	//For 3 shear stresses, sigma_xy,sigma_xz,sigma_yz
	unsigned face2[3][3] = 	{{ Face4, Face5, Face3 },
							 { Face3, Face5, Face6 },
							 { Face4, Face2, Face6 }};

	unsigned dir[3][3] 	= 	{{ 0, 1, 1 },
							 { 2, 2, 0 },
							 { 2, 2, 1 }};

	int		sign[3]		= 	{ -1, 1, -1};

	for (unsigned i=0; i<3; i++) {
		if (fabs(initStress[i+3])>Tiny) {
			for (unsigned j=0; j<3; j++) {
				load.applyGUniPressureOnBoundary(sign[j]*initStress[i+3],face2[i][j],dir[i][j]);
			}
		}
	}
}

void InitCond::applyBoundaryStressFree		( 	const Tensor6				initStress)
{
	Loading load;
	//For 3 axial stresses, sigma_xx,sigma_yy,sigma_zz
	unsigned face[3] 	= { Face3, Face4, Face6 };
	unsigned oppFace[3]	= { Face5, Face2, Face1 };
	for (unsigned i=0; i<3; i++) {
		if (fabs(initStress[i])>Tiny) {
			load.applyGUniPressureOnBoundary(-initStress[i],face[i],i);
			load.applyGUniPressureOnBoundary(initStress[i],oppFace[i],i);
		}
	}

	//For 3 shear stresses, sigma_xy,sigma_xz,sigma_yz
	unsigned face2[3][4] = 	{{ Face2, Face4, Face5, Face3 },
							 { Face1, Face6, Face5, Face3 },
							 { Face1, Face6, Face2, Face4 }};

	unsigned dir[3][4] 	= 	{{ 0, 0, 1, 1 },
							 { 0, 0, 2, 2 },
							 { 1, 1, 2, 2 }};

	int		sign[4]		= 	{ 1, -1, 1, -1};

	for (unsigned i=0; i<3; i++) {
		if (fabs(initStress[i+3])>Tiny) {
			for (unsigned j=0; j<4; j++) {
				load.applyGUniPressureOnBoundary(sign[j]*initStress[i+3],face2[i][j],dir[i][j]);
			}
		}
	}
}

RestrainMat  InitCond::getRestrainMat				( 	const Tensor6	initStress) {
	unsigned face[6] 	= { Face5, Face2, Face1, Face2, Face1, Face1 };
	unsigned dir[6]		= { 0 , 1 , 2 , 0 , 0 , 1 };
	RestrainMat outMat;
	for (unsigned i=0; i<6; i++) {
		for (unsigned j=0; j<DofPerNode; j++) {
			outMat[i][j] = false;
		}
	}
	for (unsigned i=0; i<6; i++) {
		if (fabs(initStress[i])>Tiny) {
			outMat[face[i]][dir[i]] = true;
		}
	}
	return outMat;
}

void InitCond::applyBoundaryRestrain	( 	const RestrainMat restrainMat) {
	dualOut << "Applying boundary restrain..." << std::endl;
	nodeLists.restrain.clear();
	nodeLists.free.clear();
	for (unsigned NodeID=0; NodeID < GNodeTab.size(); NodeID++) {
	    nodeLists.free.push_back(NodeID);
	}
	unsigned NodeID;
	for (unsigned i=0; i<6; i++) {
		if (isConstrain(restrainMat[i])) {
			for (unsigned j=0; j<nodeLists.boundary[i].size(); j++) {
				NodeID = nodeLists.boundary[i][j];
				std::list<Restrain>::iterator it = std::find(nodeLists.restrain.begin(),nodeLists.restrain.end(),NodeID);
				if (it==nodeLists.restrain.end()) {
					Restrain	newRe(NodeID,restrainMat[i]);
					nodeLists.restrain.push_back(newRe);
				} else {
				    it->addRestrain(restrainMat[i]);
				}
				eraseValue(nodeLists.free,NodeID);
			}
		}
	}
	dualOut << "Applying boundary restrain finished" << std::endl;
}

void InitCond::applyMinRestrain ()
{
    dualOut << "InitCond::applyMinRestrain..." << std::endl;
    std::vector<Point>      centers;
    Point center;
    center[0] = 0.5*Nx*UnitLength; center[1] = 0.5*Ny*UnitLength; center[2] = 0.0;
    centers.push_back(center);
    center[0] = 0.5*Nx*UnitLength; center[1] = 0.0; center[2] = 0.5*Nz*UnitLength;
    centers.push_back(center);
    center[0] = Nx*UnitLength; center[1] = 0.5*Ny*UnitLength; center[2] = 0.5*Nz*UnitLength;
    centers.push_back(center);
    center[0] = 0.5*Nx*UnitLength; center[1] = Ny*UnitLength; center[2] = 0.5*Nz*UnitLength;
    centers.push_back(center);
    center[0] = 0.0; center[1] = 0.5*Ny*UnitLength; center[2] = 0.5*Nz*UnitLength;
    centers.push_back(center);
    center[0] = 0.5*Nx*UnitLength; center[1] = 0.5*Ny*UnitLength; center[2] = Nz*UnitLength;
    centers.push_back(center);

    RestrainMat restrainMat;
    for (unsigned d=0; d<6; d++) {
        for (unsigned l=0; l<DofPerNode; l++) {
            restrainMat[d][l] = false;
        }
    }
    restrainMat[0][0] = true;  restrainMat[0][1] = true;
    restrainMat[1][0] = true;  restrainMat[1][2] = true;
    restrainMat[2][1] = true;  restrainMat[2][2] = true;
    restrainMat[3][0] = true;  restrainMat[3][2] = true;
    restrainMat[4][1] = true;  restrainMat[4][2] = true;
    restrainMat[5][0] = true;  restrainMat[5][1] = true;
    double minDist = Huge;
    Geometry geo;
    for (unsigned d=0; d<6; d++) {
        unsigned NodeID;
        unsigned restrainNodeID;
        for (unsigned i=0; i<nodeLists.boundary[d].size(); i++) {
            NodeID = nodeLists.boundary[d][i];
            double dist = geo.dist(GNodeTab[NodeID].coord,centers[d]);
            if (dist<minDist) {
                minDist = dist;
                restrainNodeID = NodeID;
            }
        }
        std::list<Restrain>::iterator it = std::find(nodeLists.restrain.begin(),nodeLists.restrain.end(),restrainNodeID);
        if (it==nodeLists.restrain.end()) {
            Restrain    newRe(restrainNodeID,restrainMat[d]);
            nodeLists.restrain.push_back(newRe);
        } else {
            it->addRestrain(restrainMat[d]);
        }
        eraseValue(nodeLists.free,restrainNodeID);
        std::cout << GNodeTab[restrainNodeID].coord[0] << ' ' << GNodeTab[restrainNodeID].coord[1]
                  << ' ' << GNodeTab[restrainNodeID].coord[2] << std::endl;
        minDist = Huge;
    }
}

void InitCond::outputBoundaryRestrain    () {
    dualOut << "Output boundary restrain..." << std::endl;
    char fileName[255]; //filename
    std::sprintf(fileName,"%s/%sRestrain.txt",OutputFolder,OutputFilePrefix);
    std::ofstream file(fileName, std::ios::out);
    for (std::list<Restrain>::iterator it= nodeLists.restrain.begin(); it!=nodeLists.restrain.end(); ++it) {
        unsigned NodeID = it->NodeID;
        file << GNodeTab[NodeID].coord[0] << ' ' << GNodeTab[NodeID].coord[1] << ' ' << GNodeTab[NodeID].coord[2] << ':'
                << it->Dir[0] << ' ' << it->Dir[1] << ' ' << it->Dir[2] << ' ' << it->Dir[3] << ' ' << it->Dir[4] << ' ' << it->Dir[5] << '-'
                << GNodeTab[NodeID].d[0] << ' ' << GNodeTab[NodeID].d[1] << ' ' << GNodeTab[NodeID].d[2] << ' ' <<
                GNodeTab[NodeID].d[3] << ' ' << GNodeTab[NodeID].d[4] << ' ' << GNodeTab[NodeID].d[5] << ' ' << std::endl;
    }
    dualOut << "Applying boundary restrain finished" << std::endl;
}

unsigned InitCond::genPriInjectWell	(	const std::vector<unsigned>*	p_pxCrackLatList,
                                        const double					injectRate,
                                        const double					pressure,
                                        ConnectLat*						p_conLat,
                                        GLatRemoval*                    p_glatRemoval,
                                        Vertices*						p_verti,
                                        FluidCal*						p_fCal)
{
	if (p_pxCrackLatList->empty()) {
		tripOut << "[InitCond::genPriInjectWell] the input pxCrackLatList is empty, nothing is done" << std::endl;
		return LatTab.size();
	}
	Point centre;
	centre[0] = 0.0; centre[1] = 0.0; centre[2] = 0.0;
	for (unsigned i=0; i<p_pxCrackLatList->size(); i++) {
		for (unsigned d=0; d<Dim; d++) {
			centre[d] += GNodeTab[LatTab[(*p_pxCrackLatList)[i]].nb[0]].coord[d]
			             +GNodeTab[LatTab[(*p_pxCrackLatList)[i]].nb[1]].coord[d];
		}
	}
	centre[0] /=2*p_pxCrackLatList->size();
	centre[1] /=2*p_pxCrackLatList->size();
	centre[2] /=2*p_pxCrackLatList->size();
	unsigned targetLattice = LatTab.size();
	double min = Huge;
	double diff;
	Geometry geo;
	Point mean;
	for (unsigned i=0; i<p_pxCrackLatList->size(); i++) {
		for (unsigned d=0; d<Dim; d++) {
			mean[d] = 0.5*(GNodeTab[LatTab[(*p_pxCrackLatList)[i]].nb[0]].coord[d]
			              +GNodeTab[LatTab[(*p_pxCrackLatList)[i]].nb[1]].coord[d]);
		}
		diff = geo.dist(centre,mean);
		if (min>diff) {
			min = diff;
			targetLattice = (*p_pxCrackLatList)[i];
		}
	}
	genPriInjectWellSingle (targetLattice,InjectRate,pressure,p_conLat,p_glatRemoval,p_verti,p_fCal);
	return targetLattice;
}
void InitCond::genPriInjectWell		(	Point					centre,
                                        const double			initRadius,
                                        const double			injectRate,
                                        const double			pressure,
                                        ConnectLat*				p_conLat,
                                        Vertices*				p_verti,
                                        FluidCal*				p_fCal)
{
//    GLatRemoval 	glr(1e-4);
    Geometry		geo;
    dualOut << "Generate initial small fracture..." << std::endl;
    dualOut << "Target location, (xc, yc, zc) = " << centre[0] << " , " << centre[1] << " , " << centre[2] << std::endl;

    unsigned targetNodeID = GNodeTab.size();
    double	min = initRadius;
    for (unsigned NodeID = 0; NodeID < tNodeNum; NodeID++) {
//    	std::cout << geo.dist(GNodeTab[NodeID].coord,centre) << " ";
		if (geo.dist(GNodeTab[NodeID].coord,centre)<min) {
			min = geo.dist(GNodeTab[NodeID].coord,centre);
			targetNodeID = NodeID;
		}
    }

    if (targetNodeID != GNodeTab.size()) {
    	unsigned i=0;
    	do {
    		if (GNodeTab[targetNodeID].ne==0) {
    			break;
    		} else {
				unsigned ranLat =  GNodeTab[targetNodeID].n*(std::rand()/RAND_MAX);
				unsigned LatticeID = GNodeTab[targetNodeID].nbLatticeID[ranLat];
				if (!LatTab[LatticeID].isBreak) {
					if (p_fCal->putPriInjectNode(LatticeID,injectRate,pressure,p_verti,p_conLat)) {
						dualOut << "Lattice = " << LatticeID << " is set as injection fNode" << std::endl;
					}
					break;
				}
    		}
    	} while (++i>50);
    } else {
    	dualOut << "A node that has a distance with inputted centre smaller than " << initRadius << " cannot be found" << std::endl;
    }
}

void InitCond::genPriInjectWellSingle		(	const unsigned			LatticeID,
                                        		const double			injectRate,
                                        		const double			pressure,
                                        		ConnectLat*				p_conLat,
                                        		GLatRemoval*            p_glatRemoval,
                                        		Vertices*				p_verti,
                                        		FluidCal*				p_fCal)
{
//	GLatRemoval 	glr(1e-4);
	unsigned		nbNodeID[2];
	unsigned		nNum,nbLatticeID;
	if (!LatTab[LatticeID].isBreak) {
		if (p_fCal->putPriInjectNode(LatticeID,injectRate,pressure,p_verti,p_conLat)) {
			dualOut << "Lattice = " << LatticeID << " is set as injection fNode" << std::endl;
			p_glatRemoval->setNonReConLat(LatticeID);
			//Make sure its neighbors do not fail
			for (unsigned i=0; i<2; i++) {
				nbNodeID[i] = LatTab[LatticeID].nb[i];
				nNum = GNodeTab[nbNodeID[i]].n;
				for (unsigned j=0; j<nNum; j++) {
					nbLatticeID = GNodeTab[nbNodeID[i]].nbLatticeID[j];
//					}
				}
			}
		}
	}
}

void InitCond::genGPreExRegCrack	(	const double 		xc,
                                        const double 		yc,
                                        const double 		zc,
                                        const double 		a,
                                        const double 		b,
                                        const double 		c,
                                        ConnectLat*			p_conLat,
                                        Vertices*			p_verti)
{
    double xTerm,yTerm,zTerm;

    const double		Lx = UnitLength*Nx;
    const double		Ly = UnitLength*Ny;
    const double		Lz = UnitLength*Nz;

    dualOut << "Generate Pre-existing regular crack..." << endl;
    dualOut << "Penny shape crack, (xc, yc, zc, a, b, c) = " << xc << " , " << yc << " , " << zc << " , " << a << " , " << b << " , " << c << endl;

    Clock clock;
    clock.start("genCrack");

    if ((xc+a > 1.0)||(xc-a < 0.0)||(yc+b > 1.0)||(yc-b < 0.0)||(zc + c > 1.0)||(zc - c <0.0))
    {
        cout << " The crack dimensions inputted exceed the boundary! " << endl;
    }


    vector<unsigned>		latList;
    for (unsigned NodeID=0; NodeID<tNodeNum+tbNodeNum; NodeID++)
    {
        xTerm = (GNodeTab[NodeID].coord[0]/Lx - xc)/a;
        xTerm *= xTerm;
        yTerm = (GNodeTab[NodeID].coord[1]/Ly - yc)/b;
        yTerm *= yTerm;
        zTerm = (GNodeTab[NodeID].coord[2]/Lz - zc)/c;
        zTerm *= zTerm;

        if (xTerm+yTerm+zTerm < 1.0)
        {
        	p_glatRemoval->removeUnstableNode(NodeID,&nodeLists,&latList);
        }
    }

    uniqueVector(&latList);
    p_conLat->putMultiConLat(&latList,0,0.0);
    clock.get();
}

bool	InitCond::stablizeInitModel (	Vertices*			p_verti,
										ConnectLat*			p_conLat,
										PCG_Eigen*			p_EigenPCG,
										const unsigned		maxCoordNo,
										const unsigned		maxBreakNum,
										const unsigned		maxStep)
{
	unsigned step = 0;
	PCG_Eigen	solver(maxCoordNo,100000);
	GLatForce	glatForce;
	unsigned rmLatNum;
	unsigned rmLatAcc = 0;
	unsigned dummyStep1 = 0 ,dummyStep2 = 0;
	bool isLatBroken = false;
	bool isStable = false;
	dualOut << "Force calculation and lattice removal due to initial stress..." << std::endl;
	do {

		if(!p_EigenPCG -> solve(step,true,Use_CompactStiffnessMatrix)) {
			tripOut << "Model unstable, quit calculation!" << std::endl;
			return false;
		}
	    isLatBroken = (rmLatNum = p_glatRemoval->remove (p_conLat,p_verti,p_EigenPCG,maxBreakNum,0,0.0));
	    if (Use_ReconnectLattice) {
	        stableLat (p_conLat,p_verti,p_EigenPCG,isStable,dummyStep1,dummyStep2,10,2000);
	    }
	    rmLatAcc += rmLatNum;
	    dualOut << "Step = " << step << " , Lattice removed in this step = " << rmLatNum
	    		<< " Cummulated lattice removed = " << rmLatAcc << std::endl;
	    step++;
	} while ((isLatBroken) && (step < maxStep));
	if (step>=maxStep) {
		tripOut << "Max step reached, quit the node and continue to subsequent calculation! " << std::endl;
	}
	return true;
}

unsigned InitCond::stableLat 			(	ConnectLat*						p_conLat,
    										Vertices*						p_verti,
    										PCG_Eigen*						p_EigenPCG,
    										bool&							isStable,
    										unsigned&						relaxStep,
    										unsigned&						totalStep,
    										const unsigned					maxIter,
    										const unsigned					minLatTolerance)
{
	unsigned disConCnt,reConCnt;
	unsigned oldDisConCnt = tLatticeNum, oldReConCnt = tLatticeNum;
	const unsigned minLatChangeTol = 2;
	unsigned i = 0;
	tripOut << "-------------------Start Stabilizing Lattice Modelling-------------------" << std::endl;
	do {
		oldReConCnt = reConCnt;
		reConCnt  =	p_glatRemoval 	-> reConnectLat(p_conLat,p_EigenPCG,MaxPerStep_BrokenLattice);
		oldDisConCnt = disConCnt;
		disConCnt =	p_glatRemoval	-> disConnectLat(p_verti,p_conLat,p_EigenPCG,totalStep,0.0,MaxPerStep_BrokenLattice);

		if (Use_FullPlasticModel) {
		    p_glatRemoval   -> updateReConLatShearStiffness  (p_conLat,p_EigenPCG,p_verti,ReConRatio,totalStep,0.0);
		} else if (Use_FrictionModel) {
            p_glatRemoval   -> updateReConLatShearStiffness  (p_conLat,p_EigenPCG,p_verti,totalStep,0.0);
        }
		isStable = p_EigenPCG       -> solve(totalStep,true,Use_CompactStiffnessMatrix);
		i++;
		relaxStep++;
		totalStep++;
		if (i>maxIter) {
			tripOut << "Max iteration to disconnect/reconnect lattice is reached, iteration stop and proceed" << std::endl;
			tripOut << "-------------------Finished Stabilizing Lattice Modelling-------------------" << std::endl;
			return i;
		}
		if ((disConCnt==oldReConCnt)&&(reConCnt==oldDisConCnt)) {
			tripOut << "Some lattices are oscillating between connecting and disconneting, iteration stops and proceeds" << std::endl;
			break;
		}
		unsigned maxLatChangeTol = std::max((unsigned) std::sqrt(p_glatRemoval->getLatRmNum())/minLatTolerance,minLatChangeTol);
		if (disConCnt+reConCnt< maxLatChangeTol) {
			tripOut << "Only " << disConCnt+reConCnt << " lattices connecting/disconnecting, which is smaller than tolerance "
					<< maxLatChangeTol << " iteration stops and proceeds" << std::endl;
			break;
		}
		if (!isStable) {
			tripOut << "[stableLat] - model unstable, quit calculation" << std::endl;
			return i;
		}
		oldDisConCnt = reConCnt;
		oldReConCnt = reConCnt;
	} while ((disConCnt!=0) || (reConCnt!=0));

	dualOut << "Disconnect/reconnect lattice process is stabilized in " << i << " step" << std::endl;
	tripOut << "-------------------Finished Stabilizing Lattice Modelling-------------------" << std::endl;
	return i;
}

void InitCond::applyUnbreakableNb	(	const std::vector<unsigned>*		p_pxCrackLatList)
{
	const unsigned divLine = p_pxCrackLatList->size();
	unsigned NodeID,nNum,nbLatticeID;
	for (unsigned i = 0 ; i < divLine; i++) {
		for (unsigned j=0; j<2; j++) {
			NodeID = LatTab[(*p_pxCrackLatList)[i]].nb[j];
			nNum = GNodeTab[NodeID].n;
			if (NodeID<divLine) {
				for (unsigned k=0; k<nNum; k++) {
					nbLatticeID = GNodeTab[NodeID].nbLatticeID[k];
					if ((LatTab[nbLatticeID].nb[0]<divLine)&&(LatTab[nbLatticeID].nb[1]<divLine)) {
						if (!LatTab[nbLatticeID].isBreak) {
							LatTab[nbLatticeID].et[0] = Huge;
							LatTab[nbLatticeID].et[1] = -Huge;
						} else {
							tripOut << "[InitCond::applyUnbreakableNb]: Lattice has already broken, cannot become unbreakable!" << std::endl;
						}
					}
				}
			} else if (NodeID<2*divLine) {
				for (unsigned k=0; k<nNum; k++) {
					nbLatticeID = GNodeTab[NodeID].nbLatticeID[k];
					if (((LatTab[nbLatticeID].nb[0]<2*divLine)&&(LatTab[nbLatticeID].nb[0]>=divLine))
							&&(LatTab[nbLatticeID].nb[1]<2*divLine)&&(LatTab[nbLatticeID].nb[1]>=divLine)){
						if (!LatTab[nbLatticeID].isBreak) {
							LatTab[nbLatticeID].et[0] = Huge;
							LatTab[nbLatticeID].et[1] = -Huge;
						} else {
							tripOut << "[InitCond::applyUnbreakableNb]: Lattice has already broken, cannot become unbreakable!" << std::endl;
						}
					}
				}
			}
		}
	}
}

void InitCond::applyUnbreakableNbAll	(	const std::vector<unsigned>*		p_pxCrackLatList,
											const std::array<unsigned,2>		pxCrackPerimeterNum)
{
	std::set<unsigned>	latSet;
	std::set<unsigned>	nbSet,nbSet2;
	std::vector<unsigned>	edge,edge2;
	const unsigned div1 = p_pxCrackLatList->size()-pxCrackPerimeterNum[0];
	const unsigned div2 = p_pxCrackLatList->size()-pxCrackPerimeterNum[0]-pxCrackPerimeterNum[1];
	unsigned NodeID,nNum;
	for (unsigned i=0; i<p_pxCrackLatList->size(); i++) {
		NodeID = LatTab[(*p_pxCrackLatList)[i]].nb[0];
		if (NodeID>=div1) {
			edge.push_back(i);
		} else {
			if (NodeID>=div2) {
				edge2.push_back(i);
			}
			for (unsigned j=0; j<2; j++) {
				NodeID = LatTab[(*p_pxCrackLatList)[i]].nb[j];
				nNum = GNodeTab[NodeID].n;
				for (unsigned k=0; k<nNum; k++) {
					latSet.insert(GNodeTab[NodeID].nbLatticeID[k]);
				}
			}
		}
	}
	for (unsigned i=0 ; i<edge2.size(); i++) {
		for (unsigned j=0; j<2; j++) {
			NodeID = LatTab[(*p_pxCrackLatList)[edge2[i]]].nb[j];
			nNum = GNodeTab[NodeID].n;
			for (unsigned k=0; k<nNum; k++) {
				nbSet.insert(GNodeTab[NodeID].nbNodeID[k]);
			}
		}
	}
//	const unsigned pxNodeNum = 2*p_pxCrackLatList->size();
	for (unsigned i=0 ; i<edge.size(); i++) {
		for (unsigned j=0; j<2; j++) {
			NodeID = LatTab[(*p_pxCrackLatList)[edge[i]]].nb[j];
			nbSet2.insert(NodeID);		}
	}

	unsigned LatticeID,nbNodeID,cnt=0;
	for (unsigned i=0; i<edge.size(); i++) {
		for (unsigned j=0; j<2; j++) {
			NodeID = LatTab[(*p_pxCrackLatList)[edge[i]]].nb[j];
			nbSet2.erase(LatTab[(*p_pxCrackLatList)[edge[i]]].nb[(j+1)%2]);
			nNum = GNodeTab[NodeID].n;
			for (unsigned k=0; k<nNum; k++) {
				LatticeID = GNodeTab[NodeID].nbLatticeID[k];
				nbNodeID = LatTab[LatticeID].nb[0];
				if (nbNodeID==NodeID) {
					nbNodeID = LatTab[LatticeID].nb[1];
				}
				if (std::find(nbSet.begin(),nbSet.end(),nbNodeID)!=nbSet.end()) {
					latSet.insert(LatticeID);
					cnt++;
				}
				else {
					if (std::find(nbSet2.begin(),nbSet2.end(),nbNodeID)!=nbSet2.end()) {
						latSet.insert(LatticeID);
						cnt++;
					}
				}
			}
			nbSet2.insert(LatTab[(*p_pxCrackLatList)[edge[i]]].nb[(j+1)%2]);
		}
	}
	//To remove the pre-existing crack in the set
	for (unsigned i=0; i < p_pxCrackLatList->size(); i++) {
		latSet.erase((*p_pxCrackLatList)[i]);
	}

	for (std::set<unsigned>::iterator it=latSet.begin(); it!=latSet.end(); ++it) {
		if (!LatTab[*it].isBreak) {
			LatTab[*it].et[0] = Huge;
			LatTab[*it].et[1] = -Huge;
		} else {
			tripOut << "[InitCond::applyUnbreakableNbAll]: Lattice has already broken, cannot become unbreakable!" << std::endl;
		}
	}
}

void InitCond::applyUnbreakableOuterLat	()
{
	std::array<double,2> unbreakable;
	unbreakable[0] = Huge;
	unbreakable[1] = -Huge;
	for (unsigned LatticeID=0; LatticeID<tLatticeNum; LatticeID++) {
		if ((LatTab[LatticeID].nb[0]>=tInNodeNum)||(LatTab[LatticeID].nb[1]>=tInNodeNum)) {
			LatTab[LatticeID].et = unbreakable;
		}
	}
}

void InitCond::applySmallDistortPxCrack	(	const std::vector<unsigned>*		p_pxCrackLatList)
{
	const double delta = 0.05*UnitLength;
	double sgn;
	unsigned NodeID;
	for (unsigned i=0; i<p_pxCrackLatList->size(); i++) {
		for (unsigned j=0; j<2; j++) {
			sgn = 2.0*j-1.0;
			NodeID = LatTab[(*p_pxCrackLatList)[i]].nb[j];
			GNodeTab[NodeID].coord[2] += delta*((double) std::rand()/RAND_MAX)*sgn;
		}
	}
}

void 	InitCond::modifyTensileStrain	(	std::vector<unsigned>*		p_latList,
											double						value)
{
	for (std::vector<unsigned>::iterator it=p_latList->begin(); it!=p_latList->end(); ++it) {
		LatTab[*it].et[0] = value;
	}
}

void InitCond::modifyLatStiffnessAndThreshold(const unsigned LatticeID, const double ratio)
{
    LatTab[LatticeID].k[0] *= ratio;
    LatTab[LatticeID].et[0] *=ratio;
    LatTab[LatticeID].et[1] *=ratio;
}

void InitCond::modifyLatStiffness(const unsigned LatticeID, const double ratio)
{
    LatTab[LatticeID].k[0] *= ratio;
}

void InitCond::modifyLatThreshold(const unsigned LatticeID, const double ratio)
{
    LatTab[LatticeID].et[0] *=ratio;
    LatTab[LatticeID].et[1] *=ratio;
}

void InitCond::setBoreHoleConstrain ()
{
    std::vector<unsigned> VecDir;
    VecDir.push_back(0);    VecDir.push_back(1);    VecDir.push_back(2);
    VecDir.push_back(3);    VecDir.push_back(4);    VecDir.push_back(5);
    std::vector<double> VecDisp(6,0.0);
    for (unsigned i=0; i<boreHoleNodeList.size(); i++) {
        Constrain con(boreHoleNodeList[i],VecDir,VecDisp);
        nodeLists.constrain.push_back(con);
    }
}

void InitCond::setBoreHoleNodeList (   std::vector<unsigned>*      p_latList)
{
    boreHoleNodeList = *p_latList;
}

void InitCond::setZero()
{
    for (unsigned NodeID = 0; NodeID < GNodeTab.size(); NodeID++) {
        for (unsigned k=0; k<DofPerNode; k++) {
            GNodeTab[NodeID].d0[k] = GNodeTab[NodeID].d[k];
        }
    }
    dualOut << "Set zero for displacement and external force" << std::endl;
}

