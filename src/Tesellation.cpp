/*
 * Tesellation.cpp
 *
 *  Created on: Nov 18, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "Tesellation.hpp"

using namespace voro;

std::vector			<GNode>				GNodeTab;
std::vector			<Lattice>		    LatTab;
std::vector			<Beam>		        LatBeamTab;

VoroTesell::VoroTesell() {
    // TODO Auto-generated constructor stub

}

VoroTesell::~VoroTesell() {
    // TODO Auto-generated destructor stub
    delete	p_voroInfo;
}
double avgD;

void VoroTesell::generate	                (	const double			_minDist,
												double					scaleBoundary,
												double					outerScaleLength,
												unsigned*				p_maxCoordNo,
												Vertices*				p_verti)
{
    unsigned			NodeID = 0;
    GeoID = GeometryID;
    scaleX = scaleBoundary;
    scaleY = scaleBoundary;
    scaleZ = scaleBoundary;

    if (Use_UnevenMesh) {
        if (std::max(std::max(Nx,Ny),Nz) > 20) {
            coarseness = 2.5;
        } else {
            coarseness = 1.0;
        }
    } else {
        coarseness = 1.0;
    }
    minD = _minDist;
    if (Use_const_NodeNum) {
        minD = _minDist*std::pow(PackingDensity/0.5,1.0/3.0);
    }
    dualOut << "minD = " << minD << std::endl;
    anisotropy = 1.0;
    notchLen = 0.1733;
    Clock rclock;
    rclock.start("Voro");

    bool useCompositeMesh = true;
    if (scaleBoundary < 1.0+Tiny) {
    	tripOut << "Inputed ScaleBoundary =" << scaleBoundary << " should be greater than 1.0, uniform mesh will be used" << std::endl;
    	scaleBoundary = outerScaleLength = 1.0;
    	scaleX = scaleY = scaleZ = 1.0;
        useCompositeMesh = false;
    }
    if (outerScaleLength < scaleBoundary) {
        tripOut << "Inputed scaleLength is smaller than scaleBoundary, scaleLength is set to be the same as scaleBoundary!" << std::endl;
        outerScaleLength = scaleBoundary;
    }
    unsigned    layer = 1;
    scaleX = scaleY = std::pow(scaleBoundary,layer);
    container_poly              conp    (       (-(scaleX-1.0)/2.0*Nx-scaleX*0.0-Tiny)*UnitLength,
                                                ((scaleX+1.0)/2.0*Nx+scaleX*0.0+Tiny)*UnitLength,
                                                (-(scaleY-1.0)/2.0*Ny-scaleY*0.0-Tiny)*UnitLength,
                                                ((scaleY+1.0)/2.0*Ny+scaleY*0.0+Tiny)*UnitLength,
                                                (-(scaleZ-1.0)/2.0*Nz-scaleZ*0.0-Tiny)*UnitLength,
                                                ((scaleZ+1.0)/2.0*Nz+scaleZ*0.0+Tiny)*UnitLength,
                                                Nx,Ny,Nz,
                                                false,false,false,8);

    initBoundary();

    bool isStable = true;
    unsigned	cnt_genModel=0;
    std::vector<NodeIDPair> conjNodeList;
    do {
    	conjNodeList = genRanNode 		(&conp,&NodeID);
		tInNodeNum = NodeID;

		if (useCompositeMesh) {
	            genOuterNode    (&conp,&NodeID,coarseness*minD,outerScaleLength,scaleBoundary,layer);
//		    }
		}
		if (anisotropy > 1.0) {
		    genRanNodeAnisotropic(&NodeID,&conp,0.5/anisotropy);
		}
		genRegBoundaryNodes(&conp,&NodeID);
		tNodeNum = NodeID;
		strengthAdj = std::sqrt(std::pow(Nx*Ny*Nz/tInNodeNum,0.333333333)*UnitLength);
		dualOut << "No of normal node = " << tNodeNum << std::endl;
		tbNodeNum = NodeID - tNodeNum;
		dualOut << "No of boundary node = " << tbNodeNum << std::endl;
		tvNodeNum = NodeID - tNodeNum - tbNodeNum;
		dualOut << "No of virtual node = " << tvNodeNum << std::endl;

		updateVoroInfo(&conp);
		//Return a lattice ID list for pre-existed crack
		findNeighbour(outerScaleLength,&conjNodeList);
		const double radius = CrackRadius*Nx*UnitLength;
		Point centre;
		centre[0] = UnitLength*Nx*0.5;
		centre[1] = UnitLength*Ny*0.5;
		centre[2] = UnitLength*Nz*0.5;
		avgD = std::pow(((double) Nx*Ny*Nz)/Nnum,0.3333333)*UnitLength;
		dualOut << "avgD = " << avgD << std::endl;
		if (GeoID==3) {
		    genRoughPennyCrack(centre,radius,2.0*avgD);
		}
		*p_maxCoordNo = updateNodeCoordinationNo()+1;

		unsigned		MinCoordNo = 3;
		if (KNum==4) {
			MinCoordNo = 1;
		}
		if ((testMinCoordNum(MinCoordNo))||(true)) {
			dualOut << " Every node has at least " << MinCoordNo << " neightbour, model stable" << std::endl;
			isStable = true;
		} else {
			dualOut << "*****Warning!*****" << std::endl;
			dualOut << "Model Unstable!, re-generate the model..." << std::endl;
			isStable = false;
			reset (&conp);
			NodeID = 0;
		}
		cnt_genModel++;
		if (cnt_genModel>=3) {
			tripOut << "Too many trials, stop trying to generate a stable model, please check the input" << std::endl;
			break;
		}
    } while (!isStable);
    updateVertices (p_verti);
    updateLatInfo (p_verti);
    updateNodeVoroVol();
    fillInBoundary();
    std::vector<unsigned> inBNodeList = getInBoundary(0);
    applyUnbreakableLattice(&inBNodeList);
    applyUnbreakableNode ();
    if (GeoID==4) {
        setBreakableLatticeAtBoundary ();
    }
    findPreExistingCrack(&conjNodeList);
    if (GeoID==4) {
        setUnbreakablePair ();
    }
    findWeakPlane();

    assignRandomStrength(StrengthDistModel,Sigma,0.2);
    rclock.get();
    return;
}

std::vector<NodeIDPair> VoroTesell::genRanNode	(	container_poly* 	p_conp,
                                					unsigned* 			p_NodeID)
{
    double                      n = PackingDensity;         //Percent from most dense HCF packing 0.7405
    double                      Nxyz = Nx*Ny*Nz;
    unsigned                    Nnum_HCF = sqrt(2)*Nxyz;
    if (Use_const_NodeNum) {
        n = 0.5;
    }
    Nnum = Nnum_HCF*n;
    dualOut << "Nnum = " << Nnum << std::endl;
    Mat3Dvec			NodeInPartition;
    GNodeTab.clear();
    NodeInPartition.resize(Nx*UnitLength/minD+1,Ny*UnitLength/minD+1,Nz*UnitLength/minD+1);
    dualOut << "Nx*L/D = " << Nx*UnitLength/minD << " Ny*L/D = " << Ny*UnitLength/minD << " Nz*L/D = " << Nz*UnitLength/minD << std::endl;
    dualOut << " VoroTesell::genRanNode started..." << std::endl;
    if (Use_RandomNode) {
        srand(time(NULL));  //To randomize the seed
    }
    std::cout << *p_NodeID << '\n';
    Point centre;
    centre[0] = UnitLength*Nx*0.5;
    centre[1] = UnitLength*Ny*0.5;
    centre[2] = UnitLength*Nz*0.5;
    const double radius = CrackRadius*Nx*UnitLength;
    const double outRadius = 1.5*(Nx+Ny)*0.5*0.5*UnitLength; //= 3.0*radius;
    avgD = std::pow(((double) Nx*Ny*Nz)/Nnum,0.3333333)*UnitLength;
    const double ld = avgD;
    const double lt = avgD/100.0;
    dualOut << "ld = " << ld << " lt = " << lt << std::endl;
    crackNodeNum = 0;
    Vec     normVec;
    normVec[0] = std::sin(Theta*Pi/180.0);        normVec[1] = 0.0;       normVec[2] = std::cos(Theta*Pi/180.0);
    std::vector<NodeIDPair> conjNodeList;
    if (GeoID==2) {
        genNotches (notchLen*Nx*UnitLength,ld,lt,p_NodeID,&conjNodeList,p_conp,&NodeInPartition);
    } else if ((GeoID==1)||(GeoID==4)) {
        bottomLocation = centre;
        BHdepth = Nz*UnitLength-Tiny;
        bottomLocation[2] = centre[2]-BHdepth/2.0;
        unsigned angRes = std::max(12.01,(1.1*2.0*Pi*CrackRadius*UnitLength/minD+1.0));
        angRes +=angRes%2;  // make sure angRes is divisible by 2
        genBorehole2 (bottomLocation,CrackRadius*UnitLength,BHdepth,angRes,normVec,p_NodeID,p_conp,&NodeInPartition,&conjNodeList);
    } else if (GeoID==0) {
        if (radius>Tiny) {
           genInclinedPennyCrackRan(centre,radius,normVec,ld,lt,p_NodeID,&conjNodeList,p_conp,&NodeInPartition);
        }
    }
    unsigned crackNum1 = *p_NodeID; //-regNodeNum;
    Point centre2;
    centre2[0] = UnitLength*Nx*0.35;
    centre2[1] = UnitLength*Ny*0.35;
    centre2[2] = UnitLength*Nz*0.35;
    Vec     normVec2;
    normVec2[0] = std::sqrt(0.5);        normVec2[1] = 0.0;       normVec2[2] = std::sqrt(0.5);
    if ((GeoID!=1)&&(GeoID!=4)) {
        for (unsigned i = 0; i< crackNum1; i++) {
            nodeLists.free.push_back(*p_NodeID);
            conjNodeList[i].i2 = *p_NodeID;
            putNode(p_conp,p_NodeID,GNodeTab[i].coord[0]-lt*normVec[0],
            GNodeTab[i].coord[1]-lt*normVec[1],GNodeTab[i].coord[2]-lt*normVec[2],Type_0);
        }
    }
    crackNodeNum = *p_NodeID;
    refineNode[0] = *p_NodeID;
    dualOut << "crackNodeNum = " << crackNodeNum << std::endl;
    double weakPlaneRadius = radius;
    if ((StrengthFactor>1.0)&&(radius>Tiny)) {
        if (GeoID==0) {
            genIncOuterPennyCrackRan(centre,radius,outRadius,normVec,ld,lt,p_NodeID,p_conp,&NodeInPartition);
            weakPlaneRadius = outRadius;
        }
    }
    if ((GeoID==4)&&(StrengthFactor>1.0)) {
        genNotchWeakPlane (centre,8.9,minD/10.0,p_NodeID,p_conp,&NodeInPartition);
    }
    double anisotropy=1.0;
    if (!Use_UnevenMesh) {
//    if (true) {
        genRanNode  (p_NodeID,p_conp,&NodeInPartition,centre,normVec,lt,ld,weakPlaneRadius,conjNodeList.empty());
        refineNode[1] = *p_NodeID-1;
    } else {
        double lmin=coarseness*minD;
        double rRatio;
        if (GeometryID==4) {
            rRatio = std::min(std::max(2.5*ld/(UnitLength*Ny/2.0),0.1),0.5);
            genRanNodeReduced  (p_NodeID,p_conp,&NodeInPartition,centre,normVec,lt,ld,weakPlaneRadius,1,rRatio,conjNodeList.empty());
            refineNode[1] = *p_NodeID-1;
            std::cout << "NodeID = " << *p_NodeID << std::endl;
            genRanNodeReduced  (p_NodeID,p_conp,&NodeInPartition,lmin,1,rRatio);
            std::cout << "NodeID = " << *p_NodeID << std::endl;
        } else {
            rRatio = std::min(std::max(2.5*ld/(UnitLength*Nz/2.0),0.1),0.5);
            if (std::min(Nx,Ny)<4) {
                rRatio *=2;
            }
            if (anisotropy >1.0) {
                genRanNodeReduced  (p_NodeID,p_conp,&NodeInPartition,centre,normVec,lt,ld,weakPlaneRadius,2,rRatio,
                                anisotropy,conjNodeList.empty());
            } else {
                genRanNodeReduced  (p_NodeID,p_conp,&NodeInPartition,centre,normVec,lt,ld,weakPlaneRadius,2,rRatio,conjNodeList.empty());
            }
            refineNode[1] = *p_NodeID-1;
            std::cout << "NodeID = " << *p_NodeID << std::endl;
            genRanNodeReduced  (p_NodeID,p_conp,&NodeInPartition,lmin,2,rRatio);
            std::cout << "NodeID = " << *p_NodeID << std::endl;
        }
    }


    return conjNodeList;
}

void 	VoroTesell::genRanNode					(	unsigned* 			p_NodeID,
													container_poly* 	p_conp,
													Mat3Dvec*			p_NodeInPartition,
													const Point 		centre,
													const Vec     		normVec,
													const double		lt,
													const double		ld,
													const double		weakPlaneRadius,
													const bool			isConjNodeListEmpty)
{
	unsigned nodeCnt = *p_NodeID;
	unsigned rejectCnt = 0;
	double nMaxRej = 10.0;
	if (MaxTotalStep == 0) {
	    nMaxRej = Disp_Scale; //= 0.00001;
	}

	unsigned maxRej = nMaxRej*sqrt(2)*nx*ny*nz;
	std::array<unsigned,Dim> 	p;
	Point rTest;
	if (Use_RandomNode) {
	    std::srand (std::time(NULL));
	}
	while (nodeCnt < Nnum) {
		rTest[0] = (((double) std::rand())/RAND_MAX*nx+ex)*UnitLength;
		rTest[1] = (((double) std::rand())/RAND_MAX*ny+ey)*UnitLength;
		rTest[2] = (((double) std::rand())/RAND_MAX*nz+ez)*UnitLength;
		if (Ny==1) {
		    rTest[1] = (0.5 + (((double) std::rand())/RAND_MAX-0.5)/10.0)*Ny*UnitLength;
		}
		p[0] = ( rTest[0]/minD-Tiny);
		p[1] = ( rTest[1]/minD-Tiny);
		p[2] = ( rTest[2]/minD-Tiny);


		if (testOverlap(p_NodeInPartition,minD,rTest,*p_NodeID,p)==0) {
			if ((((!isRestrictedZone(rTest,centre,normVec,0.5*lt+1.25*ld,weakPlaneRadius))
						)||(isConjNodeListEmpty))||(GeoID!=0)) {
			    if (((!isRestrictedZone (rTest,notchLen*Nx*UnitLength,1.0*ld,lt))||(GeoID!=2))&&
			        ((!isRestrictedZoneBH(rTest,bottomLocation,CrackRadius*UnitLength+minD/100.0,BHdepth))||((GeoID!=1)&&(GeoID!=4)))&&
			        (((StrengthFactor<=1.0)||(GeoID!=4))||((rTest[1]>centre[1]+1.0*ld)||(rTest[1]<centre[1]-1.0*ld)))){
                    nodeLists.free.push_back(*p_NodeID);
                    putNode(p_conp,p_NodeID,rTest,Type_0);
                    nodeCnt++;
                    rejectCnt = 0;
			    }
			}
	    } else {
	    	rejectCnt++;
	    	if (rejectCnt>maxRej) {
	    		dualOut << "Cannot place a node, force quit the loop"<< std::endl;
	    		dualOut << "Modified n (packing ratio relative to HCF packing) = "
	    				<< ((double) nodeCnt)/(std::sqrt(2.0)*Nx*Ny*Nz) << std::endl;
	    		if (MaxTotalStep==0) {
	    		    char fileName[255];
	    		    std::sprintf(fileName,"%s/%snLog.vtk",OutputSubFolder,OutputFilePrefix);
	    		    std::ofstream file(fileName, std::ios::out | std::ios::app);
	    		    file <<  nMaxRej  << '\t' << ((double) nodeCnt)/(std::sqrt(2.0)*Nx*Ny*Nz) << std::endl;
	    		}
	    		Nnum = nodeCnt;
	    		break;
	        }
	    }
	}
}

void    VoroTesell::initBoundary() {

    unsigned face[6] = {Face3, Face5, Face2, Face4, Face1, Face6};
    nx = Nx - 2;
    ny = Ny - 2;
    nz = Nz - 2;

    ex = 1;
    ey = 1;
    ez = 1;

    if (scaleX-Tiny>1.0) {
        nx = Nx;
        ex = 0;
    } else {
    }

    if (scaleY-Tiny>1.0) {
        ny = Ny;
        ey = 0;
    } else {
    }

    if (scaleZ-Tiny>1.0) {
        nz = Nz;
        ez = 0;
    } else {
    }
    if (Nz<4) {
        nz = Nz;
        ez = 0;
    }
    if ((Nx<4)||(Ny<4)){
        nx = Nx;
        ex = 0;
        ny = Ny;
        ey = 0;
    }
    dualOut << "(nx,ny,nz) = " << nx << " , " << ny << " , " << nz << std::endl;
    dualOut << "(ex,ey,ez) = " << ex << " , " << ey << " , " << ez << std::endl;
}
void    VoroTesell::genRanNodeReduced                   (   unsigned*           p_NodeID,
                                                            container_poly*     p_conp,
                                                            Mat3Dvec*           p_NodeInPartition,
                                                            const Point         centre,
                                                            const Vec           normVec,
                                                            const double        lt,
                                                            const double        ld,
                                                            const double        weakPlaneRadius,
                                                            const unsigned      dir,
                                                            double              rRatio,
                                                            const bool          isConjNodeListEmpty)
{
    dualOut << "VoroTesell::genRanNodeReduced (dense thin layer) starts...." << std::endl;
    Clock clock;
    clock.start("ranNodeRed");
    if ((dir!=0)&&(dir!=1)&&(dir!=2)) {
        tripOut << "dir inputted = " << dir << " is not 0, 1 or 2, default value of 0 will be used" << std::endl;
    }
    double                      n = PackingDensity;         //Percent from most dense HCF packing 0.7405
    if (Use_const_NodeNum) {
        n = 0.5;
    }

    double                      Nxyz = Nx*Ny*Nz;
    unsigned                    Nnum_HCF = sqrt(2)*Nxyz;
    unsigned targetNnum = Nnum_HCF*n*rRatio;
    unsigned nodeCnt = *p_NodeID;
    unsigned rejectCnt = 0;
    double nMaxRej = 10.0;
    if (MaxTotalStep == 0) {
        nMaxRej = Disp_Scale; //= 0.00001;
    }
    if (rRatio>=1.0) {
        tripOut << "[VoroTesell::genRanNodeReducedZ], rRatio = " << rRatio << " is greater than 1.0, set rRatio = 1.0" << std::endl;
        rRatio = 1.0;
    }
    unsigned maxRej = nMaxRej*sqrt(2.0)*Nx*Ny*Nz*rRatio;
    std::array<unsigned,Dim>    p;
    Point rTest;
    if (Use_RandomNode) {
        std::srand (std::time(NULL));
    }
    double mid;
    if (dir==2) {
        mid = 0.5*Nz*UnitLength;
    } else if (dir==1) {
        mid = 0.5*Ny*UnitLength;
    } else {
        mid = 0.5*Nx*UnitLength;
    }
    while (nodeCnt < targetNnum) {
        if (dir==2) {
            rTest[0] = (((double) std::rand())/RAND_MAX*nx+ex)*UnitLength;
            rTest[1] = (((double) std::rand())/RAND_MAX*ny+ey)*UnitLength;
            rTest[2] = (((double) std::rand())/RAND_MAX-0.5)*Nz*UnitLength*rRatio+mid;
        } else if (dir==1) {
            rTest[0] = (((double) std::rand())/RAND_MAX*nx+ex)*UnitLength;
            rTest[1] = (((double) std::rand())/RAND_MAX-0.5)*Ny*UnitLength*rRatio+mid;
            rTest[2] = (((double) std::rand())/RAND_MAX*nz+ez)*UnitLength;
        } else {
            rTest[0] = (((double) std::rand())/RAND_MAX-0.5)*Nx*UnitLength*rRatio+mid;
            rTest[1] = (((double) std::rand())/RAND_MAX*ny+ey)*UnitLength;
            rTest[2] = (((double) std::rand())/RAND_MAX*nz+ez)*UnitLength;
        }
/*      if (is2D==true) {
            rTest[1] = 0.5*Ny*UnitLength;
        }*/
        p[0] = ( rTest[0]/minD-Tiny);
        p[1] = ( rTest[1]/minD-Tiny);
        p[2] = ( rTest[2]/minD-Tiny);
        if (testOverlap(p_NodeInPartition,minD,rTest,*p_NodeID,p)==0) {
            if ((((!isRestrictedZone(rTest,centre,normVec,0.5*lt+1.25*ld,weakPlaneRadius))
                        )||(isConjNodeListEmpty))||(GeoID!=0)) {
                if (((!isRestrictedZone (rTest,notchLen*Nx*UnitLength,1.0*ld,lt))||(GeoID!=2))&&
                    ((!isRestrictedZoneBH(rTest,bottomLocation,CrackRadius*UnitLength+minD/100.0,BHdepth))||((GeoID!=1)&&(GeoID!=4)))&&
                    (((StrengthFactor<=1.0)||(GeoID!=4))||((rTest[1]>centre[1]+1.0*ld)||(rTest[1]<centre[1]-1.0*ld)))){
                    nodeLists.free.push_back(*p_NodeID);
                    putNode(p_conp,p_NodeID,rTest,Type_0);
                    nodeCnt++;
                    rejectCnt = 0;
                }
            }
        } else {
            rejectCnt++;
            if (rejectCnt>maxRej) {
                dualOut << "Cannot place a node, force quit the loop"<< std::endl;
                dualOut << "Modified n (packing ratio relative to HCF packing) = "
                        << PackingDensity*((double) nodeCnt)/(targetNnum) << std::endl;
                if (MaxTotalStep==0) {
                    char fileName[255];
                    std::sprintf(fileName,"%s/%snLog.vtk",OutputSubFolder,OutputFilePrefix);
                    std::ofstream file(fileName, std::ios::out | std::ios::app);
                    file <<  nMaxRej  << '\t' << ((double) nodeCnt)/(std::sqrt(2.0)*Nx*Ny*Nz) << std::endl;
                }
                break;
            }
        }
    }
    clock.get();
}

void    VoroTesell::genRanNodeReduced                   (   unsigned*           p_NodeID,
                                                            container_poly*     p_conp,
                                                            Mat3Dvec*           p_NodeInPartition,
                                                            const Point         centre,
                                                            const Vec           normVec,
                                                            const double        lt,
                                                            const double        ld,
                                                            const double        weakPlaneRadius,
                                                            const unsigned      dir,
                                                            double              rRatio,
                                                            double              anisotropy,
                                                            const bool          isConjNodeListEmpty)
{
    dualOut << "VoroTesell::genRanNodeReduced (dense thin layer, anisotropic) starts...." << std::endl;
    dualOut << " anisotropy = " << anisotropy << std::endl;
    Clock clock;
    clock.start("ranNodeRed");
    if ((dir!=0)&&(dir!=1)&&(dir!=2)) {
        tripOut << "dir inputted = " << dir << " is not 0, 1 or 2, default value of 0 will be used" << std::endl;
    }
    double                      n = PackingDensity;         //Percent from most dense HCF packing 0.7405
    if (Use_const_NodeNum) {
        n = 0.5;
    }

    double   Nxyz = Nx*Ny*Nz;
    unsigned Nnum_HCF = sqrt(2)*Nxyz;
    unsigned targetNnum = Nnum_HCF*n*rRatio*anisotropy;
    unsigned nodeCnt = *p_NodeID;
    unsigned rejectCnt = 0;
    double nMaxRej = 10.0;
    if (MaxTotalStep == 0) {
        nMaxRej = Disp_Scale; //= 0.00001;
    }
    if (rRatio>=1.0) {
        tripOut << "[VoroTesell::genRanNodeReducedZ], rRatio = " << rRatio << " is greater than 1.0, set rRatio = 1.0" << std::endl;
        rRatio = 1.0;
    }
    unsigned maxRej = nMaxRej*sqrt(2.0)*Nx*Ny*Nz*rRatio*anisotropy;
    std::array<unsigned,Dim>    p;
    Point rTest;
    if (Use_RandomNode) {
        std::srand (std::time(NULL));
    }
    double mid;
    if (dir==2) {
        mid = 0.5*Nz*UnitLength;
    } else if (dir==1) {
        mid = 0.5*Ny*UnitLength;
    } else {
        mid = 0.5*Nx*UnitLength;
    }
    while (nodeCnt < targetNnum) {
        if (dir==2) {
            rTest[0] = (((double) std::rand())/RAND_MAX*nx+ex)*UnitLength;
            rTest[1] = (((double) std::rand())/RAND_MAX*ny+ey)*UnitLength;
            rTest[2] = (((double) std::rand())/RAND_MAX-0.5)*Nz*UnitLength*rRatio+mid;
        } else if (dir==1) {
            rTest[0] = (((double) std::rand())/RAND_MAX*nx+ex)*UnitLength;
            rTest[1] = (((double) std::rand())/RAND_MAX-0.5)*Ny*UnitLength*rRatio+mid;
            rTest[2] = (((double) std::rand())/RAND_MAX*nz+ez)*UnitLength;
        } else {
            rTest[0] = (((double) std::rand())/RAND_MAX-0.5)*Nx*UnitLength*rRatio+mid;
            rTest[1] = (((double) std::rand())/RAND_MAX*ny+ey)*UnitLength;
            rTest[2] = (((double) std::rand())/RAND_MAX*nz+ez)*UnitLength;
        }
/*      if (is2D==true) {
            rTest[1] = 0.5*Ny*UnitLength;
        }*/
        p[0] = ( rTest[0]/minD-Tiny);
        p[1] = ( rTest[1]/minD-Tiny);
        p[2] = ( rTest[2]/minD-Tiny);
        if (testOverlapAni(p_NodeInPartition,minD,anisotropy,dir,rTest,*p_NodeID,p)==0) {
            if ((((!isRestrictedZone(rTest,centre,normVec,0.5*lt+1.25*ld,weakPlaneRadius))
    //              &&(!isRestrictedZone(rTest,centre2,normVec2,0.5*lt+1.25*ld,outRadius))
                        )||(isConjNodeListEmpty))||(GeoID!=0)) {
                if (((!isRestrictedZone (rTest,notchLen*Nx*UnitLength,1.0*ld,lt))||(GeoID!=2))&&
                    ((!isRestrictedZoneBH(rTest,bottomLocation,CrackRadius*UnitLength+minD/100.0,BHdepth))||((GeoID!=1)&&(GeoID!=4)))&&
                    (((StrengthFactor<=1.0)||(GeoID!=4))||((rTest[1]>centre[1]+1.0*ld)||(rTest[1]<centre[1]-1.0*ld)))){
                    nodeLists.free.push_back(*p_NodeID);
                    putNode(p_conp,p_NodeID,rTest,Type_0);
                    nodeCnt++;
                    rejectCnt = 0;
                }
            }
        } else {
            rejectCnt++;
            if (rejectCnt>maxRej) {
                dualOut << "Cannot place a node, force quit the loop"<< std::endl;
                dualOut << "Modified n (packing ratio relative to HCF packing) = "
                        << PackingDensity*((double) nodeCnt)/(targetNnum) << std::endl;
                if (MaxTotalStep==0) {
                    char fileName[255];
                    std::sprintf(fileName,"%s/%snLog.vtk",OutputSubFolder,OutputFilePrefix);
                    std::ofstream file(fileName, std::ios::out | std::ios::app);
                    file <<  nMaxRej  << '\t' << ((double) nodeCnt)/(std::sqrt(2.0)*Nx*Ny*Nz) << std::endl;
                }
                break;
            }
        }
    }
    clock.get();
}

void    VoroTesell::genRanNodeReduced                   (   unsigned*           p_NodeID,
                                                            container_poly*     p_conp,
                                                            Mat3Dvec*           p_NodeInPartition,
                                                            const double        lmin,
                                                            const unsigned      dir,
                                                            double              rRatio)
{
    dualOut << "VoroTesell::genRanNodeReduced (ramain corase mesh) starts...." << std::endl;
    Clock clock;
    clock.start("ranNodeRed");
    if ((dir!=0)&&(dir!=1)&&(dir!=2)) {
        tripOut << "dir inputted = " << dir << " is not 0, 1 or 2, default value of 0 will be used" << std::endl;
    }
    unsigned nodeCnt = 0;
    unsigned rejectCnt = 0;
    double nMaxRej = 10.0;
    double factor = (minD/lmin)*(minD/lmin)*(minD/lmin);
    if (std::min(std::min(Nx,Ny),Nz)<4) {
        factor /= (minD/lmin);
    }
    double n = PackingDensity;         //Percent from most dense HCF packing 0.7405
    if (Use_const_NodeNum) {
        n = 0.5;
    }

    unsigned face[6] = {Face3, Face5, Face2, Face4, Face1, Face6};

    double                      Nxyz = nx*ny*nz;
    unsigned                    Nnum_HCF = sqrt(2.0)*Nxyz;
    double adj,mid;

    if (dir==0) {
        adj = (nx*UnitLength*(1.0-rRatio)-lmin)/(nx*UnitLength*(1.0-rRatio));
        mid = 0.5*Nx*UnitLength;
    } else if (dir==1) {
        adj = (ny*UnitLength*(1.0-rRatio)-lmin)/(ny*UnitLength*(1.0-rRatio));
        mid = 0.5*Ny*UnitLength;
    } else if (dir==2) {
        adj = (nz*UnitLength*(1.0-rRatio)-lmin)/(nz*UnitLength*(1.0-rRatio));
        mid = 0.5*Nz*UnitLength;
    } else {
        tripOut << "[VoroTesell::genRanNodeReduced], dir = " << dir << " is neither 0, 1, 2, do nothing and leave" << std::endl;
        return;
    }
    dualOut << "adj = " << adj << " , mid " << mid << std::endl;
    dualOut << "(nx,ny,nz) = " << nx << " , " << ny << " , " << nz << std::endl;

    unsigned targetNum = Nnum_HCF*n*(1-rRatio)*factor*adj;
    dualOut << "[VoroTesell::genRanNodeReduced] - targetNum = " << targetNum << std::endl;
    dualOut << "factor = " << factor << std::endl;
    if (MaxTotalStep == 0) {
        nMaxRej = Disp_Scale; //= 0.00001;
    }
    if (rRatio>=1.0) {
        tripOut << "[VoroTesell::genRanNodeReduced], rRatio = " << rRatio << " is greater than 1.0, do nothing and leave" << std::endl;
        return;
    }
    unsigned maxRej = nMaxRej*sqrt(2.0)*nx*ny*nz*(1-rRatio)*factor;

    Point rTest;
    if (Use_RandomNode) {
        std::srand (std::time(NULL));
    }
    *p_NodeInPartition = getNodeInPartition(lmin);

    while (nodeCnt < targetNum) {
        double ran = 2.0*(((double) std::rand())/RAND_MAX-0.5);
        if (dir==2) {
            rTest[0] = (((double) std::rand())/RAND_MAX*nx+ex)*UnitLength;
            rTest[1] = (((double) std::rand())/RAND_MAX*ny+ey)*UnitLength;

            if (ran>0.0) {
                rTest[2] = mid + nz*UnitLength*(ran*(1-rRatio)+rRatio)/2.0;
            } else {
                rTest[2] = mid + nz*UnitLength*(ran*(1-rRatio)-rRatio)/2.0;
            }
        } else if (dir==1) {
            rTest[0] = (((double) std::rand())/RAND_MAX*nx+ex)*UnitLength;
            rTest[2] = (((double) std::rand())/RAND_MAX*nz+ez)*UnitLength;
            if (ran>0.0) {
                rTest[1] = mid + ny*UnitLength*(ran*(1-rRatio)+rRatio)/2.0;
            } else {
                rTest[1] = mid + ny*UnitLength*(ran*(1-rRatio)-rRatio)/2.0;
            }
        } else if (dir==0) {
            rTest[1] = (((double) std::rand())/RAND_MAX*ny+ey)*UnitLength;
            rTest[2] = (((double) std::rand())/RAND_MAX*nz+ez)*UnitLength;
            if (ran>0.0) {
                rTest[0] = mid + nx*UnitLength*(ran*(1-rRatio)+rRatio)/2.0;
            } else {
                rTest[0] = mid + nx*UnitLength*(ran*(1-rRatio)-rRatio)/2.0;
            }
        }
        std::array<unsigned,Dim>    p;
        p[0] = (std::fabs(rTest[0]/lmin-Tiny));
        p[1] = (std::fabs(rTest[1]/lmin-Tiny));
        p[2] = (std::fabs(rTest[2]/lmin-Tiny));
        if ((testOverlap(p_NodeInPartition,lmin,rTest,*p_NodeID,p)==0)) {
            nodeLists.free.push_back(*p_NodeID);
            putNode(p_conp,p_NodeID,rTest,Type_0);
            nodeCnt++;
            rejectCnt = 0;
        } else {
            rejectCnt++;
            if (rejectCnt>maxRej) {
                dualOut << "Cannot place a node, force quit the loop"<< std::endl;
                dualOut << "Modified n (packing ratio relative to HCF packing) = "
                        << PackingDensity*((double) nodeCnt)/(targetNum) << std::endl;
                if (MaxTotalStep==0) {
                    char fileName[255];
                    std::sprintf(fileName,"%s/%snLog.vtk",OutputSubFolder,OutputFilePrefix);
                    std::ofstream file(fileName, std::ios::out | std::ios::app);
                    file <<  nMaxRej  << '\t' << ((double) nodeCnt)/(std::sqrt(2.0)*Nx*Ny*Nz*(1.0-rRatio)*factor) << std::endl;
                }
                break;
            }
        }
    }
    clock.get();
}

Mat3Dvec    VoroTesell::getNodeInPartition  (    const double        lmin)
{
    Mat3Dvec            NodeInPartition;
    NodeInPartition.resize(Nx*UnitLength/lmin+1,Ny*UnitLength/lmin+1,Nz*UnitLength/lmin+1);
    for (unsigned NodeID =0; NodeID < GNodeTab.size(); NodeID++) {
        std::array<unsigned,3>    p;
        for (unsigned d=0; d<3; d++) {
            p[d] = (GNodeTab[NodeID].coord[d]/lmin-Tiny);
        }
        NodeInPartition.push_back(p,NodeID);
    }
    return NodeInPartition;
}

void 	VoroTesell::genRanNodeAnisotropic					            (	unsigned* 			p_NodeID,
																			container_poly* 	p_conp,
																			double				aniZone)
{
	unsigned nodeCnt = *p_NodeID;
	unsigned rejectCnt = 0;
	unsigned maxRej = 10.0*sqrt(2)*Nx*Ny*Nz;
	std::array<unsigned,Dim> 	p;
	Point rTest;
	dualOut << "Generating random node anisotropic , anisotropy = " << anisotropy << std::endl;
	if (anisotropy*aniZone > 1.0) {
		aniZone = 1.0 / anisotropy;
		tripOut << "[VoroTesell::genRanNodeAnisotropic] anisotropy*aniZone > 1.0, set anisotropy = 1.0/aniZone  = "
				<< aniZone << '\n';
	}
	double deltaz = Nz*UnitLength*aniZone/2.0;
	double deltaz2 = deltaz*anisotropy;
	Point centre;
	centre[0] = UnitLength*Nx*0.5;
	centre[1] = UnitLength*Ny*0.5;
	centre[2] = UnitLength*Nz*0.5;
	for (unsigned NodeID=crackNodeNum; NodeID<tInNodeNum; NodeID++) {
		double zOffset = GNodeTab[NodeID].coord[2]-centre[2];
		if (std::fabs(zOffset)<deltaz) {
		    GNodeTab[NodeID].coord[2] = centre[2]+zOffset/anisotropy;
		} else if (std::fabs(zOffset)<deltaz2) {
			nodeLists.free.push_back(*p_NodeID);
			putNode(p_conp,p_NodeID,GNodeTab[NodeID].coord[0],GNodeTab[NodeID].coord[1],centre[2]+zOffset/anisotropy,Type_0);
		}
	}
	return;
	dualOut << "Finish generating random node anisotropic" << std::endl;
}

void	VoroTesell::genSimpleCubicNodes						(	container_poly* 	p_conp,
																unsigned* 			p_NodeID)
{
	Point inPoint;
	for (unsigned i=0; i<Nx; i++) {

		for (unsigned j=0; j<Ny; j++) {

			for (unsigned k=0; k<Nz; k++) {
				inPoint[0] = (0.5+i)*UnitLength;
				inPoint[1] = (0.5+j)*UnitLength;
				inPoint[2] = (0.5+k)*UnitLength;
				nodeLists.free.push_back(*p_NodeID);
				putNode(p_conp,p_NodeID,inPoint,Type_0);
			}
		}
	}
}

void	VoroTesell::genBCC_Nodes							(	container_poly* 	p_conp,
																unsigned* 			p_NodeID)
{
	Point inPoint;
	double a = 2.0/std::sqrt(3.0);
	double a2 = a/2.0;
	unsigned nx = Nx / a;
	unsigned ny = Ny / a;
	unsigned nz = Nz / a ;
	double offsetx = (Nx - nx*a)/2.0;
	double offsety = (Ny - ny*a)/2.0;
	double offsetz = (Nz - nz*a)/2.0;
	std::cout << a << ' '<< nx << ' ' << ny << ' ' << nz << ' ' << offsetx << ' ' << offsety << ' ' << offsetz << std::endl;
	for (unsigned i=0; i<2*nx+1; i++) {
		inPoint[0] = (offsetx+a2*i)*UnitLength;
		for (unsigned j=0; j<ny-i%2+1; j++) {
			inPoint[1] = (offsety+a*j+0.5*a*(i%2))*UnitLength;
			for (unsigned k=0; k<nz-i%2+1; k++) {
				inPoint[2] = (offsetz+a*k+0.5*(i%2)*a)*UnitLength;
				nodeLists.free.push_back(*p_NodeID);
				putNode(p_conp,p_NodeID,inPoint,Type_0);
			}
		}
	}
}

void	VoroTesell::genFCC_Nodes							(	container_poly* 	p_conp,
																unsigned* 			p_NodeID)
{
	Point inPoint;
	double a = std::sqrt(2.0);
	double a2 = a/2.0;
	unsigned nx = Nx / a;
	unsigned ny = Ny / a;
	unsigned nz = Nz / a ;
	double offsetx = (Nx - nx*a)/2.0;
	double offsety = (Ny - ny*a)/2.0;
	double offsetz = (Nz - nz*a)/2.0;
	std::cout << a << ' '<< nx << ' ' << ny << ' ' << nz << ' ' << offsetx << ' ' << offsety << ' ' << offsetz << std::endl;
		for (unsigned i=0; i<2*nx+1; i++) {
			inPoint[0] = (offsetx+a2*i)*UnitLength;
			for (unsigned j=0; j<ny + 1 + ny*(i%2); j++) {
				for (unsigned k=0; k<nz-((j+1)*i)%2+1; k++) {
					inPoint[1] = (offsety+a2*(1+(i+1)%2)*j)*UnitLength;
					inPoint[2] = (offsetz+a*k+0.5*(((j+1)*i)%2)*a)*UnitLength;
					nodeLists.free.push_back(*p_NodeID);
					putNode(p_conp,p_NodeID,inPoint,Type_0);
				}
			}
		}
}

void    VoroTesell::genHCP_Nodes                            (   container_poly*     p_conp,
                                                                unsigned*           p_NodeID)
{
    Point inPoint;
    double dx = UnitLength;
    double dy = std::sqrt(3.0)/2.0*UnitLength;
    double dz = std::sqrt(6.0)/3.0*UnitLength;
    unsigned nz = Nz*UnitLength /dz ;
    unsigned ny = Ny*UnitLength /dy;
    unsigned nx = Nx-1 ;

    double offsetx = (Nx*UnitLength - nx*dx)/2.0;
    double offsety = (Ny*UnitLength - ny*dy)/2.0;
    double offsetz = (Nz*UnitLength - nz*dz)/2.0;
    std::cout << nx << ' ' << ny << ' ' << nz << ' ' << offsetx << ' ' << offsety << ' ' << offsetz << std::endl;
        for (unsigned k=0; k<nz+1; k++) {
            inPoint[2] = offsetz + dz*k;
            for (unsigned j=0; j<ny+1-k%2; j++) {
                inPoint[1] = offsety + dy*j+sqrt(3.0)/6.0*UnitLength*(k%2);
                for (unsigned i=0; i<nx+1-((j+k)%2); i++) {
                    inPoint[0] = dx*i+0.5*UnitLength*((j+k)%2);
                    nodeLists.free.push_back(*p_NodeID);
                    putNode(p_conp,p_NodeID,inPoint,Type_0);
                }
            }
        }
}

void	VoroTesell::genRegBoundaryNodes						(	container_poly* 	p_conp,
																unsigned* 			p_NodeID)
{
	Point inPoint;
	double scalex = scaleX;
	double scaley = scaleY;
	double scalez = scaleZ;

	double offsetX = -((scaleX-1.0)/2.0*Nx)*UnitLength;
	double offsetY = -((scaleY-1.0)/2.0*Ny)*UnitLength;
	double offsetZ = -((scaleZ-1.0)/2.0*Nz)*UnitLength;

	unsigned nx = Nx;
	unsigned ny = Ny;
	unsigned nz = Nz;

	if (nx>4*coarseness) {
	    nx = Nx / coarseness;
	    scalex *= ((double) Nx) / nx;
	}

	if (ny>4*coarseness) {
	   ny = Ny / coarseness;
	   scaley *= ((double) Ny) / ny;
	}

	if (nz>4*coarseness) {
	    nz = Nz / coarseness;
	    scalez *= ((double) Nz) / nz;
	}

	// Top and bottom face
	if (Nz>3) {
        for (unsigned k=0; k<2; k+=1) {
            if (k==0) {
                inPoint[2] = offsetZ;
            } else {
                inPoint[2] = ((scaleZ+1.0)/2.0*Nz)*UnitLength;
            }
            for (unsigned i=0; i<nx+1; i++) {
                inPoint[0] = offsetX + i*scalex*UnitLength;
                for (unsigned j=0; j<ny+1; j++) {
                    inPoint[1] = offsetY + j*scaley*UnitLength;
                    nodeLists.free.push_back(*p_NodeID);
                    putNode(p_conp,p_NodeID,inPoint,Type_0);
                }
            }
        }
	}

	// Left and right faces
	if ((Nx<4)&&(Ny<4)) {
	    return;
	}
        for (unsigned i=0; i<2; i+=1) {
            if (i==0) {
                inPoint[0] = offsetX;
            } else {
                inPoint[0] = ((scaleX+1.0)/2.0*Nx-0.0*scalex)*UnitLength;
            }
        for (unsigned j=0; j<ny+1; j++) {
            inPoint[1] = offsetY + j*scaley*UnitLength;
            for (unsigned k=1; k<nz; k++) {
                inPoint[2] = offsetZ + k*scalez*UnitLength;
                nodeLists.free.push_back(*p_NodeID);
                putNode(p_conp,p_NodeID,inPoint,Type_0);
            }
        }
        }


	//Front and back faces
        for (unsigned j=0; j<2; j+=1) {
            if (j==0) {
                inPoint[1] = offsetY;
            } else {
                inPoint[1] = ((scaleY+1.0)/2.0*Ny-0.0*scaley)*UnitLength;
            }
            for (unsigned i=1; i<ny; i++) {
                inPoint[0] = offsetX + i*scalex*UnitLength;
                for (unsigned k=1; k<nz; k++) {
                    inPoint[2] = offsetZ + k*scalez*UnitLength;
                    nodeLists.free.push_back(*p_NodeID);
                    putNode(p_conp,p_NodeID,inPoint,Type_0);
                }
            }
	}
}

std::vector<NodeIDPair> VoroTesell::genRegNodes     (   container_poly*         p_conp,
                                                        unsigned*               p_NodeID)
{
    dualOut << " VoroTesell::genRegNodes started..." << std::endl;
    Mat3Dvec            NodeInPartition;
    GNodeTab.clear();
    NodeInPartition.resize(Nx*UnitLength/minD,Ny*UnitLength/minD,Nz*UnitLength/minD);
    Point centre;
    centre[0] = UnitLength*Nx*0.5;
    centre[1] = UnitLength*Ny*0.5;
    centre[2] = UnitLength*Nz*0.5;
    std::vector<NodeIDPair> conjNodeList = genPennyCrackReg(centre,CrackRadius*Nx*UnitLength,0.4,0.4*minD,p_NodeID,p_conp,&NodeInPartition);


    double lx = UnitLength;
    double ly = std::sqrt(3.0)/2.0*UnitLength;
    double lz = std::sqrt(2.0/3.0)*UnitLength;
    unsigned nx = Nx;
    unsigned ny = Ny*UnitLength/ly;
    unsigned nz = (Nz-2)*UnitLength/2.0/lz;
    double offx = 0; //UnitLength/2.0;
    double offy = (Ny*UnitLength-ny*ly)/2.0;
    double offz = ((Nz-2)*UnitLength/2.0-nz*lz)/4.0;
    double upz = (Nz+2)*UnitLength/2.0;

    for (unsigned i=0; i<2; i++) {
        for (unsigned iz=0; iz< nz+1; iz++) {
            double z = offz + iz*lz + upz*(i%2);
            for (unsigned iy=0; iy< ny-iz%2+1; iy++) {
                double y = offy + iz%2*1.0/3.0*ly + iy*ly;
                for (unsigned ix=0; ix< nx-(iy+iz)%2+1; ix++) {
                    double x = offx + (iy+iz)%2*0.5*lx + ix*lx;
                    putNode(p_conp,p_NodeID,x,y,z,Type_0);
                }
            }
        }
    }
    return conjNodeList;
}

std::vector<NodeIDPair> VoroTesell::genPennyCrack	(	const Point				centre,
														const double			radius,
														unsigned* 				p_NodeID,
														container_poly* 		p_conp,
														Mat3Dvec*				p_NodeInPartition)
{
	dualOut << "VoroTesell::genPennyCrack starts...." << std::endl;
	Geometry					geo;
	std::vector<NodeIDPair>		conjNodeList;
	if (radius < 2.0/((double) std::min(std::min(Nx,Ny),Nz))) {
		dualOut << "Radius entered " << radius << " is too small, no circular crack is created " << std::endl;
		return conjNodeList;
	}
	NodeIDPair					conjugateNode(0,0);
	const double				densify = 1.5;
	double                      n = PackingDensity;
	if (Use_const_NodeNum) {
	    n = 0.5;
	}
	unsigned 					Num = (unsigned) densify*densify*2.5*(radius*radius*Pi)*
										sqrt(4.0/3.0)*n*(UnitLength/minD)*(UnitLength/minD);
	Point 						rTest;
	std::array<unsigned,Dim> 	p;
	const double 				delta = minD*UnitLength/10.0;
	std::vector<Point>			nList;
	unsigned					rejectCnt = 0;
	unsigned					maxRej = std::max(Num/4,(unsigned) 20);
	rTest[2] = centre[2] + delta/2.0;
	p[2] = ((unsigned) rTest[2]/minD-Tiny);
	unsigned i=0;
	dualOut << "Target no. of Penny Crack node = " << Num << std::endl;
	while (i < Num) {
		do {
			rTest[0] = (((double) rand())*2.0/RAND_MAX-1.0)*radius+centre[0];
			rTest[1] = (((double) rand())*2.0/RAND_MAX-1.0)*radius+centre[1];
		} while (geo.dist(rTest,centre)>radius);
		p[0] = ( rTest[0]/UnitLength/minD-Tiny);
		p[1] = ( rTest[1]/UnitLength/minD-Tiny);
		if (testOverlap(p_NodeInPartition,minD,rTest,*p_NodeID,p)==0) {
			nodeLists.free.push_back(*p_NodeID);
			conjugateNode.i1 = *p_NodeID;
			conjNodeList.push_back(conjugateNode);
			putNode(p_conp,p_NodeID,rTest,Type_0);
			nList.push_back(rTest);
			i++;
			rejectCnt = 0;
		} else {
			rejectCnt++;
			if (rejectCnt>maxRej) {
				dualOut << "Cannot place a node, force quit the loop" << std::endl;
				break;
			}
		}
	}
	for (unsigned i = 0; i< nList.size(); i++) {
		nodeLists.free.push_back(*p_NodeID);
		conjNodeList[i].i2 = *p_NodeID;
		rTest[0] = (((double) rand())*2.0/RAND_MAX-1.0)*delta/5.0;
		rTest[1] = (((double) rand())*2.0/RAND_MAX-1.0)*delta/5.0;
		rTest[2] = (((double) rand())*2.0/RAND_MAX-1.0)*delta/5.0-delta;
		putNode(p_conp,p_NodeID,nList[i][0]+rTest[0],nList[i][1]+rTest[1],nList[i][2]+rTest[2],Type_0);
	}
	dualOut << "VoroTesell::genPennyCrack finishes...." << std::endl;
	return conjNodeList;
}

std::vector<NodeIDPair> VoroTesell::genPennyCrackReg	(	const Point				centre,
															const double			radius,
															const double            ld,
															const double            lt,
															unsigned* 				p_NodeID,
															container_poly* 		p_conp,
															Mat3Dvec*				p_NodeInPartition)
{
	dualOut << "VoroTesell::genPennyCrack regular starts...." << std::endl;
	std::vector<NodeIDPair>		conjNodeList;
	if (radius < 2.0/((double) std::min(std::min(Nx,Ny),Nz))) {
		dualOut << "Radius entered " << radius << " is too small, no circular crack is created " << std::endl;
		return conjNodeList;
	}
	NodeIDPair					conjugateNode(0,0);
	std::vector<Point>			nList;
	const double				radius_m = radius-0.5*ld;
	const unsigned				rNum = (unsigned) (radius_m/ld);
	const double				deltaR = radius_m/rNum;
	const double 				delta = ld/10.0;
	//Put the centre
	Point                       inPoint,ran;
	std::array<unsigned,Dim>    p;
	inPoint[0] = centre[0];
	inPoint[1] = centre[1];
	inPoint[2] = centre[2] + lt/2.0;
	p[0] = ((unsigned) inPoint[0]/minD-Tiny);
	p[1] = ((unsigned) inPoint[1]/minD-Tiny);
	p[2] = ((unsigned) inPoint[2]/minD-Tiny);
	p_NodeInPartition -> push_back(p[0],p[1],p[2],*p_NodeID);
	nodeLists.free.push_back(*p_NodeID);
	conjugateNode.i1 = *p_NodeID;
	conjNodeList.push_back(conjugateNode);
	putNode(p_conp,p_NodeID,inPoint,Type_0);
	nList.push_back(inPoint);
	pxPerimeterNum[0] = pxPerimeterNum[1] = 0;
	std::vector<double>	ranZ;
	ranZ.push_back(lt/2.0);
	for (unsigned i=1; i<rNum+1; i++) {
		double r = deltaR*i;
		unsigned phiNum = (unsigned) (2.0*Pi*r/ld);
		double deltaPhi = 2.0*Pi/phiNum;
		double halfDelta = 0.5*deltaPhi*(i%2);
		for (unsigned j=0; j<phiNum; j++) {
			ran[0] = (((double) rand())*2.0/RAND_MAX-1.0)*delta*2.0;
			ran[1] = (((double) rand())*2.0/RAND_MAX-1.0)*delta*2.0;
			inPoint[0] = centre[0] + r*std::cos(deltaPhi*j+halfDelta)+ran[0];
			inPoint[1] = centre[1] + r*std::sin(deltaPhi*j+halfDelta)+ran[1];
			inPoint[2] = centre[2] + lt/2.0;
			ranZ.push_back(lt/2.0);
			p[0] = ((unsigned) inPoint[0]/minD-Tiny);
			p[1] = ((unsigned) inPoint[1]/minD-Tiny);
			p[2] = ((unsigned) inPoint[2]/minD-Tiny);
			p_NodeInPartition -> push_back(p[0],p[1],p[2],*p_NodeID);
			nodeLists.free.push_back(*p_NodeID);
			conjugateNode.i1 = *p_NodeID;
			conjNodeList.push_back(conjugateNode);
			putNode(p_conp,p_NodeID,inPoint,Type_0);
			nList.push_back(inPoint);
			if (i==rNum-1) {
				pxPerimeterNum[0]++;
			} else if (i==rNum-2) {
				pxPerimeterNum[1]++;
			}
		}
	}
	for (unsigned i = 0; i< nList.size(); i++) {
		nodeLists.free.push_back(*p_NodeID);
		conjNodeList[i].i2 = *p_NodeID;
		putNode(p_conp,p_NodeID,nList[i][0],nList[i][1],(nList[i][2]-2.0*ranZ[i]),Type_0);
	}
	dualOut << "VoroTesell::genPennyCrack regular finishes...." << std::endl;
	crackNodeNum = *p_NodeID;

	return conjNodeList;
}

std::vector<NodeIDPair> VoroTesell::genPennyCrackRan    (   const Point             centre,
                                                            const double            radius,
                                                            const double            ld,
                                                            const double            lt,
                                                            unsigned*               p_NodeID,
                                                            container_poly*         p_conp,
                                                            Mat3Dvec*               p_NodeInPartition)
{
    dualOut << "VoroTesell::genPennyCrack random starts...." << std::endl;
    std::vector<NodeIDPair>     conjNodeList;
    if (radius < 2.0/((double) std::min(std::min(Nx,Ny),Nz))) {
        dualOut << "Radius entered " << radius << " is too small, no circular crack is created " << std::endl;
        return conjNodeList;
    }
    NodeIDPair                  conjugateNode(0,0);
    std::vector<Point>          nList;

    //Put the centre
    Point                       inPoint;
    std::array<unsigned,Dim>    p;
    inPoint[0] = centre[0];
    inPoint[1] = centre[1];
    inPoint[2] = centre[2] + lt/2.0;
    p[0] = ((unsigned) inPoint[0]/minD-Tiny);
    p[1] = ((unsigned) inPoint[1]/minD-Tiny);
    p[2] = ((unsigned) inPoint[2]/minD-Tiny);
    p_NodeInPartition -> push_back(p[0],p[1],p[2],*p_NodeID);
    nodeLists.free.push_back(*p_NodeID);
    conjugateNode.i1 = *p_NodeID;
    conjNodeList.push_back(conjugateNode);
    putNode(p_conp,p_NodeID,inPoint,Type_0);
    nList.push_back(inPoint);
    //Put the remaining nodes
    pxPerimeterNum[0] = pxPerimeterNum[1] = 0;
    std::vector<double> ranZ;
    ranZ.push_back(lt/2.0);
    unsigned num_penny = Pi*radius*radius/(ld*ld);

    unsigned i=0,rejectCnt=0;
    double r,theta;
    double rMax = radius - ld/2.0;
    Mat2Dvec       Partition2D;
    Partition2D.resize(Nx*UnitLength/minD,Ny*UnitLength/minD);
    std::array<unsigned,2>  pp;
    unsigned maxRej = 10*num_penny;
    std::cout << num_penny << ' ' << maxRej << std::endl;
    while (i<num_penny) {
        r = ((double) rand())/RAND_MAX*rMax;
        theta = ((double) rand())/RAND_MAX*2*Pi;
        inPoint[0] = centre[0] + r*std::cos(theta);
        inPoint[1] = centre[1] + r*std::sin(theta);

        if (isInRange(inPoint)) {
            p[0] = pp[0] = ((unsigned) inPoint[0]/minD-Tiny);
            p[1] = pp[1] = ((unsigned) inPoint[1]/minD-Tiny);
            p[2] = ((unsigned) inPoint[2]/minD-Tiny);
            if (testOverlap(&Partition2D,inPoint,*p_NodeID,pp)==0) {
                ranZ.push_back(lt/2.0);
                p_NodeInPartition -> push_back(p[0],p[1],p[2],*p_NodeID);
                nodeLists.free.push_back(*p_NodeID);
                conjugateNode.i1 = *p_NodeID;
                conjNodeList.push_back(conjugateNode);
                putNode(p_conp,p_NodeID,inPoint,Type_0);
                nList.push_back(inPoint);
                i++;
                rejectCnt = 0;
            } else {
                rejectCnt++;
                if (rejectCnt>maxRej) {
                    dualOut << "Cannot place a node, force quit the loop" << std::endl;
                    break;
                }
        }
        }
    }
    for (unsigned i = 0; i< nList.size(); i++) {
        nodeLists.free.push_back(*p_NodeID);
        conjNodeList[i].i2 = *p_NodeID;
        putNode(p_conp,p_NodeID,nList[i][0],nList[i][1],(nList[i][2]-2.0*ranZ[i]),Type_0);
    }
    crackNodeNum = *p_NodeID;
    std::cout << "crackNodeNum = " << crackNodeNum << std::endl;
    nList.clear();
    ranZ.clear();
    i = 0; rejectCnt = 0;
    unsigned num_outerPenny = (Nx*Ny*UnitLength*UnitLength-Pi*radius*radius)/(ld*ld);
    maxRej = num_outerPenny/2;
    std::cout << num_outerPenny << ' ' << maxRej << std::endl;
    double dist;
    while (i<num_outerPenny) {
        inPoint[0] = ((double) rand())/RAND_MAX*Nx*UnitLength;
        inPoint[1] = ((double) rand())/RAND_MAX*Ny*UnitLength;
        dist = std::sqrt((inPoint[0]-centre[0])*(inPoint[0]-centre[0])+
                (inPoint[1]-centre[1])*(inPoint[1]-centre[1]));
        if (dist>radius) {
            p[0] = pp[0] = ((unsigned) inPoint[0]/minD-Tiny);
            p[1] = pp[1] = ((unsigned) inPoint[1]/minD-Tiny);
            p[2] = ((unsigned) inPoint[2]/minD-Tiny);
            if (testOverlap(&Partition2D,inPoint,*p_NodeID,pp)==0) {
                ranZ.push_back(lt/2.0);
                p_NodeInPartition -> push_back(p[0],p[1],p[2],*p_NodeID);
                nodeLists.free.push_back(*p_NodeID);
                conjugateNode.i1 = *p_NodeID;
                weakPlane.push_back(conjugateNode);
                putNode(p_conp,p_NodeID,inPoint,Type_0);
                nList.push_back(inPoint);
                i++;
                rejectCnt = 0;
            } else {
                rejectCnt++;
                if (rejectCnt>maxRej) {
                    dualOut << "Cannot place a node, force quit the loop" << std::endl;
                    break;
                }
            }
        }
    }
    for (unsigned i = 0; i< nList.size(); i++) {
        nodeLists.free.push_back(*p_NodeID);
        weakPlane[i].i2 = *p_NodeID;
        putNode(p_conp,p_NodeID,nList[i][0],nList[i][1],(nList[i][2]-2.0*ranZ[i]),Type_0);
    }
    outerCrackNodeNum = *p_NodeID - crackNodeNum;
    std::cout << "outerCrackNodeNum = " << outerCrackNodeNum << std::endl;
    dualOut << "VoroTesell::genPennyCrack random finishes...." << std::endl;
    return conjNodeList;
}


bool VoroTesell::genInclinedPennyCrackRan                       (   const Point                 centre,
                                                                    const double                radius,
                                                                    Vec                         normVec,
                                                                    const double                ld,
                                                                    const double                lt,
                                                                    unsigned*                   p_NodeID,
                                                                    std::vector<NodeIDPair>*    p_conjNodeList,
                                                                    container_poly*             p_conp,
                                                                    Mat3Dvec*                   p_NodeInPartition)
{
    dualOut << "VoroTesell::genPennyCrack random starts...." << std::endl;
    Geometry geo;
    double mag = geo.mag(normVec);
    normVec[0] /= mag;      normVec[1] /= mag;      normVec[2] /= mag;
    if (radius < 2.0/((double) std::min(std::min(Nx,Ny),Nz))) {
        dualOut << "Radius entered " << radius << " is too small, no circular crack is created " << std::endl;
        return false;
    }
    NodeIDPair                  conjugateNode(0,0);
    std::vector<Point>          nList;

    Point                       inPoint;
    std::array<unsigned,Dim>    p;
    inPoint[0] = centre[0] + lt*normVec[0]/2.0;
    inPoint[1] = centre[1] + lt*normVec[1]/2.0;
    inPoint[2] = centre[2] + lt*normVec[2]/2.0;
    p[0] = ( inPoint[0]/minD-Tiny);
    p[1] = ( inPoint[1]/minD-Tiny);
    p[2] = ( inPoint[2]/minD-Tiny);
    p_NodeInPartition -> push_back(p[0],p[1],p[2],*p_NodeID);
    nodeLists.free.push_back(*p_NodeID);
    conjugateNode.i1 = *p_NodeID;
    p_conjNodeList->push_back(conjugateNode);
    putNode(p_conp,p_NodeID,inPoint,Type_0);
    nList.push_back(inPoint);
    //Put the remaining nodes
    pxPerimeterNum[0] = pxPerimeterNum[1] = 0;
    unsigned num_penny = Pi*radius*radius/(ld*ld)*anisotropy*anisotropy*0.849;

    unsigned i=0,rejectCnt=0;
    double r,theta;
    double rMax = radius - ld/2.0;
    Mat2Dvec       Partition2D;
    Partition2D.resize(Nx*UnitLength/minD,Ny*UnitLength/minD);
    unsigned maxRej = 10000*num_penny;
    std::cout << num_penny << ' ' << maxRej << std::endl;
    while (i<num_penny) {
        r = ((double) rand())/RAND_MAX*rMax;
        theta = ((double) rand())/RAND_MAX*2*Pi;
        inPoint = geo.rotateVec(r*std::cos(theta),r*std::sin(theta),normVec);
        inPoint[0] += centre[0] + lt*normVec[0]/2.0;
        inPoint[1] += centre[1] + lt*normVec[1]/2.0;
        inPoint[2] += centre[2] + lt*normVec[2]/2.0;
        if (isInRange(inPoint)){
            p[0] = ((unsigned) inPoint[0]/minD-Tiny);
            p[1] = ((unsigned) inPoint[1]/minD-Tiny);
            p[2] = ((unsigned) inPoint[2]/minD-Tiny);
            if (testOverlap(p_NodeInPartition,minD/anisotropy, inPoint,*p_NodeID,p)==0) {
                nodeLists.free.push_back(*p_NodeID);
                conjugateNode.i1 = *p_NodeID;
                p_conjNodeList->push_back(conjugateNode);
                putNode(p_conp,p_NodeID,inPoint,Type_0);
                nList.push_back(inPoint);
                i++;
                rejectCnt = 0;
            } else {
                rejectCnt++;
                if (rejectCnt>maxRej) {
                    dualOut << "Cannot place a node, force quit the loop" << std::endl;
                    break;
                }
            }
        }
    }
    dualOut << "VoroTesell::genPennyCrack random finishes...." << std::endl;
    return true;
}

void VoroTesell::genIncOuterPennyCrackRan                       (   const Point             centre,
                                                                    const double            inRadius,
                                                                    const double            outRadius,
                                                                    Vec                     normVec,
                                                                    const double            ld,
                                                                    const double            lt,
                                                                    unsigned*               p_NodeID,
                                                                    container_poly*         p_conp,
                                                                    Mat3Dvec*               p_NodeInPartition)
{
    dualOut << "VoroTesell::genPennyCrack random starts...." << std::endl;
    Geometry geo;
    double mag = geo.mag(normVec);
    normVec[0] /= mag;      normVec[1] /= mag;      normVec[2] /= mag;
    NodeIDPair                  conjugateNode(0,0);
    std::vector<Point>          nList;

    unsigned i = 0;
    unsigned rejectCnt = 0;
    unsigned num_outerPenny = (Nx*Ny*UnitLength*UnitLength-Pi*inRadius*inRadius)*anisotropy*anisotropy/(ld*ld)*0.849;
    unsigned maxRej = 10*num_outerPenny;
    double r,theta;
    Point inPoint;
    std::array<unsigned,Dim>    p;
    std::cout << num_outerPenny << ' ' << maxRej << std::endl;
    while (i<num_outerPenny) {
        r = ((double) rand())/RAND_MAX*(outRadius-inRadius)+inRadius;
        theta = ((double) rand())/RAND_MAX*2*Pi;
        inPoint = geo.rotateVec(r*std::cos(theta),r*std::sin(theta),normVec);
        inPoint[0] += centre[0] + lt*normVec[0]/2.0;
        inPoint[1] += centre[1] + lt*normVec[1]/2.0;
        inPoint[2] += centre[2] + lt*normVec[2]/2.0;
        if (isInRange(inPoint)) {
            p[0] = ( inPoint[0]/minD-Tiny);
            p[1] = ( inPoint[1]/minD-Tiny);
            p[2] = ( inPoint[2]/minD-Tiny);
            if (testOverlap(p_NodeInPartition,minD,inPoint,*p_NodeID,p)==0) {
                nodeLists.free.push_back(*p_NodeID);
                conjugateNode.i1 = *p_NodeID;
                weakPlane.push_back(conjugateNode);
                putNode(p_conp,p_NodeID,inPoint,Type_0);
                nList.push_back(inPoint);
                i++;
                rejectCnt = 0;
            } else {
                rejectCnt++;
                if (rejectCnt>maxRej) {
                    dualOut << "Cannot place a node, force quit the loop" << std::endl;
                    break;
                }
            }
        }
    }
    for (unsigned i = 0; i< nList.size(); i++) {
        nodeLists.free.push_back(*p_NodeID);
        weakPlane[i].i2 = *p_NodeID;
        putNode(p_conp,p_NodeID,nList[i][0]-lt*normVec[0],
                nList[i][1]-lt*normVec[1],(nList[i][2]-lt*normVec[2]),Type_0);
    }
    outerCrackNodeNum = *p_NodeID - crackNodeNum;
    dualOut << "outerCrackNodeNum" << ' ' << outerCrackNodeNum << std::endl;
    dualOut << "VoroTesell::genPennyCrack random finishes...." << std::endl;
}

void VoroTesell::genNotchWeakPlane                              (   const Point             centre,
                                                                    const double            notchLen,
                                                                    const double            lt,
                                                                    unsigned*               p_NodeID,
                                                                    container_poly*         p_conp,
                                                                    Mat3Dvec*               p_NodeInPartition)
{
    dualOut << "VoroTesell::genNotchWeakPlane random starts...." << std::endl;
    Geometry geo;
    NodeIDPair                  conjugateNode(0,0);
    std::vector<Point>          nList;

    unsigned i = 0;
    unsigned rejectCnt = 0;
    unsigned num_outerPenny = (Nx*Nz*UnitLength*UnitLength-Nz*UnitLength*notchLen)/(avgD*avgD)*std::sqrt(2.0);
    unsigned maxRej = 10*num_outerPenny;
    double r,theta;
    double outRadius = std::sqrt(Nx*Nx+Nz*Nz)*UnitLength/2.0+avgD;
    Point inPoint;
    std::array<unsigned,Dim>    p;
    std::cout << num_outerPenny << ' ' << maxRej << std::endl;
    inPoint[1] = centre[1] + lt/2.0;
    double xMax = centre[0] + notchLen;
    double xMin = centre[0] - notchLen;
    unsigned iStart = weakPlane.size();
    while (i<num_outerPenny) {
        r = ((double) rand())/RAND_MAX*outRadius;
        theta = ((double) rand())/RAND_MAX*2*Pi;
        inPoint[0] = centre[0] + r*std::cos(theta);
        inPoint[2] = centre[2] + r*std::sin(theta);
        if ((isInRange(inPoint))&&((inPoint[0]>xMax)||(inPoint[0]<xMin))) {
            p[0] = ( inPoint[0]/minD-Tiny);
            p[1] = ( inPoint[1]/minD-Tiny);
            p[2] = ( inPoint[2]/minD-Tiny);
            if (testOverlap(p_NodeInPartition,minD,inPoint,*p_NodeID,p)==0) {
                nodeLists.free.push_back(*p_NodeID);
                conjugateNode.i1 = *p_NodeID;
                weakPlane.push_back(conjugateNode);
                putNode(p_conp,p_NodeID,inPoint,Type_0);
                nList.push_back(inPoint);
                i++;
                rejectCnt = 0;
            } else {
                rejectCnt++;
                if (rejectCnt>maxRej) {
                    dualOut << "Cannot place a node, force quit the loop" << std::endl;
                    break;
                }
            }
        }
    }
    for (unsigned i = 0; i< nList.size(); i++) {
        nodeLists.free.push_back(*p_NodeID);
        weakPlane[i+iStart].i2 = *p_NodeID;
        putNode(p_conp,p_NodeID,nList[i][0],nList[i][1]-lt,nList[i][2],Type_0);
    }
    outerCrackNodeNum = *p_NodeID - crackNodeNum;
    dualOut << "outerCrackNodeNum" << ' ' << outerCrackNodeNum << std::endl;
    dualOut << "VoroTesell::genNotchWeakPlane random finishes...." << std::endl;
}

void VoroTesell::genRoughPennyCrack						(	const Point				centre,
															const double			crackRadius,
															const double            lt)
{
	Geometry geo;
	if (crackRadius<2.0*avgD/(std::min(Nx,Ny)*UnitLength)) {
	    tripOut << "[VoroTesell::genRoughPennyCrack] The crackRadius is too small, no rough penny crack is generated!"
	            << " r in terms of avg facet width = " << crackRadius /(avgD/(std::min(Nx,Ny)*UnitLength)) << std::endl;
	    return;
	}
	dualOut << "Generating rough Penny crack..." << std::endl;
	for (unsigned NodeID=0; NodeID < GNodeTab.size(); NodeID++) {
		double dist = geo.dist(centre,GNodeTab[NodeID].coord);
		if (dist<crackRadius) {
			double t = GNodeTab[NodeID].coord[2] - centre[2];
			if ((t>0.0)&&(t<lt)) {
				for (unsigned k=0; k<GNodeTab[NodeID].nbNodeID.size(); k++) {
					unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[k];
						double t2 = GNodeTab[nbNodeID].coord[2] - centre[2];
						if (t2<0.0) {
							preExistingCrack.push_back(GNodeTab[NodeID].nbLatticeID[k]);
//						}
					}
				}
			}
		}
	}
}
void VoroTesell::genOuterPennyCrack						(	const Point				centre,
															const double			crackRadius,
															const double			outRadius,
															const double            ld,
															const double            lt,
															unsigned* 				p_NodeID,
															container_poly* 		p_conp,
															Mat3Dvec*				p_NodeInPartition)
{
	dualOut << "VoroTesell::genOuterPennyCrack starts...." << std::endl;
	NodeIDPair					conjugateNode(0,0);
	if (crackRadius < 2.0/((double) std::min(std::min(Nx,Ny),Nz))) {
		dualOut << "Radius entered " << crackRadius << " is too small, no potential crack is created " << std::endl;
		return;
	}
	std::vector<Point>			nList;
    avgD = std::pow(Nx*Ny*Nz/Nnum,0.3333333)*UnitLength;
	const double				radius_m = outRadius - crackRadius +0.5*ld - 0.5*avgD;
	const unsigned				rNum = (unsigned) (radius_m/ld);
	const double				deltaR = radius_m/rNum;
	const double 				delta = ld/10.0;

	if (outRadius < crackRadius + minD) {
		dualOut << "outRadius < crackRadius + minD, no, no potential crack is created" << std::endl;
		return;
	}
	weakPlane.clear();
	std::vector<double>	ranZ;
	for (unsigned i=1; i<rNum+1; i++) {
		double r = crackRadius - 0.5*ld + deltaR*i;
		unsigned phiNum = (unsigned) (2.0*Pi*r/ld);
		double deltaPhi = 2.0*Pi/phiNum;
		double halfDelta = 0.5*deltaPhi*(i%2);
		for (unsigned j=0; j<phiNum; j++) {
		    Point  inPoint,ran;
			ran[0] = (((double) rand())*2.0/RAND_MAX-1.0)*delta*2.0;
			ran[1] = (((double) rand())*2.0/RAND_MAX-1.0)*delta*2.0;
//			ran[2] = (((double) rand())*2.0/RAND_MAX-1.0)*deltaZ/4.0;
			inPoint[0] = centre[0] + r*std::cos(deltaPhi*j+halfDelta)+ran[0];
			inPoint[1] = centre[1] + r*std::sin(deltaPhi*j+halfDelta)+ran[1];
			inPoint[2] = centre[2] + lt/2.0;
			if (((inPoint[0]>Tiny)&&(inPoint[0]<Nx*UnitLength-Tiny))&&
			        ((inPoint[1]>Tiny)&&(inPoint[1]<Ny*UnitLength-Tiny))&&
			        ((inPoint[2]>Tiny)&&(inPoint[2]<Nz*UnitLength-Tiny)))
			{
                ranZ.push_back(lt/2.0);
                std::array<unsigned,Dim>    p;
                p[0] = ( inPoint[0]/minD-Tiny);
                p[1] = ( inPoint[1]/minD-Tiny);
                p[2] = ( inPoint[2]/minD-Tiny);
                p_NodeInPartition -> push_back(p[0],p[1],p[2],*p_NodeID);
                nodeLists.free.push_back(*p_NodeID);
                conjugateNode.i1 = *p_NodeID;
                weakPlane.push_back(conjugateNode);
                putNode(p_conp,p_NodeID,inPoint,Type_0);
                nList.push_back(inPoint);
			}
		}
	}
	for (unsigned i = 0; i< nList.size(); i++) {
		nodeLists.free.push_back(*p_NodeID);
		weakPlane[i].i2 = *p_NodeID;
		putNode(p_conp,p_NodeID,nList[i][0],nList[i][1],(nList[i][2]-2.0*ranZ[i]),Type_0);
	}
	outerCrackNodeNum = *p_NodeID - crackNodeNum;
	dualOut << "VoroTesell::genOuterPennyCrack regular finishes...." << std::endl;
	return;
}

void VoroTesell::genNotches                             (   const double                width,
                                                            const double                ld,
                                                            const double                lt,
                                                            unsigned*                   p_NodeID,
                                                            std::vector<NodeIDPair>*    p_conjNodeList,
                                                            container_poly*             p_conp,
                                                            Mat3Dvec*                   p_NodeInPartition)
{
    dualOut << "VoroTesell::genNotches starts...." << std::endl;
    double zCoord = 0.5*Nz*UnitLength - lt/2.0;
    unsigned xNum = (width-0.5*ld)/ld;
    unsigned yNum = Ny*UnitLength/ld;
    double xStep = (width-0.5*ld)/xNum;
    double yStep = Ny*UnitLength/yNum;
    Point                       inPoint;
    NodeIDPair                  conjugateNode(0,0);
    inPoint[2] = zCoord;
    for (unsigned i=0; i<yNum; i++) {
        inPoint[1] = yStep*(0.5+i);
        for (unsigned j=0; j<xNum; j++) {
            inPoint[0] = xStep*(0.5+j);
            conjugateNode.i1 = *p_NodeID;
            p_conjNodeList->push_back(conjugateNode);
            insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
        }
    }

    for (unsigned i=0; i<yNum; i++) {
        inPoint[1] = yStep*(0.5+i);
        for (unsigned j=0; j<xNum; j++) {
            inPoint[0] = Nx*UnitLength - xStep*(0.5+j);
            conjugateNode.i1 = *p_NodeID;
            p_conjNodeList->push_back(conjugateNode);
            insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
        }
    }

}

void VoroTesell::genBorehole                                    (   const Point                 centre,
                                                                    const double                radius,
                                                                    const double                t,
                                                                    const unsigned              resolution,
                                                                    Vec                         normVec,
                                                                    unsigned*                   p_NodeID,
                                                                    container_poly*             p_conp,
                                                                    Mat3Dvec*                   p_NodeInPartition,
                                                                    std::vector<NodeIDPair>*    p_conjNodeList)
{
    dualOut << "VoroTesell::genBorehole starts...." << std::endl;
    dualOut << "radius = " << radius << " t = " << t << " resolution = " << resolution << std::endl;
    Geometry geo;
    double mag = geo.mag(normVec);
    normVec[0] /= mag;      normVec[1] /= mag;      normVec[2] /= mag;
    NodeIDPair                  conjugateNode(0,0);
    NodeIDPair                  conjugateNode2(0,0);
    std::vector<Point>          nList;

    double dTheta = 2*Pi/resolution;

    Point inPoint;
    Point centreOffset;
    unsigned NodeID1 = *p_NodeID;
    unsigned interval = t/minD+1;
    double dt = t/interval;
    std::cout << t << ' ' << interval << ' ' << dt << std::endl;

    for (unsigned n=0; n<interval; n++) {
        unBreakableNodeList.push_back(NodeID1);
        for (unsigned d=0; d<Dim; d++) {
            centreOffset[d] = centre[d] + (n+0.5)*dt*normVec[d];
        }
        insertNode(centreOffset,p_NodeID,p_conp,p_NodeInPartition);
        boreHoleNodeList.push_back(NodeID1);
        for (unsigned i=0; i<resolution; i++) {
            double theta = dTheta*(i+0.5);
            inPoint = geo.rotateVec(2*radius*std::cos(theta),2*radius*std::sin(theta),normVec);
            inPoint[0] += centreOffset[0];
            inPoint[1] += centreOffset[1];
            inPoint[2] += centreOffset[2];
            conjugateNode.i1 = NodeID1;
            conjugateNode.i2 = *p_NodeID;
            conjugateNode2.i1 = *p_NodeID;
            if (i!=resolution-1) {
                conjugateNode2.i2 = *p_NodeID +1;
            } else {
                conjugateNode2.i2 = *p_NodeID -i;
            }
            if (GeoID==4) {
                if ((i!=resolution-1)&&(i!=resolution/2-1)) {
                    unBreakablePair.push_back(conjugateNode2);
                } else {
                    weakPlane.push_back(conjugateNode2);
                }
            }
            p_conjNodeList->push_back(conjugateNode);
            insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
        }
        NodeID1 = *p_NodeID;
    }

    if (GeoID==4) {
        const double notchLen = 8.9;
        double totLen = notchLen - 2.0*radius - minD/2.0;
        unsigned n = totLen/minD +1;
        double len = totLen/n;
        double sgn = 0.0;
        std::cout << "totLen, n, len " << totLen << ' ' << n << ' ' << len << std::endl;
        std::vector<NodeIDPair> tmp;
        for (unsigned i=0; i<weakPlane.size(); i++) {
            if (GNodeTab[weakPlane[i].i1].coord[0]>centre[0]) {
                sgn = 1.0;
            } else {
                sgn = -1.0;
            }
            for (unsigned k=0; k<n; k++) {
                inPoint[0] = GNodeTab[weakPlane[i].i1].coord[0] + sgn*(k+1)*len;
                inPoint[1] = GNodeTab[weakPlane[i].i1].coord[1];
                inPoint[2] = GNodeTab[weakPlane[i].i1].coord[2];
                conjugateNode2.i1 = *p_NodeID;
                insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);

                inPoint[0] = GNodeTab[weakPlane[i].i2].coord[0] + sgn*(k+1)*len;
                inPoint[1] = GNodeTab[weakPlane[i].i2].coord[1];
                inPoint[2] = GNodeTab[weakPlane[i].i2].coord[2];
                conjugateNode2.i2 = *p_NodeID;
                insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
                tmp.push_back(conjugateNode2);
            }
            inPoint[0] = GNodeTab[weakPlane[i].i1].coord[0] + sgn*(n*len+minD);
            inPoint[1] = GNodeTab[weakPlane[i].i1].coord[1];
            inPoint[2] = GNodeTab[weakPlane[i].i1].coord[2];
            insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);

            inPoint[0] = GNodeTab[weakPlane[i].i2].coord[0] + sgn*(n*len+minD);
            inPoint[1] = GNodeTab[weakPlane[i].i2].coord[1];
            inPoint[2] = GNodeTab[weakPlane[i].i2].coord[2];
            insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
        }
        weakPlane.insert(weakPlane.end(),tmp.begin(),tmp.end());
    }
//    unBreakableNodeList.push_back(NodeID1);
    for (unsigned d=0; d<Dim; d++) {
        centreOffset[d] = centre[d] + (t+0.5*dt)*normVec[d];
    }
    if ((centreOffset[0]<Nx*UnitLength)&&(centreOffset[1]<Ny*UnitLength)&&(centreOffset[2]<Nz*UnitLength)) {
        insertNode(centreOffset,p_NodeID,p_conp,p_NodeInPartition);
        for (unsigned i=0; i<resolution; i++) {
            double theta = dTheta*i;
            inPoint = geo.rotateVec(2*radius*std::cos(theta),2*radius*std::sin(theta),normVec);
            inPoint[0] += centreOffset[0];
            inPoint[1] += centreOffset[1];
            inPoint[2] += centreOffset[2];
            insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
        }
    }
    NodeID1 = *p_NodeID;
    for (unsigned d=0; d<Dim; d++) {
        centreOffset[d] = centre[d] - 0.5*dt*normVec[d];
    }
    if ((centreOffset[0]>0.0)&&(centreOffset[1]>0.0)&&(centreOffset[2]>0.0)) {
        insertNode(centreOffset,p_NodeID,p_conp,p_NodeInPartition);
        for (unsigned i=0; i<resolution; i++) {
            double theta = dTheta*i;
            inPoint = geo.rotateVec(2*radius*std::cos(theta),2*radius*std::sin(theta),normVec);
            inPoint[0] += centreOffset[0];
            inPoint[1] += centreOffset[1];
            inPoint[2] += centreOffset[2];
            insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
        }
    }
    BHNodeNum = *p_NodeID;
    dualOut << "VoroTesell::genBorehole finishes...." << std::endl;
}

void VoroTesell::genBorehole2                                   (   const Point                 centre,
                                                                    const double                radius,
                                                                    const double                t,
                                                                    const unsigned              resolution,
                                                                    Vec                         normVec,
                                                                    unsigned*                   p_NodeID,
                                                                    container_poly*             p_conp,
                                                                    Mat3Dvec*                   p_NodeInPartition,
                                                                    std::vector<NodeIDPair>*    p_conjNodeList)
{
    dualOut << "VoroTesell::genBorehole starts...." << std::endl;
    dualOut << "radius = " << radius << " t = " << t << " resolution = " << resolution << std::endl;
    Geometry geo;
    double mag = geo.mag(normVec);
    normVec[0] /= mag;      normVec[1] /= mag;      normVec[2] /= mag;
    NodeIDPair                  conjugateNode(0,0);
    NodeIDPair                  conjugateNode2(0,0);

    double dTheta = 2*Pi/resolution;
    double thk = minD/100.0;

    Point inPoint;
    Point centreOffset;
    unsigned NodeID1 = *p_NodeID;
    unsigned interval = t/minD+1;
    double dt = t/interval;
    std::cout << t << ' ' << interval << ' ' << dt << std::endl;

    for (unsigned n=0; n<interval; n++) {
        unBreakableNodeList.push_back(NodeID1);
        for (unsigned d=0; d<Dim; d++) {
            centreOffset[d] = centre[d] + (n+0.5)*dt*normVec[d];
        }
        insertNode(centreOffset,p_NodeID,p_conp,p_NodeInPartition);
        boreHoleNodeList.push_back(NodeID1);
        for (unsigned i=0; i<resolution; i++) {
            double theta = dTheta*(i+0.5);
            inPoint = geo.rotateVec((radius-thk)*std::cos(theta),(radius-thk)*std::sin(theta),normVec);
            inPoint[0] += centreOffset[0];
            inPoint[1] += centreOffset[1];
            inPoint[2] += centreOffset[2];
            conjugateNode.i1 = *p_NodeID;
            insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);

            inPoint = geo.rotateVec((radius+thk)*std::cos(theta),(radius+thk)*std::sin(theta),normVec);
            inPoint[0] += centreOffset[0];
            inPoint[1] += centreOffset[1];
            inPoint[2] += centreOffset[2];
            conjugateNode.i2 = *p_NodeID;
            conjugateNode2.i1 = *p_NodeID;
            p_conjNodeList->push_back(conjugateNode);
            insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
            if (i!=resolution-1) {
                conjugateNode2.i2 = conjugateNode2.i1 +2;
            } else {
                conjugateNode2.i2 = conjugateNode2.i1 -2*i;
            }
            if (GeoID==4) {
                if ((i!=resolution-1)&&(i!=resolution/2-1)) {
                    unBreakablePair.push_back(conjugateNode2);
                } else {
                    weakPlane.push_back(conjugateNode2);
                }
            }
        }
        NodeID1 = *p_NodeID;
    }

    if (GeoID==4) {
        const double notchLen = 8.9;
        double totLen = notchLen - thk - avgD/2.0;
        unsigned n = totLen/minD +1;
        double len = totLen/n;
        double sgn = 0.0;
        std::cout << "totLen, n, len " << totLen << ' ' << n << ' ' << len << std::endl;
        std::vector<NodeIDPair> tmp;
        for (unsigned i=0; i<weakPlane.size(); i++) {
            if (GNodeTab[weakPlane[i].i1].coord[0]>centre[0]) {
                sgn = 1.0;
            } else {
                sgn = -1.0;
            }
            for (unsigned k=0; k<n; k++) {
                inPoint[0] = GNodeTab[weakPlane[i].i1].coord[0] + sgn*(k+1)*len;
                inPoint[1] = GNodeTab[weakPlane[i].i1].coord[1];
                inPoint[2] = GNodeTab[weakPlane[i].i1].coord[2];
                conjugateNode2.i1 = *p_NodeID;
                insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);

                inPoint[0] = GNodeTab[weakPlane[i].i2].coord[0] + sgn*(k+1)*len;
                inPoint[1] = GNodeTab[weakPlane[i].i2].coord[1];
                inPoint[2] = GNodeTab[weakPlane[i].i2].coord[2];
                conjugateNode2.i2 = *p_NodeID;
                insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
                tmp.push_back(conjugateNode2);
            }
        }
        BHNodeNum = *p_NodeID;
        if (false) {
            for (unsigned i=0; i<weakPlane.size(); i++) {
                if (GNodeTab[weakPlane[i].i1].coord[0]>centre[0]) {
                    sgn = 1.0;
                } else {
                    sgn = -1.0;
                }
                inPoint[0] = GNodeTab[weakPlane[i].i1].coord[0] + sgn*(totLen+minD);
                inPoint[1] = GNodeTab[weakPlane[i].i1].coord[1];
                inPoint[2] = GNodeTab[weakPlane[i].i1].coord[2];
                insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);

                inPoint[0] = GNodeTab[weakPlane[i].i2].coord[0] + sgn*(totLen+minD);
                inPoint[1] = GNodeTab[weakPlane[i].i2].coord[1];
                inPoint[2] = GNodeTab[weakPlane[i].i2].coord[2];
                insertNode(inPoint,p_NodeID,p_conp,p_NodeInPartition);
            }
        }
        weakPlane.insert(weakPlane.end(),tmp.begin(),tmp.end());
    }
    dualOut << "VoroTesell::genBorehole finishes...." << std::endl;
}

bool VoroTesell::insertNode         (   const Point                         inPoint,
                                        unsigned*                           p_NodeID,
                                        container_poly*                     p_conp,
                                        Mat3Dvec*                           p_NodeInPartition)
{
    if (!isInRange(inPoint)) {
        tripOut << "[VoroTesell::insertNode] Point trying to be inserted is out of range!" << std::endl;
        return false;
    }
    std::array<unsigned,Dim>    p;
    p[0] = (inPoint[0]/minD-Tiny);
    p[1] = (inPoint[1]/minD-Tiny);
    p[2] = (inPoint[2]/minD-Tiny);
    p_NodeInPartition -> push_back(p[0],p[1],p[2],*p_NodeID);
    nodeLists.free.push_back(*p_NodeID);
    putNode(p_conp,p_NodeID,inPoint,Type_0);
    return true;
}


unsigned VoroTesell::testOverlap		(	Mat3Dvec* 							p_3Dvec,
                                            const double                        lmin,
                                            const Point							rTest,
                                            const unsigned						NodeID,
                                            const std::array<unsigned,Dim>		p)
{
    int istart,iend,jstart,jend,kstart,kend;
    unsigned n;
    double l;
    std::vector<unsigned>	NodeIDList;
    Geometry			geo;

    std::array<unsigned,3>  dims = p_3Dvec->getDims ();

    istart = jstart = kstart = -1;
    iend = jend = kend = 1;
    if (p[0]==0) {
        istart = 0;
    }
    if (p[0]>= dims[0]-1) {
        istart = 0;
        iend = 0;
    }
    if (p[1]==0){
        jstart = 0;
    }

    if (p[1]>= dims[1]-1) {
        jstart = 0;
        jend = 0;
    }
    if (p[2]==0) {
        kstart = 0;
    }

    if (p[2]>= dims[2]-1) {
        kstart = 0;
        kend = 0;
    }

    for (int i=istart; i<=iend; i++) {
        for (int j=jstart; j<=jend; j++) {
            for (int k=kstart; k<=kend; k++) {
                int flag = (p_3Dvec -> get(p[0]+i,p[1]+j,p[2]+k,&NodeIDList));
                if (flag==1) {
                    for (n=0; n<NodeIDList.size(); n++) {
                        l = geo.dist(rTest,GNodeTab[NodeIDList[n]].coord);
                        if (l<lmin) return 1;
                    }
                } else if (flag<0) {
                    tripOut << "[VoroTesell::testOverlap] : rTest ("
                            << rTest[0] << ',' << rTest[1] << ',' << rTest[2] << ')'
                            << " , p = (" << p[0] << ',' << p[1] << ',' << p[2] << ')' << std::endl;
                    return 2;
                }
            }
        }
    }
    p_3Dvec -> push_back(p[0],p[1],p[2],NodeID);
    return 0;
}



unsigned VoroTesell::testOverlapAni    (    Mat3Dvec*                           p_3Dvec,
                                            const double                        lmin,
                                            const double                        anistropy,
                                            const unsigned                      dir,
                                            const Point                         rTest,
                                            const unsigned                      NodeID,
                                            const std::array<unsigned,Dim>      p)
{
    if (!testInside(rTest)) {
        return false;
    }

    int istart,iend,jstart,jend,kstart,kend;
    unsigned n;
    double l;
    std::vector<unsigned>   NodeIDList;
    Geometry            geo;

    std::array<unsigned,3>  dims = p_3Dvec->getDims ();

    istart = jstart = kstart = -1;
    iend = jend = kend = 1;
    std::array<double,3> r;
    r[0] = r[1] = r[2] = lmin;
    r[dir] = lmin/anistropy;

    if (p[0]==0) {
        istart = 0;
    }
    if (p[0]== (unsigned) Nx*(UnitLength/minD)-1) {
        iend = 0;
    }
    if (p[1]==0){
        jstart = 0;
    }
    if (p[1]== (unsigned) Ny*(UnitLength/minD)-1) {
        jend = 0;
    }
    if (p[2]==0) {
        kstart = 0;
    }
    if (p[2]== (unsigned) Nz*(UnitLength/minD)-1) {
        kend = 0;
    }

    for (int i=istart; i<=iend; i++) {
        for (int j=jstart; j<=jend; j++) {
            for (int k=kstart; k<=kend; k++) {
                int flag = (p_3Dvec -> get(p[0]+i,p[1]+j,p[2]+k,&NodeIDList));
                if (flag==1) {
                    for (n=0; n<NodeIDList.size(); n++) {
//                        l = geo.dist(rTest,GNodeTab[NodeIDList[n]].coord);
                        l=0;
                        for (unsigned d=0; d<Dim; d++) {
                            l += (rTest[d]-GNodeTab[NodeIDList[n]].coord[d])*(rTest[d]-GNodeTab[NodeIDList[n]].coord[d])/(r[d]*r[d]);
                        }
                        if (l<1.0) return 1;
                    }
                } else if (flag<0) {
                    tripOut << "[VoroTesell::testOverlap] : rTest ("
                            << rTest[0] << ',' << rTest[1] << ',' << rTest[2] << ')'
                            << " , p = (" << p[0] << ',' << p[1] << ',' << p[2] << ')' << std::endl;
                    return 2;
                }
            }
        }
    }
    p_3Dvec -> push_back(p[0],p[1],p[2],NodeID);
    return 0;
}


unsigned VoroTesell::testOverlap        (   Mat2Dvec*                           p_2Dvec,
                                            const Point                         rTest,
                                            const unsigned                      NodeID,
                                            const std::array<unsigned,2>        p)
{
    int istart,iend,jstart,jend;
    unsigned n;
    double l;
    std::vector<unsigned>   NodeIDList;
    Geometry            geo;

    if ((rTest[0]<0.0)||(rTest[0]>Nx*UnitLength)||(rTest[1]<0.0)||(rTest[1]>Ny*UnitLength)||(rTest[2]<0.0)||(rTest[2]>Ny*UnitLength)) {
        return false;
    }

    std::array<unsigned,2>  dims = p_2Dvec->getDims ();


    if (p[0]==0) {
            istart = 0;
        }
    //    if (p[0]== (unsigned) Nx*(UnitLength/minD)-1) {
    if (p[0]>= dims[0]-1) {
        istart = 0;
        iend = 0;
    }

    if (p[1]==0){
        jstart = 0;
    }

    if (p[1]>= dims[1]-1) {
        jstart = 0;
        jend = 0;
    }

    for (int i=istart; i<=iend; i++) {
        for (int j=jstart; j<=jend; j++) {
            if (p_2Dvec -> get(p[0]+i,p[1]+j,&NodeIDList)==1) {
                for (n=0; n<NodeIDList.size(); n++) {
                    l = geo.dist(rTest,GNodeTab[NodeIDList[n]].coord);
                    if (l<minD) return 1;
                }
            } else {
                return 2;
            }
        }
    }
    p_2Dvec -> push_back(p[0],p[1],NodeID);
    return 0;
}

void VoroTesell::genOuterNode	(	container_poly* 	p_conp,
									unsigned* 			p_NodeID,
                                	double				scaleLength,
                                	const double		scaleBoundary,
                                	const unsigned		layer)
{
	dualOut << " VoroTesell::genOuterNode started..." << std::endl;


    unsigned oldNodeSize = GNodeTab.size();
    std::array<double,Dim>  offset;
    offset[0] = 0.5*Nx*UnitLength;
    offset[1] = 0.5*Ny*UnitLength;
    offset[2] = 0.5*Nz*UnitLength;

    double xMean = 0.5*Nx*UnitLength;
    double yMean = 0.5*Ny*UnitLength;
    double zMean = 0.5*Nz*UnitLength;

    double xInLimit = xMean/scaleLength;
    double yInLimit = yMean/scaleLength;
    double zInLimit = zMean/scaleLength;

    double xOutLimit = xMean*scaleX/scaleLength;
    double yOutLimit = yMean*scaleY/scaleLength;
    double zOutLimit = zMean*scaleZ/scaleLength;

    double scale = scaleLength;
    for (unsigned i=0; i<layer; i++) {
        if (i==layer-1) {
            xOutLimit -= scaleX/scaleLength;
            yOutLimit -= scaleY/scaleLength;
            zOutLimit -= scaleZ/scaleLength;
        }
        for (unsigned NodeID=0; NodeID<oldNodeSize; NodeID++) {
            double xd = fabs(GNodeTab[NodeID].coord[0]-xMean);
            double yd = fabs(GNodeTab[NodeID].coord[1]-yMean);
            double zd = fabs(GNodeTab[NodeID].coord[2]-zMean);
            if   (((xd>xInLimit)||(yd>yInLimit)||(zd>zInLimit))
                 &&((xd<xOutLimit)&&(yd<yOutLimit)&&(zd<zOutLimit))) {
                nodeLists.free.push_back(*p_NodeID);
                putNode(p_conp,p_NodeID,
                        (GNodeTab[NodeID].coord[0]-xMean)*scale+offset[0],
                        (GNodeTab[NodeID].coord[1]-yMean)*scale+offset[1],
                        (GNodeTab[NodeID].coord[2]-zMean)*scale+offset[2],
                        Type_Outer);
            }
        }
        scale *= scaleLength;
    }
    dualOut << "Total Outer Node Num = " << GNodeTab.size()-oldNodeSize << std::endl;
}

void VoroTesell::genOuterNode   (   container_poly*     p_conp,
                                    unsigned*           p_NodeID,
                                    const double        Lmin,
                                    double              scaleLength,
                                    const double        scaleBoundary,
                                    unsigned            layer)
{
    dualOut << " VoroTesell::genOuterNode started..." << std::endl;

    if (Use_RandomNode) {
        srand(time(NULL));  //To randomize the seed
    }
    if (layer>1) {
        tripOut << "Layer can only be 1 in current code, set layer =1" << std::endl;
        layer = 1;
    }
    double lmin = Lmin*scaleLength;
    Mat3Dvec NodeInPartition;
    NodeInPartition.resize(Nx*UnitLength/lmin+1,Ny*UnitLength/lmin+1,Nz*UnitLength/lmin+1);
    unsigned rejectCnt = 0;
    double nMaxRej = 10.0;
    if (MaxTotalStep == 0) {
        nMaxRej = Disp_Scale; //= 0.00001;
    }
    double originalAdj = (1.0+lmin/(Nx*UnitLength))*(1.0+lmin/(Ny*UnitLength))*(1.0+lmin/(Nz*UnitLength));
    double rRatio = (scaleLength*scaleLength*scaleLength-originalAdj)/(scaleLength*scaleLength*scaleLength);
    unsigned oldNodeSize = GNodeTab.size();
    std::array<double,Dim>  offset;
    offset[0] = 0.5*Nx*UnitLength;
    offset[1] = 0.5*Ny*UnitLength;
    offset[2] = 0.5*Nz*UnitLength;

    double xMean = 0.5*Nx*UnitLength;
    double yMean = 0.5*Ny*UnitLength;
    double zMean = 0.5*Nz*UnitLength;
    double xInLimit = xMean/scaleLength;
    double yInLimit = yMean/scaleLength;
    double zInLimit = zMean/scaleLength;
    double xOutLimit = xMean*scaleX/scaleLength-lmin/2;
    double yOutLimit = yMean*scaleY/scaleLength-lmin/2;
    double zOutLimit = zMean*scaleZ/scaleLength-lmin/2;


    double scale = scaleLength;
    double n = PackingDensity;         //Percent from most dense HCF packing 0.7405
    if (Use_const_NodeNum) {
        n = 0.5;
    }
    double adj = std::pow(UnitLength*scaleLength/lmin,3)*
                            ((1-2*coarseness/Nx)*(1-2*coarseness/Ny)*(1-2*coarseness/Nz));
    unsigned targetNum = sqrt(2.0)*Nx*Ny*Nz*n*adj*rRatio;
    unsigned maxRej = nMaxRej*targetNum;
    unsigned nodeCnt = 0;
    std::cout << "lmin = " << lmin << " targetNum = " << targetNum << std::endl;

    for (unsigned i=0; i<layer; i++) {
        while (nodeCnt < targetNum) {
            Point rTest;
            rTest[0] = ((double) std::rand())/RAND_MAX*Nx*UnitLength;
            rTest[1] = ((double) std::rand())/RAND_MAX*Ny*UnitLength;
            rTest[2] = ((double) std::rand())/RAND_MAX*Nz*UnitLength;
            double xd = std::fabs(rTest[0]-xMean);
            double yd = std::fabs(rTest[1]-yMean);
            double zd = std::fabs(rTest[2]-zMean);
            if   (((xd>xInLimit)||(yd>yInLimit)||(zd>zInLimit))
                    &&((xd<xOutLimit)&&(yd<yOutLimit)&&(zd<zOutLimit))) {
                std::array<unsigned,Dim>    p;
                p[0] = ( rTest[0]/lmin-Tiny);
                p[1] = ( rTest[1]/lmin-Tiny);
                p[2] = ( rTest[2]/lmin-Tiny);
                rTest[0] = (rTest[0]-xMean)*scale+offset[0];
                rTest[1] = (rTest[1]-xMean)*scale+offset[1];
                rTest[2] = (rTest[2]-xMean)*scale+offset[2];
                if (testOverlap(&NodeInPartition,lmin,rTest,*p_NodeID,p)==0) {
                    nodeLists.free.push_back(*p_NodeID);
                    putNode(p_conp,p_NodeID,rTest,Type_Outer);
                    nodeCnt++;
                    rejectCnt = 0;
                } else {
                    rejectCnt++;
                    if (rejectCnt>maxRej) {
                        dualOut << "Cannot place a node, force quit the loop"<< std::endl;
                        dualOut << "Modified n (packing ratio relative to HCF packing) = "
                                << ((double) nodeCnt)/(targetNum) << std::endl;
                        break;
                    }
                }
            }
        }
        scale *= scaleLength;
    }
    dualOut << "Total Outer Node Num = " << GNodeTab.size()-oldNodeSize << std::endl;
}

void VoroTesell::genOuterNodeAni	(	container_poly* 	p_conp,
									    unsigned* 			p_NodeID,
									    double				scaleLength,
									    const double		scaleBoundary,
									    const unsigned		layer)
{
	dualOut << " VoroTesell::genOuterNodeAni started..." << std::endl;


    unsigned oldNodeSize = GNodeTab.size();
    std::array<double,Dim>  offset;
    offset[0] = 0.5*Nx*UnitLength;
    offset[1] = 0.5*Ny*UnitLength;
    offset[2] = 0.5*Nz*UnitLength;

    double xMean = 0.5*Nx*UnitLength;
    double yMean = 0.5*Ny*UnitLength;
    double zMean = 0.5*Nz*UnitLength;
    double xInLimit = xMean/scaleLength;
    double yInLimit = yMean/scaleLength;
    double zInLimit = zMean/scaleLength;
    double xOutLimit = xMean*scaleBoundary/scaleLength;
    double yOutLimit = yMean*scaleBoundary/scaleLength;
    double zOutLimit = zMean*scaleBoundary/scaleLength;

    double scale = scaleLength;
    for (unsigned i=0; i<layer; i++) {
        for (unsigned i=0; i<tmpNodeList.size(); i++) {
            double xd = fabs(tmpNodeList[i][0]-xMean);
            double yd = fabs(tmpNodeList[i][1]-yMean);
            double zd = fabs(tmpNodeList[i][2]-zMean);

            if   (((xd>xInLimit)||(yd>yInLimit)||(zd>zInLimit))
                 &&((xd<xOutLimit)&&(yd<yOutLimit)&&(zd<zOutLimit))) {
                nodeLists.free.push_back(*p_NodeID);
                putNode(p_conp,p_NodeID,
                        (tmpNodeList[i][0]-xMean)*scale+offset[0],
                        (tmpNodeList[i][1]-yMean)*scale+offset[1],
                        (tmpNodeList[i][2]-zMean)*scale+offset[2],
                        Type_Outer);
            }
        }
        scale *= scaleLength;
    }
    dualOut << "Total Outer Node Num = " << GNodeTab.size()-oldNodeSize << std::endl;
}

void VoroTesell::inputBoundary	(	container_poly* 			p_conp,
                                    unsigned* 					p_NodeID,
                                    std::vector<double>*		p_Di)
{
    unsigned 	i,j;
    double		rx,ry,rz;
    double 		bDist = 2.5*pow(((double) Nx*Ny*Nz)/tNodeNum,0.3333);
//	double		nN = sqrt(((double) Nnum)/(Nx*Ny*Nz));
    std::cout << " Starting VoroTesell::inputBoundary..." << std::endl;
    dualOut << bDist << std::endl;

    Npx = ((double) Nx)/bDist;
    Npy = ((double) Ny)/bDist;
    Npz = ((double) Nz)/bDist;
    double 		dlx = ((double) Nx)*UnitLength/Npx;
    double 		dly = ((double) Ny)*UnitLength/Npy;
    double 		dlz = ((double) Nz)*UnitLength/Npz;
    p_Di -> resize(Dim,0.0);

    (*p_Di)[0] = dlx;
    (*p_Di)[1] = dly;
    (*p_Di)[2] = dlz;

//	*p_NodeID = tNodeNum;

    nodeLists.boundary[0].clear();
    bStartEnd[0][0] = (*p_NodeID);

    //Face 1
    rz = -0.5*dlz;
    for (i = 0; i < Npx; i++) {
    	for (j = 0; j < Npy; j++) {
            rx = (0.5 + i)*dlx;
            ry = (0.5 + j)*dly;
            nodeLists.boundary[0].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID, rx,ry,rz,Type_Boundary);
        }
    }

    bStartEnd[0][1] = (*p_NodeID)-1;
    bStartEnd[1][0] = *p_NodeID;

    //Face 2
    nodeLists.boundary[1].clear();

    ry = -0.5*dly;
    for (i = 0; i < Npx; i++) {
    	for (j = 0; j < Npz; j++) {
        	rx = (0.5 + i)*dlx;
        	rz = (0.5 + j)*dlz;
        	nodeLists.boundary[1].push_back(*p_NodeID);
        	putNode(p_conp,p_NodeID, rx,ry,rz,Type_Boundary);
    	}
    }

    bStartEnd[1][1] = (*p_NodeID)-1;
    bStartEnd[2][0] = *p_NodeID;

        //Face 3
    nodeLists.boundary[2].clear();
    rx = (Npx+0.5)*dlx;
    for (i = 0; i < Npy; i++) {
    	for (j = 0; j < Npz; j++) {
        	ry = (0.5 + i)*dly;
            rz = (0.5 + j)*dlz;
            nodeLists.boundary[2].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID,rx,ry,rz,Type_Boundary);
        }
    }

    bStartEnd[2][1] = (*p_NodeID)-1;
    bStartEnd[3][0] = *p_NodeID;

    //Face 4
    nodeLists.boundary[3].clear();
    ry = (Npy+0.5)*dly;

    for (i = 0; i < Npx; i++) {
    	for (j = 0; j < Npz; j++) {
        	rx = (0.5 + i)*dlx;
            rz = (0.5 + j)*dlz;
            nodeLists.boundary[3].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID,rx,ry,rz,Type_Boundary);
        }
    }

    bStartEnd[3][1] = (*p_NodeID)-1;
    bStartEnd[4][0] = *p_NodeID;

    //Face 5
    nodeLists.boundary[4].clear();
    rx = -0.5*dlx;
    for (i = 0; i < Npy; i++) {
    	for (j = 0; j < Npz; j++) {
        	ry = (0.5 + i)*dly;
            rz = (0.5 + j)*dlz;
            nodeLists.boundary[4].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID,rx,ry,rz,Type_Boundary);
        }
    }

    bStartEnd[4][1] = (*p_NodeID)-1;
    bStartEnd[5][0] = *p_NodeID;

    //Face 6
    nodeLists.boundary[5].clear();
    rz = (Npz+0.5)*dlz;

    for (i = 0; i < Npx; i++) {
        for (j = 0; j < Npy; j++) {
        	rx = (0.5 + i)*dlx;
            ry = (0.5 + j)*dly;
            nodeLists.boundary[5].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID,rx,ry,rz,Type_Boundary);
        }
    }
    bStartEnd[5][1] = (*p_NodeID)-1;
    std::cout << " VoroTesell::inputBoundary completed" << std::endl;
}

void VoroTesell::inputOuterBoundary	(	container_poly* 		p_conp,
                                    	unsigned* 				p_NodeID,
                                    	double					scaleLength,
                                    	const double			scaleBoundary,
										std::vector<double>*	p_Di)
{
    unsigned 		i,j;
    double		rx,ry,rz;
    double 		bDist = 2.5*pow(((double) Nx*Ny*Nz)/tNodeNum,0.3333);

    std::array<double,Dim>		offset;
    offset[0] = (1.0-scaleBoundary)/2.0*Nx*UnitLength;
    offset[1] = (1.0-scaleBoundary)/2.0*Ny*UnitLength;
    offset[2] = (1.0-scaleBoundary)/2.0*Nz*UnitLength;

/*   offset[0] = -(scaleUp-1.0)/2.0*Nx*UnitLength;
    offset[1] = -(scaleUp-1.0)/2.0*Ny*UnitLength;
    offset[2] = -(scaleUp-1.0)/2.0*Nz*UnitLength;*/

//	double		nN = sqrt(((double) Nnum)/(Nx*Ny*Nz));

    Npx = (unsigned) Nx*scaleBoundary/scaleLength/bDist;
    Npy = (unsigned) Ny*scaleBoundary/scaleLength/bDist;
    Npz = (unsigned) Nz*scaleBoundary/scaleLength/bDist;

    double 		dlx = ((double) Nx)*UnitLength/Npx*scaleBoundary;
    double 		dly = ((double) Ny)*UnitLength/Npy*scaleBoundary;
    double 		dlz = ((double) Nz)*UnitLength/Npz*scaleBoundary;

    std::cout << "Npx = " << Npx << " dlx = " << dlx << std::endl;

    p_Di -> resize(Dim,0.0);

    (*p_Di)[0] = dlx;
    (*p_Di)[1] = dly;
    (*p_Di)[2] = dlz;

//	*p_NodeID = tNodeNum;

    nodeLists.boundary[0].clear();
    bStartEnd[0][0] = (*p_NodeID);

    //Face 1
    rz = -0.5*dlz;
    for (i = 0; i < Npx; i++) {
    	for (j = 0; j < Npy; j++) {
            rx = (0.5 + i)*dlx;
            ry = (0.5 + j)*dly;
            nodeLists.boundary[0].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID, rx+offset[0],ry+offset[1],rz+offset[2],Type_Boundary);
        }
    }

    bStartEnd[0][1] = (*p_NodeID)-1;
    bStartEnd[1][0] = *p_NodeID;

    //Face 2
    nodeLists.boundary[Face2].clear();

    ry = -0.5*dly;
    for (i = 0; i < Npx; i++) {
    	for (j = 0; j < Npz; j++) {
        	rx = (0.5 + i)*dlx;
        	rz = (0.5 + j)*dlz;
        	nodeLists.boundary[Face2].push_back(*p_NodeID);
        	putNode(p_conp,p_NodeID, rx+offset[0],ry+offset[1],rz+offset[2],Type_Boundary);
    	}
    }

    bStartEnd[1][1] = (*p_NodeID)-1;
    bStartEnd[2][0] = *p_NodeID;

        //Face 3
    nodeLists.boundary[Face3].clear();
    rx = (Npx+0.5)*dlx;
    for (i = 0; i < Npy; i++) {
    	for (j = 0; j < Npz; j++) {
        	ry = (0.5 + i)*dly;
            rz = (0.5 + j)*dlz;
            nodeLists.boundary[Face3].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID, rx+offset[0],ry+offset[1],rz+offset[2],Type_Boundary);
        }
    }

    bStartEnd[2][1] = (*p_NodeID)-1;
    bStartEnd[3][0] = *p_NodeID;

    //Face 4
    nodeLists.boundary[Face4].clear();
    ry = (Npy+0.5)*dly;

    for (i = 0; i < Npx; i++) {
    	for (j = 0; j < Npz; j++) {
        	rx = (0.5 + i)*dlx;
            rz = (0.5 + j)*dlz;
            nodeLists.boundary[Face4].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID, rx+offset[0],ry+offset[1],rz+offset[2],Type_Boundary);
        }
    }

    bStartEnd[3][1] = (*p_NodeID)-1;
    bStartEnd[4][0] = *p_NodeID;

    //Face 5
    nodeLists.boundary[Face5].clear();
    rx = -0.5*dlx;
    for (i = 0; i < Npy; i++) {
    	for (j = 0; j < Npz; j++) {
        	ry = (0.5 + i)*dly;
            rz = (0.5 + j)*dlz;
            nodeLists.boundary[Face5].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID, rx+offset[0],ry+offset[1],rz+offset[2],Type_Boundary);
        }
    }

    bStartEnd[4][1] = (*p_NodeID)-1;
    bStartEnd[5][0] = *p_NodeID;

    //Face 6
    nodeLists.boundary[Face6].clear();
    rz = (Npz+0.5)*dlz;

    for (i = 0; i < Npx; i++) {
        for (j = 0; j < Npy; j++) {
        	rx = (0.5 + i)*dlx;
            ry = (0.5 + j)*dly;
            nodeLists.boundary[Face6].push_back(*p_NodeID);
            putNode(p_conp,p_NodeID, rx+offset[0],ry+offset[1],rz+offset[2],Type_Boundary);
        }
    }
    bStartEnd[5][1] = (*p_NodeID)-1;
    std::cout << " VoroTesell::inputBoundary completed" << std::endl;
}

void VoroTesell::fillInBoundary	()
{
//	unsigned face[6] = {Face1,Face2,Face3,Face4,Face5,Face6};
	unsigned bNodeID,NodeID,nbNum;
	std::vector<unsigned>	tmpList;
	for (unsigned d=0; d<6; d++) {
		unsigned n = nodeLists.boundary[d].size();
		tmpList.clear();
		for (unsigned i=0; i<n; i++) {
			bNodeID = nodeLists.boundary[d][i];
			nbNum = GNodeTab[bNodeID].n;
			for (unsigned k=0; k<nbNum; k++) {
				NodeID = GNodeTab[bNodeID].nbNodeID[k];
				if (NodeID < tNodeNum) {
					tmpList.push_back(NodeID);
				}
			}
		}
		uniqueVector(&tmpList);
		nodeLists.inBoundary[d].clear();
		nodeLists.inBoundary[d] = tmpList;
		std::cout << nodeLists.inBoundary[d].size() << " ";
	}
}

std::vector<unsigned> VoroTesell::getInBoundary	(	const unsigned		layerNum)
{
	std::vector<unsigned>	nodeList;
	std::array<bool,6> face = getFaceFixity ();
	if (layerNum==0) {
	    std::list<unsigned> nList;
	    for (unsigned d=0; d<6; d++) {
	        if (face[d]) {
	            nList.insert(nList.end(),nodeLists.boundary[d].begin(),nodeLists.boundary[d].end());
	        }
	    }
	    uniqueList(&nList);
	    nodeList.assign(nList.begin(),nList.end());
	    return nodeList;
	}
	std::vector<std::list<unsigned> >   layer;
	layer.resize(layerNum);
	for (unsigned d=0; d<6; d++) {
	    if (face[d]) {
	        layer[0].insert(layer[0].end(),nodeLists.inBoundary[d].begin(),nodeLists.inBoundary[d].end());
	    }
	}
	uniqueList(&layer[0]);
	nodeList.assign(layer[0].begin(),layer[0].end());
	if (layerNum==1) {
	    return nodeList;
	} else {
	    unsigned bNodeID,NodeID;
	    for (std::list<unsigned>::iterator it = layer[0].begin(); it!=layer[0].end(); it++) {
	        bNodeID = *it;
	        unsigned nbNum = GNodeTab[bNodeID].n;
	        for (unsigned k=0; k<nbNum; k++) {
	            NodeID = GNodeTab[bNodeID].nbNodeID[k];
	            if (NodeID < tNodeNum) {
	                layer[1].push_back(NodeID);
	            }
	        }
	    }
	    uniqueList(&layer[1]);
	    nodeList.insert(nodeList.end(),layer[1].begin(),layer[1].end());
	    if (layerNum==2) {
	        uniqueVector(&nodeList);
	        return nodeList;
	    } else {
	        for (unsigned l=2; l<layerNum; l++) {
	            for (std::list<unsigned>::iterator it = layer[l-1].begin(); it!=layer[l-1].end(); it++) {
	                bNodeID = *it;
	                unsigned nbNum = GNodeTab[bNodeID].n;
	                for (unsigned k=0; k<nbNum; k++) {
	                    layer[l].push_back(NodeID);
	                }
	            }
	            uniqueList(&layer[l]);
	            for (std::list<unsigned>::iterator it = layer[l-1].begin(); it!=layer[l-1].end(); it++) {
	                layer[l].remove(*it);
	            }
	            nodeList.insert(nodeList.end(),layer[l].begin(),layer[l].end());
	        }
	        uniqueVector(&nodeList);
	        return nodeList;
	    }
	}
}

void VoroTesell::applyUnbreakableLattice	(	std::vector<unsigned>*	p_nodeList)
{
	unsigned NodeID;
	for (unsigned i=0; i<p_nodeList->size(); i++) {
		NodeID = (*p_nodeList)[i];
		for (unsigned k=0; k<GNodeTab[NodeID].n; k++) {
			LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[0] = 1000.0*Huge;
			LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[1] = -1000.0*Huge;
		}
	}
}

void VoroTesell::applyUnbreakableNode ()
{
    if (Use_UnevenMesh) {
        for (unsigned NodeID = refineNode[1]+1; NodeID<GNodeTab.size(); NodeID++) {
            unBreakableNodeList.push_back(NodeID);
        }
    }
    for (unsigned i=0; i<unBreakableNodeList.size(); i++) {
        for (unsigned k=0; k<GNodeTab[unBreakableNodeList[i]].nbLatticeID.size(); k++) {
            LatTab[GNodeTab[unBreakableNodeList[i]].nbLatticeID[k]].et[0] = 1000.0*Huge;
            LatTab[GNodeTab[unBreakableNodeList[i]].nbLatticeID[k]].et[1] = -1000.0*Huge;
        }
    }
}

void VoroTesell::applyUnbreakableInBoundary	()
{
	unsigned face[6] = {Face1,Face2,Face3,Face4,Face5,Face6};
	unsigned NodeID;
	for (unsigned d=0; d<6; d++) {
		for (unsigned i=0; i<nodeLists.inBoundary[face[d]].size(); i++) {
			NodeID = nodeLists.inBoundary[face[d]][i];
			for (unsigned k=0; k<GNodeTab[NodeID].n; k++) {
				LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[0] = 1000.0*Huge;
				LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[1] = -1000.0*Huge;
			}
		}
	}
}

void VoroTesell::applyUnbreakableCoarseMesh	()
{
	for (unsigned NodeID=tInNodeNum; NodeID<tNodeNum; NodeID++) {
		for (unsigned k=0; k<GNodeTab[NodeID].n; k++) {
			LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[0] = 1000.0*Huge;
			LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[1] = -1000.0*Huge;
		}
	}
}

bool VoroTesell::testSameBoundary		(	const unsigned					Node1,
        									const unsigned					Node2)
{
	unsigned face;
    for (face=0; face<6; face++) {
        if ((Node1>=bStartEnd[face][0])&&(Node1<=bStartEnd[face][1])) {
            break;
        }
    }
    if ((Node2>=bStartEnd[face][0])&&(Node2<=bStartEnd[face][1])) {
        return true;
    } else {
        return false;
    }
}

bool VoroTesell::testMinCoordNum		(	const unsigned					MinCoordNum)
{
    unsigned min = 100;

    for (unsigned NodeID = 0; NodeID < tNodeNum+tbNodeNum; NodeID++) {
        if (GNodeTab[NodeID].ne < min) {
            min = GNodeTab[NodeID].ne;
        }
        if (GNodeTab[NodeID].ne <MinCoordNum) {
            std::cout << GNodeTab[NodeID].coord[0] << ' ' << GNodeTab[NodeID].coord[1] << ' ' << GNodeTab[NodeID].coord[2] << std::endl;
        }
    }
    dualOut << "Min coordination number = " << min << std::endl;
    if (min < MinCoordNum) {
        return false;
    } else {
        return true;
    }
}

void VoroTesell::findNeighbour	(	const double					scaleLength,
						            std::vector<NodeIDPair>*		p_conjNodeList)
{
    int 		NodeID,nbNodeID;
    unsigned 	LatticeID = 0;
    double		minArea = MinArea*UnitLength*UnitLength;
    std::array<double,2>		weakMicroAxialStrength;
    weakMicroAxialStrength[0] 	= MicroAxialStrength[0]/10;
    weakMicroAxialStrength[1] 	= MicroAxialStrength[1]/10;

    std::vector<unsigned>		b_NodeList;
    std::vector<unsigned>		b_nbNodeList;
    std::vector<double>			b_AreaList;
    std::vector<unsigned>		b_nList;

    std::vector<unsigned>		v_NodeList;
    std::vector<unsigned>		v_nbNodeList;
    std::vector<double>			v_AreaList;
    std::vector<unsigned>		v_nList;

    std::vector<unsigned>		bb_NodeList;
    bool						is_bb = false;

    std::cout << "VoroTesell::findNeighbour..." << std::endl;
//    std::vector<NodeIDPair>::iterator it_cNodeList;
//    std::vector<int>::iterator it_nbNodeID;
//    unsigned conNodeNum = 2*p_conjNodeList->size();
//    unsigned weakNodeNum = 2*p_conjNodeList->size() + 2*weakPlane.size();
//    unsigned halfConNodeNum = p_conjNodeList->size();
//    unsigned halfweakNodeNum = 2*p_conjNodeList->size()+weakPlane.size();
//    unsigned kk;
    unsigned rejCnt = 0;
    for (unsigned i = 0; i < p_voroInfo-> IDList.size(); i++) {
        NodeID = p_voroInfo-> IDList[i];
/*        if ((unsigned) NodeID < p_conjNodeList->size()) {
        	it_cNodeList = std::find(p_conjNodeList->begin(),p_conjNodeList->end(),((unsigned) NodeID));
        	if (it_cNodeList!=p_conjNodeList->end()) {
        		it_nbNodeID = std::find(p_voroInfo->nbList[i].begin(),p_voroInfo->nbList[i].end(),((int) it_cNodeList->i2));
        		kk = it_nbNodeID - p_voroInfo->nbList[i].begin();
        		if (it_nbNodeID!=p_voroInfo->nbList[i].end()) {
        		    preExistingCrack.push_back(LatticeID);
        			putLattice(NodeID,*it_nbNodeID,NormLatStiffness,MicroAxialStrength,StrengthFactor,p_voroInfo -> area[i][kk],&LatticeID,kk,Use_LEFM);
        		} else {
        			tripOut << "[VoroTesell::findNeighbour] Finding pre-existing crack Lattice error! cannot find nbNodeID" << std::endl;
        		}
        	} else {
        		tripOut << "[VoroTesell::findNeighbour] Finding pre-existing crack Lattice error! cannot find NodeID" << std::endl;
        	}
        } else if (((unsigned) NodeID < 2*p_conjNodeList->size() + weakPlane.size())&&((unsigned) NodeID >= 2*p_conjNodeList->size())) {
        	it_cNodeList = std::find(weakPlane.begin(),weakPlane.end(),((unsigned) NodeID));
        	if (it_cNodeList!=weakPlane.end()) {
        		it_nbNodeID = std::find(p_voroInfo->nbList[i].begin(),p_voroInfo->nbList[i].end(),((int) it_cNodeList->i2));
        	    kk = it_nbNodeID - p_voroInfo->nbList[i].begin();
        	    if (it_nbNodeID!=p_voroInfo->nbList[i].end()) {
//        	    	pxCrackLatID.push_back(LatticeID);
        	        weakPlaneLatID.push_back(LatticeID);
        	        putLattice(NodeID,*it_nbNodeID,NormLatStiffness,MicroAxialStrength,1.0,p_voroInfo -> area[i][kk],&LatticeID,kk,Use_LEFM);
        	    } else {
        	    	tripOut << "[VoroTesell::findNeighbour] Finding weak plane Lattice error! cannot find nbNodeID" << std::endl;
        	    	std::cout << "Weak plane: " << NodeID << std::endl;
        	    }
        	} else {
        		tripOut << "[VoroTesell::findNeighbour] Finding weak plane Lattice error! cannot find NodeID" << std::endl;
        		std::cout << "NodeID = " << NodeID << std::endl;
        	}
        }*/
        if (NodeID < tNodeNum) {
            for (unsigned k=0; k < p_voroInfo->nbList[i].size(); k++) {
                nbNodeID = p_voroInfo->nbList[i][k];
                if (((nbNodeID > NodeID)&&(nbNodeID >= 0))
//                    &&((NodeID >= halfConNodeNum)||((nbNodeID>conNodeNum)||(nbNodeID<halfConNodeNum)))&&
//					(((NodeID >= halfweakNodeNum)||(NodeID < conNodeNum))||((nbNodeID>weakNodeNum)||(nbNodeID<halfweakNodeNum)))
                ) {
                    if (nbNodeID < tNodeNum) {
                    	if ((NodeID<tInNodeNum)||(nbNodeID<tInNodeNum)) {
                    		if (p_voroInfo->area[i][k] > minArea) {
                    			putLattice(NodeID,nbNodeID,NormLatStiffness,MicroAxialStrength,StrengthFactor,p_voroInfo -> area[i][k],&LatticeID,k,Use_LEFM);
                    		} else {
								rejCnt++;
							}
                    	} else {
                    		if (p_voroInfo->area[i][k] > minArea*scaleLength*scaleLength) {
                    			putLattice(NodeID,nbNodeID,NormLatStiffness,MicroAxialStrength,StrengthFactor,p_voroInfo -> area[i][k],&LatticeID,k,Use_LEFM);
                    		} else {
                    			rejCnt++;
                    		}
                    	}
                    } else if (nbNodeID < tNodeNum+tbNodeNum) {
                        b_NodeList.push_back(NodeID);
                        b_nbNodeList.push_back(nbNodeID);
                        b_nList.push_back(k);
                        b_AreaList.push_back(p_voroInfo -> area[i][k]);
                        is_bb = true;
                    }
                } else if (nbNodeID < 0) {
                    updateBoundaryList(NodeID,nbNodeID,p_voroInfo -> area[i][k]);
/*                    std::cout << "nbNodeID = " << nbNodeID << ' ' << GNodeTab[NodeID].coord[0] << ' ' << GNodeTab[NodeID].coord[1]
                              << ' ' << GNodeTab[NodeID].coord[2]<< ' ' << p_voroInfo -> area[i][k] << '\n';*/
                }
            }
            if (is_bb) {
                bb_NodeList.push_back(NodeID);
                is_bb = false;
            }
        } else {
            for (unsigned k=0; k < p_voroInfo->nbList[i].size(); k++) {
                nbNodeID = p_voroInfo->nbList[i][k];
                if (nbNodeID>NodeID) {
                    if (testSameBoundary(NodeID,nbNodeID)) {
                        v_NodeList.push_back(NodeID);
                        v_nbNodeList.push_back(nbNodeID);
                        v_nList.push_back(k);
                        v_AreaList.push_back(p_voroInfo -> area[i][k]);
                    }
                } else if (nbNodeID < 0) {
                    updateBoundaryList(NodeID,nbNodeID,p_voroInfo -> area[i][k]);
                    std::cout << "nbNodeID = " << nbNodeID << ' ' << GNodeTab[NodeID].coord[0] << ' ' << GNodeTab[NodeID].coord[1]
                              << ' ' << GNodeTab[NodeID].coord[2]<< ' ' << p_voroInfo -> area[i][k] << '\n';
                }
            }
        }
    }
    //Generate lattice for pre-existing crack
/*    for (unsigned i = 0; i<p_conjNodeList->size(); i++) {
    	preExistingCrack.push_back(LatticeID);
    	putLattice((*p_conjNodeList)[i].first,(*p_conjNodeList)[i].second,NormLatStiffness,LatticeBreakDisp,p_voroInfo -> area[i][k],&LatticeID,k);
    }*/
    dualOut << " No of lattice rejected as plane area is smaller than " << minArea/(UnitLength*UnitLength) << "*D*D = " << rejCnt
            << " Reject percentage = " << ((double) rejCnt)/LatticeID*100. << std::endl;

    tLatticeNum = LatticeID;
    std::array<double,2>		Unbreakable;
    Unbreakable[0] = 1000.*Huge;
    Unbreakable[1] = -1000.*Huge;

    for (unsigned i=0; i<b_NodeList.size(); i++) {
        putLattice(b_NodeList[i],b_nbNodeList[i],NormLatStiffness,Unbreakable,1.0,b_AreaList[i],&LatticeID,b_nList[i],Use_LEFM);
    }
    tbLatticeNum = LatticeID - tLatticeNum;

    for (unsigned i=0; i<v_NodeList.size(); i++) {
        putLattice(v_NodeList[i],v_nbNodeList[i],NormLatStiffness*Tiny,Unbreakable,1.0,v_AreaList[i],&LatticeID,v_nList[i],Use_LEFM);		// Tiny stiffness
//    	putLattice(v_NodeList[i],v_nbNodeList[i],NormLatStiffness,Unbreakable,v_AreaList[i],&LatticeID,v_nList[i]);
    }
    tvLatticeNum = LatticeID - tLatticeNum - tbLatticeNum;

    adjCellArea (&bb_NodeList);
    double tmpArea = getMeanCellArea(0,tLatticeNum);

    dualOut << "No of inside lattice = " << tLatticeNum << " no of boundary lattice = " << tbLatticeNum
    		<< " no of virtual lattice = " << tvLatticeNum << std::endl;
    dualOut << "Mean Area becomes = " << tmpArea << std::endl;
    dualOut << "Area SD becomes = " << getCellAreaSD(0,tLatticeNum,tmpArea) << std::endl;
    return;
}

void VoroTesell::findPreExistingCrack (     std::vector<NodeIDPair>*     p_conjNodeList)
{
    dualOut << "VoroTesell::findPreExistingCrack..." << std::endl;
    for (unsigned i=0; i< p_conjNodeList->size(); i++) {
        unsigned NodeID = (*p_conjNodeList)[i].i1;
        for (unsigned k=0; k<GNodeTab[NodeID].nbNodeID.size(); k++) {
            if (GNodeTab[NodeID].nbNodeID[k]==(*p_conjNodeList)[i].i2) {
                unsigned LatticeID = GNodeTab[NodeID].nbLatticeID[k];
                preExistingCrack.push_back(LatticeID);
                LatTab[LatticeID].et[0] *= StrengthFactor;
                LatTab[LatticeID].et[1] *= StrengthFactor;
//                if (GeoID==1) {
//                    LatTab[LatticeID].length /= 2.0;
//                }
            }
        }
    }
}

void VoroTesell::findWeakPlane ()
{
    dualOut << "VoroTesell::findWeakPlane..." << std::endl;
    for (unsigned i=0; i< weakPlane.size(); i++) {
        unsigned NodeID = weakPlane[i].i1;
        for (unsigned k=0; k<GNodeTab[NodeID].nbNodeID.size(); k++) {
            if (GNodeTab[NodeID].nbNodeID[k]==weakPlane[i].i2) {
                weakPlaneLatID.push_back(GNodeTab[NodeID].nbLatticeID[k]);
                if (GeoID!=4) {
                    LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[0] /= StrengthFactor;
                } else {
                    if ((NodeID<BHNodeNum) || (weakPlane[i].i2<BHNodeNum)) {
                        LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[0] = Tiny;
                        LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[1] = -Tiny;
                        LatTab[GNodeTab[NodeID].nbLatticeID[k]].initOpening = 0.25;
                    } else {
                        LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[0] /= StrengthFactor;
                    }
                }
            }
        }
    }
}

void VoroTesell::setUnbreakablePair ()
{
    dualOut << "VoroTesell::setUnbreakablePair..." << std::endl;
    for (unsigned NodeID=0; NodeID<BHNodeNum; NodeID++) {
        for (unsigned k=0; k<GNodeTab[NodeID].nbNodeID.size(); k++) {
            if (GNodeTab[NodeID].nbNodeID[k]<BHNodeNum) {
                LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[0] = Huge;
            }
        }
    }
    /*
    for (unsigned i=0; i< unBreakablePair.size(); i++) {
        unsigned NodeID = unBreakablePair[i].i1;
        for (unsigned k=0; k<GNodeTab[NodeID].nbNodeID.size(); k++) {
            if (GNodeTab[NodeID].nbNodeID[k]==unBreakablePair[i].i2) {
                LatTab[GNodeTab[NodeID].nbLatticeID[k]].et[0] = Huge;
            }
        }
    }
    */
}

void VoroTesell::setBreakableLatticeAtBoundary ()
{
    dualOut << "VoroTesell::setBreakableLatticeAtBoundary..." << std::endl;
    std::vector<unsigned>    face;
    face.push_back(Face1);
    face.push_back(Face6);
    double yMid = Ny*UnitLength*0.5;
    for (unsigned d=0; d<face.size(); d++) {
        for (unsigned i=0; i<nodeLists.boundary[face[d]].size(); i++) {
            unsigned NodeID = nodeLists.boundary[face[d]][i];
            for (unsigned k=0; k<GNodeTab[NodeID].nbNodeID.size(); k++) {
                unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[k];
                if ((GNodeTab[NodeID].coord[1]-yMid)*(GNodeTab[nbNodeID].coord[1]-yMid)<0.0) {
                    unsigned LatticeID = GNodeTab[NodeID].nbLatticeID[k];
                    LatTab[LatticeID].et[0] = LatTab[LatticeID].area*MicroAxialStrength[0]/LatTab[LatticeID].k[0];
                    break;
                }
            }
        }
    }
}

void VoroTesell::putNode	(	container_poly* 		p_conp,
                                unsigned*				p_NodeID,
                                const double 			rx,
                                const double 			ry,
                                const double 			rz,
                                const unsigned 			type)
{
    p_conp -> put(*p_NodeID,rx,ry,rz,UnitLength);
    GNode	newGNode(type,rx,ry,rz);
    GNodeTab.push_back(newGNode);
    (*p_NodeID)++;
}

void VoroTesell::putNode	(	container_poly* 		p_conp,
                                unsigned*				p_NodeID,
                                const Point				p,
                                const unsigned 			type)
{
    p_conp -> put(*p_NodeID,p[0],p[1],p[2],UnitLength);
    GNode	newGNode(type,p);
    GNodeTab.push_back(newGNode);
    (*p_NodeID)++;
}

void VoroTesell::putLattice	(	const unsigned 				NodeID,
                                const unsigned				nbNodeID,
                                const double				latStiffness,
                                const std::array<double,2>	latBreakStress,
								const double				factor,
                                const double				latArea,
                                unsigned*					p_LatticeID,
                                const unsigned				vIndex,
                                const bool					useLEFM)
{
    Geometry	geo;
    double len = geo.dist(GNodeTab[NodeID].coord,GNodeTab[nbNodeID].coord);
    double kn = latStiffness*latArea/len;
    std::array<double,2> et;
    if (useLEFM) {
    	et[0] = sqrt(2.0*ReConRatio*factor*latArea/kn);
    } else {
    	et[0] = latArea*latBreakStress[0]*factor/kn;
    }
    et[1] = latBreakStress[1];
    if (KNum==1) {
        Lattice newLat(NodeID,nbNodeID,len,latArea,kn,et);
        LatTab.push_back(newLat);
    } else if (KNum==2) {
        double ks = ShearToAxialStiffnessRatio*kn;
        Lattice newLat(NodeID,nbNodeID,len,latArea,kn,ks,et);
        LatTab.push_back(newLat);
    } else if (KNum==4){
        double ks = ShearToAxialStiffnessRatio*kn;
        double k_theta = MomentToAxialStiffnessRatio*latStiffness;
        Lattice newLat(NodeID,nbNodeID,len,latArea,kn,ks,k_theta,et,MicroShearStrength);
        LatTab.push_back(newLat);
    } else {
        tripOut << "[VoroTesell::putLattice] KNum = " << KNum << " is not defined in code" << std::endl;
    }
    GNodeTab[NodeID].setNeighbour(nbNodeID,*p_LatticeID);
    GNodeTab[nbNodeID].setNeighbour(NodeID,*p_LatticeID);
    p_voroInfo->vListIndex[NodeID].push_back(vIndex);
    p_voroInfo->vListIndex[nbNodeID].push_back(-1);
    (*p_LatticeID)++;
}

void VoroTesell::assignRandomStrength   (   const unsigned      type,
                                            const double        sigma,
                                            const double        min)
{
    std::default_random_engine generator;
    if (Use_RandomNode) {
        generator.seed(time(NULL));
    }
    std::normal_distribution<double> normDist(1.0,sigma);
    std::lognormal_distribution<double> logNormDist(0.0,sigma);
    std::weibull_distribution<double> weibullDist(sigma,1.0);
    std::vector<double>     nList;
    std::vector<unsigned>   latList;
    double ranRatio;
    double sum = 0.0;
    double alpha,a,b;

    switch (type) {
    case 0:
        dualOut << "No distribution is assigned to lattice strength. " << std::endl;
        return;
    case 1:
        dualOut << "Assign normal distribution to strength, mean = 1.0, sd = " << sigma
                << " min = " << min << std::endl;

        for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
            if (LatTab[LatticeID].et[0]<Huge) {
                do {
                    ranRatio = normDist(generator);
                } while (ranRatio < min);
                nList.push_back(ranRatio);
                latList.push_back(LatticeID);
            }

        }

        break;
    case 2:
        dualOut << "Assign lognormal distribution to strength, mu = 1.0, sigma = " << sigma
                << " min = " << min << std::endl;
        for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
            if (LatTab[LatticeID].et[0]<Huge) {
                do {
                    ranRatio = logNormDist(generator);
                } while (ranRatio < min);
                nList.push_back(ranRatio);
                latList.push_back(LatticeID);
            }
        }

        break;
    case 3:
        dualOut << "Assign Weibull distribution to strength, mu = 1.0, beta = " << sigma
                << " min = " << min << std::endl;
        for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
            if (LatTab[LatticeID].et[0]<Huge) {
                do {
                    ranRatio = weibullDist(generator);
                } while (ranRatio < min);
                nList.push_back(ranRatio);
                latList.push_back(LatticeID);
            }
        }
        break;
    }

    for (unsigned i=0; i<nList.size(); i++) {
        sum += nList[i];
    }
    alpha = sum / nList.size();

    dualOut << "Mean = " << alpha;
    a = (1.0 - min)/(alpha - min);
    b = min*(1.0 - alpha)/(min - alpha);
    dualOut << " (a,b) = " << a << ' ' << b;
    for (unsigned i=0; i<nList.size(); i++) {
        nList[i] *= a;
        nList[i] += b;
    }
    sum = 0.0;
    std::array<unsigned,2> tmp;
    tmp[0] = tmp[1] = 0;
    Statistics stat(OutputSubFolder,OutputFilePrefix,0,2*pxPerimeterNum.size(),tmp);
    for (unsigned i=0; i<nList.size(); i++) {
        if (nList[i]<min) {
            dualOut << nList[i] << ' ';
        }
        sum += nList[i];
        ranRatio = nList[i];
        unsigned LatticeID = latList[i];
        if (Use_RanLatStrength) {
            LatTab[LatticeID].et[0] *= ranRatio;
            LatTab[LatticeID].shearStr *= ranRatio;
        }
        if (Use_RanLatStiffness) {
            for (unsigned i=0; i<KNum; i++) {
                LatTab[LatticeID].k[i] *= ranRatio;
            }
        }
    }
    stat.printPDF(&nList,"ranRatio",0,200,0.0);
    dualOut << " Modified mean = " << sum / nList.size() << std::endl;
}

void VoroTesell::putBeam    (   const unsigned              NodeID,
                                const unsigned              nbNodeID,
                                const double                latStiffness,
                                const std::array<double,2>  latBreakStress,
                                const double                factor,
                                const double                alpha,
                                const double                latArea,
                                unsigned*                   p_LatticeID,
                                const unsigned              vIndex,
                                const bool                  useLEFM)
{
    Geometry    geo;
    double len = geo.dist(GNodeTab[NodeID].coord,GNodeTab[nbNodeID].coord);
    double kn = latStiffness*latArea/len;
    std::array<double,2> et;
    if (useLEFM) {
        et[0] = std::sqrt(2.0*ReConRatio*factor*latArea/kn);
    } else {
        et[0] = latArea*latBreakStress[0]*factor/kn;
    }
    et[1] = latBreakStress[1];
    double ks = alpha*kn;
    Beam newBeam(NodeID,nbNodeID,len,latArea,kn,ks,0.0,0.0,0.0,et);
    LatBeamTab.push_back(newBeam);
    GNodeTab[NodeID].setNeighbour(nbNodeID,*p_LatticeID);
    GNodeTab[nbNodeID].setNeighbour(NodeID,*p_LatticeID);
    p_voroInfo->vListIndex[NodeID].push_back(vIndex);
    p_voroInfo->vListIndex[nbNodeID].push_back(-1);
    (*p_LatticeID)++;
}

void VoroTesell::calLatParameters ()
{
    for (unsigned LatticeID = 0; LatticeID < LatTab.size(); LatticeID++) {
        LatTab[LatticeID].k[0] *= LatTab[LatticeID].area / (LatTab[LatticeID].length);
        LatTab[LatticeID].et[0] *= sqrt(2.0*ReConRatio*LatTab[LatticeID].area/LatTab[LatticeID].k[0]);
    }
}

void VoroTesell::updateBoundaryList     (   const unsigned          NodeID,
                                            const int               vFace,
                                            const double            area)
{
    switch (vFace)
    {
    case (-1):
        nodeLists.boundary[Face5].push_back(NodeID);
        nodeLists.boundaryArea[Face5].push_back(area);
        break;

    case (-2):
        nodeLists.boundary[Face3].push_back(NodeID);
        nodeLists.boundaryArea[Face3].push_back(area);
        break;

    case (-3):
        nodeLists.boundary[Face2].push_back(NodeID);
        nodeLists.boundaryArea[Face2].push_back(area);
        break;

    case (-4):
        nodeLists.boundary[Face4].push_back(NodeID);
        nodeLists.boundaryArea[Face4].push_back(area);
        break;

    case (-5):
        nodeLists.boundary[Face1].push_back(NodeID);
        nodeLists.boundaryArea[Face1].push_back(area);
        break;

    case (-6):
        nodeLists.boundary[Face6].push_back(NodeID);
        nodeLists.boundaryArea[Face6].push_back(area);
        break;
    }
}

void VoroTesell::updateVertices		(	Vertices*					p_verti)
{
    Clock vClock;

    vClock.start("Vertices");
    p_verti->initialize(tLatticeNum+tbLatticeNum+tvLatticeNum,Tiny);

    for (unsigned NodeID=0; NodeID < tNodeNum+tbNodeNum; NodeID++)
    {
        p_verti->inputNode(NodeID,&p_voroInfo->vCoord[NodeID],&p_voroInfo->vList[NodeID],&p_voroInfo->vListIndex[NodeID]);
    }
    vClock.get();
}

void    VoroTesell::updateLatInfo       (   Vertices*                       p_verti)
{
    Geometry geo;
    unsigned NodeID, nbNodeID;
    UniVec      e1,e2;
    for (unsigned LatticeID=0; LatticeID < LatTab.size(); LatticeID++) {
        Surface v = p_verti-> getVertiCoordArray (LatticeID);
        Point centroid = geo.getCentroid(LatTab[LatticeID].area,v);
        LatTab[LatticeID].centroid = centroid;
        NodeID = LatTab[LatticeID].nb[0];
        nbNodeID = LatTab[LatticeID].nb[1];
        // set the centroid to be origin
        for (unsigned i=0; i<v.size(); i++) {
            for (unsigned d=0; d<Dim; d++) {
                v[i][d] -= centroid[d];
            }
        }
        double e1Len = 0.0;
        for (unsigned d=0; d<Dim; d++) {
            LatTab[LatticeID].axes[0][d] = (GNodeTab[nbNodeID].coord[d] - GNodeTab[NodeID].coord[d])
                    / LatTab[LatticeID].length;
            e1[d] = v[0][d];
            e1Len += e1[d]*e1[d];
        }
        e1Len = std::sqrt(e1Len);
        for (unsigned d=0; d<Dim; d++) {
            e1[d] /=e1Len;
        }
        e2 = geo.crossProd(LatTab[LatticeID].axes[0],e1);

        if ((e1!=e1)||(e2!=e2)) {
//            if (LatTab[LatticeID].k[2]!=LatTab[LatticeID].k[2]) {
            unsigned NodeID1 = LatTab[LatticeID].nb[0];
            unsigned NodeID2 = LatTab[LatticeID].nb[1];
            tripOut << "LatticeID = " << LatticeID
                    << " NodeIDs coord - "
                    << GNodeTab[NodeID1].coord[0] << ' ' << GNodeTab[NodeID1].coord[1] << ' ' << GNodeTab[NodeID1].coord[2] << '|'
                    << GNodeTab[NodeID2].coord[0] << ' ' << GNodeTab[NodeID2].coord[1] << ' ' << GNodeTab[NodeID1].coord[2] << '|'
                    << "Centroid = " << centroid[0] << ' ' << centroid[1] << ' ' << centroid[2] << std::endl;
//                        }
        }

/*        double e1e1=geo.dot(e1,e1);
        double e2e2=geo.dot(e2,e2);
        double e1e2=geo.dot(e1,e2);
        if (std::fabs(e1e1)-1.0>Tiny) {
            std::cout << "Unit vector e1.e1 = 1.0 check error, e1e1 = " << e1e1 << std::endl;
        }
        if (std::fabs(e2e2)-1.0>Tiny) {
            std::cout << "Unit vector e2.e2 = 1.0 check error, e2e2 = " << e2e2 << std::endl;
        }
        if (std::fabs(e1e2)>Tiny) {
            std::cout << "Unit vector e1.e2 = 0.0 check error, e1e2 = " << e1e2 << std::endl;
        }*/
        std::vector<std::array<double,2> >  v2D;
        std::array<double,2>    coord2D;
        for (unsigned i=0; i<v.size(); i++) {
            coord2D[0] = geo.dot(v[i],e1);
            coord2D[1] = geo.dot(v[i],e2);
            v2D.push_back(coord2D);
        }
        coord2D[0] = geo.dot(v[0],e1);
        coord2D[1] = geo.dot(v[0],e2);
        v2D.push_back(coord2D);

        std::array<double,3> I=geo.getSecMomentOfArea2D(v2D);
//        std::cout << I[0] << ' ' << I[1] << ' ' << I[2] << '\n';
        double theta;
        if (I[0]==I[1]) {
        	theta = Pi/4.0;
/*        	if ((I[0]-I[1])>0.0) {
        		theta = Pi/4.0;
        	} else {
        		theta = -Pi/4.0;
        	}*/
        } else {
        	theta = 0.5*std::atan((-2.0*I[2])/(I[0]-I[1]));
        }
//        std::cout << theta << '\n';
//        if (I[1]>I[0]) {
            for (unsigned d=0; d<Dim; d++) {
                LatTab[LatticeID].axes[1][d] = std::cos(theta)*e1[d] + std::sin(theta)*e2[d];
//                LatTab[LatticeID].axes[2][d] = -std::sin(theta)*e1[d] + std::cos(theta)*e2[d];
            }
            LatTab[LatticeID].axes[2] =  geo.crossProd(LatTab[LatticeID].axes[0],LatTab[LatticeID].axes[1]);
/*        } else {
            for (unsigned d=0; d<Dim; d++) {
                LatTab[LatticeID].axes[2][d] = std::cos(theta)*e1[d] + std::sin(theta)*e2[d];
                LatTab[LatticeID].axes[1][d] = -std::sin(theta)*e1[d] + std::cos(theta)*e2[d];
            }
        }*/
        if (KNum==4) {
        	double Iu = I[0]*std::cos(theta)*std::cos(theta) + I[1]*std::sin(theta)*std::sin(theta) - I[2]*std::sin(2.0*theta);
        	double Iv = I[0]*std::sin(theta)*std::sin(theta) + I[1]*std::cos(theta)*std::cos(theta) + I[2]*std::sin(2.0*theta);
            LatTab[LatticeID].k[2] *= Iu/LatTab[LatticeID].length;
            LatTab[LatticeID].k[3] *= Iv/LatTab[LatticeID].length;

        }
        Point intersect,intersect_local;
        Point d;
        for (unsigned i=0; i<Dim; i++) {
            d[i] = centroid[i]-GNodeTab[NodeID].coord[i];
        }
        double k = geo.dot(d,LatTab[LatticeID].axes[0]);
        for (unsigned i=0; i<Dim; i++) {
            intersect[i] =GNodeTab[NodeID].coord[i] + k*LatTab[LatticeID].axes[0][i];
            intersect_local[i] = centroid[i] - intersect[i];
        }
        LatTab[LatticeID].intersect = intersect;
        LatTab[LatticeID].offset[0] = geo.dot(intersect_local,LatTab[LatticeID].axes[1]);
        LatTab[LatticeID].offset[1] = geo.dot(intersect_local,LatTab[LatticeID].axes[2]);
    }
}


//Transfer data from Voro library to p_voroInfo
void VoroTesell::updateVoroInfo		(	container_poly* 			p_con)
{
    unsigned i = 0;
    unsigned v,n,k,NodeID,vID;
    unsigned fOrder;
    std::array<double,3> coord;
   std::vector<double> vertices;
   std::vector<int> vList;

    voronoicell_neighbor c;
    c_loop_all vl(*p_con);
    std::cout << "Starting VoroTesell::updateVoroInfo..." << std::endl;

    p_voroInfo	= new VoroInfo(GNodeTab.size());

    if(vl.start()) do {
            p_con -> compute_cell(c,vl);
            NodeID = vl.id[vl.ijk][vl.q];
            p_voroInfo -> IDList.push_back(NodeID);
            p_voroInfo -> vol.push_back(c.volume());
            c.centroid(coord[0],coord[1],coord[2]);
            p_voroInfo -> centroid.push_back(coord);
            c.vertices(GNodeTab[NodeID].coord[0],GNodeTab[NodeID].coord[1],GNodeTab[NodeID].coord[2],vertices);
            p_voroInfo -> vCoord[NodeID] = vertices;
            c.face_vertices(vList);
            c.neighbors ( p_voroInfo -> nbList[i]);
            c.face_areas( p_voroInfo -> area[i]);
            v = n = 0;
            while (v<vList.size()) {
                fOrder = vList[v];
                p_voroInfo -> vList[NodeID].resize(n+1);
                v++;
                for (k=0; k<fOrder; k++) {
                    vID = vList[v];
                    p_voroInfo -> vList[NodeID][n].push_back(vID);
                    v++;
                }
                n++;
            }
            i++;
        }
        while (vl.inc());
}

bool VoroTesell::isRestrictedZone           (   const Point         testPoint,
                                                const Point         centre,
                                                const Vec           normVec,
                                                const double        offsetTol,
                                                const double        radius)
{
    const double    D = -(centre[0]*normVec[0]+centre[1]*normVec[1]+centre[2]*normVec[2]);
    double offset = testPoint[0]*normVec[0]+testPoint[1]*normVec[1]+testPoint[2]*normVec[2]+D;
    if (std::fabs(offset)>offsetTol) {
        return false;
    }
    Point   offCentre;
    offCentre[0] = centre[0] + normVec[0]*offset;
    offCentre[1] = centre[1] + normVec[1]*offset;
    offCentre[2] = centre[2] + normVec[2]*offset;
    Geometry geo;
    double dist = geo.dist(offCentre,testPoint);
    if (dist > radius) {
        return false;
    }
    return true;
}

bool VoroTesell::isRestrictedZone           (   const Point         testPoint,
                                                const double        width,
                                                const double        ld,
                                                const double        t)
{
    if ((testPoint[2] < 0.5*Nz*UnitLength - t/2.0 - ld) || (testPoint[2] > 0.5*Nz*UnitLength + t/2.0 + ld)) {
        return false;
    } else {
        if ((testPoint[0] > width- ld/2.0) && (testPoint[0] < Nx*UnitLength - (width- ld/2.0))) {
            return false;
        }
    }
    return true;
}

bool VoroTesell::isRestrictedZoneBH           (     const Point         testPoint,
                                                    const Point         bottom,
                                                    const double        radius,
                                                    const double        depth)
{
    if ((testPoint[2]>bottom[2])&&(testPoint[2]<bottom[2]+depth)) {
        double dist = std::sqrt((testPoint[0]-bottom[0])*(testPoint[0]-bottom[0])+(testPoint[1]-bottom[1])*(testPoint[1]-bottom[1]));
//        std::cout << dist << ' ';
        if (dist < radius) {
//            std::cout << testPoint[0] << ' ' << testPoint[1] << ' ' << testPoint[2] << std::endl;
            return true;
        }
    }
    return false;
}

bool VoroTesell::isInRange                  (   const Point         point)
{
    return ((point[0]>0.0)&&(point[0]<Nx*UnitLength)&&(point[1]>0.0)&&(point[1]<Ny*UnitLength)
            &&(point[2]>0.0)&&(point[2]<Nz*UnitLength));
}

unsigned VoroTesell::updateNodeCoordinationNo()
{
    unsigned nMax = 0;
    for (unsigned i=0; i<GNodeTab.size(); i++) {
        GNodeTab[i].n = GNodeTab[i].nbNodeID.size();
        if (GNodeTab[i].n > nMax) { nMax = GNodeTab[i].n; }
        for (unsigned k=0; k<GNodeTab[i].n; k++) {
            if (LatTab[GNodeTab[i].nbLatticeID[k]].k[0] > 200*Tiny) {
                GNodeTab[i].ne++;
            }
        }
    }

    dualOut << "Max coordination No = " << nMax << std::endl;
    return nMax;
}

void VoroTesell::updateNodeVoroVol		()
{
    unsigned NodeID;
    for (unsigned i = 0; i < p_voroInfo->IDList.size(); i++) {
        NodeID = p_voroInfo->IDList[i];
        GNodeTab[NodeID].v = p_voroInfo->vol[i];
        for (unsigned d=0; d<3; d++) {
            GNodeTab[NodeID].centroid[d] = GNodeTab[NodeID].coord[d]+p_voroInfo->centroid[i][d];
        }
    }
}

void VoroTesell::adjCellArea			(  std::vector<unsigned>*			p_NodeList)
{
    unsigned 	i,NodeID,LatticeID,k;
    double 		meanArea,areaSD;

    meanArea = getMeanCellArea (tLatticeNum/2,tLatticeNum);
    dualOut << "Mean Area = " << meanArea << std::endl;
    areaSD = getCellAreaSD (tLatticeNum/2,tLatticeNum,meanArea);
    dualOut << "Area SD = " << areaSD << std::endl;

    for (i = 0; i < p_NodeList -> size(); i++)
    {
        NodeID = (*p_NodeList)[i];
        for (k = 0; k < GNodeTab[NodeID].n ; k++ )
        {
            LatticeID = GNodeTab[NodeID].nbLatticeID[k];
            if (LatTab[LatticeID].area > meanArea + 2.0*areaSD)
            {
                LatTab[LatticeID].area = meanArea;
            }
        }
    }
}

double VoroTesell::getMeanCellArea (	const unsigned		start,
										const unsigned		end)
{
    unsigned LatticeID;
    double	 sum = 0.0;

    for (LatticeID = start; LatticeID < end; LatticeID++)
    {
        sum += LatTab[LatticeID].area;
    }

    return sum/tLatticeNum;
}

double VoroTesell::getCellAreaSD (	const unsigned		start,
									const unsigned		end,
									const double		mean)
{
    unsigned LatticeID;
    double	 sqSum = 0.0;

    for (LatticeID = start; LatticeID < end; LatticeID++)
    {
        sqSum += (LatTab[LatticeID].area -mean)*(LatTab[LatticeID].area -mean);
    }

    return sqrt(sqSum/tLatticeNum);
}

std::vector<unsigned> VoroTesell::getPreExCrackLatID ()
{
	return preExistingCrack;
}

std::vector<unsigned> VoroTesell::getWeakPlaneLatID ()
{
    return weakPlaneLatID;
}

std::array<unsigned,2> VoroTesell::getPxCrackPerimeterNum ()
{
	return pxPerimeterNum;
}

std::vector<unsigned>  VoroTesell::getBoreHoleNodeList ()
{
    return boreHoleNodeList;
}
void VoroTesell::reset (	container_poly* 	p_conp)
{
	GNodeTab.clear();
	LatTab.clear();
	nodeLists.free.clear();
	unsigned face[6] = {Face1,Face2,Face3,Face4,Face5,Face6};
	for (unsigned d=0; d<6; d++) {
		nodeLists.boundary[face[d]].clear();
		nodeLists.inBoundary[face[d]].clear();
	}
	p_conp->clear();
}

void VoroTesell::writeCrackNodeIDPair (     std::vector<NodeIDPair>*        p_conjNodeList)
{
    char fileName[255]; //filename
    std::sprintf(fileName,"%s/%sCrackNodeIDPair.txt",OutputFolder,OutputFilePrefix);
    std::ofstream file(fileName, std::ios::out);
    std::cout << "Generating CrackNodeIDPair file..."<< std::endl;
    for (unsigned i=0; i< p_conjNodeList->size(); i++) {
        file << (*p_conjNodeList)[i].i1 << ' ' <<  (*p_conjNodeList)[i].i2 << '\n';
    }
    std::cout << "CrackNodeIDPair file is generated"<< std::endl;
}

void VoroTesell::setRegFace ()
{
    unsigned face[6] = {Face3, Face5, Face2, Face4, Face1, Face6};
    for (unsigned d=0; d<3; d++) {
        isRegFace[face[d]] = false;
        isRegFace[face[d+3]] = false;
        if (std::fabs(InitStress[d])>Tiny) {
            isRegFace[face[d]] = true;
            isRegFace[face[d+3]] = true;
        }
    }

    for (unsigned f=0; f<6; f++) {
        bool isReg = false;
        for (unsigned d=0; d<6; d++) {
            isReg = (isReg || BoundaryMat[f][d]);
        }
        isRegFace[f] = isRegFace[f] || isReg ;
    }

    isRegFace[LoadFace] = true;

    if (Nx < 4) {
        isRegFace[face[0]] = false;
        isRegFace[face[1]] = false;
    }
    if (Ny < 4) {
        isRegFace[face[2]] = false;
        isRegFace[face[3]] = false;
    }
    if (Nz < 4) {
        isRegFace[face[4]] = false;
        isRegFace[face[5]] = false;
    }
}

void VoroTesell::writeWeakPlaneNodeIDPair ()
{
    char fileName[255]; //filename
    std::sprintf(fileName,"%s/%sWeakPlaneNodeIDPair.txt",OutputFolder,OutputFilePrefix);
    std::ofstream file(fileName, std::ios::out);
    std::cout << "Generating CrackNodeIDPair file..."<< std::endl;
    for (unsigned i=0; i< weakPlane.size(); i++) {
        file << weakPlane[i].i1 << ' ' <<  weakPlane[i].i2 << '\n';
    }
    std::cout << "WeakPlaneNodeIDPair file is generated"<< std::endl;
}
