/*
 * GLatRemoval.cpp
 *
 *  Created on: Nov 29, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "GLatRemoval.hpp"

using namespace std;

GLatRemoval::GLatRemoval()  {
	cntNodeRm = 0;
	cntLatRm = 0;
	cntLatShearRm = 0;
	isRecord = false;
	accSE = 0.0;
	double nominalVol = (Nx*Ny*Nz*UnitLength*UnitLength*UnitLength/tInNodeNum)/6.0;
	chkCoplanarTol = 1e-4;
}

GLatRemoval::GLatRemoval (	const double	_chkCoplanarTol) {
    	cntNodeRm = 0;
    	cntLatRm = 0;
    	cntLatShearRm = 0;
    	isRecord = false;
    	accSE = 0.0;
    	double nominalVol = (Nx*Ny*Nz*UnitLength*UnitLength*UnitLength/tInNodeNum)/6.0;
    	chkCoplanarTol = nominalVol*_chkCoplanarTol;
    	dualOut << "Nominal Volume = " << nominalVol << " chkCoplanarTol is set as " << chkCoplanarTol << std::endl;
    }

GLatRemoval::~GLatRemoval() {
    // TODO Auto-generated destructor stub
}

//Changed
unsigned GLatRemoval::remove	(	ConnectLat*				p_conLat,
                                    Vertices*				p_verti,
                                    FailStats*				p_fStats,
                                    PCG_Eigen*				p_EigenPCG,
                                    EnergyCal*              p_eCal,
                                    const unsigned			step,
                                    const double			time,
                                    double* 				maxLatNormForce ,
                                    const unsigned 			maxBreak)
{
    std::vector<double>	latNormForceList;
    std::vector<iDouble> exceed;
    std::vector<iTriple> exceedTri;
    Clock clock_lr;
    char fileName[255]; //filename
    std::sprintf(fileName,"%s/%sFailLat.txt",OutputSubFolder,OutputFilePrefix);
    std::ofstream file(fileName, std::ios::out | std::ios::app);
    static bool isFirst = true;
    if (isFirst){
        file << "step" << '\t' << "cntTensileStep" << '\t' << "cntShearStep" << '\t' << "cntLatRm" << '\t' <<  "cntLatShearRm" << '\t'
             << "cntLatReCon" << std::endl;
        isFirst = false;
    }

    clock_lr.start("LR");
    exceed.clear();
    if (KNum==1) {
        findFailLatSimple (&exceed,&latNormForceList);
    } else {
        if (LatFailureCriterion==1) {
            findFailLatMC(&exceedTri,&latNormForceList);
        } else if (LatFailureCriterion==2) {
            findFailLatTensionAndShear(&exceedTri,&latNormForceList);
        } else {
            findFailLatTensile (&exceed,&latNormForceList);
        }
    }
    unsigned maxNo;
    if (maxBreak<latNormForceList.size()/100) {
    	maxNo = maxBreak;
    } else {
    	maxNo = latNormForceList.size()/100;
    }
    partial_sort(latNormForceList.begin(),latNormForceList.begin()+maxNo,latNormForceList.end(),reverseSortDir);
    double sum = 0.;
    double lmax = 0.0;
    unsigned cntTensile = 0;
    unsigned cntShear = 0;
    for (unsigned i=0 ; i < maxNo ; i++) {
    	sum += latNormForceList[i];
    	if (latNormForceList[i]>lmax) {
    	    lmax = latNormForceList[i];
    	}
    }
    *maxLatNormForce = sum / maxNo;
    dualOut << "maxLatNormForce = " << lmax << '\n';
    dualOut << "maxLatNormForce avg among " << maxNo << " lattices = " << *maxLatNormForce << '\n';
    if (LatFailureCriterion==0) {
        if (exceed.size()==0) {
            dualOut << " No Lattice to be break !" << std::endl;
            return 0;
        } else {
            unsigned LatticeID;
            if (exceed.size() > maxBreak) {
                //For sorting from greatest to smallest
                partial_sort(exceed.begin(),exceed.begin()+maxBreak, exceed.end());
                    for (unsigned i = 0; i < maxBreak; i++) {
                        LatticeID = exceed[i].index;
                        recordFailLatStressState(LatticeID);
                        removeLattice (LatticeID,p_conLat,p_verti,p_fStats,p_EigenPCG,p_eCal,step,time);
                        cntTensile++;
                    }
                dualOut << " No. of Lattice exceed threshold = " << exceed.size()
                        << " No. of Lattice removed = " << maxBreak << std::endl;
            } else {
                dualOut << " No. of Lattice exceed threshold is less than or equal to the max allowed value, all " << exceed.size() << " lattices exceed the threshold are removed"
                        << std::endl;
                    for (unsigned i = 0; i < exceed.size(); i++) {
                        LatticeID = exceed[i].index;
                        recordFailLatStressState(LatticeID);
                        removeLattice (LatticeID,p_conLat,p_verti,p_fStats,p_EigenPCG,p_eCal,step,time);
                        cntTensile++;
                    }
            }
            dualOut << " Cumulative no. of Lattice removed = " << cntLatRm << std::endl;
            dualOut << " Cumulative no. of unstable node removed = " << cntNodeRm << std::endl;
            file << step << '\t' << cntTensile << '\t' << cntShear << '\t' << cntLatRm << '\t' <<  cntLatShearRm << '\t'
                    << p_conLat->reConList.size() << std::endl;
            clock_lr.get();
            return std::min((int) exceed.size(),(int) maxBreak);
        }
    } else {
        if (exceedTri.size()==0) {
            dualOut << " No Lattice to be break !" << std::endl;
            return 0;
        } else {
            unsigned LatticeID;
            if (exceedTri.size() > maxBreak) {
                partial_sort(exceedTri.begin(),exceedTri.begin()+maxBreak, exceedTri.end());
                for (unsigned i = 0; i < maxBreak; i++) {
                    LatticeID = exceedTri[i].index;
                    recordFailLatStressState(LatticeID);
                    if (exceedTri[i].data2<Tiny) {
                        removeLattice (LatticeID,p_conLat,p_verti,p_fStats,p_EigenPCG,p_eCal,step,time);
                        cntTensile++;
                    } else {
                        double factorS = ReConRatio;
                        if (Use_FrictionModel) {
                            factorS = 1.0/exceedTri[i].data1;
                        }
                        removeLatShear (LatticeID,factorS,exceedTri[i].data2,p_conLat,p_verti,p_fStats,p_EigenPCG,p_eCal,step,time);
                        cntShear++;
                    }
                }
                dualOut << " No. of Lattice exceed threshold = " << exceedTri.size()
                        << " No. of Lattice removed = " << maxBreak << std::endl;
            } else {
                dualOut << " No. of Lattice exceed threshold is less than or equal to the max allowed value, all "
                        << exceedTri.size() << " lattices exceed the threshold are removed"
                        << std::endl;
                for (unsigned i = 0; i < exceedTri.size(); i++) {
                    LatticeID = exceedTri[i].index;
                    recordFailLatStressState(LatticeID);
                    if (exceedTri[i].data2<Tiny) {
                        removeLattice (LatticeID,p_conLat,p_verti,p_fStats,p_EigenPCG,p_eCal,step,time);
                        cntTensile++;
                    } else {
                        double factorS = ReConRatio;
                        if (Use_FrictionModel) {
                            factorS = 1.0/exceedTri[i].data1;
                        }
                        removeLatShear (LatticeID,factorS,exceedTri[i].data2,p_conLat,p_verti,p_fStats,p_EigenPCG,p_eCal,step,time);
                        cntShear++;
                    }
                }
            }
            dualOut << " Cumulative no. of Lattice removed = " << cntLatRm << endl;
            dualOut << " Cumulative no. of Lattice removed (Shear only) = " << cntLatShearRm << std::endl;
            dualOut << " Cumulative no. of unstable node removed = " << cntNodeRm << endl;
            file << step << '\t' << cntTensile << '\t' << cntShear << '\t' << cntLatRm << '\t' <<  cntLatShearRm << '\t'
                 << p_conLat->reConList.size() << std::endl;
            clock_lr.get();
            return std::min((int) exceedTri.size(),(int) maxBreak);
        }
    }
}


unsigned GLatRemoval::remove	(	ConnectLat*				p_conLat,
                                    Vertices*				p_verti,
                                    PCG_Eigen*				p_EigenPCG,
                                    const unsigned 			maxBreak,
                                    const unsigned			step,
                                    const double			time)
{
    std::vector<double>	latNormForceList;
    std::vector<iDouble> exceed;
    Clock clock_lr;

    clock_lr.start("LR");
    exceed.clear();
    if (KNum==1) {
        findFailLatSimple(&exceed,&latNormForceList);
    } else {
    	findFailLatResultant(&exceed,&latNormForceList);
    }
    unsigned maxNo;
    if (maxBreak<latNormForceList.size()/100) {
    	maxNo = maxBreak;
    } else {
    	maxNo = latNormForceList.size()/100;
    }
    partial_sort(latNormForceList.begin(),latNormForceList.begin()+maxNo,latNormForceList.end(),reverseSortDir);

    if (exceed.size()==0) {
        dualOut << " No Lattice to be break !" << endl;
        return 0;
    } else {
		unsigned LatticeID;
		if (exceed.size() > maxBreak) {
			//For sorting from greatest to smallest
			partial_sort(exceed.begin(),exceed.begin()+maxBreak, exceed.end());
			for (unsigned i = 0; i < maxBreak; i++) {
				LatticeID = exceed[i].index;
				removeLattice (LatticeID,p_conLat,p_verti,p_EigenPCG,step,time);
			}
			dualOut << " No. of Lattice exceed threshold = " << exceed.size()
					<< " No. of Lattice removed = " << maxBreak << endl;
		} else {
			dualOut << " No. of Lattice exceed threshold is less than or equal to the max allowed value, all " << exceed.size() << " lattices exceed the threshold are removed"
					<< std::endl;
			for (unsigned i = 0; i < exceed.size(); i++) {
				LatticeID = exceed[i].index;
				removeLattice (LatticeID,p_conLat,p_verti,p_EigenPCG,step,time);
			}
		}
		dualOut << " Cumulative no. of Lattice removed = " << cntLatRm << endl;
		dualOut << " Cumulative no. of unstable node removed = " << cntNodeRm << endl;
        clock_lr.get();
        return std::min((int) exceed.size(),(int) maxBreak);
    }
}

void GLatRemoval::findFailLatSimple	(	std::vector<iDouble>*	p_exceed,
                                        double* 				maxLatNorm)
{
    double 		elongation,norm;
    GLatForce l;
    for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
        elongation = l.getLatExtension(LatticeID);
        norm = elongation/LatTab[LatticeID].et[0];
        if (norm > *maxLatNorm) {
            *maxLatNorm = norm;
        }
        if (!LatTab[LatticeID].isBreak) {
			if (elongation > LatTab[LatticeID].et[0]) {
				iDouble newitem(elongation/LatTab[LatticeID].et[0],LatticeID);
				p_exceed -> push_back(newitem);
			}
			else if (elongation < LatTab[LatticeID].et[1]) {
				iDouble newitem(elongation / LatTab[LatticeID].et[1],LatticeID);
				p_exceed -> push_back(newitem);
			}
        }
    }
}

void GLatRemoval::findFailLatSimple	(	std::vector<iDouble>*	p_exceed,
                                        std::vector<double>* 	p_latNormForceList)
{
    unsigned 	LatticeID;
    double 		elongation,norm;

    GLatForce l;
    for (LatticeID = 0; LatticeID < tLatticeNum; LatticeID++)
    {
    	if (LatTab[LatticeID].isBreak==false) {
			elongation = l.getLatExtension(LatticeID);
			norm = max(elongation/LatTab[LatticeID].et[0],elongation/LatTab[LatticeID].et[1]);

			p_latNormForceList -> push_back(norm);

			if (!LatTab[LatticeID].isBreak) {
			    if (norm>1.0) {
			        iDouble newitem(elongation/LatTab[LatticeID].et[0],LatticeID);
			        p_exceed -> push_back(newitem);
			    }
			}
    	}
    }
}


void GLatRemoval::findFailLatTensile   (    std::vector<iDouble>*   p_exceed,
                                        	std::vector<double>*    p_latNormForceList)
{
    double      norm;

    PCG_Eigen pcg;
    GLatForce l;
    for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++)
    {
        if (LatTab[LatticeID].isBreak==false) {
        	Tensor memForce = pcg.getMemForce(LatticeID);
        	double e = -memForce[0]/LatTab[LatticeID].k[0];
        	norm = max(e/LatTab[LatticeID].et[0],e/LatTab[LatticeID].et[1]);
        	p_latNormForceList -> push_back(norm);
        	if (!LatTab[LatticeID].isBreak) {
        		if (e > LatTab[LatticeID].et[0]) {
        			iDouble newitem(e/LatTab[LatticeID].et[0],LatticeID);
        			p_exceed -> push_back(newitem);
        		} else if (e < LatTab[LatticeID].et[1]) {
        			iDouble newitem(e / LatTab[LatticeID].et[1],LatticeID);
        			p_exceed -> push_back(newitem);
        		}
        	}
        }
    }
}

void GLatRemoval::findFailLatMC         (    std::vector<iTriple>*   p_exceedTri,
                                             std::vector<double>*    p_latNormForceList)
{
    double      norm;
    PCG_Eigen pcg;
    GLatForce l;
//    const double phi = 30.0;
    const double tanPhi = std::tan(Phi*Pi/180.0);
    for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
        if ((LatTab[LatticeID].isBreak==false)&&(LatTab[LatticeID].et[0]<Huge)) {
            double tensileStr = LatTab[LatticeID].et[0]*LatTab[LatticeID].k[0]/LatTab[LatticeID].area;
            double c = LatTab[LatticeID].shearStr;
            Tensor memForce = pcg.getMemForce(LatticeID);
            double sigma_n = memForce[0]/LatTab[LatticeID].area;
            double sigma_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
            double shearStr = c+std::max(sigma_n*tanPhi,0.0);
            norm = std::max(-sigma_n/tensileStr,sigma_s/shearStr);
            p_latNormForceList -> push_back(norm);
            if (!LatTab[LatticeID].isBreak) {
                if (norm > 1.0) {
                    double fShear;
                    if ((-sigma_n/tensileStr) > 1.0) {
                        fShear = 0.0;
                    } else {
                        fShear = shearStr;
                    }
                    iTriple newitem(norm,fShear,LatticeID);
                    p_exceedTri -> push_back(newitem);

                }
            }
        }
    }
}

void GLatRemoval::findFailLatTensionAndShear         (      std::vector<iTriple>*   p_exceedTri,
                                                            std::vector<double>*    p_latNormForceList)
{
    double      norm;
    PCG_Eigen pcg;
    GLatForce l;
    for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
        if ((LatTab[LatticeID].isBreak==false)&&(LatTab[LatticeID].et[0]<Huge)) {
            double tensileStr = LatTab[LatticeID].et[0]*LatTab[LatticeID].k[0]/LatTab[LatticeID].area;
            Tensor memForce = pcg.getMemForce(LatticeID);
            double sigma_n = memForce[0]/LatTab[LatticeID].area;
            double sigma_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
            norm = std::max(-sigma_n/tensileStr,sigma_s/LatTab[LatticeID].shearStr);
            p_latNormForceList -> push_back(norm);
            if (!LatTab[LatticeID].isBreak) {
                if (norm > 1.0) {
                    double fShear;
                    if ( (-sigma_n/tensileStr) > 1.0) {
                        fShear = 0.0;
                    } else {
                        fShear = LatTab[LatticeID].shearStr;
                    }
                    iTriple newitem(norm,fShear,LatticeID);
                    p_exceedTri -> push_back(newitem);
                }
            }
        }
    }
}

void GLatRemoval::recordFailLatStressState  (   const unsigned      LatticeID)
{
    char fileName[255];
    std::sprintf(fileName,"%s/%sfLatStress.txt",OutputSubFolder,OutputFilePrefix);
    std::ofstream file(fileName, std::ios::out | std::ios::app);
    PCG_Eigen pcg;
    Tensor memForce = pcg.getMemForce(LatticeID);
    double sigma_n = memForce[0]/LatTab[LatticeID].area;
    double sigma_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
    file << sigma_n << '\t' << sigma_s << std::endl;
}

void GLatRemoval::findFailLatResultant   (    	std::vector<iDouble>*   p_exceed,
                                        		std::vector<double>*    p_latNormForceList)
{
    double      norm;
    PCG_Eigen pcg;
    GLatForce l;
    for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
        if (LatTab[LatticeID].isBreak==false) {
        	Tensor memForce = pcg.getMemForce(LatticeID);
        	double e = std::sqrt(memForce[0]*memForce[0]+memForce[1]*memForce[1]
        	            +memForce[2]*memForce[2])/LatTab[LatticeID].k[0];
        	norm = e/LatTab[LatticeID].et[0];
        	p_latNormForceList -> push_back(norm);
        	if (!LatTab[LatticeID].isBreak) {
        		if ((e > LatTab[LatticeID].et[0])&&(memForce[0]<0.0)) {
        			iDouble newitem(e/LatTab[LatticeID].et[0],LatticeID);
        			p_exceed -> push_back(newitem);
        		} else if ((e < LatTab[LatticeID].et[1])&&(memForce[0]>0.0)) {
        			iDouble newitem(e / -LatTab[LatticeID].et[1],LatticeID);
        			p_exceed -> push_back(newitem);
        		}
        	}
        }
    }
}

//Changed
void GLatRemoval::checkNbNode	(	const unsigned 			LatticeID,
									const unsigned			step,
									const double			time,
                                    ConnectLat*				p_conLat,
                                    Vertices*				p_verti,
                                    PCG_Eigen*				p_EigenPCG)
{
    unsigned nbNodeID,nbLatticeID;
    std::vector<unsigned> rmDirList;
    std::vector<Point>		chkPoint;
    Geometry geo;
    unsigned minNb = 3;
    if (DofPerNode==6) {
    	minNb = 1;
    }
    for (unsigned i=0; i<2; i++)
    {
        nbNodeID = LatTab[LatticeID].nb[i];
        if (GNodeTab[nbNodeID].type!= Type_Unstable) {
        	if (GNodeTab[nbNodeID].ne<minNb) {
        		for (unsigned k = 0 ; k < GNodeTab[nbNodeID].n; k++) {
        			nbLatticeID = GNodeTab[nbNodeID].nbLatticeID[k];
        		    if (!(LatTab[nbLatticeID].isBreak)) {
        		    	rmDirList.push_back(k);
        		    }
        		}
        		removeUnstableNode (nbNodeID,step,time,p_conLat,p_verti,p_EigenPCG,rmDirList);
        	} else if (minNb==3) {
        		chkPoint.push_back(GNodeTab[nbNodeID].coord);
        		for (unsigned k = 0 ; k < GNodeTab[nbNodeID].n; k++) {
        		    nbLatticeID = GNodeTab[nbNodeID].nbLatticeID[k];
        		    if (LatTab[LatticeID].factor*LatTab[LatticeID].k[0]>10*Tiny) {
        		    	rmDirList.push_back(k);
        		    	if (LatTab[nbLatticeID].nb[0]!=nbNodeID) {
        		    		chkPoint.push_back(GNodeTab[LatTab[nbLatticeID].nb[0]].coord);
        		    	} else {
        		    		chkPoint.push_back(GNodeTab[LatTab[nbLatticeID].nb[1]].coord);
        		    	}
        		    }
        		}
        		if ((geo.isCoplanar(&chkPoint,chkCoplanarTol))&&(DofPerNode==3)&&(minNb!=1)) {
        			tripOut << "Co-planar nodes is detected and Node is unstable, ne = " << GNodeTab[nbNodeID].ne << std::endl;
        			removeUnstableNode (nbNodeID,step,time,p_conLat,p_verti,p_EigenPCG,rmDirList);
        		}
        	}
        }
        rmDirList.clear();
    }
}

// Changed
void GLatRemoval::removeUnstableNode	(	const unsigned 					NodeID,
											const unsigned					step,
											const double					time,
        									ConnectLat*						p_conLat,
        									Vertices*						p_verti,
        									PCG_Eigen*						p_EigenPCG,
        									const std::vector<unsigned> 	rmDirList)
{
    unsigned dir,LatticeID;
    std::vector<unsigned> rmLatList;
    if ((GNodeTab[NodeID].type!= Type_Unstable)&&(GNodeTab[NodeID].type!= Type_Virt)) {
        GNodeTab[NodeID].type = Type_Unstable;
        for (unsigned i = 0; i < rmDirList.size(); i++) {
            dir = rmDirList[i];
            LatticeID = GNodeTab[NodeID].nbLatticeID[dir];
            if (LatticeID>=LatTab.size()) {
            	tripOut << "[GLatRemoval::removeUnstableNode] - LatticeID = " << LatticeID << " is out of range" << std::endl;
            	tripOut << "NodeID = " << NodeID << " dir = " << rmDirList[i] << " GNodeTab[NodeID].nbLatticeID.size() = "
            			<< GNodeTab[NodeID].nbLatticeID.size() << std::endl;
            }
            if (LatTab[LatticeID].isBreak == false)
            {
                removeLattice(LatticeID,p_conLat,p_verti,p_EigenPCG,step,time);
                rmLatList.push_back(LatticeID);
            }
        }
        nodeLists.unstable.push_back(NodeID);
        eraseValue(nodeLists.free,NodeID);
        nodeLists.removeRestrain(NodeID);
        //Remove external force on unstable node
        for (unsigned i=0; i<DofPerNode; i++)
        {
            GNodeTab[NodeID].d[i] = 0.0;
            GNodeTab[NodeID].extF[i] = 0.0;
            GNodeTab[NodeID].extdF[i] = 0.0;

        }
        p_EigenPCG	-> updateDOF(true);
        tripOut << "Unstable node, ID = " << NodeID << " coordinates = ( "
        		<< GNodeTab[NodeID].coord[0] << " , " << GNodeTab[NodeID].coord[1] << " , " << GNodeTab[NodeID].coord[2] << " ) " << std::endl;
        tripOut << "cntNodeRm = " << cntNodeRm << std::endl;
        cntNodeRm++;
    }
}

void GLatRemoval::removeLattice	(	    const unsigned 			LatticeID,
                                        ConnectLat*				p_conLat,
                                        Vertices*				p_verti,
                                        FailStats*				p_fStats,
                                        PCG_Eigen*				p_EigenPCG,
                                        EnergyCal*              p_eCal,
                                        const unsigned			step,
                                        const double			time)
{
    GLatForce l;
    if ((LatticeID < tLatticeNum+tbLatticeNum)&& !(LatTab[LatticeID].isBreak))			//May have problem when introduce split nodes
    {
        LatTab[LatticeID].isBreak = true;
        GNodeTab[LatTab[LatticeID].nb[0]].ne--;
        GNodeTab[LatTab[LatticeID].nb[1]].ne--;
        p_EigenPCG		-> updateSM_full(LatticeID,Tiny,Tiny);
        LatTab[LatticeID].factor = Tiny;
        LatTab[LatticeID].factorShear = Tiny;
        p_conLat		-> putLattice(LatticeID,step,time,0,p_verti,p_conLat);
        p_fStats 		-> pushLat (LatticeID,step);
        p_eCal          -> accArea(LatTab[LatticeID].area);
        cntLatRm++;
        checkNbNode (LatticeID,step,time,p_conLat,p_verti,p_EigenPCG);
    }
}

void GLatRemoval::removeLatShear     (   const unsigned          LatticeID,
                                         const double            factorS,
                                         const double            fShear,
                                         ConnectLat*             p_conLat,
                                         Vertices*               p_verti,
                                         FailStats*              p_fStats,
                                         PCG_Eigen*              p_EigenPCG,
                                         EnergyCal*              p_eCal,
                                         const unsigned          step,
                                         const double            time)
{
    if ((LatticeID < tLatticeNum+tbLatticeNum)&& !(LatTab[LatticeID].isBreak)) {
        if (true) {
            LatTab[LatticeID].isBreak = true;
            GNodeTab[LatTab[LatticeID].nb[0]].ne--;
            GNodeTab[LatTab[LatticeID].nb[1]].ne--;
            p_EigenPCG      -> updateSM_full(LatticeID,ReConRatio,factorS);
            LatTab[LatticeID].factor = ReConRatio;
            LatTab[LatticeID].factorShear = factorS;
            LatTab[LatticeID].fShear = fShear;
        }
        p_conLat        -> putLattice(LatticeID,step,time,1,p_verti,p_conLat);
        p_fStats        -> pushLat (LatticeID,step);
        p_eCal          -> accArea(LatTab[LatticeID].area);
        cntLatRm++;
        cntLatShearRm++;
    }
}

void GLatRemoval::removeLattice	(	const unsigned 			LatticeID,
                                    ConnectLat*				p_conLat,
                                    Vertices*				p_verti,
                                    PCG_Eigen*				p_EigenPCG,
                                    const unsigned			step,
                                    const double			time)
{
    GLatForce l;
    double e;
    if ((LatticeID < tLatticeNum+tbLatticeNum)&& !(LatTab[LatticeID].isBreak) ) {
        e = l.getLatElongation(LatticeID);
        accSE += 0.5*LatTab[LatticeID].k[0]*e*e;
        LatTab[LatticeID].isBreak = true;

        GNodeTab[LatTab[LatticeID].nb[0]].ne--;
        GNodeTab[LatTab[LatticeID].nb[1]].ne--;
        p_EigenPCG		-> updateSM_full (LatticeID,Tiny,Tiny);
        LatTab[LatticeID].factor = Tiny;
        LatTab[LatticeID].factorShear = Tiny;
        p_conLat		-> putLattice(LatticeID,step,time,0,p_verti,p_conLat);
        cntLatRm++;
        checkNbNode (LatticeID,step,time,p_conLat,p_verti,p_EigenPCG);
    }
}

//For creating pre-existing fracture only
void GLatRemoval::removeUnstableNode	(	const unsigned 				NodeID,
        									NodeLists*					p_nodeLists,
        									std::vector<unsigned>*		p_latList)
{
    if ((GNodeTab[NodeID].type!= Type_Unstable)&&(GNodeTab[NodeID].type!= Type_Virt)) {
        GNodeTab[NodeID].type = Type_Unstable;
        for (unsigned i = 0; i < GNodeTab[NodeID].n; i++) {
            unsigned LatticeID = GNodeTab[NodeID].nbLatticeID[i];
            if (LatTab[LatticeID].isBreak == false) {
                removeLattice(LatticeID);
                p_latList->push_back(LatticeID);
            }
        }
        p_nodeLists	->unstable.push_back(NodeID);
        eraseValue(p_nodeLists->free,NodeID);
        p_nodeLists	->removeRestrain(NodeID);
        tripOut << "cntNodeRm = " << cntNodeRm << std::endl;
        cntNodeRm++;
    }
}

//For creating pre-existing fracture only
void GLatRemoval::removeLattice	(	const unsigned 			LatticeID)
{
    LatTab[LatticeID].isBreak = true;

    GNodeTab[LatTab[LatticeID].nb[0]].ne--;
    GNodeTab[LatTab[LatticeID].nb[1]].ne--;
    LatTab[LatticeID].k[0] = Tiny;
}



void GLatRemoval::disConnectLat	(	Vertices*						p_verti,
									ConnectLat*						p_conLat,
									PCG_Eigen*						p_EigenPCG,
									const unsigned					step,
									const double					time,
									const std::vector<unsigned>*	p_latList)
{

	for (unsigned i=0; i<p_latList->size(); i++) {
		if (std::find(p_conLat->reConList.begin(),p_conLat->reConList.end(),
				(*p_latList)[i])!=p_conLat->reConList.end()) {
			disConnectLatSingle((*p_latList)[i],step,time,p_verti,p_conLat,p_EigenPCG);
			p_conLat->reConList.erase((*p_latList)[i]);
		}
	}
}

bool GLatRemoval::updateReConLatShearStiffness  (    ConnectLat*                 p_conLat,
                                                     PCG_Eigen*                  p_EigenPCG,
                                                     Vertices*                   p_verti,
                                                     const unsigned              step,
                                                     const double                time)
{
    dualOut << "GLatRemoval::updateReConLatShearStiffness..." << std::endl;
    Clock   clock;
    clock.start("updateReConLatShearStiffness");
    bool convergence = true;
    unsigned disConLatCnt = 0;
    double tolerence = 0.05;
    for (std::set<unsigned>::iterator it=p_conLat->reConList.begin(); it!=p_conLat->reConList.end(); ++it) {
        unsigned LatticeID = *it;
        Tensor memForce = p_EigenPCG -> getMemForce(LatticeID);
        double sigma_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
        double factorS = std::max(std::min(1.0,LatTab[LatticeID].factorShear*LatTab[LatticeID].shearStr/sigma_s),0.0);
        p_EigenPCG -> updateSM_full (LatticeID,ReConRatio,factorS);
        LatTab[LatticeID].factor = ReConRatio;
        LatTab[LatticeID].factorShear = factorS;
        if ((factorS/LatTab[LatticeID].factorShear<1-tolerence)||(factorS/LatTab[LatticeID].factorShear>1+tolerence)) {
            convergence = false;
        }
    }
    clock.get();
    dualOut << "[GLatRemoval::updateReConLatShearStiffness] = No of reConLat disconnected = " << disConLatCnt << std::endl;
    return convergence;
}

bool GLatRemoval::updateReConLatShearStiffness  (    ConnectLat*                 p_conLat,
                                                     PCG_Eigen*                  p_EigenPCG,
                                                     Vertices*                   p_verti,
                                                     const double                gamma,
                                                     const unsigned              step,
                                                     const double                time)
{
    dualOut << "GLatRemoval::updateReConLatShearStiffnessCap..." << std::endl;
    bool convergence = true;
    unsigned disConLatCnt = 0;
    const double tol = 0.005;
    Clock   clock;
    if (false) {
        clock.start("Serial code");
        for (std::set<unsigned>::iterator it=p_conLat->reConList.begin(); it!=p_conLat->reConList.end(); ++it) {
            unsigned LatticeID = *it;
            Tensor memForce = p_EigenPCG -> getMemForce(LatticeID);
            double sigma_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
            double factorS = gamma*std::max(std::min(1.0,gamma*LatTab[LatticeID].shearStr/sigma_s),0.0);
            if (std::fabs((factorS-LatTab[LatticeID].factorShear)/LatTab[LatticeID].factorShear)>tol) {
                p_EigenPCG -> updateSM_full (LatticeID,ReConRatio,factorS);
                LatTab[LatticeID].factor = ReConRatio;
                LatTab[LatticeID].factorShear = factorS;
            }
        }
        clock.get();
    } else {
        clock.start("Parallel code");
        std::vector<unsigned> reConList{std::begin(p_conLat->reConList), std::end(p_conLat->reConList)};
        p_EigenPCG -> updateSM_full(reConList,ReConRatio,gamma);
        clock.get();
    }
    dualOut << "[GLatRemoval::updateReConLatShearStiffnessCap] = No of reConLat disconnected = " << disConLatCnt << std::endl;
    return convergence;
}

unsigned GLatRemoval::disConnectLat	(	Vertices*					p_verti,
										ConnectLat*					p_conLat,
										PCG_Eigen*					p_EigenPCG,
										const unsigned				step,
										const double				time,
										const unsigned              maxNo)
{
	GLatForce		glforce;
	static unsigned cntAcc = 0;
	unsigned		cnt = 0;
	std::vector<iDouble> exceed;
	std::vector<unsigned>   latList;
	Clock clock;
	clock.start("disConnectLat");
	std::cout << "Disconnecting Lattice..." << std::endl;
	for (std::set<unsigned>::iterator it=p_conLat->reConList.begin(); it!=p_conLat->reConList.end();++it) {
	    double apecture = glforce.getCrackOpening(*it);
		if (apecture > p_conLat->minApecture) {
			    iDouble newitem(apecture/LatTab[*it].et[0],*it);
			    exceed.push_back(newitem);
		}
	}
	if (exceed.size() > maxNo) {
	    partial_sort(exceed.begin(),exceed.begin()+maxNo, exceed.end());
	    for (unsigned i = 0; i < maxNo; i++) {
	        unsigned LatticeID = exceed[i].index;
	        disConnectLatSingle(LatticeID,step,time,p_verti,p_conLat,p_EigenPCG);
	        latList.push_back(LatticeID);
	        cnt++;
	    }
	    p_EigenPCG -> updateSM_full(latList,Tiny,Tiny);
	    dualOut << " No. of Lattice exceed threshold = " << exceed.size()
	            << " No. of Lattice reconnected = " << maxNo << std::endl;
	} else {
	    for (unsigned i = 0; i < exceed.size(); i++) {
	        unsigned LatticeID = exceed[i].index;
	        disConnectLatSingle(LatticeID,step,time,p_verti,p_conLat,p_EigenPCG);
	        latList.push_back(LatticeID);
	        cnt++;
	    }
	    p_EigenPCG -> updateSM_full(latList,Tiny,Tiny);
	    dualOut << " No. of Lattice exceed threshold is less than or equal to the max allowed value, all " << exceed.size() << " lattices exceed the threshold are reconnected"
	            << std::endl;
	}
	cntAcc +=cnt;
	if (cnt!=0) {
	    dualOut << "No of re-connected lattice disconnected in this step = " << cnt << std::endl
	    		<< "Total re-connected lattice disconnected = " << p_conLat->disConList.size() << std::endl;
	} else {
		dualOut << "No re-connected lattice disconnected" << std::endl;
	}
	clock.get();
	return exceed.size();
}
void GLatRemoval::disConnectLatSingle	(	const unsigned			LatticeID,
											const unsigned			step,
											const double			time,
											Vertices*				p_verti,
											ConnectLat*				p_conLat,
											PCG_Eigen*				p_EigenPCG)
{
	unsigned nbNodeID;
	for (unsigned i=0; i<2; i++) {
		nbNodeID = LatTab[LatticeID].nb[i];
		if (GNodeTab[nbNodeID].ne<1) {
			tripOut << "[GLatRemoval::disConnectLatSingle] - LatticeID = " << LatticeID << " nbNodeID = "
					<< nbNodeID << " ne < 1, please check!!" << std::endl;
		} else {
			GNodeTab[nbNodeID].ne--;
			checkNbNode	(LatticeID,step,time,p_conLat,p_verti,p_EigenPCG);
		}
	}
	p_conLat->reConList.erase(LatticeID);
	p_conLat->disConList.insert(LatticeID);
}

unsigned GLatRemoval::reConnectLat	(	ConnectLat*				p_conLat,
										PCG_Eigen*				p_EigenPCG,
										const unsigned          maxNo)
{
	static unsigned cntAcc = 0;
	unsigned		cnt = 0;
    GLatForce		glforce;
    std::vector<iDouble> exceed;
    std::vector<unsigned> latList;
    std::vector<double> shearRatioList;
    std::cout << "Reconnecting Lattice..." << std::endl;
    Clock clock;
    clock.start("reConnectLat");
    for (std::set<unsigned>::iterator it=p_conLat->disConList.begin(); it!=p_conLat->disConList.end();++it) {
        double apecture = glforce.getCrackOpening(*it);
        if ((apecture<p_conLat->minApecture)&&(!p_conLat->isOneEndUnstable(*it))) {
            if (!std::binary_search(nonReConList.begin(),nonReConList.end(),*it)) {
                iDouble newitem(-apecture/LatTab[*it].et[0],*it);
                exceed.push_back(newitem);
            }
        }
    }
    if (exceed.size() > maxNo) {
        partial_sort(exceed.begin(),exceed.begin()+maxNo, exceed.end());
        for (unsigned i = 0; i < maxNo; i++) {
            unsigned LatticeID = exceed[i].index;
            shearRatioList.push_back(reConnectLatSingle(LatticeID,p_EigenPCG,p_conLat));
            latList.push_back(LatticeID);
            cnt++;
        }
        p_EigenPCG -> updateSM_full(latList,ReConRatio,shearRatioList);
        dualOut << " No. of Lattice exceed threshold = " << exceed.size()
                << " No. of Lattice disconnected = " << maxNo << std::endl;
    } else {
        for (unsigned i = 0; i < exceed.size(); i++) {
            unsigned LatticeID = exceed[i].index;
            shearRatioList.push_back(reConnectLatSingle(LatticeID,p_EigenPCG,p_conLat));
            latList.push_back(LatticeID);
            cnt++;
        }
        p_EigenPCG -> updateSM_full(latList,ReConRatio,shearRatioList);
        dualOut << " No. of Lattice exceed threshold is less than or equal to the max allowed value, all " << exceed.size() << " lattices exceed the threshold are disconnected"
                << std::endl;
    }
    cntAcc +=cnt;
    if (cnt!=0) {
    	dualOut << "No of broken lattice re-connected in this step = " << cnt
    			<< " Total re-connected lattice = " << p_conLat->reConList.size() << std::endl;
    } else {
    	dualOut << "NO broken lattice re-connected in this step" << std::endl;
    }
    clock.get();
    return exceed.size();
}

double GLatRemoval::reConnectLatSingle	(	const unsigned		LatticeID,
											PCG_Eigen*			p_EigenPCG,
											ConnectLat*         p_conLat)
{
    if (LatticeID==PriInjLattice) {
        dualOut << "The primary injection point cannot be reconnected" << std::endl;
        return 0.0;
    }
    double ratioShear = ReConRatio;
    if (Use_FrictionModel) {
        Tensor memForce = p_EigenPCG -> getMemForce(LatticeID);
        double sigma_n = memForce[0]/LatTab[LatticeID].area;
        double shearStr = LatTab[LatticeID].shearStr+sigma_n*std::tan(Phi*Pi/180.0);
        double sigma_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
        ratioShear = LatTab[LatticeID].factorShear*std::min(1.0,shearStr/sigma_s);
    }
	p_conLat->reConList.insert(LatticeID);
	p_conLat->disConList.erase(LatticeID);
	for (unsigned i=0; i<2; i++) {
		unsigned nbNodeID = LatTab[LatticeID].nb[i];
		if (GNodeTab[nbNodeID].n-GNodeTab[nbNodeID].ne>=1) {
			GNodeTab[nbNodeID].ne++;
		} else {
			dualOut << "[GLatRemoval::reConnectLatSingle] - LatticeID = " << LatticeID << " nbNodeID = "
					<< nbNodeID << " n - ne < 1, please check!!" << std::endl;
		}
	}
	return ratioShear;
}

unsigned GLatRemoval::checkAllLatticeSM     (   PCG_Eigen*          p_EigenPCG,
                                                ConnectLat*         p_conLat)
{
    std::vector<unsigned>   abNormalLatList = p_EigenPCG->  testSM();
    static std::vector<unsigned> oldList;
    static std::vector<unsigned> oldType;
    for (unsigned i=0; i<oldList.size(); i++) {
        unsigned LatticeID = oldList[i];
        p_conLat -> setType(LatticeID,2);
        GNodeTab[LatTab[LatticeID].nb[0]].type = 0;
        GNodeTab[LatTab[LatticeID].nb[1]].type = 0;
    }
    oldType.clear();
    for (unsigned i=0; i<abNormalLatList.size(); i++) {
        unsigned LatticeID = abNormalLatList[i];
        oldType.push_back(p_conLat -> setType(LatticeID,3));
        GNodeTab[LatTab[LatticeID].nb[0]].type = 999;
        GNodeTab[LatTab[LatticeID].nb[1]].type = 999;
    }
    oldList = abNormalLatList;
    return abNormalLatList.size();
}

bool GLatRemoval::isLatAtBoundary 	(	const 	std::vector<iDouble> rmList)
{
    unsigned i,k,NodeID[2];
    list<Restrain>::iterator it;

    for (i=0; i<rmList.size(); i++)
    {
        for (k=0; k<2; k++)
        {
            NodeID[k] = LatTab[rmList[i].index].nb[k];
            it = find(nodeLists.restrain.begin(),nodeLists.restrain.end(),Restrain(NodeID[k]));
            if (it!=nodeLists.restrain.end())
            {
                return true;
            }
        }
    }
    return false;
}

double GLatRemoval::getMaxLatNormForce (	unsigned 					maxBreak)
{
	double 		elongation,norm;
	std::vector<double> latNormForceList;
	GLatForce l;
	for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
		elongation = l.getLatElongation(LatticeID);
	    norm = max(elongation/LatTab[LatticeID].et[0],elongation/LatTab[LatticeID].et[1]);
	    latNormForceList.push_back(norm);
	}
	unsigned maxNo;
	if (maxBreak<latNormForceList.size()/100) {
	  	maxNo = maxBreak;
	} else {
	  	maxNo = latNormForceList.size()/100;
	}

	partial_sort(latNormForceList.begin(),latNormForceList.begin()+maxNo,latNormForceList.end(),reverseSortDir);

	double sum = 0;
	for (unsigned i=0 ; i < maxNo ; i++) {
	    sum += latNormForceList[i];
	}
	return sum / maxNo;
}

void GLatRemoval::setNonReConLatList      (   std::vector<unsigned>*       p_latList)
{
    nonReConList.insert(nonReConList.end(),p_latList->begin(),p_latList->end());
    std::sort(nonReConList.begin(),nonReConList.end());
}

void GLatRemoval::setNonReConLat			(	unsigned 			LatticeID)
{
    nonReConList.push_back(LatticeID);
}

