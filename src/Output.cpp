/*
 * Output.cpp
 *
 *  Created on: Oct 2, 2013
 *      Author: john
 */
#include "Output.hpp"

using namespace std;

Output::Output (	char* 	inPath,
                    char* 	inPrefix)
{
    path = inPath;
    prefix = inPrefix;

    if (inPrefix[0]!='\0')
    {
        strcat(prefix,"-");
    }
}

void Output::writeAll					(	const unsigned 					step,
                                            const double                    extForce,
                                			FluidCal*						p_fCal,
                                			ConnectLat*						p_conLat,
                                			Vertices*						p_verti,
                                			PCG_Eigen*                      p_EigenPCG)
{
	writeNodeAndLattice(step,p_EigenPCG);
	writeFluidNetwork(step,p_fCal);
	writeFailNetwork(step,p_conLat,p_verti,p_EigenPCG);
	writeForceChain(step,extForce,p_EigenPCG);
	double minR = avgD/(std::min(std::min(UnitLength*Nx,UnitLength*Ny),UnitLength*Nz));
	if (((GeometryID==0)||(GeometryID==3))&&(CrackRadius>minR)) {
	    Lattice_fPerimeter (p_conLat,step);
	}
}

void Output::writeLatAndNetwork					(	const unsigned 					step,
                                                    const double                    extForce,
                                					ConnectLat*						p_conLat,
                                					Vertices*						p_verti,
                                					PCG_Eigen*                      p_EigenPCG)
{
	writeNodeAndLattice(step,p_EigenPCG);
	writeFailNetwork(step,p_conLat,p_verti,p_EigenPCG);
	writeForceChain(step,extForce,p_EigenPCG);
	double minR = avgD/(std::min(std::min(UnitLength*Nx,UnitLength*Ny),UnitLength*Nz));
	if (((GeometryID==0)||(GeometryID==3))&&(CrackRadius>minR)) {
	    Lattice_fPerimeter (p_conLat,step);
	}
}


void Output::writeNodeAndLattice		(	const unsigned 					step,
                                            PCG_Eigen*                      p_EigenPCG)
{
    char fileName[255]; //filename

    sprintf(fileName,"%s/%scombined%04d.vtk",path,prefix,step);

    ofstream file(fileName, ios::out);

    tDispNode = tNodeNum + (tbNodeNum+tvNodeNum)*Use_RanLatStrength;
    tDispLattice = tLatticeNum + (tbLatticeNum + tvLatticeNum) * isTestCase;

    dualOut << "Generating combined node and lattice file..."<<endl;

    if (file) {
        file.precision(5);
        file << scientific;
        file << "# vtk DataFile Version 3.0" << endl;
        file << "My lattices" << endl;
        file << "ASCII" << endl;
        file << "DATASET POLYDATA" << endl;
        file << "POINTS "<< tDispNode << " float" << endl;

        if (Disp_Scale>0.0) {
			for (unsigned NodeID=0; NodeID<tDispNode; NodeID++) {
			    file << GNodeTab[NodeID].coord[0]  + GNodeTab[NodeID].d[0]*Disp_Scale	<< " "
					 << GNodeTab[NodeID].coord[1]  + GNodeTab[NodeID].d[1]*Disp_Scale	<< " "
					 << GNodeTab[NodeID].coord[2]  + GNodeTab[NodeID].d[2]*Disp_Scale	<< endl;
			}
        } else {
        	for (unsigned NodeID=0; NodeID<tDispNode; NodeID++) {
        	file << GNodeTab[NodeID].coord[0]  << " "
        	     << GNodeTab[NodeID].coord[1] 	<< " "
        	     << GNodeTab[NodeID].coord[2] 	<< endl;
        	}
        }

        file << "LINES "<< tDispLattice << " " <<3*tDispLattice <<endl;
        for (unsigned k=0; k<tDispLattice; k++) {
        	file << "2 " << LatTab[k].nb[0] << " " << LatTab[k].nb[1] <<endl;
        }
        file << "POINT_DATA " << tDispNode << endl;
        Node_d (&file);
        Node_stress (&file);
        file << "CELL_DATA "<< tDispLattice <<endl;
        Lattice_fStress(&file,p_EigenPCG);
    }
    file.close();
    dualOut << "A combined node and lattice file is generated" << endl;
}

void Output::writeFluidNetwork		(	const unsigned 					step,
                                		FluidCal*						p_fCal)
{
	std::vector<Point>							fNode = p_fCal->getAllFNode();
	std::vector<std::pair<unsigned,unsigned> >	pipe = p_fCal->getAllPipe();
	tDispfNode 	= fNode.size();
	tDispPipe	= pipe.size();
	if ((tDispfNode==0)||(tDispPipe==0)) {
		return;
	}
	char fileName[255]; //filename
	sprintf(fileName,"%s/%sfluid%04d.vtk",path,prefix,step);
	ofstream file(fileName, ios::out);
    dualOut << "Generating Fluid Network file..."<<endl;
	if ((file))
	{
		file.precision(5);
	    file << scientific;
	    file << "# vtk DataFile Version 3.0" << endl;
	    file << "Fluid Network" << endl;
	    file << "ASCII" << endl;
	    file << "DATASET POLYDATA" << endl;
	    file << "POINTS "<< tDispfNode << " float" << endl;

	    for (unsigned fNodeID=0; fNodeID<tDispfNode; fNodeID++) {
	        file << fNode[fNodeID][0] << " "
	             << fNode[fNodeID][1] << " "
	             << fNode[fNodeID][2] << std::endl;
	    }

	    file << "LINES "<< tDispPipe << " " <<3*tDispPipe <<endl;
	    for (unsigned k=0; k<tDispPipe; k++) {
	    	file << "2 " << pipe[k].first << " " << pipe[k].second << std::endl;
	    }

	    file << "POINT_DATA " << tDispfNode << endl;
	    Fluid_pressure (&file,p_fCal);
	    Fluid_flowRate (&file,p_fCal);
	    Fluid_leakOff (&file,p_fCal);
	}
	file.close();
	dualOut << "A Fluid Network file is generated" << endl;
}

void Output::writeForceChain        (   const unsigned                  step,
                                        const double                    extForce,
                                        PCG_Eigen*                      p_EigenPCG)
{
    char fileName[255]; //filename
    sprintf(fileName,"%s/%sforceChain%04d.vtk",path,prefix,step);
    ofstream file(fileName, ios::out);
    dualOut << "Generating force chain file..."<<endl;
    if ((file))
    {
        file.precision(5);
        file << scientific;
        file << "# vtk DataFile Version 3.0" << endl;
        file << "Fluid Network" << endl;
        file << "ASCII" << endl;
        file << "DATASET POLYDATA" << endl;
        file << "POINTS "<< 2*tDispLattice << " float" << endl;

        for (unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            for (unsigned n=0; n<2; n++) {
                unsigned NodeID = LatTab[LatticeID].nb[n];
                file << GNodeTab[NodeID].coord[0]  + GNodeTab[NodeID].d[0]*Disp_Scale << " "
                     << GNodeTab[NodeID].coord[1]  + GNodeTab[NodeID].d[1]*Disp_Scale << " "
                     << GNodeTab[NodeID].coord[2]  + GNodeTab[NodeID].d[2]*Disp_Scale << std::endl;
            }
        }

        file << "LINES "<< tDispLattice << " " <<3*tDispLattice <<endl;
        for (unsigned k=0; k<tDispLattice; k++) {
            file << "2 " << 2*k << " " << 2*k+1 << std::endl;
        }

        file << "POINT_DATA " << 2*tDispLattice << std::endl;
        file << "SCALARS F_n float 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;

        std::vector<Tensor>     memForceList;
        for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            memForceList.push_back(p_EigenPCG -> getMemForce(LatticeID));
        }

        for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            file << memForceList[LatticeID][0] << std::endl
                 << memForceList[LatticeID][0] << std::endl;
        }

        file << "SCALARS F_s float 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;

        for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            double Fs = std::sqrt(memForceList[LatticeID][1]*memForceList[LatticeID][1]+
                            memForceList[LatticeID][2]*memForceList[LatticeID][2]);
            file << Fs << std::endl
                 << Fs << std::endl;
        }

        if (std::fabs(extForce-1.0)>Tiny) {
            file << "SCALARS F_n1 float 1" << std::endl;
            file << "LOOKUP_TABLE default" << std::endl;

            double aAvg = 0.0;
            for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
                aAvg += LatTab[LatticeID].area;
            }
            aAvg /= tDispLattice;

            for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
                file << memForceList[LatticeID][0]/(extForce*aAvg) << std::endl
                     << memForceList[LatticeID][0]/(extForce*aAvg) << std::endl;
            }

            file << "SCALARS F_s1 float 1" << std::endl;
            file << "LOOKUP_TABLE default" << std::endl;

            double absExtForce = std::fabs(extForce);
            for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
                double Fs = std::sqrt(memForceList[LatticeID][1]*memForceList[LatticeID][1]+
                            memForceList[LatticeID][2]*memForceList[LatticeID][2]);
                file << Fs/absExtForce << std::endl
                     << Fs/absExtForce << std::endl;
            }
        }
    }
    file.close();
    dualOut << "A force chain file is generated" << endl;
}

void Output::writeGNode		(	const unsigned 					step)
{
    char fileName[255]; //filename
    sprintf(fileName,"%s/%snodes%04d.vtk",path,prefix,step);
    ofstream nodeFile(fileName, ios::out);
    tDispNode = tNodeNum + (tbNodeNum+tvNodeNum*0)*Use_RanLatStrength;
    dualOut << "Generating GNode file..."<<endl;
    if (nodeFile) {
        nodeFile.precision(5);
        nodeFile << scientific;
        nodeFile << "# vtk DataFile Version 3.0" << endl;
        nodeFile << "My lattices" << endl;
        nodeFile << "ASCII" << endl;
        nodeFile << "DATASET UNSTRUCTURED_GRID" << endl;
        nodeFile << "POINTS "<< tDispNode << " float" << endl;

        for (unsigned NodeID=0; NodeID<tDispNode; NodeID++) {
            nodeFile << GNodeTab[NodeID].coord[0] + GNodeTab[NodeID].d[0]*Disp_Scale	<< " "
                     << GNodeTab[NodeID].coord[1] + GNodeTab[NodeID].d[1]*Disp_Scale	<< " "
                     << GNodeTab[NodeID].coord[2] + GNodeTab[NodeID].d[2]*Disp_Scale	<< endl;
        }
        nodeFile << "POINT_DATA " << tDispNode << endl;
        Node_d (&nodeFile);
    }
    nodeFile.close();
    dualOut << "A GNode file is generated" << endl;
}

void Output::writePriDir 		(	const unsigned 				step,
		                            ConnectLat*					p_conLat,
		                            Vertices*					p_verti)
{
	char fileName[255]; //filename
	dualOut << "Generating PriDir file..."<<endl;
	tDispLattice = tLatticeNum;
	sprintf(fileName,"%s/%sPriDir%04d.vtk",path,prefix,step);
	std::ofstream file(fileName, ios::out);
	std::vector<unsigned> latLists = p_conLat->getAllGroupVec();
	if (file) {
		file.precision(5);
		file << scientific;
		file << "# vtk DataFile Version 3.0" << endl;
		file << "My lattices" << endl;
		file << "ASCII" << endl;
		file << "DATASET UNSTRUCTURED_GRID" << endl;
		file << "POINTS "<< latLists.size() << " float" << endl;
		std::vector<Point>	cList;
		std::vector<BasisVecs> bvList;
		Geometry geo;
		unsigned LatticeID;
		for (unsigned i=0; i<latLists.size(); i++) {
			LatticeID = latLists[i];
			Surface v = p_verti-> getVertiCoordArray (LatticeID);
			Point centroid = geo.getCentroid(LatTab[LatticeID].area,v);
			cList.push_back(centroid);
			file << centroid[0] << ' ' << centroid[1] << ' ' << centroid[2] << '\n';
			BasisVecs basisVec = geo.getBasisVec(GNodeTab[LatTab[LatticeID].nb[0]].coord,
											 GNodeTab[LatTab[LatticeID].nb[1]].coord,v,centroid);
			bvList.push_back(basisVec);
		}
		file << "POINT_DATA " << latLists.size() << endl;
		file << "VECTORS PDir_0 float" << endl;
		for (unsigned i=0; i<latLists.size(); i++) {
			file << bvList[i][0][0] << " " << bvList[i][0][1] << " " << bvList[i][0][2] << endl;
		}
		file << "VECTORS PDir_1 float" << endl;
		for (unsigned i=0; i<latLists.size(); i++) {
			file << bvList[i][1][0] << " " << bvList[i][1][1] << " " << bvList[i][1][2] << endl;
		}
		file << "VECTORS PDir_2 float" << endl;
		for (unsigned i=0; i<latLists.size(); i++) {
			file << bvList[i][2][0] << " " << bvList[i][2][1] << " " << bvList[i][2][2] << endl;
		}
	}

}

void Output::writeCentroid ()
{
    unsigned NodeID;
    char fileName[255]; //filename

    sprintf(fileName,"%s/%scentroid.vtk",path,prefix);

    ofstream centFile(fileName, ios::out);

    tDispNode = tNodeNum + (tbNodeNum+tvNodeNum*0)*Use_RanLatStrength;

    dualOut << "Generating Centroid file..."<<endl;

    if (centFile)
    {
        centFile.precision(5);
        centFile << scientific;
        centFile << "# vtk DataFile Version 3.0" << endl;
        centFile << "My lattices" << endl;
        centFile << "ASCII" << endl;
        centFile << "DATASET UNSTRUCTURED_GRID" << endl;
        centFile << "POINTS "<< tDispNode << " float" << endl;

        for (NodeID=0; NodeID<tDispNode; NodeID++)
        {
            centFile << GNodeTab[NodeID].centroid[0] 	<< " "
                     << GNodeTab[NodeID].centroid[1] 	<< " "
                     << GNodeTab[NodeID].centroid[2] 	<< endl;
        }

        centFile << "POINT_DATA " << tDispNode << endl;
        centFile << "SCALARS NodeID int 1" << endl;
        centFile << "LOOKUP_TABLE default" << endl;

        for (NodeID=0; NodeID<tDispNode; NodeID++)
        {
            centFile << NodeID << endl;
        }
    }
    centFile.close();
    dualOut << "A Centroid file is generated" << endl;
}

void Output::writeLatAxes       (   const unsigned				step,
									ConnectLat*                 p_conLat)
{
    char fileName[255]; //filename

    dualOut << "Generating Lattice local axes file..."<<endl;
    tDispLattice = tLatticeNum;

    sprintf(fileName,"%s/%sfLatAxes%04d.vtk",path,prefix,step);
    ofstream file(fileName, ios::out);

    std::vector<std::vector<unsigned> > latLists = p_conLat->getAllGroup();
    unsigned nNum = p_conLat->getLatNum();
    if (file){
        file.precision(5);
        file << scientific;
        file << "# vtk DataFile Version 3.0" << std::endl;
        file << "Local Axes of failed lattices" << std::endl;
        file << "ASCII" << endl;
        file << "DATASET UNSTRUCTURED_GRID" << std::endl;
        file << "POINTS "<< nNum << " float" << std::endl;
        for (unsigned i=0; i<latLists.size(); i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
                file << LatTab[latLists[i][j]].centroid[0] << ' '
                     << LatTab[latLists[i][j]].centroid[1] << ' '
                     << LatTab[latLists[i][j]].centroid[2] << std::endl;
            }
        }

        file << "POINT_DATA " << nNum << std::endl;
        file << "VECTORS PDir_0 float" << endl;
        for (unsigned i=0; i<latLists.size(); i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
                file << LatTab[latLists[i][j]].axes[0][0] << ' '
                     << LatTab[latLists[i][j]].axes[0][1] << ' '
                     << LatTab[latLists[i][j]].axes[0][2] << std::endl;
            }
        }
        file << "VECTORS PDir_1 float" << endl;
        for (unsigned i=0; i<latLists.size(); i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
                file << LatTab[latLists[i][j]].axes[1][0] << ' '
                     << LatTab[latLists[i][j]].axes[1][1] << ' '
                     << LatTab[latLists[i][j]].axes[1][2] << std::endl;
            }
        }
        file << "VECTORS PDir_2 float" << endl;
        for (unsigned i=0; i<latLists.size(); i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
                file << LatTab[latLists[i][j]].axes[2][0] << ' '
                     << LatTab[latLists[i][j]].axes[2][1] << ' '
                     << LatTab[latLists[i][j]].axes[2][2] << std::endl;
            }
        }
    }
}

void Output::writeIntersect       (     const unsigned              step,
                                        ConnectLat*                 p_conLat)
{
    char fileName[255]; //filename

    dualOut << "Generating Intersect file..."<<endl;
    tDispLattice = tLatticeNum;

    sprintf(fileName,"%s/%sfLatIntersect%04d.vtk",path,prefix,step);
    ofstream file(fileName, ios::out);

    std::vector<std::vector<unsigned> > latLists = p_conLat->getAllGroup();
    unsigned nNum = p_conLat->getLatNum();
    if (file){
        file.precision(5);
        file << scientific;
        file << "# vtk DataFile Version 3.0" << std::endl;
        file << "Local Axes of failed lattices" << std::endl;
        file << "ASCII" << endl;
        file << "DATASET UNSTRUCTURED_GRID" << std::endl;
        file << "POINTS "<< nNum << " float" << std::endl;
        for (unsigned i=0; i<latLists.size(); i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
                file << LatTab[latLists[i][j]].intersect[0] << ' '
                     << LatTab[latLists[i][j]].intersect[1] << ' '
                     << LatTab[latLists[i][j]].intersect[2] << std::endl;
            }
        }

        file << "POINT_DATA " << nNum << std::endl;
        file << "SCALARS yc float" << endl;
        file << "LOOKUP_TABLE default" << endl;
        for (unsigned i=0; i<latLists.size(); i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
                file << LatTab[latLists[i][j]].offset[0] <<  std::endl;
            }
        }
        file << "SCALARS zc float" << endl;
        file << "LOOKUP_TABLE default" << endl;
        for (unsigned i=0; i<latLists.size(); i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
                file << LatTab[latLists[i][j]].offset[1] << std::endl;
            }
        }
    }
}

void Output::writeLatAxesAll       ()
{
    char fileName[255]; //filename

    dualOut << "Generating Lattice local axes all file..."<<endl;
    tDispLattice = tLatticeNum;

    sprintf(fileName,"%s/%sfLatAxesAll.vtk",path,prefix);
    ofstream file(fileName, ios::out);

    if (file){
        file.precision(5);
        file << scientific;
        file << "# vtk DataFile Version 3.0" << std::endl;
        file << "Local Axes of failed lattices" << std::endl;
        file << "ASCII" << endl;
        file << "DATASET UNSTRUCTURED_GRID" << std::endl;
        file << "POINTS "<< LatTab.size() << " float" << std::endl;
        for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
            file << LatTab[LatticeID].centroid[0] << ' '
                    << LatTab[LatticeID].centroid[1] << ' '
                    << LatTab[LatticeID].centroid[2] << std::endl;
        }

        file << "POINT_DATA " << LatTab.size() << std::endl;
        file << "VECTORS PDir_0 float" << endl;
        for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
            file << LatTab[LatticeID].axes[0][0] << ' '
                    << LatTab[LatticeID].axes[0][1] << ' '
                    << LatTab[LatticeID].axes[0][2] << std::endl;
        }
        file << "VECTORS PDir_1 float" << endl;
        for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
            file << LatTab[LatticeID].axes[1][0] << ' '
                    << LatTab[LatticeID].axes[1][1] << ' '
                    << LatTab[LatticeID].axes[1][2] << std::endl;
        }
        file << "VECTORS PDir_2 float" << endl;
        for (unsigned LatticeID=0; LatticeID<LatTab.size(); LatticeID++) {
            file << LatTab[LatticeID].axes[2][0] << ' '
                    << LatTab[LatticeID].axes[2][1] << ' '
                    << LatTab[LatticeID].axes[2][2] << std::endl;
        }
    }
}

void Output::writeFailNetwork	(	const unsigned 				step,
                                    ConnectLat*					p_conLat,
                                    Vertices*					p_verti,
                                    PCG_Eigen*                  p_EigenPCG)
{
    unsigned i,j,k,latNum,listNum,vNum,tvNum,vCnt;
    std::vector<float> vCoord;

    char fileName[255]; //filename

    dualOut << "Generating FailNetwork file..."<<endl;
    tDispLattice = tLatticeNum;

    sprintf(fileName,"%s/%sfNetwork%04d.vtk",path,prefix,step);
    ofstream fSurfFile(fileName, ios::out);

    std::vector<std::vector<unsigned> > latLists = p_conLat->getAllGroup();
    listNum = p_conLat->getListSize();

    latNum = 0;
    for (i=0; i<latLists.size(); i++)
    {
        latNum += latLists[i].size();
    }

    tvNum = 0;
    for (i=0; i<latLists.size(); i++)
    {
        for (j=0; j<latLists[i].size(); j++)
        {
            tvNum += p_verti->getVertiNum(latLists[i][j]);
        }
    }

    if (fSurfFile)
    {
        fSurfFile.precision(5);
        fSurfFile << scientific;
        fSurfFile << "# vtk DataFile Version 3.0" << endl;
        fSurfFile << "My lattices" << endl;
        fSurfFile << "ASCII" << endl;
        fSurfFile << "DATASET POLYDATA" << endl;
        fSurfFile << "POINTS "<< tvNum << " float" << endl;
        for (i=0; i<listNum; i++)
        {
            for (j=0; j<latLists[i].size(); j++)
            {
                vCoord = p_verti-> getLatVertiCoord(latLists[i][j]);
                for (k=0; k<vCoord.size(); k+=3)
                {
                    fSurfFile 	<< vCoord[k] 	<< " "
                                << vCoord[k+1] 	<< " "
                                << vCoord[k+2] 	<< '\n';
                }
            }
        }

        fSurfFile << "POLYGONS " << latNum << " " << latNum + tvNum << endl;
        vCnt = 0;

        for (i=0; i<listNum; i++)
        {
            for (j=0; j<latLists[i].size(); j++)
            {
                vNum = p_verti->getVertiNum(latLists[i][j]);
                fSurfFile << vNum << " ";
                for (k=0; k<vNum; k++)
                {
                    fSurfFile << vCnt++ << " ";
                }
                fSurfFile << '\n';
            }
        }

        fSurfFile << "CELL_DATA " << latNum << endl;
        fSurfFile << "SCALARS Group int" << endl;
        fSurfFile << "LOOKUP_TABLE default" << endl;

        for (i=0; i<listNum; i++)
        {
            for (j=0; j<latLists[i].size(); j++)
            {
                fSurfFile << i << '\n';
            }
        }

        if (ProblemID==0) {
            fSurfFile << "SCALARS Step int" << endl;
            fSurfFile << "LOOKUP_TABLE default" << endl;
            std::vector<std::vector<unsigned> >	latStepLists = p_conLat->getAllStep();
            for (i=0; i<listNum; i++) {
                for (j=0; j<latLists[i].size(); j++) {
                    fSurfFile << latStepLists[i][j]  << '\n';
                }
            }
        } else {
            fSurfFile << "SCALARS Time float" << endl;
            fSurfFile << "LOOKUP_TABLE default" << endl;
            std::vector<std::vector<double> >	latTimeLists = p_conLat->getAllTime();
            for (i=0; i<listNum; i++) {
                for (j=0; j<latLists[i].size(); j++) {
                    fSurfFile << latTimeLists[i][j]  << '\n';
                }
            }
        }
        fSurfFile << "SCALARS Crack_Opening float" << endl;
        fSurfFile << "LOOKUP_TABLE default" << endl;

        GLatForce lforce;
        for (i=0; i<listNum; i++) {
        	for (j=0; j<latLists[i].size(); j++) {
            	fSurfFile << lforce.getCrackOpening(latLists[i][j])-LatTab[latLists[i][j]].initOpening << '\n';
            }
        }

        if ((Use_FrictionModel)||(Use_FullPlasticModel)) {
            fSurfFile << "SCALARS ShearFactor float" << endl;
            fSurfFile << "LOOKUP_TABLE default" << endl;
            for (i=0; i<listNum; i++) {
                for (j=0; j<latLists[i].size(); j++) {
                    fSurfFile << LatTab[latLists[i][j]].factorShear << '\n';
                }
            }

            fSurfFile << "SCALARS ShearStress float" << endl;
            fSurfFile << "LOOKUP_TABLE default" << endl;
            for (i=0; i<listNum; i++) {
                for (j=0; j<latLists[i].size(); j++) {
                    Tensor memForce = p_EigenPCG -> getMemForce(latLists[i][j]);
                    fSurfFile << std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[latLists[i][j]].area << '\n';
                }
            }
        }

        std::vector<std::vector<unsigned> > latTypeLists = p_conLat->getAllType();
        fSurfFile << "SCALARS Type int" << endl;
        fSurfFile << "LOOKUP_TABLE default" << endl;

        for (i=0; i<listNum; i++) {
            for (j=0; j<latLists[i].size(); j++) {
                fSurfFile << latTypeLists[i][j]  << '\n';
            }
        }
    }
    dualOut << "A FailNetwork file is generated..."<<endl;
}

void Output::writeMSplot		(	const unsigned 				step,
                                    ConnectLat*					p_conLat,
                                    Vertices*					p_verti)
{
    std::vector<float> vCoord;
    char fileName[255]; //filename

    dualOut << "Generating Microseismic file..."<<endl;
    tDispLattice = tLatticeNum;

    std::sprintf(fileName,"%s/%sMSplot%04d.vtk",path,prefix,step);
    std::ofstream msPlotFile(fileName, ios::out);
    std::vector<std::vector<unsigned> > latLists = p_conLat->getAllGroup();
    unsigned listNum = p_conLat->getListSize();

    unsigned latNum = 0;
    for (unsigned i=0; i<latLists.size(); i++) {
        latNum += latLists[i].size();
    }

    if (msPlotFile) {
        msPlotFile.precision(5);
        msPlotFile << scientific;
        msPlotFile << "# vtk DataFile Version 3.0" << endl;
        msPlotFile << "My lattices" << endl;
        msPlotFile << "ASCII" << endl;
        msPlotFile << "DATASET UNSTRUCTURED_GRID" << endl;
        msPlotFile << "POINTS "<< latNum << " float" << endl;
        Geometry geo;
        for (unsigned i=0; i<listNum; i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
            	Surface v = p_verti-> getVertiCoordArray (latLists[i][j]);
            	Point centroid = geo.getCentroid(LatTab[latLists[i][j]].area,v);
            	msPlotFile << centroid[0] << ' ' << centroid[1] << ' ' << centroid[2] << '\n';
            }
        }
        GLatForce lforce;
        msPlotFile << "POINT_DATA " << latNum << endl;
        msPlotFile << "SCALARS Energy_release float" << endl;
        msPlotFile << "LOOKUP_TABLE default" << endl;
        for (unsigned i=0; i<listNum; i++) {
        	for (unsigned j=0; j<latLists[i].size(); j++) {
             	msPlotFile << lforce.getEnergyRelease(latLists[i][j]) << '\n';
            }
        }

        msPlotFile << "SCALARS Group int" << endl;
        msPlotFile << "LOOKUP_TABLE default" << endl;

        for (unsigned i=0; i<listNum; i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
                msPlotFile << i << '\n';
            }
        }

        msPlotFile << "SCALARS Step int" << endl;
        msPlotFile << "LOOKUP_TABLE default" << endl;
        std::vector<std::vector<unsigned> >	latStepLists = p_conLat->getAllStep();

        for (unsigned i=0; i<listNum; i++) {
        	for (unsigned j=0; j<latLists[i].size(); j++) {
                msPlotFile << latStepLists[i][j]  << '\n';
            }
        }
        msPlotFile << "SCALARS Time float" << endl;
        msPlotFile << "LOOKUP_TABLE default" << endl;
        std::vector<std::vector<double> >	latTimeLists = p_conLat->getAllTime();
        for (unsigned i=0; i<listNum; i++) {
            for (unsigned j=0; j<latLists[i].size(); j++) {
                msPlotFile << latTimeLists[i][j]  << '\n';
            }
        }

        msPlotFile << "SCALARS Crack_Opening float" << endl;
        msPlotFile << "LOOKUP_TABLE default" << endl;
        for (unsigned i=0; i<listNum; i++) {
        	for (unsigned j=0; j<latLists[i].size(); j++) {
            	msPlotFile << lforce.getCrackOpening(latLists[i][j]) << '\n';
            }
        }
    }
}

void Output::writeFailSurface 	(	const unsigned 				step,
                                    const vector<unsigned>*		p_failLatList,
                                    Vertices*					p_verti)
{
    unsigned i,k,vCnt,LatticeID;
    unsigned VertNum = 0;
    vector<unsigned>	printLatList,nList;
    vector<float>		vCoord;
    char fileName[255]; //filename

    sprintf(fileName,"%s/%sfSurf%04d.vtk",path,prefix,step);
    ofstream fSurfFile(fileName, ios::out);

    dualOut << "Generating FailSurface file..."<<endl;
    tDispLattice = tLatticeNum; //+tvLatticeNum * Disp_VirtLattice;

    if ((p_failLatList->empty())||(step==0))
    {
        dualOut << " No failure lattice, an empty file is generated!" << endl;
        fSurfFile << "# vtk DataFile Version 3.0" << endl;
        fSurfFile << "My lattices" << endl;
        fSurfFile << "ASCII" << endl;
        fSurfFile << "DATASET POLYDATA" << endl;
        fSurfFile << "POINTS "<< 0 << " float" << endl;
        return;
    }

    for (i=0; i<p_failLatList->size(); i++)
    {
        LatticeID = (*p_failLatList)[i];
        if (LatticeID < tDispLattice)
        {
            VertNum += p_verti->getVertiNum(LatticeID);
            printLatList.push_back(LatticeID);
        }
    }

    if (fSurfFile)
    {
        fSurfFile.precision(5);
        fSurfFile << scientific;
        fSurfFile << "# vtk DataFile Version 3.0" << endl;
        fSurfFile << "My lattices" << endl;
        fSurfFile << "ASCII" << endl;
        fSurfFile << "DATASET POLYDATA" << endl;
        fSurfFile << "POINTS "<< VertNum << " float" << endl;
        for (i=0; i<printLatList.size(); i++)
        {
            LatticeID = printLatList[i];
            vCoord = p_verti->getLatVertiCoord(LatticeID);

            for (k=0; k<vCoord.size(); k+=3)
            {
                fSurfFile 	<< vCoord[k] 	<< " "
                            << vCoord[k+1] 	<< " "
                            << vCoord[k+2] 	<< '\n';
            }
            nList.push_back(vCoord.size()/3);
        }

        fSurfFile << "POLYGONS " << printLatList.size() << " " << VertNum + printLatList.size() << endl;
        vCnt = 0;

        for (i=0; i<printLatList.size(); i++)
        {
            fSurfFile << nList[i] << " ";
            for (k=0; k<nList[i]; k++)
            {
                fSurfFile << vCnt++ << " ";
            }
            fSurfFile << '\n';
        }

        fSurfFile << "CELL_DATA " << printLatList.size() << endl;
        fSurfFile << "SCALARS area float" << endl;
        fSurfFile << "LOOKUP_TABLE default" << endl;

        for (i=0; i<printLatList.size(); i++)
        {
            LatticeID = printLatList[i];
            fSurfFile << LatTab[LatticeID].area << "\n";
        }

        fSurfFile << "SCALARS et0 float" << endl;
        fSurfFile << "LOOKUP_TABLE default" << endl;

        for (i=0; i<printLatList.size(); i++)
        {
            LatticeID = printLatList[i];
            fSurfFile << LatTab[LatticeID].et[0] << "\n";
        }
    }
    fSurfFile.close();
}

void Output::writeAllFailSurface   ( Vertices*                   p_verti)
{
    unsigned i,k,vCnt;
    unsigned VertNum = 0;
    vector<unsigned>    printLatList,nList;
    vector<float>       vCoord;
    char fileName[255]; //filename

    sprintf(fileName,"%s/%sfSurfAll.vtk",path,prefix);
    ofstream fSurfFile(fileName, ios::out);

    dualOut << "Generating FailSurface file..."<<endl;
    tDispLattice = tLatticeNum; //+tvLatticeNum * Disp_VirtLattice;

    for (unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++)
    {
        VertNum += p_verti->getVertiNum(LatticeID);
        printLatList.push_back(LatticeID);
    }

    if (fSurfFile)
    {
        fSurfFile.precision(5);
        fSurfFile << scientific;
        fSurfFile << "# vtk DataFile Version 3.0" << endl;
        fSurfFile << "My lattices" << endl;
        fSurfFile << "ASCII" << endl;
        fSurfFile << "DATASET POLYDATA" << endl;
        fSurfFile << "POINTS "<< VertNum << " float" << endl;
        for (i=0; i<printLatList.size(); i++)
        {
            vCoord = p_verti->getLatVertiCoord(printLatList[i]);

            for (k=0; k<vCoord.size(); k+=3)
            {
                fSurfFile   << vCoord[k]    << " "
                            << vCoord[k+1]  << " "
                            << vCoord[k+2]  << '\n';
            }
            nList.push_back(vCoord.size()/3);
        }

        fSurfFile << "POLYGONS " << printLatList.size() << " " << VertNum + printLatList.size() << endl;
        vCnt = 0;

        for (i=0; i<printLatList.size(); i++)
        {
            fSurfFile << nList[i] << " ";
            for (k=0; k<nList[i]; k++)
            {
                fSurfFile << vCnt++ << " ";
            }
            fSurfFile << '\n';
        }

        fSurfFile << "CELL_DATA " << printLatList.size() << endl;
        fSurfFile << "SCALARS area float" << endl;
        fSurfFile << "LOOKUP_TABLE default" << endl;

        for (i=0; i<printLatList.size(); i++)
        {
            fSurfFile << LatTab[printLatList[i]].area << "\n";
        }

        fSurfFile << "SCALARS et0 float" << endl;
        fSurfFile << "LOOKUP_TABLE default" << endl;

        for (i=0; i<printLatList.size(); i++)
        {
            fSurfFile << LatTab[printLatList[i]].et[0] << "\n";
        }
    }
    fSurfFile.close();
}

void Output::writettFailSurface 	(Vertices*					p_verti)
{
    unsigned k,n,vCnt,LatticeID;
    unsigned VertNum = 0;
    vector<float>		vCoord;
    vector<unsigned>	printLatList,nList;
    char fileName[255]; //filename

    dualOut << "Generating test FailSurface file..."<<endl;

    sprintf(fileName,"%s/%stest_ffSurf.vtk",path,prefix);

    ofstream testFile(fileName, ios::out);

    unsigned NodeID = tNodeNum/2;

    for (n=0; n<GNodeTab[NodeID].n; n++)
    {
        LatticeID = GNodeTab[NodeID].nbLatticeID[n];
        VertNum += p_verti->getVertiNum(LatticeID);
    }

    testFile.precision(5);
    testFile << scientific;
    testFile << "# vtk DataFile Version 3.0" << endl;
    testFile << "My lattices" << endl;
    testFile << "ASCII" << endl;
    testFile << "DATASET POLYDATA" << endl;
    testFile << "POINTS "<< VertNum << " float" << endl;

    for (n=0; n<GNodeTab[NodeID].n; n++)
    {
        LatticeID = GNodeTab[NodeID].nbLatticeID[n];
        vCoord = p_verti->getLatVertiCoord(LatticeID);
        for (k=0; k<vCoord.size(); k+=3)
        {
            testFile << vCoord[k] 		<< " "
                     << vCoord[k+1] 	<< " "
                     << vCoord[k+2] 	<< '\n';
        }
        nList.push_back(vCoord.size()/3);
    }

    testFile << "POLYGONS " << GNodeTab[NodeID].n << " " << VertNum + GNodeTab[NodeID].n << endl;
    vCnt = 0;
    for (n=0; n<GNodeTab[NodeID].n; n++)
    {
        testFile << nList[n] << " ";
        for (k=0; k<nList[n]; k++)
        {
            testFile << vCnt++ << " ";
        }
        testFile << '\n';
    }
}

void Output::writeGLattice		(	const unsigned 			step)
{
    char fileName[40]; //filename
    unsigned k,nb1,nb2;
    unsigned LatticeID;

    sprintf(fileName,"%s/%slattice%04d.vtk",path,prefix,step);
    ofstream latticeFile(fileName, ios::out);

    tDispLattice = tLatticeNum + (tbLatticeNum + tvLatticeNum) * isTestCase;

    dualOut << "Generating lattice file..."<<endl;
    if (latticeFile)
    {
        latticeFile.precision(5);
        latticeFile << scientific;
        latticeFile << "# vtk DataFile Version 3.0" << endl;
        latticeFile << "Lattice" << endl;
        latticeFile << "ASCII" << endl;
        latticeFile << "DATASET POLYDATA" << endl;
        latticeFile << "POINTS " << tDispNode << " float" << endl;
        for (unsigned NodeID=0; NodeID<tDispNode; NodeID++)
        {
        	latticeFile << GNodeTab[NodeID].coord[0] + GNodeTab[NodeID].d[0]*Disp_Scale	<< " "
                        << GNodeTab[NodeID].coord[1] + GNodeTab[NodeID].d[1]*Disp_Scale	<< " "
                        << GNodeTab[NodeID].coord[2] + GNodeTab[NodeID].d[2]*Disp_Scale	<< endl;
        }

        latticeFile << "LINES "<< tDispLattice << " " <<3*tDispLattice <<endl;
        for (k=0; k<tDispLattice; k++) {
        	latticeFile << "2 " << LatTab[k].nb[0] << " " << LatTab[k].nb[1] <<endl;
        }

        latticeFile << "CELL_DATA "<< tDispLattice <<endl;

        GLattice_ke(&latticeFile, "GLat_ke");
        GLattice_normF(&latticeFile,"Lattice_normF");
    }

    latticeFile.close();
    dualOut << "A lattice file is generated" << endl;
}

void Output::Node_d 	(	std::ostream*			p_file)
{
    cout << "Writing Node_d..." <<endl;

    *p_file << "SCALARS Disp_0 float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    bool isNormalCase = (ProblemID!=999) ;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << GNodeTab[NodeID].d[0]-GNodeTab[NodeID].d0[0]*isNormalCase <<endl;
    }

    *p_file << "SCALARS Disp_1 float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << GNodeTab[NodeID].d[1] - GNodeTab[NodeID].d0[1]*isNormalCase <<endl;
    }

    *p_file << "SCALARS Disp_2 float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << GNodeTab[NodeID].d[2] - GNodeTab[NodeID].d0[2]*isNormalCase <<endl;
    }
    if (DofPerNode==6) {
        *p_file << "SCALARS Disp_3 float 1" << endl;
        *p_file << "LOOKUP_TABLE default" << endl;

        for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
            *p_file << GNodeTab[NodeID].d[3] - GNodeTab[NodeID].d0[3]*isNormalCase <<endl;
        }
        *p_file << "SCALARS Disp_4 float 1" << endl;
        *p_file << "LOOKUP_TABLE default" << endl;

        for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
            *p_file << GNodeTab[NodeID].d[4] - GNodeTab[NodeID].d0[4]*isNormalCase <<endl;
        }
        *p_file << "SCALARS Disp_5 float 1" << endl;
        *p_file << "LOOKUP_TABLE default" << endl;

        for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
            *p_file << GNodeTab[NodeID].d[5] - GNodeTab[NodeID].d0[5]*isNormalCase <<endl;
        }
    }
}

void Output::Node_isFree     (   std::ostream*           p_file)
{
    cout << "Writing Node_isFree..." <<endl;

    *p_file << "SCALARS isFree int 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
        if (std::binary_search(nodeLists.free.begin(),nodeLists.free.end(),NodeID)) {
            *p_file << 1 <<endl;
        } else if (std::find(nodeLists.restrain.begin(),nodeLists.restrain.end(),NodeID)!=nodeLists.restrain.end()){
            *p_file << 0 <<endl;
        } else {
            tripOut << "[Output::Node_isFree]: a node is neither in nodeLists.free nor in nodeLists.constrain list" << std::endl;
        }
    }
}

void Output::Node_residual 	(	std::ostream*			p_file)
{
    cout << "Writing Node_residual..." <<endl;

    *p_file << "SCALARS Res_0 float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;
    PCG_Eigen pcg;
    std::array<double,DofPerNode>	residual;
    double maxForce;
    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	residual = pcg.getResidual(NodeID);
    	*p_file << residual[0] << std::endl;
    }

    *p_file << "SCALARS Res_1 float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	residual = pcg.getResidual(NodeID);
//    	maxForce = pcg.getMaxLatForce(NodeID);
    	*p_file << residual[1] << std::endl;
    }

    *p_file << "SCALARS Res_2 float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	residual = pcg.getResidual(NodeID);
    	*p_file << residual[2] << std::endl;
    }
    if (DofPerNode==6) {
    	double maxMoment;
        *p_file << "SCALARS Res_3 float 1" << endl;
        *p_file << "LOOKUP_TABLE default" << endl;
        for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
        	residual = pcg.getResidual(NodeID);
        	*p_file << residual[3] << std::endl;
        }

        *p_file << "SCALARS Res_4 float 1" << endl;
        *p_file << "LOOKUP_TABLE default" << endl;
        for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
        	residual = pcg.getResidual(NodeID);
        	*p_file << residual[4] << std::endl;
        }

        *p_file << "SCALARS Res_5 float 1" << endl;
        *p_file << "LOOKUP_TABLE default" << endl;
        for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
        	residual = pcg.getResidual(NodeID);
            *p_file << residual[5] <<endl;
        }
    }
}

void Output::Node_ID 	(	std::ostream*           p_file)
{
    std::cout << "Writing NodeID..." <<endl;

    *p_file << "SCALARS NodeID int 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
        *p_file << NodeID <<endl;
    }
}

void Output::Node_type    (   std::ostream*           p_file)
{
    std::cout << "Writing Node type..." <<endl;

    *p_file << "SCALARS NodeType int 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
        *p_file << GNodeTab[NodeID].type <<endl;
    }
}

void Output::Node_extF 		(	std::ostream*			p_file)
{
    cout << "Writing Node_extF..." << std::endl;
    *p_file << "SCALARS Ext_F0 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++)
    {
    	*p_file << GNodeTab[NodeID].extF[0] << std::endl;
    }

    *p_file << "SCALARS Ext_F1 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++)
    {
    	*p_file << GNodeTab[NodeID].extF[1] << std::endl;
    }

    *p_file << "SCALARS Ext_F2 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++)
    {
    	*p_file << GNodeTab[NodeID].extF[2] << std::endl;
    }
}

void Output::Node_extF_Mag 		(	std::ostream*			p_file)
{
    std::cout << "Writing Node_extF..." << std::endl;
    *p_file << "SCALARS Ext_F float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    Geometry geo;
    double mag;
    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << geo.mag2(GNodeTab[NodeID].extF) << std::endl;
    }
}

void Output::Node_stress 	(	std::ostream*				p_file)
{
	Stress stress;
    stress.calBasic();
    cout << "Writing Node_stress..." <<endl;

    *p_file << "SCALARS Stress_xx float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++)
    {
    	*p_file << stress.getStress(NodeID,0) << endl;
    }

    *p_file << "SCALARS Stress_yy float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++)
    {
    	*p_file << stress.getStress(NodeID,1) << endl;
    }

    *p_file << "SCALARS Stress_zz float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++)
    {
    	*p_file << stress.getStress(NodeID,2) << endl;
    }

    *p_file << "SCALARS Stress_xy float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++)
    {
    	*p_file << stress.getStress(NodeID,3) << endl;
    }

    *p_file << "SCALARS Stress_xz float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++)
    {
    	*p_file << stress.getStress(NodeID,4) << endl;
    }

    *p_file << "SCALARS Stress_yz float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++)
    {
    	*p_file << stress.getStress(NodeID,5) << endl;
    }
}

void Output::Node_PriStress (	std::ostream*			p_file)
{
    Stress stress;
    stress.calPriStress();
    cout << "Writing Node_Pristress..." << std::endl;

    *p_file << "SCALARS PStress_0 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;
    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << stress.getPriStress(NodeID,0) << std::endl;
    }
    *p_file << "VECTORS PDir_0 float" << std::endl;
    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << stress.getPriDir(NodeID,0,0) << " " << stress.getPriDir(NodeID,0,1) << " " << stress.getPriDir(NodeID,0,2) << std::endl;
    }
    *p_file << "SCALARS PStress_1 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << stress.getPriStress(NodeID,1) << std::endl;
    }

    *p_file << "VECTORS PDir_1 float" << std::endl;
    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << stress.getPriDir(NodeID,1,0) << " " << stress.getPriDir(NodeID,1,1) << " " << stress.getPriDir(NodeID,1,2) << std::endl;
    }

    *p_file << "SCALARS PStress_2 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << stress.getPriStress(NodeID,2) << std::endl;
    }

    *p_file << "VECTORS PDir_2 float" << std::endl;
    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << stress.getPriDir(NodeID,2,0) << " " << stress.getPriDir(NodeID,2,1) << " " << stress.getPriDir(NodeID,2,2) << std::endl;
    }

    *p_file << "SCALARS Stress_p float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
        *p_file << stress.getStress_p(NodeID) << std::endl;
    }

    *p_file << "SCALARS Stress_q float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
        *p_file << stress.getStress_q(NodeID) << std::endl;
    }

}

void Output::writeStress		(	const unsigned 			step)
{
    unsigned NodeID;
    unsigned i,start,end,pNodeNum;
    char fileName[255]; //filename

    tDispNode = tNodeNum + (tbNodeNum+tvNodeNum*0)*Use_RanLatStrength;

    sprintf(fileName,"%s/%sstress%04d.vtk",path,prefix,step);

    ofstream nodeFile(fileName, ios::out);

    Stress stress;

    stress.calPriStress();

    dualOut << "Generating Stress file..."<<endl;

    if (nodeFile)
    {
        nodeFile.precision(5);
        nodeFile << scientific;
        nodeFile << "# vtk DataFile Version 3.0" << endl;
        nodeFile << "My lattices" << endl;
        nodeFile << "ASCII" << endl;
        nodeFile << "DATASET UNSTRUCTURED_GRID" << endl;
        nodeFile << "POINTS "<< tDispNode << " float" << endl;

        for (NodeID = 0 ; NodeID < tDispNode; NodeID++)
        {
            nodeFile << GNodeTab[NodeID].coord[0] + GNodeTab[NodeID].d[0]*Disp_Scale	<< " "
                     << GNodeTab[NodeID].coord[1] + GNodeTab[NodeID].d[1]*Disp_Scale	<< " "
                     << GNodeTab[NodeID].coord[2] + GNodeTab[NodeID].d[2]*Disp_Scale	<< endl;
        }

        nodeFile << "POINT_DATA " << tDispNode << endl;

        Node_stress(&nodeFile);

        nodeFile << "SCALARS PStress_0 float 1" << endl;
        nodeFile << "LOOKUP_TABLE default" << endl;

        for (NodeID = 0; NodeID < tDispNode; NodeID++)
        {
            nodeFile << stress.getPriStress(NodeID,0) << endl;
        }

        nodeFile << "VECTORS PDir_0 float" << endl;
        for (NodeID = 0 ; NodeID < tDispNode; NodeID++)
        {
            nodeFile << stress.getPriDir(NodeID,0,0)*stress.getPriStress(NodeID,0) << " "
                     << stress.getPriDir(NodeID,0,1)*stress.getPriStress(NodeID,0) << " "
                     << stress.getPriDir(NodeID,0,2)*stress.getPriStress(NodeID,0) << endl;
        }

        nodeFile << "SCALARS PStress_1 float 1" << endl;
        nodeFile << "LOOKUP_TABLE default" << endl;

        for (NodeID=0; NodeID < tDispNode; NodeID++)
        {
            nodeFile << stress.getPriStress(NodeID,1) << endl;
        }

        nodeFile << "VECTORS PDir_1 float" << endl;
        for (NodeID=0; NodeID < tDispNode; NodeID++)
        {
            nodeFile << stress.getPriDir(NodeID,1,0)*stress.getPriStress(NodeID,1) << " "
                     << stress.getPriDir(NodeID,1,1)*stress.getPriStress(NodeID,1) << " "
                     << stress.getPriDir(NodeID,1,2)*stress.getPriStress(NodeID,1) << endl;
        }

        nodeFile << "SCALARS PStress_2 float 1" << endl;
        nodeFile << "LOOKUP_TABLE default" << endl;

        for (NodeID=0; NodeID < tDispNode; NodeID++)
        {
            nodeFile << stress.getPriStress(NodeID,2) << endl;
        }

        nodeFile << "VECTORS PDir_2 float" << endl;
        for (NodeID=0; NodeID < tDispNode; NodeID++)
        {
            nodeFile << stress.getPriDir(NodeID,2,0)*stress.getPriStress(NodeID,2) << " "
                     << stress.getPriDir(NodeID,2,1)*stress.getPriStress(NodeID,2) << " "
                     << stress.getPriDir(NodeID,2,2)*stress.getPriStress(NodeID,2) << endl;
        }
    }
    nodeFile.close();

    dualOut << "A Stress file is generated" << endl;
}

void Output::GNode_ne 		(	std::ostream*			p_file)
{
    cout << "Writing GNode_ne..." << endl;
    *p_file << "SCALARS ne int 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;
    for (unsigned NodeID=0; NodeID < tDispNode; NodeID++) {
    	*p_file << GNodeTab[NodeID].ne <<endl;
    }
}

void Output::GNode_v 	(	const char* 				fileName)
{
    unsigned NodeID;
    cout << "Writing GNode_vol..." <<endl;
    ofstream nodeFile(fileName, ios::out | ios::app);

    nodeFile << "SCALARS v float 1" << endl;
    nodeFile << "LOOKUP_TABLE default" << endl;

    for (NodeID=0; NodeID < tDispNode; NodeID++)
    {
        nodeFile << GNodeTab[NodeID].v <<endl;
    }
}

void Output::GLattice_normF	(	std::ostream*				p_file,
                                const char* 				varName)
{
    cout << "Writing " << varName << "..." <<endl;

    *p_file << "SCALARS "<< varName << "  float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;
    GLatForce l;
    Geometry geo;
    for(unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++)
    {
    	if (LatTab[LatticeID].isBreak) {
    		*p_file << 0.0 << std::endl;
    	} else {
    		if (true) {
    			*p_file << std::max(l.getLatElongation(LatticeID)/LatTab[LatticeID].et[0],
    					l.getLatElongation(LatticeID)/LatTab[LatticeID].et[1]) << endl;
    		} else {
    			unsigned NodeID = LatTab[LatticeID].nb[0];
    			unsigned nbNodeID = LatTab[LatticeID].nb[1];
    			std::array<double,Dim>  delta,e;
    			for (unsigned d=0; d<Dim; d++) {
    				delta[d] = GNodeTab[nbNodeID].d[d] + GNodeTab[nbNodeID].d0[d]
    			                 - GNodeTab[NodeID].d[d] - GNodeTab[NodeID].d0[d];
    			}
    			for (unsigned d=0; d<Dim; d++) {
    				e[d] = geo.dot(delta,LatTab[LatticeID].axes[d]);
    			}
    			double ratio = LatTab[LatticeID].k[1]/LatTab[LatticeID].k[0];
    			double e_tot = std::sqrt(e[0]*e[0]+ratio*e[1]*e[1]+ratio*e[2]*e[2]);

    			if (e[0]<0.0) {
    				*p_file << e_tot / LatTab[LatticeID].et[1];
    			} else {
    				*p_file << e_tot / LatTab[LatticeID].et[0];
    			}
    		}
    	}
    }
}

void Output::GLattice_area			(	const char* 				fileName,
                                        const char* 				varName)
{
    cout << "Writing " << varName << "..." <<endl;
    ofstream latticeFile(fileName, ios::out | ios::app);

    latticeFile << "SCALARS "<< varName << "  float 1" << endl;
    latticeFile << "LOOKUP_TABLE default" << endl;

    for(unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++)
    {
        latticeFile << LatTab[LatticeID].area << endl;
    }
}

void Output::GLattice_Ft			(	std::ostream*				p_file,
										const char* 				varName)
{
	*p_file << "SCALARS "<< varName << "  float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;
    for(unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        *p_file << LatTab[LatticeID].et[0] << endl;
    }
}

void Output::GLattice_e				(	std::ostream*				p_file,
										const char* 				varName)
{
	*p_file << "SCALARS "<< varName << "  float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;
    GLatForce l;
    for(unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        *p_file << l.getLatElongation(LatticeID) << endl;
    }
}

void Output::Lattice_F 				(	std::ostream*			p_file)
{
    std::cout << "Writing LatticeF ..." <<std::endl;
    *p_file << "SCALARS F0 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;
    Geometry geo;
    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
    	std::array<double,Dim> delta;
    	unsigned NodeID = LatTab[LatticeID].nb[0];
    	unsigned nbNodeID = LatTab[LatticeID].nb[1];
    	for (unsigned d=0; d<Dim; d++) {
    		delta[d] = GNodeTab[nbNodeID].d[d] + GNodeTab[nbNodeID].d0[d]
    	               - GNodeTab[NodeID].d[d] - GNodeTab[NodeID].d0[d];
    	}
    	*p_file << geo.dot(delta,LatTab[LatticeID].axes[0])*LatTab[LatticeID].k[0] << std::endl;
    }
    if (KNum>1) {
    	*p_file << "SCALARS F1 float 1" << std::endl;
    	*p_file << "LOOKUP_TABLE default" << std::endl;
    	for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
    	    std::array<double,Dim> delta;
    	    unsigned NodeID = LatTab[LatticeID].nb[0];
    	    unsigned nbNodeID = LatTab[LatticeID].nb[1];
    	    for (unsigned d=0; d<Dim; d++) {
    	    	delta[d] = GNodeTab[nbNodeID].d[d] + GNodeTab[nbNodeID].d0[d]
    	                   - GNodeTab[NodeID].d[d] - GNodeTab[NodeID].d0[d];
    	    }
    	    *p_file << geo.dot(delta,LatTab[LatticeID].axes[1])*LatTab[LatticeID].k[1] << std::endl;
    	}

    	*p_file << "SCALARS F2 float 1" << std::endl;
    	*p_file << "LOOKUP_TABLE default" << std::endl;
    	for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
    		std::array<double,Dim> delta;
    	    unsigned NodeID = LatTab[LatticeID].nb[0];
    	    unsigned nbNodeID = LatTab[LatticeID].nb[1];
    	    for (unsigned d=0; d<Dim; d++) {
    	    	delta[d] = GNodeTab[nbNodeID].d[d] + GNodeTab[nbNodeID].d0[d]
    	                    - GNodeTab[NodeID].d[d] - GNodeTab[NodeID].d0[d];
    	    }
    	    *p_file << geo.dot(delta,LatTab[LatticeID].axes[2])*LatTab[LatticeID].k[1] << std::endl;
    	}
    }
    if (KNum>2) {
    	*p_file << "SCALARS F3 float 1" << std::endl;
    	*p_file << "LOOKUP_TABLE default" << std::endl;
    	for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
    		std::array<double,Dim> delta;
    		unsigned NodeID = LatTab[LatticeID].nb[0];
    		unsigned nbNodeID = LatTab[LatticeID].nb[1];
    		for (unsigned d=0; d<Dim; d++) {
    			delta[d] = GNodeTab[nbNodeID].d[d+Dim] + GNodeTab[nbNodeID].d0[d+Dim]
    	    	            - GNodeTab[NodeID].d[d+Dim] - GNodeTab[NodeID].d0[d+Dim];
    		}
    		*p_file << geo.dot(delta,LatTab[LatticeID].axes[0])*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]) << std::endl;
    	}
    	*p_file << "SCALARS F4 float 1" << std::endl;
    	*p_file << "LOOKUP_TABLE default" << std::endl;
    	for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
    		std::array<double,Dim> delta;
    		unsigned NodeID = LatTab[LatticeID].nb[0];
    		unsigned nbNodeID = LatTab[LatticeID].nb[1];
    		for (unsigned d=0; d<Dim; d++) {
    			delta[d] = GNodeTab[nbNodeID].d[d+Dim] + GNodeTab[nbNodeID].d0[d+Dim]
    	    	    	   - GNodeTab[NodeID].d[d+Dim] - GNodeTab[NodeID].d0[d+Dim];
    	    }
    		*p_file << geo.dot(delta,LatTab[LatticeID].axes[1])*(LatTab[LatticeID].k[2]) << std::endl;
    	}
    	*p_file << "SCALARS F5 float 1" << std::endl;
    	*p_file << "LOOKUP_TABLE default" << std::endl;
    	for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
    		std::array<double,Dim> delta;
    		unsigned NodeID = LatTab[LatticeID].nb[0];
    		unsigned nbNodeID = LatTab[LatticeID].nb[1];
    		for (unsigned d=0; d<Dim; d++) {
    			delta[d] = GNodeTab[nbNodeID].d[d+Dim] + GNodeTab[nbNodeID].d0[d+Dim]
    	      	    	   - GNodeTab[NodeID].d[d+Dim] - GNodeTab[NodeID].d0[d+Dim];
    		}
    		*p_file << geo.dot(delta,LatTab[LatticeID].axes[2])*(LatTab[LatticeID].k[3]) << std::endl;
    	}
    }
}

void Output::Lattice_fLocal              (      std::ostream*           p_file,
                                                PCG_Eigen*              p_EigenPCG)
{
    std::cout << "Writing LatticeF ..." <<std::endl;
    *p_file << "SCALARS Fm0 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;
    std::vector<Tensor>     memForceList;
    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        memForceList.push_back(p_EigenPCG -> getMemForce(LatticeID));
    }

    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        *p_file << memForceList[LatticeID][0] << std::endl;
    }
    *p_file << "SCALARS Fm1 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        *p_file << memForceList[LatticeID][1] << std::endl;
    }
    *p_file << "SCALARS Fm2 float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        *p_file << memForceList[LatticeID][2] << std::endl;
    }

    if (KNum==4) {
        *p_file << "SCALARS Fm3 float 1" << std::endl;
        *p_file << "LOOKUP_TABLE default" << std::endl;
        for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            *p_file << memForceList[LatticeID][3] << std::endl;
        }

        *p_file << "SCALARS Fm4 float 1" << std::endl;
        *p_file << "LOOKUP_TABLE default" << std::endl;
        for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            *p_file << memForceList[LatticeID][4] << std::endl;
        }

        *p_file << "SCALARS Fm5 float 1" << std::endl;
        *p_file << "LOOKUP_TABLE default" << std::endl;
        for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            *p_file << memForceList[LatticeID][5] << std::endl;
        }
    }
}

void Output::Lattice_fStress             (      std::ostream*           p_file,
                                                PCG_Eigen*              p_EigenPCG)
{
    std::cout << "Writing LatticeStress ..." <<std::endl;
    *p_file << "SCALARS sigma_n float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;
    std::vector<Tensor>     memForceList;
    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        memForceList.push_back(p_EigenPCG -> getMemForce(LatticeID));
    }

    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        *p_file << memForceList[LatticeID][0]/LatTab[LatticeID].area << std::endl;
    }
    *p_file << "SCALARS sigma_s float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        double Fs = std::sqrt(memForceList[LatticeID][1]*memForceList[LatticeID][1]+
                memForceList[LatticeID][2]*memForceList[LatticeID][2]);
        *p_file << Fs/LatTab[LatticeID].area << std::endl;
    }

    *p_file << "SCALARS F_n float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        *p_file << memForceList[LatticeID][0] << std::endl;
    }

    *p_file << "SCALARS F_s float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;

    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        double Fs = std::sqrt(memForceList[LatticeID][1]*memForceList[LatticeID][1]+
                    memForceList[LatticeID][2]*memForceList[LatticeID][2]);
        *p_file << Fs << std::endl;
    }
}

void Output::Lattice_fPerimeter             (   ConnectLat*             p_conLat,
                                                const unsigned          step)
{
    std::cout << "Writing LatticePeri ..." <<std::endl;
    char fileName[255]; //filename
    std::sprintf(fileName,"%s/%scombined%04d.vtk",path,prefix,step);
    std::ofstream nodeFile(fileName, std::ios::out|std::ios::app);
    nodeFile << "SCALARS fLatPeri int 1" << std::endl;
    nodeFile << "LOOKUP_TABLE default" << std::endl;
    std::array<unsigned,3>  xyz;
    xyz[0] = 0;
    xyz[1] = 1;
    xyz[2] = 2;
    if ((GeometryID==4)) {
        xyz[0] = 0;
        xyz[1] = 2;
        xyz[2] = 1;
    }
    std::vector<std::vector<unsigned> > perList = p_conLat->getPerimeter(p_conLat->getMainCrackCenter(),avgD,xyz);
    std::vector<bool>   isPer(tDispLattice,false);
    for (unsigned i=0; i<perList.size(); i++) {
        for (unsigned j=0; j<perList[i].size(); j++)
        isPer[perList[i][j]] = true;
    }
    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        nodeFile << isPer[LatticeID] << std::endl;
    }
}

void Output::Lattice_fLocal_simp              (      std::ostream*           p_file,
                                                     PCG_Eigen*              p_EigenPCG)
{
    std::cout << "Writing LatticeF ..." <<std::endl;
    *p_file << "SCALARS Fn float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;
    std::vector<Tensor>     memForceList;
    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        memForceList.push_back(p_EigenPCG -> getMemForce(LatticeID));
    }

    for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
        *p_file << memForceList[LatticeID][0]/LatTab[LatticeID].area << std::endl;
    }
    *p_file << "SCALARS Fs float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << std::endl;
    if (!Use_FrictionModel) {
        for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            *p_file << std::sqrt(memForceList[LatticeID][1]*memForceList[LatticeID][1]+
                    memForceList[LatticeID][2]*memForceList[LatticeID][2])/LatTab[LatticeID].area << std::endl;
        }
    } else {
        std::vector<double> sigma_s_list;
        for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            double sigma_s = std::sqrt(memForceList[LatticeID][1]*memForceList[LatticeID][1]+
                    memForceList[LatticeID][2]*memForceList[LatticeID][2])/LatTab[LatticeID].area;
            *p_file << sigma_s << std::endl;
            sigma_s_list.push_back(sigma_s);
        }
        *p_file << "SCALARS f float 1" << std::endl;
        *p_file << "LOOKUP_TABLE default" << std::endl;

        for ( unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
            double sigma_n = memForceList[LatticeID][0]/LatTab[LatticeID].area;
            double shearStr = LatTab[LatticeID].shearStr+sigma_n*std::tan(Phi*Pi/180.0);
            *p_file << sigma_s_list[LatticeID]/shearStr << std::endl;
        }
    }
}

void Output::Lattice_offset 				(	std::ostream*			p_file)
{
	if (KNum<4) {
		tripOut << "No lattice offset for KNum < 4" << '\n';
		return;
	}
	std::cout << "Writing Lattice_offset ..." <<std::endl;
	*p_file << "SCALARS offset_0 float 1" << std::endl;
	*p_file << "LOOKUP_TABLE default" << std::endl;
	for (unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
		*p_file << LatTab[LatticeID].offset[0] << std::endl;
	}
	*p_file << "SCALARS offset_1 float 1" << std::endl;
	*p_file << "LOOKUP_TABLE default" << std::endl;
	for (unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
		*p_file << LatTab[LatticeID].offset[1] << std::endl;
	}
}

void Output::Lattice_isBreak 		(	const char* 			fileName,
                                        const char* 			varName)
{
    unsigned LatticeID;

    cout << "Writing " << varName << "..." <<endl;
    ofstream latticeFile(fileName, ios::out | ios::app);

    latticeFile << "SCALARS "<< varName << "  int 1" << std::endl;
    latticeFile << "LOOKUP_TABLE default" << endl;

    for(LatticeID=0; LatticeID<tDispLattice; LatticeID++)
    {
        latticeFile << LatTab[LatticeID].isBreak << std::endl;
    }
}

void Output::Lattice_ShearToAxialRatio 		(	std::ostream*			p_file)
{
    std::cout << "Writing Lattice_shearToAxialRatio..." << std::endl;

    *p_file << "SCALARS ShearToAxialRatio float 1" << std::endl;
    *p_file << "LOOKUP_TABLE default" << endl;
    PCG_Eigen pcg;
    for(unsigned LatticeID=0; LatticeID<tDispLattice; LatticeID++) {
    	Tensor memForce = pcg.getMemForce(LatticeID);
    	double d0 = memForce[0];
    	double d1 = memForce[1];
    	double d2 = memForce[2];
    	if (std::fabs(d0)<Tiny) {
    		*p_file << 0.0 << std::endl;
    	} else {
    		*p_file << std::sqrt(d1*d1+d2*d2)/std::fabs(d0) << std::endl;
    	}
    }
}

void Output::GLattice_ke			(	std::ostream*			p_file,
                                        const char* 			varName)
{
    unsigned LatticeID;

    cout << "Writing " << varName << "..." <<endl;

    *p_file << "SCALARS "<< varName << "  float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    for(LatticeID=0; LatticeID<tDispLattice; LatticeID++)
    {
    	*p_file << LatTab[LatticeID].factor*LatTab[LatticeID].k[0] << endl;
    }
}

void Output::Fluid_width 	(	std::ostream*			p_file,
                            	FluidCal* 				p_fCal)
{
    cout << "Writing Fluid_width..." <<endl;

    *p_file << "SCALARS Width float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    std::vector<double>		pipeWidth = p_fCal->getAllPipeWidth();
    for (unsigned k=0; k < tDispPipe; k++) {
    	*p_file << pipeWidth[k] <<endl;
    }
}

void Output::Fluid_pressure 	(	std::ostream*			p_file,
                            		FluidCal*	 			p_fCal)
{
    cout << "Writing Fluid_pressure..." <<endl;

    *p_file << "SCALARS Pressure float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    std::vector<double>		fPressure = p_fCal->getPressure(0.0);
    for (unsigned i=0; i < tDispfNode; i++)
    {
    	*p_file << fPressure[i] <<endl;
    }
}

void Output::Fluid_flowRate 	(	std::ostream*			p_file,
                            		FluidCal*	 			p_fCal)
{
    cout << "Writing Fluid_flow rate..." <<endl;

    *p_file << "SCALARS Q float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    std::vector<double>		q = p_fCal->getFlowRate();
    for (unsigned i=0; i < tDispfNode; i++)
    {
    	*p_file << q[i] <<endl;
    }
}

void Output::Fluid_leakOff      (   std::ostream*           p_file,
                                    FluidCal*               p_fCal)
{
    cout << "Writing Total Leak off..." <<endl;

    *p_file << "SCALARS LeakOff float 1" << endl;
    *p_file << "LOOKUP_TABLE default" << endl;

    std::vector<double>     leakRate = p_fCal->getLeakOff();
    for (unsigned i=0; i < tDispfNode; i++) {
        *p_file << leakRate[i] << std::endl;
    }
}

void Output::test		(	const char* 			fileName)
{
    unsigned dir,NodeID;
    cout << "Writing Test file..." <<endl;
    ofstream testFile(fileName, ios::out);

    testFile << "---------------Start----------------" << endl << endl;
    for (NodeID = 0; NodeID < tNodeNum+tvNodeNum; NodeID++)
    {
        testFile << " NodeID : " << NodeID << endl;
        testFile << "Nb Node List    : ";
        for (dir = 0; dir < NeighbourNum; dir++)
        {
            testFile << GNodeTab[NodeID].nbNodeID[dir] << " ";
        }
        testFile << endl;
        testFile << "Nb Lattice List : ";
        for (dir = 0; dir < NeighbourNum; dir++)
        {
            testFile << GNodeTab[NodeID].nbLatticeID[dir] << " ";
        }
        testFile << endl;
    }
    testFile << "---------------End----------------" << endl << endl;
}


void Output::printSingleCell	()
{
	char fileName[255];
	sprintf(fileName,"%s/%stest_nodes.vtk",path,prefix);
	ofstream nodeFile(fileName, ios::out);

	unsigned n = nodeLists.unstable.size();
	unsigned nn=0;
	for (unsigned i=0; i<n; i++)
	{
		nn += GNodeTab[nodeLists.unstable[i]].n;
		nn++;
	}

	if (nodeFile)
	{
	    nodeFile.precision(5);
	    nodeFile << scientific;
	    nodeFile << "# vtk DataFile Version 3.0" << endl;
	    nodeFile << "My lattices" << endl;
	    nodeFile << "ASCII" << endl;
	    nodeFile << "DATASET UNSTRUCTURED_GRID" << endl;
	    nodeFile << "POINTS "<< nn << " float" << endl;

	    for (unsigned j=0; j<n; j++)
	    {
			unsigned NodeID = nodeLists.unstable[j];
					nodeFile << GNodeTab[NodeID].coord[0] 	<< " "
							 << GNodeTab[NodeID].coord[1]	<< " "
							 << GNodeTab[NodeID].coord[2]	<< endl;

			unsigned nbNodeID;
			for (unsigned i=0; i < GNodeTab[NodeID].n; i++)
			{
				nbNodeID = GNodeTab[NodeID].nbNodeID[i];
				nodeFile << GNodeTab[nbNodeID].coord[0]  << " "
							<< GNodeTab[nbNodeID].coord[1]  << " "
							<< GNodeTab[nbNodeID].coord[2]  << endl;
			}
	    }
	}

	sprintf(fileName,"%s/%stest_lattices.vtk",path,prefix);
	ofstream latticeFile(fileName, ios::out);

	unsigned ln = nn - n;

	if (latticeFile)
	    {
	        latticeFile.precision(5);
	        latticeFile << scientific;
	        latticeFile << "# vtk DataFile Version 3.0" << endl;
	        latticeFile << "Lattice" << endl;
	        latticeFile << "ASCII" << endl;
	        latticeFile << "DATASET POLYDATA" << endl;
	        latticeFile << "POINTS " << 2*ln << " float" << endl;

	        for (unsigned j=0; j<n; j++)
	        {
	        unsigned nb1,nb2,LatticeID;
	        unsigned NodeID = nodeLists.unstable[j];

	        for (unsigned i=0; i < GNodeTab[NodeID].n; i++)
	        {
	        	LatticeID = GNodeTab[NodeID].nbLatticeID[i];
	            nb1 = LatTab[LatticeID].nb[0];
	            nb2 = LatTab[LatticeID].nb[1];

	            latticeFile << GNodeTab[nb1].coord[0]  << " "
	                        << GNodeTab[nb1].coord[1]  << " "
	                        << GNodeTab[nb1].coord[2]  << endl;

	            latticeFile << GNodeTab[nb2].coord[0]  << " "
	                        << GNodeTab[nb2].coord[1]  << " "
	                        << GNodeTab[nb2].coord[2]  << endl;

	        }
	        }
	        latticeFile << "LINES "<< ln << " " <<3*ln <<endl;
	                for (unsigned i=0; i<ln; i++)
	                {
	                    latticeFile << "2 " << 2*i << " " << 2*i+1 <<endl;
	                }
	    }
}

void Output::writeOuterBox			()
{
	char fileName[255]; //filename

	sprintf(fileName,"%s/%soutBox.vtk",path,prefix);
	ofstream file(fileName, ios::out);
    dualOut << "Generating Outer box file..."<<endl;
	if (file)
	{
		double Dx = Nx*UnitLength;
		double Dy = Ny*UnitLength;
		double Dz = Nz*UnitLength;
		double neg = -(ScaleBoundary-1.0)/2.0;
		double pos = (ScaleBoundary+1.0)/2.0;
		file.precision(5);
	    file << scientific;
	    file << "# vtk DataFile Version 3.0" << endl;
	    file << "Outer Box" << endl;
	    file << "ASCII" << endl;
	    file << "DATASET POLYDATA" << endl;
	    file << "POINTS "<< 8 << " float" << endl;

	    file << neg*Dx << " " << neg*Dy << " " << neg*Dz << endl;
	    file << pos*Dx << " " << neg*Dy << " " << neg*Dz << endl;
	    file << pos*Dx << " " << pos*Dy << " " << neg*Dz << endl;
	    file << neg*Dx << " " << pos*Dy << " " << neg*Dz << endl;
	    file << neg*Dx << " " << neg*Dy << " " << pos*Dz << endl;
	    file << pos*Dx << " " << neg*Dy << " " << pos*Dz << endl;
	    file << pos*Dx << " " << pos*Dy << " " << pos*Dz << endl;
	    file << neg*Dx << " " << pos*Dy << " " << pos*Dz << endl;

	    file << "POLYGONS "<< 6 << " " <<30 <<endl;

	    file << "4 0 1 2 3"<< endl;
	    file << "4 4 5 6 7"<< endl;
	    file << "4 0 1 5 4"<< endl;
	    file << "4 2 3 7 6"<< endl;
	    file << "4 0 4 7 3"<< endl;
	    file << "4 1 2 6 5"<< endl;

	}
	file.close();
	dualOut << "A Outer box file is generated" << endl;
}

void Output::writeInnerBox			()
{
	char fileName[255]; //filename

	sprintf(fileName,"%s/%sinBox.vtk",path,prefix);
	ofstream file(fileName, ios::out);
    dualOut << "Generating Inner box file..."<<endl;
	if (file)
	{
		double Dx = Nx*UnitLength;
		double Dy = Ny*UnitLength;
		double Dz = Nz*UnitLength;
		double neg = 0.0;
		double pos = 1.0;
		file.precision(5);
	    file << scientific;
	    file << "# vtk DataFile Version 3.0" << endl;
	    file << "Outer Box" << endl;
	    file << "ASCII" << endl;
	    file << "DATASET POLYDATA" << endl;
	    file << "POINTS "<< 8 << " float" << endl;

	    file << neg*Dx << " " << neg*Dy << " " << neg*Dz << endl;
	    file << pos*Dx << " " << neg*Dy << " " << neg*Dz << endl;
	    file << pos*Dx << " " << pos*Dy << " " << neg*Dz << endl;
	    file << neg*Dx << " " << pos*Dy << " " << neg*Dz << endl;
	    file << neg*Dx << " " << neg*Dy << " " << pos*Dz << endl;
	    file << pos*Dx << " " << neg*Dy << " " << pos*Dz << endl;
	    file << pos*Dx << " " << pos*Dy << " " << pos*Dz << endl;
	    file << neg*Dx << " " << pos*Dy << " " << pos*Dz << endl;

	    file << "POLYGONS "<< 6 << " " << 30 <<endl;

	    file << "4 0 1 2 3"<< endl;
	    file << "4 4 5 6 7"<< endl;
	    file << "4 0 1 5 4"<< endl;
	    file << "4 2 3 7 6"<< endl;
	    file << "4 0 4 7 3"<< endl;
	    file << "4 1 2 6 5"<< endl;

	}
	file.close();
	dualOut << "A Inner box file is generated" << endl;
}

void Output::writeRestrain          ()
{
    char fileName[255]; //filename
    std::sprintf(fileName,"%s/%sRestrain.txt",path,prefix);
    std::ofstream file(fileName, ios::out);
    std::cout << "Generating restrain file..."<< std::endl;
    file << "Format: NodeID, restrains (3 translations, 3 rotations (optional), 1 means restrained, 0 means free)" << std::endl;
    for (std::list<Restrain>::iterator it=nodeLists.restrain.begin(); it!=nodeLists.restrain.end(); ++it) {
        file << it->NodeID << ' ';
        for (unsigned d=0; d<DofPerNode; d++) {
            file << it->Dir[d] << ' ';
        }
        file << '\n';
    }
    std::cout << "Finish restrain file..."<< std::endl;
}

void Output::writeRawNode   (       const unsigned      totalStep)
{
    char fileName[255]; //filename
    std::sprintf(fileName,"%s/%sRawNodeData-%04d.txt",path,prefix,totalStep);
    std::ofstream file(fileName, ios::out);
    std::cout << "Generating raw node data..."<< std::endl;
    file << "Node coord (x)" << '\t' << "Node coord (y)" << '\t' << "Node coord (z)" << '\t'
            << "Node disp (x)" << '\t' << "Node disp (y)" << '\t' << "Node disp (z) " << '\t' << "Node coord no." << '\t'
            << "Sigma_xx" << '\t' << "Sigma_yy" << '\t' << "Sigma_zz" << '\t' << "Stress_xy" << '\t' << "Stress_xz" << '\t' << "Stress_yz" << std::endl;
    Stress stress;
    stress.calBasic();
    for (unsigned NodeID = 0; NodeID < tNodeNum; NodeID++) {
        file << GNodeTab[NodeID].coord[0] << '\t' << GNodeTab[NodeID].coord[1] << '\t' << GNodeTab[NodeID].coord[2]
             << '\t' << GNodeTab[NodeID].d[0] + GNodeTab[NodeID].d0[0] << '\t' << GNodeTab[NodeID].d[1] + GNodeTab[NodeID].d0[1] << '\t'
             << GNodeTab[NodeID].d[2] + GNodeTab[NodeID].d0[2] << '\t' << GNodeTab[NodeID].n << '\t'
             << stress.getStress(NodeID,0) << '\t' << stress.getStress(NodeID,1) << '\t' << stress.getStress(NodeID,2) << '\t'
             << stress.getStress(NodeID,3) << '\t' << stress.getStress(NodeID,4) << '\t' << stress.getStress(NodeID,5) << std::endl;
    }
}

void Output::writeRawLattice  (    PCG_Eigen*               p_EigenPCG,
                                   const unsigned           totalStep)
{
    char fileName[255]; //filename
    std::sprintf(fileName,"%s/%sRawLatticeData-%04d.txt",path,prefix,totalStep);
    std::ofstream file(fileName, ios::out);
    std::cout << "Generating raw lattice data..."<< std::endl;
    file << "Node1 coord (x)" << '\t' << "Node1 coord (y)" << '\t' << "Node1 coord (z)" << '\t'
         << "Node2 coord (x)" << '\t' << "Node2 coord (y)" << '\t' << "Node2 coord (z)" << '\t'
         << "Node1 disp (x)" << '\t' << "Node1 disp (y)" << '\t' << "Node1 disp (z) " << '\t'
         << "Node2 disp (x)" << '\t' << "Node2 disp (y)" << '\t' << "Node2 disp (z) " << '\t'
         << "Local axis 1 (x)" << '\t' << "Local axis 1 (y)" << '\t' << "Local axis 1 (z)" << '\t'
         << "Local axis 2 (x)" << '\t' << "Local axis 2 (y)" << '\t' << "Local axis 2 (z)" << '\t'
         << "Local axis 3 (x)" << '\t' << "Local axis 3 (y)" << '\t' << "Local axis 3 (z)" << '\t'
         << "Length" << '\t' << "Area" << '\t' << "kn" << '\t' << "ks" << '\t' << "kr1" << '\t' << "kr2" << '\t'
         << "Axial force" << '\t' << "Shear force1" << '\t' << "Shear force2" << '\t' << "isBreak" << '\t' << "fShear" << std::endl;
    for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
        unsigned NodeID1 = LatTab[LatticeID].nb[0];
        unsigned NodeID2 = LatTab[LatticeID].nb[1];
        Tensor     memForce = p_EigenPCG -> getMemForce(LatticeID);
        file << GNodeTab[NodeID1].coord[0] << '\t' << GNodeTab[NodeID1].coord[1] << '\t' << GNodeTab[NodeID1].coord[2] << '\t'
             << GNodeTab[NodeID2].coord[0] << '\t' << GNodeTab[NodeID2].coord[1] << '\t' << GNodeTab[NodeID2].coord[2] << '\t'
             << GNodeTab[NodeID1].d[0] + GNodeTab[NodeID1].d0[0] << '\t' << GNodeTab[NodeID1].d[1] + GNodeTab[NodeID1].d0[1] << '\t'
             << GNodeTab[NodeID1].d[2] + GNodeTab[NodeID1].d0[2] << '\t' << GNodeTab[NodeID2].d[0] + GNodeTab[NodeID2].d0[0] << '\t'
             << GNodeTab[NodeID2].d[1] + GNodeTab[NodeID2].d0[1] << '\t' << GNodeTab[NodeID2].d[2] + GNodeTab[NodeID2].d0[2] << '\t'
             << LatTab[LatticeID].axes[0][0] << '\t' << LatTab[LatticeID].axes[0][1] << '\t' << LatTab[LatticeID].axes[0][2] << '\t'
             << LatTab[LatticeID].axes[1][0] << '\t' << LatTab[LatticeID].axes[1][1] << '\t' << LatTab[LatticeID].axes[1][2] << '\t'
             << LatTab[LatticeID].axes[2][0] << '\t' << LatTab[LatticeID].axes[2][1] << '\t' << LatTab[LatticeID].axes[2][2] << '\t'
             << LatTab[LatticeID].length << '\t' << LatTab[LatticeID].area << '\t'
             << LatTab[LatticeID].k[0]*LatTab[LatticeID].factor <<  '\t' << LatTab[LatticeID].k[1]*LatTab[LatticeID].factorShear << '\t' << LatTab[LatticeID].k[2]*LatTab[LatticeID].factor << '\t' << LatTab[LatticeID].k[3]*LatTab[LatticeID].factor << '\t'
             << memForce[0] << '\t' << memForce[1] << '\t' << memForce[2] << '\t' << LatTab[LatticeID].isBreak << '\t' << LatTab[LatticeID].fShear << std::endl;
    }
}
