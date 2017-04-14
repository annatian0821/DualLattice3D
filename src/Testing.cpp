/*
 * Testing.cpp
 *
 *  Created on: Dec 31, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "Testing.hpp"

using namespace std;

void Testing::check()
{
    testMinCoordNum(3);
}

void Testing::check							(	const double			accSE,
                                                const unsigned			totalStep,
                                                const bool				isRelax)
{
    checkEnergy (accSE,totalStep,isRelax);
}

double	Testing::checkEnergy				(	const double			accSE,
                                                const unsigned			totalStep,
                                                const bool				isRelax)
{
    char	fileName[255];
    double 	tExtWorkDone,tStrainEnergy,residual;

    sprintf(fileName,"%s/%shist-Energy.txt",path,prefix);

    tExtWorkDone = getTotalExtWorkDoneSimple (Face6,2,isRelax);
    tStrainEnergy = getTotalStrainEnergy();
    residual = (tExtWorkDone - tStrainEnergy - accSE)/tExtWorkDone;

    dualOut << "Total external work done simple = " << tExtWorkDone << endl;

    dualOut << "Total Strain Energy = "	<< tStrainEnergy << endl;

    dualOut << "Accumlated surface energy to crack = " << accSE << endl;
    dualOut << "Residual = (Ext work done - total strain energy - acc. surface energy)/Ext work done = "
            << residual << endl;
    dualOut << "Residual = (Ext work done - total strain energy)/Ext work done = "
            << (tExtWorkDone - tStrainEnergy )/tExtWorkDone << endl;

    if (totalStep==0)
    {
        ofstream eFile(fileName, ios::out);
        eFile	<< totalStep << '\t' << tExtWorkDone << '\t' << tStrainEnergy << '\t' << accSE << '\t' << residual << endl;
    }
    else
    {
        ofstream eFile(fileName, ios::out|ios::app);
        eFile	<< totalStep << '\t' << tExtWorkDone << '\t' << tStrainEnergy << '\t' << accSE << '\t' << residual << endl;
    }

    return residual;
}

bool Testing::testMinCoordNum		(	const unsigned					MinCoordNum)
{
    unsigned min = 100;

    for (std::vector<unsigned>::iterator it = nodeLists.free.begin(); it != nodeLists.free.end(); ++it) {
        unsigned NodeID = *it;
        if (GNodeTab[NodeID].ne < min) {
            min = GNodeTab[NodeID].ne;
        }
    }

    for (std::list<Restrain>::iterator rit = nodeLists.restrain.begin() ; rit != nodeLists.restrain.end(); ++rit) {
        unsigned NodeID = rit -> NodeID;
        if (GNodeTab[NodeID].ne < min) {
            min = GNodeTab[NodeID].ne;
        }
    }

    dualOut << "Min coordination number = " << min;

    if (min < MinCoordNum) {
        dualOut << " Model unstable!!!!!" << endl;
        return false;
    } else {
        dualOut << " Model stable" << endl;
        return true;
    }
}

double Testing::getTotalStrainEnergy		()
{
    double	totalSE = 0.0;
    double	inSE = 0.0;
    GLatForce   l;
    double e;
    for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++)
    {
        e = l.getLatElongation(LatticeID);
        totalSE += 0.5*LatTab[LatticeID].k[0]*e*e;
    }

    inSE = totalSE ;

    for (unsigned LatticeID = tLatticeNum; LatticeID < tLatticeNum + tbLatticeNum; LatticeID++) {
        e = l.getLatElongation(LatticeID);
        totalSE += 0.5*LatTab[LatticeID].k[0]*e*e;
    }


    for (unsigned LatticeID = tLatticeNum + tbLatticeNum;
            LatticeID < tLatticeNum + tbLatticeNum + tvLatticeNum; LatticeID++) {
        e = l.getLatElongation(LatticeID);
        totalSE += 0.5*LatTab[LatticeID].k[0]*e*e;
    }

    cout << " Percentage of strain energy stored in boundary lattice and virtual lattice = " << (totalSE - inSE)/totalSE << endl;

    return totalSE;
}

double Testing::getTotalExtWorkDoneSimple	(	const unsigned 			face,
                                                const unsigned 			dir,
                                                const bool				isRelax)
{
    unsigned i,NodeID;
    double sum = 0.0;
    double totalWD  = 0.0;
    double totalWD2 = 0.0;

    if (isRelax)
    {
        for (i=0; i<nodeLists.boundary[face].size(); i++)
        {
            NodeID = nodeLists.boundary[face][i];
            totalWD += 0.5*GNodeTab[NodeID].extF[dir] * (GNodeTab[NodeID].d[dir]);
            sum += (GNodeTab[NodeID].d[dir] - old_bDisp[i]);
            old_bDisp[i] = GNodeTab[NodeID].d[dir];
        }
        cout << sum/nodeLists.boundary[face].size();
        accWD = totalWD;
    }
    else
    {
        for (i=0; i<nodeLists.boundary[face].size(); i++)
        {
            NodeID = nodeLists.boundary[face][i];
            totalWD += (GNodeTab[NodeID].extF[dir] - 0.5*GNodeTab[NodeID].extdF[dir]) * (GNodeTab[NodeID].d[dir] - old_bDisp[i]);
            old_bDisp[i] = GNodeTab[NodeID].d[dir];
        }
        accWD += totalWD;
    }


    for (i=0; i<nodeLists.boundary[face].size(); i++)
    {
        totalWD2 += GNodeTab[nodeLists.boundary[face][i]].extF[dir]*GNodeTab[nodeLists.boundary[face][i]].d[dir];
    }

    cout << "accWD = " << accWD << " totalWD2 = " << totalWD2/2.0 << endl;

    return accWD;
}

void Testing::checkForce					(	vector<vector<double> >*	p_bForce,
                                                const unsigned 				face,
                                                const unsigned				step)
{
    unsigned i,d,NodeID;
    unsigned f = 0;

    double	diff;

    char	fileName[255];
    const char* varName1 = "chk";

    dualOut << "Checking forces at boundaries " << endl;

    sprintf(fileName,"%s/%sstat-%s-%04d.txt",path,prefix,varName1,step);

    ofstream chkFile(fileName); //, ios::out | ios::app);

    if (face > 2)
    {
        f = face-2;
    }
    else
    {
        f = face;
    }

    chkFile << "Face = " << face << '\n';

    for (i=0; i<nodeLists.boundary[face].size(); i++)
    {
        NodeID = nodeLists.boundary[face][i];
        chkFile << NodeID << '\t';
        for (d=0; d < Dim ; d++)
        {
            diff = (*p_bForce)[f*Dim+d][i] / GNodeTab[NodeID].extF[d];
            chkFile <<  diff << '\t';
        }
        chkFile << '\n';
    }
}

double Testing::getTotalExtWorkDoneComplete	(	const unsigned			step)
{
    unsigned i,b,d,dir;
    unsigned NodeID;
    double sumWD = 0.0;
    double inWD = 0.0;
    vector<double>	totalWD(6*Dim,0.0);
    unsigned face[6] = {Face1, Face2, Face3, Face4, Face5, Face6};
    vector<vector<double> >	bForce;

    bForce.resize(6*Dim);

    getBoundaryReaction(&bForce,step);

    checkForce(&bForce,Face6,step);

    inWD = insideReaction(step);

    for (b=0; b<6; b++)
    {
        for (i=0; i<nodeLists.boundary[face[b]].size(); i++)
        {
            NodeID = nodeLists.boundary[face[b]][i];
            for (dir=0; dir<Dim; dir++)
            {
                totalWD[b*Dim+dir] += 	0.5* (bForce[b*Dim+dir][i] + old_AllbForce[b*Dim+dir][i])*
                                        (GNodeTab[NodeID].d[dir] - old_AllbDisp[b*Dim+dir][i]);

                old_AllbDisp[b*Dim+dir][i] = GNodeTab[NodeID].d[dir];
                old_AllbForce[b*Dim+dir][i] = bForce[b*Dim+dir][i];
            }
        }
    }

    for (b=0; b<6; b++)
    {
        cout << "Incremental Work done in face " << b << " are = ";
        for (d=0; d<Dim; d++)
        {
            cout << totalWD[b*Dim+d] << " ";
            sumWD += totalWD[b*Dim+d];
        }
        cout << endl;
    }

    accAllWD +=sumWD+inWD;

    return accAllWD;
}

void Testing::getBoundaryReaction	(	std::vector<std::vector<double> >*		p_bForce,
                                        const unsigned					        step)
{
    unsigned i,b,d,k,NodeID,nbNodeID,LatticeID;
    unsigned face[6] = {Face1, Face2, Face3, Face4, Face5, Face6};
    double	 force[Dim] = {0.0};

    char	fileName[255];
    const char* varName1 = "bb";

    dualOut << "Writing reactions at boundaries " << endl;

    sprintf(fileName,"%s/%sstat-%s-%04d.txt",path,prefix,varName1,step);

    ofstream bFile(fileName); //, ios::out | ios::app);

    dualOut << "Calculating reaction at boundaries " << endl;

    for (b=0; b<6; b++)
    {
        for (d=0; d<Dim; d++)
        {
            (*p_bForce)[b*Dim+d].clear();
        }
    }
    GLatForce l;

    for (b=0; b<6; b++)
    {
        bFile << "Face = " << b << '\n';
        for (i=0; i<nodeLists.boundary[face[b]].size(); i++)
        {
            NodeID = nodeLists.boundary[face[b]][i];
            for (k=0; k<GNodeTab[NodeID].n ; k++)
            {
                LatticeID = GNodeTab[NodeID].nbLatticeID[k];
                nbNodeID = GNodeTab[NodeID].nbNodeID[k];
                double forceMag = LatTab[LatticeID].k[0]*l.getLatElongation(LatticeID);
                for (d=0; d<Dim; d++)
                {
                    force[d] += forceMag *(GNodeTab[NodeID].coord[d] - GNodeTab[nbNodeID].coord[d])/LatTab[LatticeID].length;
                }
            }

            bFile << NodeID << '\t' ;
            for (d=0; d<Dim; d++)
            {
                (*p_bForce)[b*Dim+d].push_back(force[d]);
                bFile << force[d] << '\t';
                force [d] = 0.0;
            }
            bFile << '\n';
        }
        bFile << "----------------------------------------------------------" << endl;
    }
}

void Testing::writeBoundaryReaction		(	const unsigned 			step)
{
    unsigned i,b,d,k,NodeID,nbNodeID,LatticeID;
    unsigned face[6] = {Face1, Face2, Face3, Face4, Face5, Face6};
    double	 force[Dim] = {0.0};

    char	fileName[255];
    const char* path = OutputSubFolder;
    const char* prefix = OutputFilePrefix;
    const char* varName1 = "b";

    dualOut << "Writing reactions at boundaries " << endl;

    sprintf(fileName,"%s/%sstat-%s-%04d.txt",path,prefix,varName1,step);

    ofstream bFile(fileName); //, ios::out | ios::app);
    GLatForce l;
    for (b=0; b<6; b++)
    {
        bFile << "Face = " << b << '\n';
        for (i=0; i<nodeLists.boundary[face[b]].size(); i++)
        {
            NodeID = nodeLists.boundary[face[b]][i];
            for (k=0; k<GNodeTab[NodeID].n ; k++)
            {
                LatticeID = GNodeTab[NodeID].nbLatticeID[k];
                nbNodeID = GNodeTab[NodeID].nbNodeID[k];
                double forceMag = LatTab[LatticeID].k[0]*l.getLatElongation(LatticeID);
                for (d=0; d<Dim; d++)
                {
                    force[d] += forceMag*(GNodeTab[NodeID].coord[d] - GNodeTab[nbNodeID].coord[d])/LatTab[LatticeID].length;
                }
            }
            bFile << NodeID << '\t' ;
            for (d=0; d<3; d++)
            {
                bFile << force[d] << '\t';
                force[d] =0.0;
            }
            bFile << '\n';
        }
        bFile << "----------------------------------------------------------" << endl;
    }
}

double Testing::insideReaction			(	const unsigned 		step)
{
    unsigned nbNodeID,LatticeID;
    double	 force[Dim] = {0.0};
    double	 WD = 0.0;

    char		fileName[255];
    const char* varName1 = "in";

    dualOut << "Writing reactions inside " << endl;

    sprintf(fileName,"%s/%sstat-%s-%04d.txt",path,prefix,varName1,step);

    ofstream inFile(fileName); //, ios::out | ios::app);
    GLatForce l;
    for (unsigned NodeID=0; NodeID<tNodeNum; NodeID++) {
        for (unsigned k=0; k<GNodeTab[NodeID].n ; k++)
        {
            LatticeID = GNodeTab[NodeID].nbLatticeID[k];
            nbNodeID = GNodeTab[NodeID].nbNodeID[k];
            double forceMag = LatTab[LatticeID].k[0]*l.getLatElongation(LatticeID);
            for (unsigned d=0; d<Dim; d++)
            {
                force[d] += forceMag*(GNodeTab[NodeID].coord[d] - GNodeTab[nbNodeID].coord[d])/LatTab[LatticeID].length;
            }
        }
        inFile << NodeID << '\t' ;
        for (unsigned d=0; d<3; d++)
        {
            inFile << force[d] << '\t';
            WD += 0.5*force[d]*GNodeTab[NodeID].d[d];
            force[d] =0.0;
        }
        inFile << '\n';
    }

    dualOut << "Total work done by internal nodes (should be zero) = " << WD << endl;
    return WD;
}


double	Testing::calAllResidual				()
{
    GLatForce lForce;

    double residual = 0.0;
    double resInside;
    unsigned face[6] = {Face1, Face2, Face3, Face4, Face5, Face6};
    std::vector<unsigned>			calDir;

    lForce.calAllLatForce(true,true,true);

    for (unsigned d=0; d<Dim; d++) {
        calDir.push_back(d);
    }

    for (std::vector<unsigned>::iterator   it = nodeLists.free.begin(); it != nodeLists.free.end(); ++it) {
        unsigned NodeID = *it;
        residual += calSingleResidual (NodeID,calDir);
    }

    std::cout << residual << " ";
    resInside = residual ;

    for (unsigned b=0; b<6; b++) {
        calDir.clear();
        for (unsigned d=0; d<Dim; d++) {
            if (!BoundaryMat[face[b]][d]) {
                calDir.push_back(d);
            }
        }

        for (unsigned i = 0; i < nodeLists.boundary[face[b]].size() ; i++) {
            unsigned NodeID = nodeLists.boundary[face[b]][i];
            residual += calSingleResidual (NodeID,calDir);
        }
    }

    std::cout << residual - resInside << " ";

    return residual;
}

double	Testing::calSingleResidual			(	const unsigned			    NodeID,
                                                std::vector<unsigned>		calDir)
{
    unsigned nbNodeID,LatticeID;
    unsigned d;

    unsigned dn = calDir.size();
    double	 force[Dim] = {0.0};
    double 	 residual = 0.0;

    GLatForce l;
    for (unsigned k=0; k<GNodeTab[NodeID].n ; k++)
    {
        LatticeID = GNodeTab[NodeID].nbLatticeID[k];
        nbNodeID = GNodeTab[NodeID].nbNodeID[k];
        double forceMag = LatTab[LatticeID].k[0]*l.getLatElongation(LatticeID);
        for (d=0; d<dn; d++)
        {
            force[calDir[d]] += forceMag*(GNodeTab[NodeID].coord[calDir[d]] - GNodeTab[nbNodeID].coord[calDir[d]])/LatTab[LatticeID].length;
        }
    }

    for (d = 0; d < dn ; d++)
    {
        residual += fabs(GNodeTab[NodeID].extF[calDir[d]] - force[calDir[d]]);
    }

    return residual;
}
