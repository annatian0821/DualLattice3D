/*
 * Vertices.cpp
 *
 *  Created on: Jan 28, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "Vertices.hpp"
using namespace std;

Vertices::Vertices() {
    // TODO Auto-generated constructor stub

}

Vertices::~Vertices() {
    // TODO Auto-generated destructor stub
    delete p_vFile;
}

void Vertices::initialize	(	const unsigned			LatNum,
                                const double			tol)
{
    char fileName[255];
    maxSize = LatNum;
    latList.resize(LatNum);
    tolerance = tol;

    sprintf(fileName,"%s/%svertices.vtk",OutputFolder,OutputFilePrefix);
    p_vFile = new ofstream (fileName,  ios::out | ios::app);
};

void Vertices::inputNode			(	const unsigned					NodeID,
                                        const vector<double>*			p_vCoord,
                                        const vector<vector<int> >*		p_vList,
                                        const vector<int>*				p_vListIndex)
{
    unsigned i;
    unsigned nbNodeID,LatticeID;
    vector<unsigned> latList,nbLatList,vIndex,vIDList,chkVIDList;

    latList = getInputList (&(GNodeTab[NodeID].nbLatticeID),p_vListIndex,&vIndex);

    for (i=0; i<latList.size(); i++)
    {
        LatticeID = latList[i];
        nbNodeID = LatTab[LatticeID].nb[1];
        if (nbNodeID==NodeID)
        {
            nbNodeID = LatTab[LatticeID].nb[0];
        }
        nbLatList = getCheckLatList (LatticeID,&GNodeTab[NodeID].nbLatticeID,&GNodeTab[nbNodeID].nbLatticeID);
        chkVIDList = getCheckVertices (&nbLatList);

        vIDList = getLatList (&chkVIDList,p_vCoord,&((*p_vList)[vIndex[i]]));
        setLatList (LatticeID,&vIDList);
    }
}

vector<unsigned> Vertices::getInputList				(	const std::vector<unsigned>*	p_LatList,
        												const std::vector<int>*			p_vListIndex,
        												std::vector<unsigned>*			p_vIndex)
{
    unsigned i;
    vector<unsigned> inputList;

    p_vIndex->clear();

    for (i=0; i<p_vListIndex->size(); i++)
    {
        if ((*p_LatList)[i]<maxSize)
        {
            if ((*p_vListIndex)[i]!=-1)
            {
                inputList.push_back((*p_LatList)[i]);
                p_vIndex->push_back((*p_vListIndex)[i]);
            }
        }
    }
    return inputList;
}

vector<unsigned>	Vertices::getCheckLatList		(	const unsigned					LatticeID,
        const std::vector<unsigned>*	p_selfLatList,
        const std::vector<unsigned>*	p_nbLatList)
{
    unsigned i;
    vector<unsigned> checkLatList;
    vector<unsigned> mergeList = *p_selfLatList;
    mergeList.insert(mergeList.end(),p_nbLatList->begin(),p_nbLatList->end());

    for (i=0; i<mergeList.size(); i++)
    {
        if (isLatticeInputted(mergeList[i]))
        {
            checkLatList.push_back(mergeList[i]);
        }
    }
    uniqueVector (&checkLatList);
    return checkLatList;
}

vector<unsigned>	Vertices::getCheckVertices		(	const vector<unsigned>*			p_nbLatList)
{
    unsigned i,j;
    vector<unsigned>	chkVList;

    for (i=0; i<p_nbLatList->size(); i++)
    {
        for (j=0; j<latList[(*p_nbLatList)[i]].size(); j++)
        {
            chkVList.push_back(latList[(*p_nbLatList)[i]][j]);
        }
    }
    uniqueVector (&chkVList);
    return chkVList;
}

vector<unsigned>	Vertices::getLatList		(	const vector<unsigned>*		p_chkVList,
        											const vector<double>*		p_vCoord,
        											const vector<int>*			p_vList)
{
    unsigned i,j;
    double	len;
    bool 			isRepeat = false;
    array<double,3>	x,a;
    vector<unsigned>	outVID;

    Geometry geo;

    for (i=0; i<p_vList->size(); i++)
    {
        x[0] = (*p_vCoord)[3*(*p_vList)[i]];
        x[1] = (*p_vCoord)[3*(*p_vList)[i]+1];
        x[2] = (*p_vCoord)[3*(*p_vList)[i]+2];

        isRepeat = false;

        for (j=0; j<p_chkVList->size(); j++)
        {
            a[0] = coord[3*(*p_chkVList)[j]];
            a[1] = coord[3*(*p_chkVList)[j]+1];
            a[2] = coord[3*(*p_chkVList)[j]+2];

            len = geo.dist(x,a);
            if (len < tolerance)
            {
                outVID.push_back((*p_chkVList)[j]);
                isRepeat = true;
                break;
            }
        }
        if (!isRepeat)
        {
            outVID.push_back(putVertice(x[0],x[1],x[2]));
        }
    }
    return outVID;
}

Line			Vertices::getCommonEdge			(	const unsigned			LatID1,
        											const unsigned			LatID2)
{
    unsigned i,j;
    unsigned cntEq = 0;
    Point	x,a;
    Geometry geo;
    std::vector<Point> commonVerti;
    Point dummy;
    dummy[0]=0.0; dummy[1]=0.0; dummy[2]=0.0;

    Line edge;
    if (LatID1==LatID2) {
    	dualOut << "Two lattice are the same, return a dummy line" << std::endl;
    	edge[0] = edge[1] = dummy;
    	return edge;
    }

    for (i = 0 ; i<latList[LatID1].size(); i++)
    {
        x[0] = coord[3*latList[LatID1][i]];
        x[1] = coord[3*latList[LatID1][i]+1];
        x[2] = coord[3*latList[LatID1][i]+2];
        for (j = 0; j<latList[LatID2].size(); j++)
        {
            a[0] = coord[3*latList[LatID2][j]];
            a[1] = coord[3*latList[LatID2][j]+1];
            a[2] = coord[3*latList[LatID2][j]+2];

            if (geo.dist(x,a) < tolerance)
            {
            	commonVerti.push_back(x);
                break;
            }
        }
    }

    if (commonVerti.size()==2) {
    	edge[0] = commonVerti[0];
    	edge[1] = commonVerti[1];
        return edge;
    } else if (commonVerti.size()>2) {
    	dualOut << "Vertices::getCommonEdge - the common verti found is more than 2, only the first 2 vertices are used" << std::endl;
    	edge[0] = commonVerti[0];
    	edge[1] = commonVerti[1];
    	return edge;
    } else {
    	dualOut << "Vertices::getCommonEdge - the common verti found is less than 2, return zero line and please check!!!" << std::endl;
    	edge[0] = edge[1] = dummy;
        return edge;
    }
}


void	Vertices::setLatList	(	const unsigned					LatticeID,
                                    const vector<unsigned>*			p_vID_List)
{
    latList[LatticeID] = *p_vID_List;
}

void	Vertices::printAll	()
{
    unsigned i,LatticeID;
    *p_vFile << "*************************************************"<< endl;
    *p_vFile << "*************************************************"<< endl;
    *p_vFile << "*************************************************"<< endl;
    for (LatticeID =0; LatticeID < maxSize; LatticeID++)
    {
        *p_vFile << LatticeID << " - ";
        for (i=0; i<latList[LatticeID].size(); i++)
        {
            *p_vFile << latList[LatticeID][i] << " ";
        }
        *p_vFile << endl;
    }
}
