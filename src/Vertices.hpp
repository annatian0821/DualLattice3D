/*
 * Vertices.hpp
 *
 *  Created on: Jan 28, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef VERTICES_HPP_
#define VERTICES_HPP_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "MiscTools.hpp"

typedef std::vector<std::vector<unsigned> >	vList;
//typedef std::vector<std::array<float,Dim> > vaList;

class Vertices
{
    friend class ConnectLat;
public:
    Vertices();
    virtual ~Vertices();
    void initialize			(	const unsigned			LatNum,
                                const double			tol);

    void inputNode						(	const unsigned							NodeID,
                                            const std::vector<double>*				p_vCoord,
                                            const std::vector<std::vector<int> >*	p_vList,
                                            const std::vector<int>*					p_vListIndex);

    inline std::vector<float>		getLatVertiCoord	(	const unsigned			LatticeID);

    inline Surface					getVertiCoordArray	(	const unsigned			LatticeID);

    inline unsigned 	getVertiNum			(	const unsigned			LatticeID);

    inline bool			isConnected			(	const unsigned			LatID1,
            const unsigned			LatID2);

    inline bool			isConnectedExact	(	const unsigned			LatID1,
            									const unsigned			LatID2);

    Line				getCommonEdge		(	const unsigned			LatID1,
            									const unsigned			LatID2);

    void	printAll	();
    std::ofstream* 						p_vFile;

private:
    std::vector<double>					coord;
    std::vector<std::vector<unsigned> >	latList;
    double								tolerance;
    unsigned							maxSize;



    inline	unsigned putVertice	(	const double	x,
                                    const double	y,
                                    const double	z);

    void	setLatList			(	const unsigned					LatticeID,
                                    const std::vector<unsigned>*	p_vID_List);

    inline bool	isLatticeNotInputted			(	const unsigned		LatticeID);

    inline bool isLatticeInputted				(	const unsigned		LatticeID);

    std::vector<unsigned> 	getInputList		(	const std::vector<unsigned>*	p_LatList,
            const std::vector<int>*			p_vListIndex,
            std::vector<unsigned>*			p_vIndex);

    std::vector<unsigned>	getCheckLatList		(	const unsigned					LatticeID,
            const std::vector<unsigned>*	p_selfLatList,
            const std::vector<unsigned>*	p_nbLatList);

    std::vector<unsigned>	getCheckVertices	(	const std::vector<unsigned>*	p_nbLatList);

    std::vector<unsigned>	getLatList			(	const std::vector<unsigned>*	p_chkVList,
            const std::vector<double>*		p_vCoord,
            const std::vector<int>*			p_vList);

};

unsigned 				Vertices::putVertice	(	const double	x,
        const double	y,
        const double	z)
{
    coord.push_back(x);
    coord.push_back(y);
    coord.push_back(z);
    return coord.size()/3 - 1;
}

std::vector<float>		Vertices::getLatVertiCoord	(	const unsigned			LatticeID)
{
    unsigned i;
    std::vector<float>		vCoord;
    for (i=0; i<latList[LatticeID].size(); i++)
    {
        vCoord.push_back(coord[3*latList[LatticeID][i]]);
        vCoord.push_back(coord[3*latList[LatticeID][i]+1]);
        vCoord.push_back(coord[3*latList[LatticeID][i]+2]);
    }
    return vCoord;
}

Surface			Vertices::getVertiCoordArray	(	const unsigned			LatticeID)
{
	Surface out;
	Point coordArr;
	for (unsigned i=0; i<latList[LatticeID].size(); i++) {
		for (unsigned j=0; j<Dim; j++) {
			coordArr[j] = (coord[3*latList[LatticeID][i]+j]);
		}
	    out.push_back(coordArr);
	}
	return out;
};

unsigned 		Vertices::getVertiNum			(	const unsigned			LatticeID)
{
    return latList[LatticeID].size();
}

bool			Vertices::isConnected			(	const unsigned			LatID1,
        											const unsigned			LatID2)
{
    unsigned sum1,sum2;
    std::vector<unsigned>	vList12;

    sum1 = latList[LatID1].size() + latList[LatID2].size();

    vList12.reserve(sum1);

    vList12.insert(vList12.end(),latList[LatID1].begin(), latList[LatID1].end());
    vList12.insert(vList12.end(),latList[LatID2].begin(), latList[LatID2].end());

    uniqueVector(&vList12);
    sum2 = vList12.size();

    if (sum1-sum2>=2)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool			Vertices::isConnectedExact		(	const unsigned			LatID1,
        											const unsigned			LatID2)
{
    unsigned i,j;
    unsigned cntEq = 0;
    Point	x,a;
    Geometry geo;

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
                cntEq++;
                break;
            }
        }
    }

    if (cntEq>=2)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool			Vertices::isLatticeNotInputted	(	const unsigned			LatticeID)
{
    if (LatticeID<maxSize)
    {
        return latList[LatticeID].empty();
    }
    else
    {
        return false;
    }
}

bool Vertices::isLatticeInputted				(	const unsigned			LatticeID)
{
    if (LatticeID<maxSize)
    {
        return !latList[LatticeID].empty();
    }
    else
    {
        return false;
    }
}

#endif /* VERTICES_HPP_ */
