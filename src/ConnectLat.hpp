/*
 * ConnectLat.hpp
 *
 *  Created on: Jan 30, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef CONNECTLAT_HPP_
#define CONNECTLAT_HPP_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "Vertices.hpp"
#include "ProProcess.hpp"
struct FLatData {
   FLatData (	unsigned _LatticeID, unsigned _step, double _time, bool _isReconnect) {
	   LatticeID = _LatticeID;
	   time = _time;
       step = _step;
       isReconnect = _isReconnect;
       type = 0;
   }
   FLatData (   unsigned _LatticeID, unsigned _step, double _time, unsigned _type, bool _isReconnect) {
       LatticeID = _LatticeID;
       time = _time;
       step = _step;
       isReconnect = _isReconnect;
       type = _type;
   }
   unsigned LatticeID;
   double time;
   unsigned step;
   unsigned type;
   bool	 isReconnect;
   bool operator< (const FLatData& rhs) const { return LatticeID < rhs.LatticeID; }
   bool operator== (const FLatData& rhs) const { return LatticeID == rhs.LatticeID; }
   bool operator== (const unsigned _LatticeID) const { return LatticeID == _LatticeID; }
};


typedef std::list<FLatData>					LatGroup;
typedef std::list<LatGroup>					LatNetwork;


class ConnectLat {
    friend class Vertcies;
    friend class GLatRemoval;

    struct compareLatGroupReverse {
        bool operator() (const LatGroup &a,const LatGroup &b) {
            return a.size()> b.size();
        }
    };
public:
    ConnectLat();
    ConnectLat	(	const double		_minApecture) {
    	minApecture = _minApecture;
    	centre[0] = Nx*UnitLength*0.5; centre[1] = Ny*UnitLength*0.5; centre[2] = Nz*UnitLength*0.5;
    	cntStep = 0 ;cntTotal = 0;
    	sumDzStep = 0.0 ;sumDzTotal = 0.0;
    	sumSqDzStep = 0.0 ;sumSqDzTotal = 0.0;
    	isLogOffset = false;
    }
    virtual ~ConnectLat();

    void		putLattice				(	const unsigned 			LatticeID,
    										const unsigned			step,
    										const double			time,
    										const unsigned          type,
                                        	Vertices*				p_verti,
                                        	ConnectLat*				p_self);

    void		putMultiConLat			(	std::vector<unsigned>* 	p_latList,
    										const unsigned			step,
    										const double			time);

    unsigned	getListSize				();

    unsigned	getLatNum				();

    unsigned 	reConnectLat			();

    void 		disConnectLat			(   const std::vector<unsigned>*	p_latList);

    double		getTotalFluidVol		();

    std::vector<std::vector<unsigned> >		getAllStep		();

    std::vector<std::vector<double> >		getAllTime		();

    std::vector<std::vector<unsigned> >     getAllType      ();

    std::vector<unsigned>  					getGroup			(	const unsigned 			groupID);
    std::vector<unsigned>  					getGroupActive		(	const unsigned 			groupID,
																	const double			minApecture);

    std::vector<std::vector<unsigned> > 	getAllGroup		    ();

    std::vector<unsigned> 					getAllGroupVec	    ();

    std::vector<unsigned>                   getAllLat           (   const unsigned          typeID);

    std::vector<unsigned>                   getClusterList      ();

    void		setLogOffset			(	const bool			_isLogOffset) { isLogOffset = _isLogOffset;};
    void		putOffset				(	const unsigned		LatticeID,
    										Vertices*			p_verti);

    void		logOffset				(	const double		time);

    void		logOffset				(	const unsigned		step);

    void		logConOffset			(	const double		time,
    										Vertices*			p_verti);

    void		logConOffset			(	const unsigned		step,
        									Vertices*			p_verti);

    void        logFractureVol          (   const double        sumF,
                                            const unsigned      totalStep);

    void		logFLat					(	const double		newTime,
    										const unsigned 		rmLatAcc);

    void		swapGroupFrontBack		();

    void		printAll				(	const unsigned			step);
    unsigned    setType                 (   const unsigned          _LatticeID,
                                            const unsigned          _type);
    void		reset					();

    std::vector<unsigned>  getReConLat  ();

    std::vector<unsigned>  getDisConLat ();

    Point       getMainCrackCenter      ();

    std::vector<std::vector<unsigned> >  getPerimeter (   const Point         center,
                                                          const double        lavg,
                                                          const std::array<unsigned,3>  xyz);

private:
    LatNetwork									latNetwork;
    std::set<unsigned>							reConList;
    std::set<unsigned>                          disConList;
    double										minApecture;
    Point										centre;
    bool										isLogOffset;
    unsigned									cntStep,cntTotal;
    double										sumDzStep,sumDzTotal;
    double										sumSqDzStep,sumSqDzTotal;

    void	mergeLatList		(	std::vector<std::list<std::list<FLatData> >::iterator>* 	p_itListToMerge,
                                    const Vertices*												p_verti);
    void	rmInnerCrack			();


    void 	reConnectLatSingle	(	const unsigned			LatticeID);

    bool	setIsReconnect		(	const unsigned			LatticeID,
    								const bool				inBool);

    bool	isListEmpty			(	const std::list<FLatData>& listFLatData)
    {
    	return listFLatData.empty();
    }

    bool	isOneEndUnstable	(	const FLatData& fLatData) {
    	return (GNodeTab[LatTab[fLatData.LatticeID].nb[0]].type==Type_Unstable)
    	       	||(GNodeTab[LatTab[fLatData.LatticeID].nb[1]].type==Type_Unstable);
    }

    bool    isOneEndUnstable    (   const unsigned LatticeID) {
        return (GNodeTab[LatTab[LatticeID].nb[0]].type==Type_Unstable)
                ||(GNodeTab[LatTab[LatticeID].nb[1]].type==Type_Unstable);
    }

    bool	isBothEndStable		(	const FLatData& fLatData) {
    	return (GNodeTab[LatTab[fLatData.LatticeID].nb[0]].type!=Type_Unstable)
    	       	&&(GNodeTab[LatTab[fLatData.LatticeID].nb[1]].type!=Type_Unstable);
    }

    void    sortLatGroup ();


//    template <class T>
    struct 	isEmptyList
    {
        bool operator() (const std::list<FLatData>&	testList)
        {
            return testList.empty();
        }
    };
    struct isBothEndUnstable
    {
    	bool	operator()	(	const FLatData& fLatData)
       {
       	return (GNodeTab[LatTab[fLatData.LatticeID].nb[0]].type==Type_Unstable)
       			&&(GNodeTab[LatTab[fLatData.LatticeID].nb[1]].type==Type_Unstable);
       }
    };
};

#endif /* CONNECTLAT_HPP_ */
