/*
 * Statistics.hpp
 *
 *  Created on: Nov 12, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "MiscTools.hpp"
#include "ConnectLat.hpp"
#include "ProProcess.hpp"
#include "FluidCal.hpp"
#include "Stress.hpp"
#include "PCG_Eigen.hpp"
#include <eigen3/Eigen/Dense>



typedef std::array<std::array<std::vector<std::vector<double> >,3>,3>   LatPriDirList;
class Statistics {
public:
    Statistics();

    Statistics	(	char* 		                inPath,
                    char* 		                inPrefix,
                    unsigned	                _maxCoordNo,
                    unsigned                    _pxNodeNum,
                    std::array<unsigned,2>      _refineNode)
    {
        path = inPath;
        prefix = inPrefix;
        lastRecordStep = 0;
        lastLatNum = 0;
        if (inPrefix[0]!='\0')
        {
            strcat(prefix,"-");
        }
        pxNodeNum = _pxNodeNum;
        maxCoordNo = _maxCoordNo;
        tailLen = 0.05;
        refineNode = _refineNode;

        xyz[0] = 0;
        xyz[1] = 1;
        xyz[2] = 2;
        if ((GeometryID==4)) {
            xyz[0] = 0;
            xyz[1] = 2;
            xyz[2] = 1;
        }
        getInnerBoundaryNodeList ();
        fillRefineLatList ();
    };
    virtual ~Statistics();

    void setCrackMonthLatList   (   const double                width,
                                    std::vector<unsigned>*      p_latList);

    void Fluid_pressure         (   FluidCal*                   p_fCal,
                                    unsigned                    step);

    void compute		(	const unsigned 		resolution);

    void writeCoordNo 	(	const unsigned		                step,
                            const std::vector<unsigned>&        nodeList);

    void writeInLatAll 	(	const unsigned			pxNodeNum,
            				const unsigned			step,
                            const unsigned			resolution,
                            ConnectLat*             p_conLat);

    void writeStress	(	const std::vector<unsigned>*    p_weakPlaneLatID,
    						const double                    sumF,
    						const unsigned				    timeStep,
    						const double				    time);

    void writeStress    (   const std::vector<unsigned>*    p_weakPlaneLatID,
                            ConnectLat*                     p_conLat,
                            const double                    sumF,
                            const unsigned                  timeStep,
                            const double                    time);

    void writeWeakPlaneStress    (   ConnectLat*                     p_conLat,
                                     PCG_Eigen*                      p_EigenPCG,
                                     std::vector<unsigned>*          p_weakPlaneLatID,
                                     const double                    sumF,
                                     const unsigned                  step);

    void writeLatStress (   const double                sumF,
                            const unsigned              timeStep,
                            ConnectLat*                 p_conLat,
                            PCG_Eigen*                  p_EigenPCG);

    void writeStressStrain	(	const unsigned			faceDir,
    							const unsigned			dir,
    							const double			sumF,
    							const unsigned          totalStep);

    void writeStressStrain  (   const double            sumF,
                                const unsigned          totalStep);

    void writeCrackMouthDispHist      (     const double                sumF);

    void writeMaxNorm	(	const double				maxNorm,
    						const double				sumF,
    						const unsigned				step);

    void writeTotalFluidVol	(	ConnectLat*					p_conLat,
    							const unsigned				step);

    void writelatPolar      (   const unsigned              resolution,
                                const unsigned              step,
                                const double                sumF,
                                ConnectLat*                 p_conLat,
                                PCG_Eigen*                  p_EigenPCG);

    void writelatPolar      (   const unsigned              resolution,
                                const unsigned              coarseness,
                                const unsigned              step,
                                const double                sumF,
                                ConnectLat*                 p_conLat,
                                PCG_Eigen*                  p_EigenPCG);

    bool writeFailLatPolar  (   const unsigned             resolution,
                                std::vector<unsigned>      latList,
                                const unsigned             step,
                                const double               sumF,
                                const char*                vName,
                                PCG_Eigen*                 p_EigenPCG);

    void writeLatCluster    (   ConnectLat*                 p_conLat,
                                const unsigned              resolution,
                                const unsigned              step);

    void writeRLatStats     (   ConnectLat*                 p_conLat,
                                PCG_Eigen*                  p_EigenPCG,
                                const double                sumF,
                                const unsigned              step);

    void writeLatClusterAgg (   ConnectLat*                 p_conLat,
                                const unsigned              step);

    void printPDF           (   std::vector<double>*        p_vList,
                                const char*                 vName,
                                const unsigned              step,
                                const unsigned              resolution,
                                const double                tailSize);

    LatPriDirList fillLatPriDirList (  const unsigned      resolution,
                                       const unsigned      step);

    bool fillInnerLatList   (   const double                  minAreaFactor,
                                const std::array<double,6>    limits,
                                ConnectLat*                   p_conLat);

    std::array<double,4> getStatInfo  (   const std::vector<double>  &vList);

protected:

    struct RLat {
        RLat ( double    _r,
               double    _a,
               double    _dz,
               unsigned  _LatticeID)
            {
                r = _r;
                a = _a;
                dz = _dz;
                LatticeID = _LatticeID;
            }

        unsigned                        LatticeID;
        double                          r;
        double                          a;
        double                          dz;

        bool operator<  (   const RLat& rhs) const { return r < rhs.r; }
        bool operator== (   const RLat& rhs) const { return r == rhs.r; }
        bool operator== (   const unsigned rhs) const { return LatticeID == rhs; }
    };

    unsigned	maxCoordNo;
    unsigned    pxNodeNum;
    double 		unitAng;
    double 		unitAngInDeg;
    double      tailLen;
    std::array<unsigned,2>      refineNode;
    std::vector<unsigned>       refineLatList;
    std::array<unsigned,3>      xyz;

    unsigned	lastRecordStep;
    unsigned	lastLatNum;
    std::array<std::vector<std::vector<unsigned> >,3>     latAngList;
    std::vector<unsigned>                   innerNodeList;
    std::vector<unsigned>                   innerLatList;

    char*		path;
    char*		prefix;

    std::vector<unsigned>	                eLatIDList;
    std::vector<double>		                eLatAngList[Dim];
    std::vector<unsigned>	                statLatAng[Dim];
    std::vector<unsigned> 	                statLatStiffness[Dim];
    std::vector<unsigned> 	                cnt_dir[Dim];
    std::vector<unsigned>	                inBoundary[6];
    std::array<std::vector<unsigned>,2>     crackMonthLatList;

    double mean_dir[Dim];
    double sd_dir[Dim];

    void fill_eLatIDList();

    void fill_eLatIDListSimple();

    void updateStat (const double var, const unsigned resolution,const double min,const double max, std::vector<unsigned>* p_Stat);

    void fillRefineLatList ();

    void getInnerBoundaryNodeList ();

    std::vector<RLat> getRLatList             ( std::vector<std::vector<unsigned> >      &periList,
                                                std::vector<unsigned>*                   p_latList,
                                                const Point                              center);

    std::vector<RLat> getRNodeList            ( std::vector<std::vector<unsigned> >      &periList,
                                                std::vector<unsigned>                    &nodeList,
                                                const Point                              center);
    void getPlaneInfo 		(	std::string 			plane,
                                unsigned* 				xyz0,
                                unsigned* 				xyz1,
                                unsigned* 				ID);


    void latDir				(	const char* 			plane,
                                const unsigned 			resolution);

    void glatDir			(	const char* 			plane,
                                const unsigned 			resolution);

    void glatDir			(	const char* 			plane,
                                const unsigned 			resolution,
                                const unsigned			step);

    void latStiffness		(	const char* 			plane,
                                const unsigned 			resolution);

    void writeLatLen 		(	const unsigned			step,
                                const unsigned			resolution);

    void writeLatArea 		(	const unsigned			step,
                                const unsigned			resolution,
                                const bool				isLogScale);

    void writeInLatArea 	(	const unsigned			step,
                                const unsigned			resolution,
                                const bool				isLogScale);

    void writeNodeVol       (   const unsigned                      step,
                                const unsigned                      resolution,
                                const std::vector<unsigned>&        nodeList);

    void writeLatBreakDisp (	const unsigned			step,
                                const unsigned			resolution,
                                const bool				isLogScale);

    double avgFaceDisp		(	const unsigned			face,
    							const unsigned			dir,
    							const bool				isInner);

    double avgFaceCoord		(	const unsigned			face,
    							const unsigned			dir,
    							const bool				isInner);

    void print              (   const std::vector<unsigned>*    p_stat,
                                const unsigned                  resolution,
                                const double                    min,
                                const double                    max,
                                const unsigned                  pID,
                                const char*                     varName1,
                                const char*                     varName2);

    void print2 			(	const std::vector<double>*	p_vList,
    							const char*					vName,
    							const unsigned				step,
                                const unsigned				resolution,
                                const bool					isLogScale);

    void printPDF           (   std::vector<std::vector<double> >        &vList,
                                const char*                              vName,
                                const unsigned                           step,
                                const unsigned                           resolution,
                                const double                             tailSize);

    void printData2D        (  std::vector<std::pair<double,std::pair<double,double> > >   &pList,
                               const char*                                                 vName,
                               const double                                                sumF,
                               const unsigned                                              step,
                               const unsigned                                              resolution1,
                               const unsigned                                              resolution2);

    void printCluster       (   std::vector<unsigned>*                  p_vList,
                                const char*                             vName,
                                const unsigned                          step,
                                const unsigned                          resolution);

    void getStatInfo 		(	const std::vector<double>*	p_vList,
    							double*						p_mean,
    							double*						p_SD);

    std::vector<RLat> getRLatList ( std::vector<std::vector<unsigned> >      &periList,
                                    const Point                              center);

    bool fillInnerNodeList  ();

    bool fillInnerNodeList  (   const std::array<double,6>  limits);

    bool fillInnerLatList   (   const double                minAreaFactor);

    bool fillInnerLatList   (   const double                minAreaFactor,
                                ConnectLat*                 p_conLat);

    void rmBoundaryLat      (   std::vector<unsigned>       &latList);

    void limitInnerLatList  (   const double                minAreaFactor,
                                const std::array<double,6>  limits);

    unsigned fillAngList    (   const unsigned              resolution,
                                ConnectLat*                 p_conLat);

    bool isInside           (   const unsigned                  LatticeID,
                                const std::array<double,6>      limits);

    bool isInside           (   const Point                     centroid,
                                const std::array<double,6>      limits);

    std::array<double,6>    getLimits ( const Point     centre);
};

class FailStats : public Statistics {
public:
    FailStats() {};
    FailStats		(	char* 			inPath,
                        char* 			inPrefix,
                        const unsigned	_maxCoordNo)
    {
        path = inPath;
        prefix = inPrefix;
        lastRecordStep = 0;
        lastLatNum = 0;
        cnt = 0;
        if (inPrefix[0]!='\0')
        {
            strcat(prefix,"-");
        }
        maxCoordNo = _maxCoordNo;
    };

    virtual ~FailStats() {};

    void pushNode		(	const unsigned 		NodeID,
                            const unsigned 		step)
    {
        fNodePerStep.resize(step+1);
        fNodePerStep[step].push_back(NodeID);
    };

    void pushLat		(	const unsigned 		LatticeID,
                            const unsigned 		step)
    {
        fLatPerStep.resize(step+1);
        fLatPerStep[step].push_back(LatticeID);
    };

    void write			(	const unsigned 			step,
                            const unsigned 			resolution);
//							const std::string 	fileName);
private:
    unsigned							cnt;
    std::vector<std::vector<unsigned> >	fNodePerStep;
    std::vector<std::vector<unsigned> >	fLatPerStep;
};

class EnergyCal {
	friend class PCG_Eigen;
public:
	void    fillForceNodeList       (   const unsigned  faceID,
	                                    const unsigned  dir);

    void    record      (   FluidCal*                   p_fCal);

    void    record      ();

    void    compute     (   double                      time);

    void    compute     (   unsigned                    step,
                            double                      sumF);

    void    accArea     (   double                      fracArea);

    void    accumulate  ();

    void    setAreaZero ();

    void    accBrokenLatEnergy  (   double              energy);

    void    setBrokenLatEnergyZero  ();

    void    fillForceNodeList       (   const unsigned  faceID);
private:
    double 	 getTotalStrainEnergy	(	const unsigned		LatticeID);

    Eigen::MatrixXd    getLocalSM   (   const unsigned      LatticeID);

    std::vector<double>     oldApecture;
    std::vector<double>     oldPressure;
    std::vector<unsigned>   forceNodeIDList;
    std::vector<double>     oldDisp;
    std::vector<unsigned>   fLatID;
    unsigned                loadDir;
    double                  accER;
    double                  ER;
    double                  oldTotalStrainEnergy;
    double                  totalFracArea;
    double                  accFracArea;
    double                  totalBrokenLatEnergy;
};


#endif /* STATISTICS_H_ */
