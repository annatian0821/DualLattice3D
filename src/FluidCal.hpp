/*
 * FluidCal.hpp
 *
 *  Created on: May 20, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef FLUIDCAL_HPP_
#define FLUIDCAL_HPP_

#include "FluidNetwork.hpp"
#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include <eigen3/Eigen/Sparse>

struct PRpair {
	PRpair () {};
	PRpair (	const double &_p, const double &_r) {
		p = _p;
		r = _r;
	};
	double p;
	double r;

	bool operator< 	(	const PRpair& rhs) const { return fabs(r) < fabs(rhs.r); }

	bool operator== (	const PRpair& rhs) const { return r == rhs.r; }

	bool operator== (	const double rhs)  const { return r == rhs; }
};

class FluidCal: public FluidNetwork {
public:
	FluidCal();
	FluidCal (	const double 		_viscosity,
				const double		_initTime,
				const double        _timeStep,
                const double        _p_tolerance,
				const double 		_q_tolerance,
				const double		_minApecture,
				const double		_maxIncr,
				const double		_maxDecr,
				const double        _baseSensitiveFactor,
				const double		_pMinAbs,
				const double		_pMaxAbs,
				const unsigned		_DOFLimit);

	virtual ~FluidCal();
	bool solve							    ();

/*	bool calculate_static					(	std::pair<bool,bool>& 			qStatus,
												const unsigned					step);*/

	bool checkGlobalConvergence				(   bool                            isMassBalanceOK);

//	void updateUniformP                     (   const double                    pNew);

	std::array<bool,2> fillQ                ();

	void setInitialArea                     (   ConnectLat*                     p_conLat);

	void setBHInitStorage                   ();

	void setNotchList                       (   const std::vector<unsigned>*    p_latList);

	bool interpolate_p                      ();

	double  updateP_Old                     ();

	std::vector<unsigned>   getLatID        ();

	std::vector<IDAndInfo>	getIDAndPressure();

	std::vector<double>		getPressure		(	const double					dummyValue);

	std::vector<double>     getStorageRate  ();

	std::vector<double>     getLeakOff      ();

	void                    interpolatingP0();

    void                    updateRes       ();

    void updatingPressureProfile            (    bool       isUniformPressure);

	double                  updateP0        ();

    bool                    interpolate_pi  ();

    void                    setP0           (   double      pNew);

    double                  getP0           ();

	void                    fill_pList      (   double                          P);

	std::vector<double>		getFlowRate		();

	void                    initizeNextStep ();

	double  getTime							() 	{ return time; };

	void    setTime                         (   double      _time)
	                                            { time = _time; };

	double  getTimeStep                     ()  { return timeStep;};

	double  getSensitiveFactor              ()  { return sensitiveFactor; };

	bool    isPChange                       ()  { return isInjPChange; };
protected:
private:
	double						time;
	double                      timeStep;
	unsigned                    totalTimeStep;
//	double                      initTimeStep;
	double						maxIncrInUse;
	double						maxDecrInUse;
	double                      sensitiveFactor;
	double                      p_tolerance;
	double                      baseSensitiveFactor;
	double                      maxIncrBase;
	double                      maxDecrBase;
	double						viscosity;
	double						q_tolerance;
	double						normRes;
	double                      oldBHstorage;
	double						pMin;
	double						pMinAbs;
	double						pMax;
	double						pMaxAbs;
	double						oldInjP;
	double                      oldTotalArea;
	double						fluidEfficieny;
	double                      initBHstorage;
	unsigned					DOFLimit;
	unsigned					nbPriNum;
	unsigned					DOFseg1,DOFseg2;
	bool						isInjPChange;
	bool                        isFirstStep;
	bool                        checkpMin;
	bool                        checkpMax;
	bool                        isBounded;
	bool                        isFirstP0Iter;
	std::array<PRpair,2> 		iterHist;
	std::vector<unsigned>		DOFtoFNode;
	std::vector<CustomPair<unsigned,double> >       notchList;
	std::vector<unsigned>		fNodeToDOF;
	std::vector<unsigned>       wetFNodeList;
	std::set<unsigned>			activeDOF;
	Eigen::VectorXd 			p;				//pressure at fluid node
	Eigen::VectorXd 			q;				//discharge at fluid node, may be modified when filling H
	Eigen::VectorXd             l;              //leak off rate
	Eigen::VectorXd 			qOld;			//original discharge at fluid node
	unsigned					calNodeNum;
	SpMat						H;
	std::vector<SpMat>			bH;

	void initialize			                ();

	void updateP0                           (   double                          pNew);

	std::pair<bool,double> findPlimits	    ();

	double getLeakOff                       (   const unsigned                  fNodeID);

	double getDeltaLeakOff                  (   const unsigned                  fNodeID);

	double getTotalLeakOff                  ();

	double getDeltaLeakOff 					();

	double getBHstorage                     ();

	void updateLeakOff                      ();

	void setInitLeakOff                     ();

	void fillQInitial						(	const std::vector<IDpair>		qinList);

	void enforceGlobalContinuity 			(	const double					residual,
												const double					sum_s,
												const std::vector<double>*		p_storage);

	bool applyGlobalContinuity				();

	void fillP 								();

	void updateOldp                         ();

	void updateOldP_static                  (   const double                    pNew);

	void set_newP_lastStep                  ();

	bool calNewp_static                     ();

	double getHij							(	unsigned						pipeID);

	bool fillH								();

	bool serialSolve 						(	const unsigned					minIter,
		                                    	const unsigned					maxIter,
		                                    	const double 					precision);
	bool useEigenCGSolver					(	unsigned						maxIteration,
	                                        	double							precision);

	bool sparseCholeskySolve 				();

	Eigen::VectorXd getPreCondJacobi		();

/*	double priInjPressureAdjust				(	const double					pNow,
												const double					pNew);*/

/*	double UpdateInjPressureHist            (   const double                    pNow);*/
	double getP0SecondOrderInterpolate      (   const double                    pNow,
                                                const double                    rNow);
	double  update_p_bounds                 (   const double                    pNow,
                                                const double                    pNew);
/*	double searchP0					        (	const double					pNow,
												const double					rNow);*/

	double linearPredict					(	const double					pNow,
												const double					rNow);


	bool updateActiveConFNode				();

	void printVectors						(	const unsigned 					step,
												const bool 						isPrintp,
												const bool 						isPrintq,
												const bool 						isPriH);
	template <typename T>
	void log								(	std::ostream*				p_file,
												const T						data);

	double  getResidual 					();

	double  getActiveTotalArea              ();

	void logPQ								(	const unsigned				step,
												const std::vector<double>	P,
												const std::vector<double>	Q);

	void logRes								(	const unsigned				step,
												const double				nowP,
												const double				newP,
												const double				normRes);
	void 	test();
};

#endif /* FLUIDCAL_HPP_ */
