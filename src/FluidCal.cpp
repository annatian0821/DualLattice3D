/*
 * FluidCal.cpp
 *
 *  Created on: May 20, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "FluidCal.hpp"
/*
FluidCal::FluidCal() {
	// TODO Auto-generated constructor stub

}
*/
/*
FluidCal::~FluidCal() {
	// TODO Auto-generated destructor stub
}
*/
FluidCal::~FluidCal() {}
FluidCal::FluidCal() {
//	maxCvDist = 0.0;
//	maxCoordNo = 0;
	q_tolerance = 1e-6;
	minApecture = 1e-5;
	DOFLimit = 10000;
	time = 0.0;
	oldTotalArea = 0.0;
}

FluidCal::FluidCal			(	const double 		_viscosity,
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
								const unsigned		_DOFLimit) {
	viscosity = _viscosity;
	p_tolerance = _p_tolerance;
	q_tolerance = _q_tolerance;
	minApecture = _minApecture;
	maxIncrBase = maxIncrInUse = _maxIncr;
	maxDecrBase = maxDecrInUse = _maxDecr;
    sensitiveFactor = baseSensitiveFactor = _baseSensitiveFactor;

	oldInjP = InitForce;
//	p0 = InitForce;
	pMinAbs = _pMinAbs;
	pMin = InitForce*LoadDecrFactor;
	pMax = InitForce*LoadIncrFactor;
	pMaxAbs = _pMaxAbs;
	DOFLimit = _DOFLimit;
	time = _initTime;
	initBHstorage = 0.0;
	oldBHstorage = 0.0;
	timeStep = _timeStep;
	isInjPChange = true;
	isFirstStep = true;
	isBounded = false;
	checkpMin = checkpMax = false;
}


bool FluidCal::solve		    ()
{
	std::pair<bool,bool>		status(false,false);
//	totalTimeStep = _totalTimeStep;
	const double    ratioEnhance = 10.0;
	if (pipeList.empty()) {
		calNodeNum = 0;
		tripOut << "There is no pipe available for fluid flow, no calculation is done!" << '\n';
		return false;
	}
	dualOut << "Start calculation of fluid pressure..." << " time = " << time << '\n';
	Clock clock;
	clock.start("Fluid");
//	bool isFluidNetworkChanged = true;
//	unsigned fCalStep = 0;
//	dualOut << "--------fCalStep = " << fCalStep << " starts--------" << '\n';
//	qStatus = fillQ(initialize());
	maxIncrInUse = maxIncrBase;
	maxDecrInUse = maxDecrBase;
	double p0 = this->getPriInjP();
	if ((calNodeNum==0)||(!fillH())) {
	    dualOut << "Only injection node is active, solving is not required, pressure is set to 0.0" << '\n';
	    p.fill(0.0);
	    std::cout << " maxDecrBase  = " << maxDecrBase << " maxDecrInUse = " << maxDecrInUse << std::endl;
	} else {
	    unsigned minIter = (unsigned) (2.5*std::sqrt(fNodeList.size()));
	    unsigned maxIter = 15*minIter;
	    if (10*calNodeNum<DOFLimit) {
	        status.second = sparseCholeskySolve ();
	        dualOut << "Cholesky residual = " << getResidual() << '\n';
	        if (!status.second) {
	            tripOut << "Cholesky is not stable, try to solve by PCG..." << '\n';
	            status.second = serialSolve(minIter,maxIter,PCG_Precision);
	            dualOut << "PCG residual = " << getResidual() << '\n';
	            if (!status.second) {
	                tripOut << "Calculation unstable, reset to hydrostatic pressure and start again";
	                return false;
	        }
	    }
/*	        else {
	        status.second = serialSolve(minIter,maxIter,PCG_Precision);
	        dualOut << "PCG residual = " << getResidual() << '\n';
	        if (!status.second) {
	            tripOut << "PCG is not stable, try to solve by Cholesky..." << '\n';
	            status.second = sparseCholeskySolve ();
	            dualOut << "Cholesky residual = " << getResidual() << '\n';
	            if (!status.second) {
	                return false;
	            }
	        }*/
	    }
	}
//	fCalStep++;
	status.second = !(!(p0!=p0)||(p0>Huge));
//	bool is_pConverge = calNewp();
//	qStatus.second = (is_pConverge);
//	isInjPChange = false;
	dualOut << "--------fluid calculation finished" << '\n';
//	dualOut << "Fluid calculation completed..." << '\n';
	clock.get();
	return true;
}

/*
bool FluidCal::calculate_static		(	std::pair<bool,bool>& 	qStatus,
										const unsigned			totalTimeStep)
{
	std::pair<bool,bool>		status(false,false);
	if (pipeList.empty()) {
		calNodeNum = 0;
		tripOut << "There is no pipe available for fluid flow, no calculation is done!" << '\n';
		return false;
	}
	dualOut << "Start calculation of fluid pressure..." << " time = " << time << '\n';
	qStatus = fillQ(initialize());
	qStatus.second = true;
//	p.resize(calNodeNum);
//	p.fill(p0);
//	updateOldp();
//	updateActiveConFNode();
	isInjPChange = false;
	dualOut << "Fluid static calculation completed..." << '\n';
	return true;
}
*/

/*
bool FluidCal::updateCalInfo	(	std::pair<bool,bool>& 		qStatus)
{
	static unsigned convergeCnt = 0;
	static bool		isQbounded = false;
	static unsigned relaxStep = 0;
	if (qStatus.second) {
        dualOut << "p0 = " << p0 << " pMax-pMin = " << pMax-pMin << " pMax = " << pMax << " pMin = " << pMin << std::endl;
        if ((qStatus.first)) {
            dualOut << "Global flow rate convergence achieved, proceed to next time step" << '\n';
            initizeNextStep(convergeCnt,relaxStep,isQbounded);
            isInjPChange = true;
            relaxStep = 0;
            return true;
        } else if (((pMax-pMin)/p0<1e-3)) {
            dualOut << "Pressure difference is small, proceed to next time step" << '\n';
            initizeNextStep(convergeCnt,relaxStep,isQbounded);
            isInjPChange = true;
            relaxStep = 0;
            return true;
        }
        this ->updateOldP_static(updatePriInjPressure());
        dualOut << "Injection pressure updated" << std::endl;
	}
	relaxStep++;
	return false;
}
*/

/*
void FluidCal::set_newP_lastStep () {
    for (unsigned fNodeID=0; fNodeID < fNodeList.size(); fNodeID++) {
        fNodeList[fNodeID].new_p = fNodeList[fNodeID].current_p;
    }
}
*/
void FluidCal::initizeNextStep	()
{
	char fileName[255];
	static unsigned old_fNodeNum = 0;
	static int old_activefNodeNum = 0;
	static bool isFirst = true;
	dualOut << "FluidCal::initizeNextStep ... " << std::endl;
	std::sprintf(fileName,"%s/%spHist.txt",OutputSubFolder,OutputFilePrefix);
	std::ofstream file(fileName, std::ios::out | std::ios::app);
	updateLeakOff ();
	if (isFirst) {
	    file << "time" << '\t' << "P" << '\t' << "increase in total fluid node" << '\t'
	         << "increase in active fluid node" << '\t' << "calNodeNum" << '\t' << "Fluid node num" << '\t'
	         << "active node ratio" << '\t'<< "fluidEfficieny" << '\n';
	    file << 0 << '\t' << 0 << '\t' << 0 << '\t'
	             << 0 << '\t' << calNodeNum << '\t' << fNodeList.size() << '\t'
	             << 1.0 << '\t'<< 1.0 << 0.0 << 1.0 << '\n';
	    isFirst = false;
	}
	file << time << '\t' << getPriInjP() << '\t' << fNodeList.size() - old_fNodeNum << '\t'
	        << ((int) calNodeNum) - old_activefNodeNum << '\t' << calNodeNum << '\t' << fNodeList.size() << '\t'
	        << ((double) calNodeNum+1)/fNodeList.size()<< '\t'<< fluidEfficieny  << '\n';
	dualOut << " Current time = " << time << " timeStep = " << timeStep << std::endl;
	time += timeStep;
	dualOut << " Current time new = " << time << " timeStep = " << timeStep << std::endl;
	double p0 = this->getPriInjP();
	pMin = p0;
	pMax = p0*(time+timeStep)/time;
	isBounded = false;
//	updateUniformP(p0);
	old_fNodeNum = fNodeList.size();
	old_activefNodeNum = calNodeNum;
	checkpMin = checkpMax = false;
	double totalArea = getActiveTotalArea ();
	if (totalArea > 2.0*oldTotalArea) {
		dualOut << "totalArea, oldArea = " << totalArea << ',' << oldTotalArea << std::endl;
	    timeStep *=2.0;
	    dualOut << "************TimeStep is now increased to " << timeStep << "************" << std::endl;
	    oldTotalArea *= 2.0;
	}
//	convergeCnt++;
//	p0 = getPriInjP();
//	set_newP_lastStep();
//	isQbounded = false;
	isFirstStep = false;
	isFirstIter = true;
//	relaxStep = 0;
}

/*
double FluidCal::UpdateInjPressureHist       (   const double        pNow)
{
    if (normRes>0.0) {
        pMin = pNow;
        iterHist[0] = PRpair(pMin,normRes);
        tripOut << "Set pMin = pNow = " << pMin << '\n';
    } else {
        pMax = pNow;
        iterHist[1] = PRpair(pMax,normRes);
        tripOut << "Set pMax = pNew = " << pMax << '\n';
    }
}*/

double FluidCal::getActiveTotalArea ()
{
    double totalArea = 0.0;
    for (unsigned i=0; i<fNodeList.size(); i++) {
        if (fNodeList[i].isActive()) {
            double LatticeID = fNodeList[i].LatticeID;
            if (ProblemID==2) {
                if (!std::binary_search(initCrackList.begin(),initCrackList.end(),LatticeID)) {
                    totalArea += LatTab[fNodeList[i].LatticeID].area;
                }
            } else {
                totalArea += LatTab[fNodeList[i].LatticeID].area;
            }
        }
    }
    return totalArea;
}

void FluidCal::setNotchList          (   const std::vector<unsigned>*    p_latList)
{
    dualOut << "FluidCal::setNotchList..." << std::endl;
    notchList.clear();
    GLatForce lforce;
    CustomPair<unsigned,double> cp;
    for (unsigned i=0 ; i< p_latList->size(); i++ ) {
        unsigned LatticeID = (*p_latList)[i];
        cp.i1 = LatticeID;
        cp.i2 = lforce.getCrackOpening(LatticeID)*LatTab[LatticeID].area;
        notchList.push_back(cp);
    }
    std::sort(notchList.begin(),notchList.end());
}

void FluidCal::setBHInitStorage ()
{
    GLatForce lforce;
    initBHstorage = 0.0;
    if ((GeometryID==1)||(GeometryID==4)) {
        for (unsigned i=0; i<initCrackList.size(); i++) {
            unsigned LatticeID = initCrackList[i];
            initBHstorage += lforce.getCrackOpening(LatticeID)*LatTab[LatticeID].area;
        }
    }
    dualOut << "Initial borehole storage = " << initBHstorage << std::endl;
}

double FluidCal::getBHstorage ()
{
    double storage = 0.0;
    GLatForce lforce;
    for (unsigned i=0; i<initCrackList.size(); i++) {
        unsigned LatticeID = initCrackList[i];
        storage += lforce.getCrackOpening(LatticeID)*LatTab[LatticeID].area;
    }
    std::cout << "Storage , initBHstorage = " << storage << ',' << initBHstorage << std::endl;
    return std::max(0.0,storage - initBHstorage);
}
void FluidCal::setInitLeakOff ()
{
    dualOut << "Set initial leak off..." << std::endl;
    double totalInitLeakOff = 0.0;
    for (unsigned i=0; i<fNodeList.size(); i++) {
        fNodeList[i].t_exp = InitTime-timeStep;
        fNodeList[i].leakVol = 2.0*fNodeList[i].Cl*std::sqrt(InitTime-timeStep);
        totalInitLeakOff += fNodeList[i].leakVol;
    }
    dualOut << "Initial leak off volume = " << totalInitLeakOff
            << " initial fluid efficiency = " << (InjectRate*InitTime-totalInitLeakOff)/(InjectRate*InitTime) << std::endl;
}


void FluidCal::setInitialArea (     ConnectLat*         p_conLat)
{
    std::vector<unsigned>   group0 = p_conLat->getGroup(0);
    oldTotalArea = 0;
    for (unsigned i=0; i<group0.size(); i++) {
        oldTotalArea += LatTab[group0[i]].area;
    }
    if (ProblemID==2) {
        for (unsigned i=0; i<initCrackList.size(); i++) {
            oldTotalArea -= LatTab[initCrackList[i]].area;
        }
    }
    dualOut << "Set oldTotalArea = " << oldTotalArea << std::endl;
}

std::pair<bool,double> FluidCal::findPlimits	()
{
    const double adjRatio = 1.000;
    std::pair<bool,double>  out;
    if (!checkpMin) {
        if (normRes>0.0) {
            checkpMin = true;

            iterHist[0] = PRpair(pMin,normRes);
            tripOut << "pMin is captured, now try to capture pMax" << '\n';
            pMin = getPriInjP();
            if (normRes<1.0) {
                pMax = std::min(std::min(pMin/(1.0-normRes)*adjRatio,maxIncrInUse*pMin),pMaxAbs);
            } else {
                pMax = maxIncrInUse*pMin;
            }
            out.first = false;
            out.second = pMin;
        } else {
            pMin = getPriInjP();
            if (std::fabs(pMin-pMinAbs)<Tiny) {
                time +=timeStep;
                tripOut << "pMin is too large but cannot be further reduced as it reaches pMinAbs already,increase the time by one timeStep "
                        << ", time = " << time << '\n';
            }
            pMin = std::max(std::max(pMin/(1.0-normRes)/adjRatio,maxDecrInUse*pMin),pMinAbs);
            tripOut << "pMin is too large, try again with pMin = pMin/(1.0-normRes)/adjRatio "
                    << "or pMinAbs, whichever is larger, pMin = " << pMin << '\n';
            out.second = pMin;
            out.first = false;
        }
    } else if (!checkpMax) {
        if (normRes<0.0) {
            checkpMax = true;
            iterHist[1] = PRpair(pMax,normRes);
            tripOut << "pMin and pMax are captured, now solution is bounded" << '\n';
            out.first = true;
            out.second = pMax;
        } else {
            pMax = std::min(pMax/(1.0-normRes)*adjRatio,pMaxAbs);
            tripOut << "pMax is too small, try again with pMax = pMin-normRes/(normRes+1.0)*pMin*adjRatio "
                    <<    "or pMaxAbs, whichever is smaller, pMax = " << pMax << '\n';
            out.second = pMax;
            out.first = false;
        }
    }
    dualOut << " out.First = " << out.first << " out.second = " << out.second << "checkMax = "
            << checkpMax << " checkMin" << checkpMin << std::endl;
    dualOut << pMin << ' ' << pMax << ' ' << std::endl;
    return out;
}
    /*
    dualOut << "Searching pMax and pMin..." << std::endl;
	const double adjRatio = 1.005;
	if (!checkpMin) {
		if (normRes>0.0) {
			checkpMin = true;
			iterHist[0] = PRpair(pMin,normRes);
			tripOut << "pMin is captured" << '\n';
			pMin = p0;
			dualOut << "Update pMin = " << pMin << " , pMax = " << pMax << std::endl;
		}   else {
			if ((pMin-pMinAbs)<Tiny) {
				time +=timeStep;
				tripOut << "pMin is too large but cannot be further reduced as it reaches pMinAbs already,increase the time by one timeStep "
						<< ", time = " << time << '\n';
			}
			pMin = std::max(std::min(std::max(p0/(1.0-normRes)/adjRatio,maxDecrInUse*maxDecrInUse*p0)
			        ,maxDecrInUse*maxDecrInUse*p0),pMinAbs);
			tripOut << "pMin is too large, try again with pMin = pMin/(1.0-normRes)/adjRatio "
					<< "or pMinAbs, whichever is larger, pMin = " << pMin << '\n';
            dualOut << "Update pMin = " << pMin << " , pMax = " << pMax << " , pNow = " << p0 << std::endl;
		}
	}
	if (!checkpMax) {
		if (normRes<0.0) {
			checkpMax = true;
			iterHist[1] = PRpair(pMax,normRes);
			tripOut << "pMin and pMax are captured, now solution is bounded" << '\n';
			pMax = p0;
		}
		else {
		    dualOut << "normRes = " << normRes << std::endl;
		    pMax = std::min(p0/(std::fabs(1.0-normRes))*adjRatio,maxIncrInUse*p0);
			dualOut << "pMax is too small, try again with pMax = pMax/(1.0-normRes)*pMin*adjRatio "
					<<	"or pMaxAbs, whichever is smaller, pMax = " << pMax << '\n';
			dualOut << "Update pMax = " << pMax << " , pMin = " << pMin << " , pNow = " << p0 << std::endl;
		}
	}
	dualOut << "FluidCal::logConvgHist ... " << std::endl;
	char  fileName[255];
	std::sprintf(fileName,"%s/%sconvHist%04d.txt",OutputSubFolder,OutputFilePrefix,totalTimeStep);
	std::ofstream file(fileName, std::ios::out | std::ios::app);
	bool isBounded = (checkpMax&&checkpMin);
	if (isBounded) {
	    dualOut << "Solution bounded, pMax = " << pMax << " pMin = " << pMin << std::endl;
	}
	return isBounded;
}*/
/*
double FluidCal::update_p_bounds                       (   const double        pNow,
                                                           const double        pNew)
{
    if (normRes>0.0) {
        pMin = pNow;
        iterHist[0] = PRpair(pMin,normRes);
        tripOut << "Set pMin = pNow = " << pMin << '\n';
    } else {
        pMax = pNow;
        iterHist[1] = PRpair(pMax,normRes);
        tripOut << "Set pMax = pNew = " << pMax << '\n';
    }
    if ((pNew>pMin)&&(pNew<pMax)) {
            return pNew;
        }
    dualOut << "pNew is adjusted to (minP+maxP)/2 = " << (pMin+pMax)/2.0 << '\n';
    return (pMin+pMax)/2.0;
}*/


std::array<bool,2> FluidCal::updateQ					()
{
	//Storage of primary injection nodes which are excluded in solving system of equations
	this->initialize();
	if (priINodeList.empty()) {
		tripOut << "---WARNING---No injection node!!!" << '\n';
	}
	const std::vector<IDpair> qinList = this->initialize();
	q.resize(calNodeNum);
	q.fill(0.0);
	l.resize(calNodeNum);
	l.fill(0.0);
	unsigned LatticeID;
	double sum_qin = 0.0;
	double sum_s = 0.0;
	double sum_l = 0.0;
//	double p0 = getPriInjP();
	static double oldDeltaQ = Huge;

	GLatForce lforce;
	for (std::vector<PriINode>::iterator it=priINodeList.begin(); it!=priINodeList.end(); ++it) {
		LatticeID = fNodeList[it->fNodeID].LatticeID;
		sum_qin += it->qin;
		double s;
		if (!std::binary_search(initCrackList.begin(),initCrackList.end(),LatticeID)) {
		    s = std::max(lforce.getCrackOpening(LatticeID),0.0)*LatTab[LatticeID].area;
		} else {
		    s = 0.0;
		}
		sum_s += -s/time;
	}
	std::vector<double>		apertureList;
	dualOut << "----------print q vector--------" << std::endl;
	double avgAperture = 0.0;
	for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
	    LatticeID = fNodeList[DOFtoFNode[DOF]].LatticeID;
	    double aperture = std::max(lforce.getCrackOpening(LatticeID),0.0);
	    apertureList.push_back(aperture);
	    avgAperture += aperture;
	}

	avgAperture /= calNodeNum;
	for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
		LatticeID = fNodeList[DOFtoFNode[DOF]].LatticeID;
		double s=0.0,leak=0.0;
		    if (fNodeList[DOFtoFNode[DOF]].old_p>Tiny) {
		            s = apertureList[DOF]*LatTab[LatticeID].area;
		            leak = getDeltaLeakOff(DOFtoFNode[DOF]);
		    }
		    if (GeometryID ==4) {
		        unsigned i = std::lower_bound(notchList.begin(),notchList.end(),LatticeID)-notchList.begin();
		        if (notchList[i]==LatticeID) {
		            sum_s -= -notchList[i].i2/time;
		        }
		    }
		q[DOF] = -(s+leak)/time;
		l[DOF] = -leak/time;
		sum_s += -s/time;
		sum_l += -leak/time;
	}

	double bHstorage = -getBHstorage();
	double newDeltaLeakOff = -getDeltaLeakOff();
	double leakOff = -getTotalLeakOff();
	double rate_leakOff = (newDeltaLeakOff+leakOff)/time;
	double rate_BHstorage = bHstorage/time;
	for (unsigned i=0; i< qinList.size(); i++) {
		sum_qin +=iNodeList[qinList[i].first].qin;
	}
	double deltaQ	= sum_s+rate_BHstorage+sum_qin+sum_l;
	oldDeltaQ = deltaQ;
	std::array<bool,2>    qStatus;
	qStatus[1] = ((std::fabs(deltaQ)<1e-3)&&(!isInjPChange));
	static double oldSum_s = 0.0;
	static double oldSum_l = 0.0;
//	oldSum_s = sum_s;
//	oldSum_l = sum_l;
	for (unsigned i=0; i< qinList.size(); i++) {
		q[fNodeToDOF[qinList[i].second]] += iNodeList[qinList[i].first].qin;
		sum_qin +=iNodeList[qinList[i].first].qin;
	}
	dualOut << "sum_s, sum_qin, rate_leakOff, rate_BHstorage = "
	        << sum_s << "," << sum_qin << ' ' << rate_leakOff << ' ' << rate_BHstorage << '\n';
	if (sum_qin+rate_leakOff<0.0) {
	    dualOut << "Leak off is too large, increase time by one timeStep" << std::endl;
	    time += timeStep;
	    qStatus[0] = true;
	    return qStatus;
	}
	normRes = (sum_s+rate_BHstorage+sum_qin+sum_l)/
	        (std::fabs(sum_qin)+std::fabs(sum_l)+std::abs(sum_s)+std::abs(rate_BHstorage));
	dualOut << "normRes = (sum_s+sum_qin+rate_leakOff)/(|sum_qin|+|rate_leakOff|+|sum_s|) = " << normRes << '\n';
	fluidEfficieny = -(sum_s+rate_BHstorage)/(sum_qin);
	dualOut << "Fluid efficiency = " << fluidEfficieny << std::endl;
	oldBHstorage = bHstorage;
	if (std::fabs(normRes)<q_tolerance) {
		qStatus[0] = true;
		return qStatus;
	} else {
		qStatus[0] = false;
		return qStatus;
	}
}

void FluidCal::updateRes                    ()
{
    std::vector<IDpair>      qinList = initialize();
    dualOut << " Updating normRes..." << std::endl;
    if (priINodeList.empty()) {
        tripOut << "---WARNING---No injection node!!!" << '\n';
    }
    unsigned LatticeID;
    double sum_qin = 0.0;
    double sum_s = 0.0;
    double sum_l = 0.0;

    GLatForce lforce;
    for (std::vector<PriINode>::iterator it=priINodeList.begin(); it!=priINodeList.end(); ++it) {
        LatticeID = fNodeList[it->fNodeID].LatticeID;
        sum_qin += it->qin;
        double s;
        if (!std::binary_search(initCrackList.begin(),initCrackList.end(),LatticeID)) {
            s = std::max(lforce.getCrackOpening(LatticeID),0.0)*LatTab[LatticeID].area;
        } else {
            s = 0.0;
        }
        sum_s += -s/time;
    }
    std::vector<double>     apertureList;

    double avgAperture = 0.0;
    for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
        LatticeID = fNodeList[DOFtoFNode[DOF]].LatticeID;
        double aperture = std::max(lforce.getCrackOpening(LatticeID),0.0);
        apertureList.push_back(aperture);
        avgAperture += aperture;
    }
    avgAperture /= calNodeNum;
    double   p0 = this->getPriInjP();
    static unsigned step = 0;
    char fileName[255];
    sprintf(fileName,"%s/%spProf%04d.vtk",OutputSubFolder,OutputFilePrefix,step);
    std::ofstream file(fileName, std::ios::out | std::ios::app);
    for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
        LatticeID = fNodeList[DOFtoFNode[DOF]].LatticeID;
        double s=0.0,leak=0.0;
        if (!std::binary_search(initCrackList.begin(),initCrackList.end(),LatticeID)) {
            if (fNodeList[DOFtoFNode[DOF]].old_p/p0>-p_tolerance) {
            if (apertureList[DOF]>50*avgAperture) {
                s = 0.0;
                leak = 0.0;
                file << "1" << ' ';
            } else {
                s = apertureList[DOF]*LatTab[LatticeID].area;
                leak = getDeltaLeakOff(DOFtoFNode[DOF]);
                file << "2" << ' ';
            }
        } else {
            s = 0.0;
            leak = 0.0;
            file << "3" << ' ';
        }
        file << fNodeList[DOFtoFNode[DOF]].old_p/p0 << ' ' << s << ' ' << apertureList[DOF] << std::endl;
            if (GeometryID ==4) {
                unsigned i = std::lower_bound(notchList.begin(),notchList.end(),LatticeID)-notchList.begin();
                if (notchList[i]==LatticeID) {
                    sum_s -= -notchList[i].i2/time;
                }
            }
        }
        sum_s += -s/time;
        sum_l += -leak/time;
    }
    double bHstorage = -getBHstorage();
    double newDeltaLeakOff = -getDeltaLeakOff();
    double leakOff = -getTotalLeakOff();
    double rate_leakOff = (newDeltaLeakOff+leakOff)/time;
    double rate_BHstorage = bHstorage/time;
    for (unsigned i=0; i< qinList.size(); i++) {
        sum_qin +=iNodeList[qinList[i].first].qin;
    }
    double deltaQ = Huge;
    static double       oldSum_s = 0.0;
    static double       oldSum_l = 0.0;
    deltaQ      = (sum_s - oldSum_s + (sum_l-oldSum_l)+(bHstorage-oldBHstorage)/time) /(sum_s+sum_l);
    if ((oldSum_s+oldSum_l)<Tiny) {
        deltaQ = Huge;
    }

    oldSum_s = sum_s;
    oldSum_l = sum_l;
    for (unsigned i=0; i< qinList.size(); i++) {
        q[fNodeToDOF[qinList[i].second]] += iNodeList[qinList[i].first].qin;
        sum_qin +=iNodeList[qinList[i].first].qin;
    }
    dualOut << "sum_s, sum_qin, rate_leakOff, rate_BHstorage = "
            << sum_s << "," << sum_qin << ' ' << rate_leakOff << ' ' << rate_BHstorage << '\n';
    normRes = (sum_s+rate_BHstorage+sum_qin+rate_leakOff)/
                (std::fabs(sum_qin)+std::fabs(rate_leakOff)+std::abs(sum_s)+std::abs(rate_BHstorage));
    dualOut << "normRes = (sum_s+sum_qin+rate_leakOff)/(|sum_qin|+|rate_leakOff|+|sum_s|) = " << normRes << '\n';
    fluidEfficieny = -(sum_s+rate_BHstorage)/(sum_qin);
    dualOut << "Fluid efficiency = " << fluidEfficieny << std::endl;
//    oldDeltaLeakOff = newDeltaLeakOff;
    oldBHstorage = bHstorage;
}

double FluidCal::getTotalLeakOff ()
{
    double totalLeakOff = 0.0;
    for (unsigned fNodeID=0; fNodeID<fNodeList.size(); fNodeID++) {
        totalLeakOff += fNodeList[fNodeID].leakVol;
    }
    return totalLeakOff;
}


double FluidCal::getDeltaLeakOff () {
	double delta = 0.0;
	for (unsigned i=0; i<wetFNodeList.size(); i++) {
	    unsigned fNodeID = wetFNodeList[i];
		delta += 2.0*fNodeList[fNodeID].Cl*LatTab[fNodeList[fNodeID].LatticeID].area
                *(std::sqrt(fNodeList[fNodeID].t_exp+timeStep)-std::sqrt(fNodeList[fNodeID].t_exp));
	}
	return delta;
}

double FluidCal::getDeltaLeakOff (   const unsigned      fNodeID) {
    return 2.0*fNodeList[fNodeID].Cl*LatTab[fNodeList[fNodeID].LatticeID].area
                *(std::sqrt(fNodeList[fNodeID].t_exp+timeStep)-std::sqrt(fNodeList[fNodeID].t_exp));
}

double FluidCal::getLeakOff (   const unsigned      fNodeID)
{
    return fNodeList[fNodeID].leakVol+2.0*fNodeList[fNodeID].Cl*LatTab[fNodeList[fNodeID].LatticeID].area
            *(std::sqrt(fNodeList[fNodeID].t_exp+timeStep));
}

void FluidCal::updateLeakOff ()
{
    for (unsigned i=0; i<wetFNodeList.size(); i++) {
        unsigned fNodeID = wetFNodeList[i];
        double currentLeakOff = 2.0*fNodeList[fNodeID].Cl*LatTab[fNodeList[fNodeID].LatticeID].area
                *(std::sqrt(fNodeList[fNodeID].t_exp+timeStep));
        fNodeList[fNodeID].leakVol = currentLeakOff;
        fNodeList[fNodeID].t_exp += timeStep;
    }
    for (std::vector<PriINode>::iterator it=priINodeList.begin(); it!=priINodeList.end(); ++it) {
        double currentLeakOff = 2.0*fNodeList[it->fNodeID].Cl*LatTab[fNodeList[it->fNodeID].LatticeID].area
                        *(std::sqrt(fNodeList[it->fNodeID].t_exp+timeStep));
        fNodeList[it->fNodeID].leakVol = currentLeakOff;
        fNodeList[it->fNodeID].t_exp += timeStep;
    }
}

void FluidCal::enforceGlobalContinuity 	(	const double				residual,
											const double				sum_s,
											const std::vector<double>*	p_storage)
{
	for (unsigned DOF = 0; DOF < calNodeNum; DOF++) {
		q[DOF] -= residual*((*p_storage)[DOF])/sum_s;
	}
}

void FluidCal::fillP ()
{
	p.resize(calNodeNum);
	p.fill(0.0);
	for (unsigned DOF = 0; DOF < calNodeNum; DOF++) {
		p[DOF] = fNodeList[DOFtoFNode[DOF]].old_p;
	}
}

void FluidCal::setP0 (double P)
{
    for (unsigned fNodeID=0; fNodeID < fNodeList.size(); fNodeID++) {
        fNodeList[fNodeID].new_p = fNodeList[fNodeID].old_p = P;
    }
}

void FluidCal::updateOldp () {
    for (unsigned i=0; i<priINodeList.size(); i++) {
        fNodeList[priINodeList[i].fNodeID].new_p = this->getPriInjP();
    }
    for (unsigned DOF = 0; DOF < calNodeNum; DOF++) {
        fNodeList[DOFtoFNode[DOF]].old_p = p[DOF];
        wetFNodeList.clear();
        if(p[DOF]>0.0) {
            wetFNodeList.push_back(DOFtoFNode[DOF]);
        }
        if (DOF<10) {
            std::cout << p[DOF] << ' ';
        }
    }
	std::cout << std::endl;
}

void FluidCal::updateOldP_static (const double        pNew)
{
    for (unsigned fNodeID=0; fNodeID < fNodeList.size(); fNodeID++) {
        fNodeList[fNodeID].new_p = pNew;
    }
    for (unsigned i=0; i<priINodeList.size(); i++) {
        fNodeList[priINodeList[i].fNodeID].new_p = pNew;
    }
    for (unsigned i=0; i<priINodeList.size(); i++) {
        priINodeList[i].p = pNew;
    }
    std::cout << std::endl;
}

bool FluidCal::interpolate_pi ()
{
    static double old_global_alpha = 0.5;
    int cnt = 0;
    dualOut << "Interpolating new pressure profile..." << std::endl;
    double min_gamma = 0.1;
    unsigned minCnt = (std::sqrt(fNodeList.size()));
    for (unsigned fNodeID = 0; fNodeID < fNodeList.size(); fNodeID++) {
        fNodeList[fNodeID].current_p = 0.0;
    }
    unsigned activeNodeNum = 0;
    for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
        unsigned fNodeID = DOFtoFNode[DOF];
        fNodeList[fNodeID].current_p = std::max(p[DOF],0.0);
        if (fNodeList[fNodeID].current_p>Tiny) {
            activeNodeNum++;
        }
    }
    dualOut << " activeNodeNum = " << activeNodeNum << std::endl;
    dualOut << std::endl;
    for (std::vector<PriINode>::iterator it=priINodeList.begin(); it!=priINodeList.end(); ++it) {
        fNodeList[it->fNodeID].new_p = fNodeList[it->fNodeID].current_p
                = fNodeList[it->fNodeID].old_p = getPriInjP();
    }
    double global_alpha = 0.0;
    double p0 = this->getPriInjP();
    for (unsigned fNodeID=0; fNodeID<fNodeList.size(); fNodeID++) {
        double deltaP = std::fabs(fNodeList[fNodeID].current_p - fNodeList[fNodeID].old_p)/p0;
            global_alpha += deltaP;
    }
    if (cnt<minCnt) {
        cnt = minCnt;
    }
    global_alpha /= (cnt+1);
    sensitiveFactor = baseSensitiveFactor;
    if (global_alpha<100*p_tolerance) {
        if (sensitiveFactor>=1.0) {
            sensitiveFactor = std::sqrt(baseSensitiveFactor);
        } else {
            sensitiveFactor *= sensitiveFactor;
        }
        if (global_alpha<10*p_tolerance) {
            sensitiveFactor = 0.0;
        }
    }
    dualOut << " global alpha = " << global_alpha ;
    dualOut     << " old global alpha " << old_global_alpha << std::endl
            << " sensitive factor = " << sensitiveFactor << std::endl;
    dualOut << " display pressure profile" << std::endl;
    dualOut << "current_p , old_p , alpha , beta" << std::endl;
    for (unsigned fNodeID = 0; fNodeID < fNodeList.size(); fNodeID++) {
        double deltaP = std::fabs(fNodeList[fNodeID].current_p - fNodeList[fNodeID].old_p)/
                (fNodeList[fNodeID].current_p + fNodeList[fNodeID].old_p);
        double alpha = 0;
        double beta = 1;
        if (deltaP < 10*p_tolerance) {
            fNodeList[fNodeID].new_p = fNodeList[fNodeID].current_p;
        } else {
            double alpha = deltaP;
            double beta = (1-min_gamma)*std::pow((1-alpha),sensitiveFactor)+min_gamma;
          fNodeList[fNodeID].new_p = std::max(fNodeList[fNodeID].current_p*beta + (1-beta)*fNodeList[fNodeID].old_p,0.0);
          if (false) {
              dualOut << fNodeList[fNodeID].new_p << ' ' << fNodeList[fNodeID].current_p << ' ' << fNodeList[fNodeID].old_p << ' '
                      << alpha << ' ' << beta << std::endl;
          }
       }
       dualOut << fNodeList[fNodeID].current_p << ' ' << fNodeList[fNodeID].old_p << ' ' << alpha << ' ' << beta << std::endl;
    }
    for (unsigned fNodeID = 0; fNodeID < fNodeList.size(); fNodeID++) {
        fNodeList[fNodeID].old_p = fNodeList[fNodeID].new_p;
    }
    old_global_alpha = global_alpha;
    bool isConverge = (global_alpha < p_tolerance);

    if (isConverge) {
        dualOut << "Global converges " << std::endl;
    }
    if ((cnt<minCnt)&&(!isConverge)) {
        dualOut << "Local converges " << std::endl;
        for (unsigned fNodeID = 0; fNodeID < fNodeList.size(); fNodeID++) {
            fNodeList[fNodeID].old_p = fNodeList[fNodeID].new_p;
            isConverge = true;
        }
    }
    return (isConverge);
}

bool FluidCal::calNewp_static()
{
    static double old_alpha = 0.5;
    double current_alpha;
    int cnt = 0;
    for (unsigned fNodeID = 0; fNodeID < fNodeList.size(); fNodeID++) {
        fNodeList[fNodeID].new_p = 0.0;
    }
    for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
        unsigned fNodeID = DOFtoFNode[DOF];
        fNodeList[fNodeID].new_p = std::max(p[DOF],0.0);
    }
    return true;
}

/*
void FluidCal::updateUniformP 	(	const double		pNew)
{
	dualOut << "FluidCal::updateUniformP... p = " << pNew << std::endl;
	for (unsigned i=0; i<priINodeList.size(); i++) {
		priINodeList[i].p = pNew;
	}

    for (std::vector<PriINode>::iterator it=priINodeList.begin(); it!=priINodeList.end(); ++it) {
        fNodeList[it->fNodeID].new_p = fNodeList[it->fNodeID].old_p = pNew;
    }
    for (unsigned fNodeID =0; fNodeID<fNodeList.size(); fNodeID++) {
        fNodeList[fNodeID].new_p = fNodeList[fNodeID].old_p = pNew;
    }
}
*/
std::vector<IDpair> FluidCal::initialize	()
{
	DOFtoFNode.clear();
	dualFNode.clear();
	std::vector<unsigned>	nbPriINodeList = getPriINodeNbList();
	nbPriNum = nbPriINodeList.size();
	unsigned fNodeID;
	for (unsigned i=0; i<nbPriNum; i++) {
		fNodeID = nbPriINodeList[i];
		if (conFNode.find(fNodeID)!=conFNode.end()) {
			DOFtoFNode.push_back(fNodeID);
		}
	}
	DOFseg1 = DOFtoFNode.size();
	std::vector<unsigned>	chkList = DOFtoFNode;
	for (unsigned i=0; i<priINodeList.size(); i++) {
		chkList.push_back(priINodeList[i].fNodeID);
	}
	std::vector<IDpair>	qinList;
	for (unsigned i=0; i<iNodeList.size(); i++) {
		fNodeID = iNodeList[i].fNodeID;
		if (conFNode.find(fNodeID)!=conFNode.end()) {
			IDpair	iNodeAndfNodeID (i,fNodeID);
			qinList.push_back(iNodeAndfNodeID);
			if (std::find(chkList.begin(),chkList.end(),fNodeID)==chkList.end()) {
				DOFtoFNode.push_back(fNodeID);
			} else {
				dualFNode.push_back(fNodeID);
			}
		}
	}
	DOFseg2 = DOFtoFNode.size();
	chkList.insert(chkList.end(),DOFtoFNode.begin()+DOFseg1,DOFtoFNode.end());
	for (std::set<unsigned>::iterator it = conFNode.begin(); it != conFNode.end(); ++it) {
		fNodeID = *it;
			if(std::find(chkList.begin(),chkList.end(),*it)==chkList.end()) {
			DOFtoFNode.push_back(fNodeID);
		}
	}
	calNodeNum = DOFtoFNode.size();
	fNodeToDOF.resize(fNodeList.size(),calNodeNum);
	for (unsigned DOF=0; DOF<DOFtoFNode.size();DOF++) {
		fNodeToDOF[DOFtoFNode[DOF]]=DOF;
	}
	return qinList;
}

void FluidCal::interpolatingP0					()
{
	dualOut << "Updating primary injection pressure ..." << '\n';
	double pNew;
	if (!isBounded) {
	    std::pair<bool,double>  out;
	    out = findPlimits ();
	    isBounded = out.first;
	    pNew = out.second;
	    dualOut << " Bounds updated " << "isBounded " << isBounded << std::endl;
	} else {
	    double pNew = getP0SecondOrderInterpolate(this->getPriInjP(),normRes);
//	    p0 = pOutput
//	    oldInjP = pNew;
	    dualOut << " Interpolation complete" << std::endl;
	}
	this->updatePriInjP(pNew);
}

void FluidCal::updatingPressureProfile (    bool       isHydrostatic)
{
    double p0 = this->getPriInjP();
    for (unsigned fNodeID = 0; fNodeID < fNodeList.size(); fNodeID++) {
        fNodeList[fNodeID].new_p = p0;
    }
    if (isHydrostatic) {
        return;
    } else {
        for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
            fNodeList[DOFtoFNode[DOF]].new_p = p[DOF];
        }
    }
}

double FluidCal::getP0SecondOrderInterpolate		        (	const double	pNow,
											                    const double	rNow)
{
    double old_p1 = iterHist[0].p;
    double old_p2 = iterHist[1].p;
    double old_r1 = iterHist[0].r;
    double old_r2 = iterHist[1].r;

    dualOut << "-----Befor Cal------"<< std::endl;
    dualOut << "oldp_1    oldp_r2     oldr1     old_r2" << std::endl;
    dualOut << old_p1 << ' ' << old_p2 << ' ' << old_r1 << ' ' << old_r2 << std::endl;

	unsigned method;
	bool isFail = false;
	dualOut << " Interpolating new injection pressure, r(n) = " << rNow << " r(n-1) " << old_p1 << " r(n-2)"<< std::endl;
	double pNew = 0.0;
	bool is_rNow_large = std::fabs(rNow)>1.0;
	if (!is_rNow_large) {
        if (!isFirstIter) {
            double p1p2 = old_p1 - old_p2;
            double p0p2 = pNow - old_p2;
            double p0p1 = pNow - old_p1;
            double d1 	= (rNow*(2*p0p1+p1p2)*(p1p2)-old_r1*(p0p2*p0p2)+old_r2*(p0p1*p0p1))/(p0p1*p1p2*p0p2);
            double d2	= 2.0*old_r2/(p1p2*p0p2)-2.0*old_r1/(p0p1*p1p2)+2.0*rNow/(p0p1*p0p2);
            old_p2 = old_p1;
            old_p1 = pNow;
            old_r2 = old_r1;
            old_r1 = rNow;
            dualOut << "High order predict, " << " d1 = " << d1 << " d2 = " << d2 << std::endl;
            isFail = ((d1!=d1)||(d2!=d2));
            if (p0p1<Tiny) {
                dualOut << " p0p1 < Tiny, pNow = " << p0p1;
                dualOut << " p0, p1, p2 = " << pNow << ' ' << old_p1 << ' ' << old_p2 << std::endl;
            }
            if (p0p2<Tiny) {
                dualOut << " p0p2 < Tiny, old_p1 = " << p0p2;
                dualOut << " p0, p1, p2 = " << pNow << ' ' << old_p1 << ' ' << old_p2 << std::endl;
            }
            if (p1p2<Tiny) {
                dualOut << " p1p2 < Tiny, p1p2 = " << old_p2;
                dualOut << " p0, p1, p2 = " << pNow << ' ' << old_p1 << ' ' << old_p2 << std::endl;
            }
            pNew = pNow - rNow/d1*(1+rNow*d2/(2.0*d1*d1));
            dualOut << "High order predict pNew = " << pNew << " d1 = " << d1 << " d2 = " << d2 << std::endl;
            method = 3;
        }

        if (((pNew>pMax)||(pNew<pMin))||(!isFail)) {
            tripOut << "High order predict pNew = " << pNew << " is out of range, "
                    << "try secant method, pMin = " << pMin << " pMax = " << pMax;
            pNew = linearPredict(pNow,rNow);
            dualOut << "Now, pNew = " << pNew << '\n';
            method = 2;
        }
	}
	if (is_rNow_large) {
	    dualOut << " NormRes = " << rNow << " is too large, use method of bisection " << std::endl;
	}
	if ((pNew>pMax)||(pNew<pMin)) {
	    dualOut << "secant method pNew = " << pNew << " is out of range, "
	            << "type method of bisection" << '\n';
	    if (rNow>0.0) {
	        pNew = (pNow+pMax)/2.0;
	    } else {
	        pNew = (pNow+pMin)/2.0;
	    }
	    dualOut << "Now, pNew = " << pNew << '\n';
	    method = 1;
	}
	if (rNow>0.0) {
	    pMin = pNow;
	} else {
	    pMax = pNow;
	}
	char fileName[255];
	dualOut << "FluidCal::logConvgHist ... " << std::endl;
	std::sprintf(fileName,"%s/%sconvHist%04d.txt",OutputSubFolder,OutputFilePrefix,totalTimeStep);
	std::ofstream file(fileName, std::ios::out | std::ios::app);
	if (isFirstIter) {
	    file << "Time ="  <<'\t' << time << '\n'
	         << "pNew"  << '\t' << "res" << '\t' << " pMax" << '\t' << " pMin" << '\t' << "method" << std::endl;
	}
	file << pNew  << '\t' << rNow << '\t' << pMax << '\t' << pMin << '\t' << method << std::endl;

	dualOut << "-----After Cal------"<< std::endl;
	dualOut << "oldp_1    oldp_r2     oldr1     old_r2" << std::endl;
	dualOut << old_p1 << ' ' << old_p2 << ' ' << old_r1 << ' ' << old_r2 << std::endl;

    iterHist[0].p = iterHist[1].p;
    iterHist[0].r = iterHist[1].r;
    iterHist[1].p = pNow;
    iterHist[1].r = rNow;

	isFirstIter = false;
	return pNew;
}

double FluidCal::linearPredict			(	const double	pNow,
											const double	rNow)
{
	static double old_p1 = 0.0;
	static double old_r1 = 1.0;
	double d1 	= (rNow - old_r1)/(pNow - old_p1);
	if ((d1>Huge)||(d1!=d1)) {
	    return Huge;
	}
	old_p1 = pNow;
	old_r1 = rNow;
	return pNow - rNow/d1;
}

/*
double FluidCal::linearPredict			()
{
	double secant = (iterHist[1].r - iterHist[0].r)/(iterHist[1].p - iterHist[0].p);
	return iterHist[1].p - iterHist[1].r/secant;
}*/

void FluidCal::logRes			(	const unsigned		step,
									const double		nowP,
									const double		newP,
									const double		normRes)
{
	char fileName[255];
	sprintf(fileName,"%s/%sfluidRes%04d.vtk",OutputSubFolder,OutputFilePrefix,step);
	std::ofstream file(fileName, std::ios::out | std::ios::app);
	file << step << ' ' << nowP << ' ' << newP << ' ' << normRes << '\n';
}

template <typename T>
void FluidCal::log				(	std::ostream*	p_file,
										const T			data)
{
	*p_file << ' ' << data;
}

bool FluidCal::serialSolve (	const unsigned	minIter,
                                const unsigned	maxIter,
                                const double 	precision)
{
    Clock	sClock("PCG Serial");
    Eigen::VectorXd r = q-H*p;
    Eigen::VectorXd m = getPreCondJacobi();
    Eigen::VectorXd d = m.cwiseProduct(r);
    double newDelta;
    double init_delta = newDelta = d.dot(r);
    dualOut << "gIteration starts... precision*init_delta =" << std::max(precision*init_delta,precision) <<'\n';
    double 	alpha,beta;
    double	oldDelta;
    Eigen::VectorXd qq(q.size()), s(q.size());
    unsigned i = 0;
    while (((newDelta > fabs(precision*precision*init_delta))||(i<minIter))&&(i<maxIter))
    {
        qq = H*d;
        alpha = newDelta/d.dot(qq);
        p += d*alpha;
        r -= qq*alpha;
        s = m.cwiseProduct(r);
        oldDelta = newDelta;
        newDelta = r.dot(s);
        beta = newDelta / oldDelta;
        d = s + beta*d;
        i++;
    }
    sClock.get();
    double convergence = newDelta / precision/init_delta;
    dualOut << "Iteration cycle = " << i << " new_delta / precision*init_delta = " << convergence << '\n';
    if (convergence > 50.0) {
    	errOut << "Sereve divergence, calculation stop!"<< '\n';
        return false;
    }
    return true;
}

bool FluidCal::useEigenCGSolver	(	unsigned		maxIteration,
                                    double			precision)
{
    using namespace Eigen;
    Clock	sClock;
    ConjugateGradient<SpMat,Lower,DiagonalPreconditioner<double> >	ecg;
    ecg.setMaxIterations(maxIteration);
    ecg.setTolerance(precision);
    sClock.start("EigenCGSolver");
    ecg.compute(H);
    p = ecg.solve(q);
    sClock.get();

    dualOut << "Iteration No = " 		<< ecg.iterations() << '\n';
    dualOut << "Estimated error = " 	<< ecg.error() << '\n';

    if (ecg.error()>50*sqrt(precision))
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool FluidCal::sparseCholeskySolve ()
{
	Clock clock("sparseCholeskySolve");
	Eigen::SimplicialCholesky<SpMat>		sCholesky(H);
	p = sCholesky.solve(q);
	clock.get();
	return (sCholesky.info()==Eigen::Success);
}

Eigen::VectorXd FluidCal::getPreCondJacobi	()
{
	Eigen::VectorXd		m;
    m.resize(calNodeNum);
    for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
    	m[DOF] = 1.0/H.coeff(DOF,DOF);
    }
    return m;
}

bool FluidCal::fillH	()
{
	using Eigen::VectorXd;
	std::cout << "Filling fluid matrix..."<< '\n';
	H.resize(0,0);
	H.resize(calNodeNum,calNodeNum);
	H.reserve(VectorXd::Constant(calNodeNum,maxCoordNo));
//	tripOut << H << '\n';
	unsigned fNodeID,nbDOF,nbFNodeID,nbNum;
	double	hij,sum_hij=0.0;
	for (unsigned DOF=0; DOF<DOFseg1; DOF++) {
		fNodeID 	= DOFtoFNode[DOF];
		nbNum 		= fNodeList[fNodeID].nbActive.size();
		for (unsigned k=0; k<nbNum; k++) {
//			tripOut << H << '\n';
			nbFNodeID = fNodeList[fNodeID].nbActive[k].fNodeID;
			if (conFNode.find(fNodeID)!=conFNode.end()) {
				hij = getHij(fNodeList[fNodeID].nbActive[k].pipeID);
				H.coeffRef(DOF,DOF) 		+= hij;
				sum_hij +=hij;
				//To handle the case where the nbFNode is a primary injection node
				std::vector<PriINode>::iterator it = std::find(priINodeList.begin(),priINodeList.end(),nbFNodeID);
				if (it==priINodeList.end()) {
					nbDOF = fNodeToDOF[nbFNodeID];
					if (nbDOF<calNodeNum) {
						H.coeffRef(DOF,nbDOF)		-= hij;
						sum_hij -=hij;
					}
				} else {
					q[DOF] += hij*priINodeList[it-priINodeList.begin()].p;
				}
			}
		}
	}
	for (unsigned DOF=DOFseg1; DOF<calNodeNum; DOF++) {
		fNodeID 	= DOFtoFNode[DOF];
		nbNum 		= fNodeList[fNodeID].nbActive.size();
		for (unsigned k=0; k<nbNum; k++) {
			nbDOF = fNodeToDOF[fNodeList[fNodeID].nbActive[k].fNodeID];
			if (nbDOF<calNodeNum) {
				hij = getHij(fNodeList[fNodeID].nbActive[k].pipeID);
				H.coeffRef(DOF,DOF) 		+= hij;
				H.coeffRef(DOF,nbDOF)		-= hij;
				sum_hij +=hij;
			}
		}
	}
	if (sum_hij<Tiny*Tiny*Tiny) {
		tripOut << "The permeability of pipes are too low, impossible to carry out fluid calculation!" << '\n';
		return false;
	}
	H.makeCompressed();
	std::cout << "Finish filling fluid matrix..."<< '\n';
	return true;
}

double FluidCal::getHij			(	unsigned						pipeID)
{
	GLatForce lforce;
	unsigned LatID1 = fNodeList[pipeList[pipeID].nbFNodeID[0]].LatticeID;
	unsigned LatID2 = fNodeList[pipeList[pipeID].nbFNodeID[1]].LatticeID;
	if (ProblemID==2) {
        if ((std::binary_search(initCrackList.begin(),initCrackList.end(),LatID1))||
                (std::binary_search(initCrackList.begin(),initCrackList.end(),LatID2))) {
            return Huge;
        }
	}
	double apecture1 = lforce.getCrackOpening(LatID1);
	double apecture2 = lforce.getCrackOpening(LatID2);
	return pipeList[pipeID].width/(12.0*viscosity*
	        (	pipeList[pipeID].halfLength[0]/(apecture1*apecture1*apecture1)+
	            pipeList[pipeID].halfLength[1]/(apecture2*apecture2*apecture2)));
	return 0.0;
}

std::vector<IDAndInfo>	FluidCal::getIDAndPressure	()
{
	std::vector<IDAndInfo>		out;
	out.resize(fNodeList.size());
	unsigned fNodeID;
	double p0 = getPriInjP();
	for (fNodeID=0; fNodeID<fNodeList.size(); fNodeID++) {
		out[fNodeID].first 	= fNodeList[fNodeID].LatticeID;
		out[fNodeID].second = 0.0;
	}
	for (unsigned fNodeID = 0; fNodeID < fNodeList.size(); fNodeID++) {
	    out[fNodeID].second = fNodeList[fNodeID].new_p;
	}
	return out;
}


std::vector<unsigned>   FluidCal::getLatID  ()
{
    std::vector<unsigned>   latList;

    for (unsigned fNodeID = 0; fNodeID<fNodeList.size(); fNodeID++) {
        latList.push_back(fNodeList[fNodeID].LatticeID);
    }
    return latList;
}

bool	FluidCal::updateActiveConFNode		()
{
	bool isChanged = false;
	unsigned i = 0;
	static unsigned acc = 0;
	for (unsigned DOF=0; DOF<DOFtoFNode.size(); DOF++) {
		if (p[DOF]<0.0) {
			setInactiveFNode(fNodeList[DOFtoFNode[DOF]].LatticeID);
			isChanged = true;
			i++;
		}
	}
	acc+=i;
	dualOut << "No of conFNode become inactive in this step = " << i << " accumlated inactive conFNode = " << acc << '\n';
	if (i==0) {
		acc = 0;
	} else {
		connectedFNodeSet();
	}
	return isChanged;
}

std::vector<double>	FluidCal::getPressure	(	const double	dummyValue)
{
	std::vector<double>		out;
	out.resize(fNodeList.size(),dummyValue);
	for (unsigned fNodeID = 0; fNodeID < fNodeList.size(); fNodeID++) {
	    out[fNodeID] = fNodeList[fNodeID].new_p;
	}
	return out;
}

std::vector<double> FluidCal::getStorageRate   ()
{
    std::vector<double>     out;
    out.resize(fNodeList.size(),0.0);
    unsigned fNodeID;
    if (!DOFtoFNode.empty()) {
        for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
            fNodeID = DOFtoFNode[DOF];
            out[fNodeID] = q[DOF];
        }
    }
    GLatForce lforce;
    for (unsigned i=0; i<priINodeList.size(); i++) {
        fNodeID = priINodeList[i].fNodeID;
        unsigned LatticeID = fNodeList[fNodeID].LatticeID;
        double currentStorage = std::max(lforce.getCrackOpening(LatticeID),0.0)*LatTab[LatticeID].area;
        out[fNodeID] = currentStorage/time;
    }
    return out;
}

std::vector<double> FluidCal::getLeakOff   ()
{
    std::vector<double>     out;
    out.resize(fNodeList.size(),0.0);
    if (!DOFtoFNode.empty()) {
        for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
            out[DOFtoFNode[DOF]] = l[DOF];
        }
    }
    return out;
}

double FluidCal::getPriInjP()
{
	return fNodeList[priINodeList[0].fNodeID].new_p;
}

void FluidCal::updatePriInjP(   double      pNew)
{
    fNodeList[priINodeList[0].fNodeID].new_p = pNew;
}

std::vector<double>	FluidCal::getFlowRate	()
{
	std::vector<double>		out;
	out.resize(fNodeList.size(),0.0);
	if (!DOFtoFNode.empty()) {
		for (unsigned DOF=0; DOF<calNodeNum; DOF++) {
			out[DOFtoFNode[DOF]] = q[DOF];
		}
	}
	return out;
}

double  FluidCal::getResidual ()
{
	Eigen::VectorXd res = H*p-q;
	return res.dot(res)/p.size();
}

void FluidCal::logPQ				(	const unsigned				step,
										const std::vector<double>	P,
										const std::vector<double>	Q)
{
	char fileName[255];
	sprintf(fileName,"%s/%sfPQ%04d.vtk",OutputSubFolder,OutputFilePrefix,step);
	std::ofstream file(fileName, std::ios::out | std::ios::app);
	for (unsigned i=0; i<p.size(); i++) {
		file << P[i] << ' ' << Q[i] << '\n';
	}
}

void FluidCal::printVectors	(	    const unsigned 		step,
									const bool 			isPrintp,
									const bool 			isPrintq,
									const bool 			isPrintH)
{
	char fileName[255];
	sprintf(fileName,"%s/%sFluidVecs%04d.vtk",OutputSubFolder,OutputFilePrefix,step);
	std::ofstream file(fileName, std::ios::out | std::ios::app);
	if (isPrintp) {
		file << "p = " << '\n' << p << '\n';
	}
	if (isPrintq) {
		file << "q = " << '\n' << q << '\n';
	}
	if (isPrintH) {
		file << "H = " << '\n' << H << '\n';
	}
}
