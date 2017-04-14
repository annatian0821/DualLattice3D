/*
 * PCG_Eigen.cpp
 *
 *  Created on: Feb 15, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "PCG_Eigen.hpp"

using namespace std;

PCG_Eigen::PCG_Eigen() {
    // TODO Auto-generated constructor stub

}

PCG_Eigen::~PCG_Eigen() {
    // TODO Auto-generated destructor stub
}

bool PCG_Eigen::solve	(			const unsigned 					step,
                                    const bool 						isRelax,
                                    const bool						useCompSM)
{
    unsigned i_min = (unsigned) (1.5*max(max(Nx,Ny),Nz));
    unsigned i_max = 15*i_min;
    Clock							lClock;
    GLatForce						tLatForce;
    omp_set_num_threads(ThreadNumUsed);
    lClock.start("CG");
    initialize(isRelax,useCompSM);
    lClock.get();
    bool isConverged;
    if ((!Use_OpenMP)||(ThreadNumUsed<2)) {
    	isConverged = integralSolveSerial(i_min,i_max,PCG_Precision);
    } else {
    	isConverged = integralSolveParallel2(i_min,i_max,PCG_Precision,useCompSM);
    }
    if (useCompSM) {
        updateNodeStateCompact (step,isRelax);
    } else {
        updateNodeStateSimple (step,isRelax);
    }
    updateConstrainForce();
    isUpdateDOF = false;

    lClock.get();
    cnt_PCG++;
    if (isConverged) {
        return true;
    } else {
        dualOut << " Severe divergence, calculation stop!!!" << endl;
        return false;
    }
}


bool PCG_Eigen::integralSolve2 (	const unsigned	minIter,
                                    const unsigned	maxIter,
                                    const double 	precision,
                                    const bool		useCompSM)
{
    unsigned i;

    double	init_delta,newDelta,convergence;
    double 	alpha,beta;
    double	oldDelta;

    Eigen::VectorXd		vs(vr.size()),vd(vr.size()),vq(vr.size()),vm(vr.size());

    unsigned 	chunkSize = vr.size()/chunkNum;
    unsigned	chunkSizeEnd = vr.size() - chunkSize*(chunkNum-1);

    Clock	sClock;

    sClock.start("EigenIntegralSolver");
    i = 0;
    omp_set_num_threads(ThreadNumUsed);
    cout << "Eigen Thread No = " << Eigen::nbThreads() << endl;

    getPreCondJacobi(&vm);

    Clock	tClock;
    tClock.start("Parallel");

#pragma omp parallel for firstprivate (chunkSize,chunkSizeEnd)
   for (unsigned i=0; i<chunkNum; i++)
    {
    	if (i==chunkNum-1)
    	{
    		vd.tail(chunkSizeEnd) = vm.tail(chunkSizeEnd).cwiseProduct(vr.tail(chunkSizeEnd));
    	}
    	else
    	{
    		vd.segment(i*chunkSize,chunkSize) = vm.segment(i*chunkSize,chunkSize).cwiseProduct(vr.segment(i*chunkSize,chunkSize));
    	}
    }
    tClock.get();

    tClock.start("Serial");
    vd = vm.cwiseProduct(vr);
    tClock.get();


    init_delta = newDelta = vd.dot(vr);
    dualOut << "gIteration starts... precision*init_delta =" << precision*init_delta <<endl;

    cout << minIter << " " << maxIter << endl;

    while (((newDelta > fabs(precision*precision*init_delta))||(i<minIter))&&(i<maxIter))
    {
        vq = SM*vd;

        if (!useCompSM)
        {
            setRestrainZero	(&vq);
        }

        alpha = newDelta/vd.dot(vq);

        vx += vd*alpha;
        vr -= vq*alpha;
        vs = vm.cwiseProduct(vr);

        oldDelta = newDelta;
        newDelta = vr.dot(vs);
        beta = newDelta / oldDelta;

        vd = vs + beta*vd;
        i++;
    }
    sClock.get();

    convergence = newDelta / precision/init_delta;

    dualOut << "Iteration cycle = " << i << " new_delta / precision*init_delta = " << convergence << endl;

    if (convergence > 50.0)
    {
        dualOut << "Sereve divergence, calculation stop!"<< endl;
        return false;
    }
    return true;
}

bool PCG_Eigen::integralSolveSerial (		const unsigned	minIter,
                                    		const unsigned	maxIter,
                                    		const double 	precision)
{
    unsigned i;

    double	init_delta,newDelta,convergence;
    double 	alpha,beta;
    double	oldDelta;

    Eigen::VectorXd		vs(vr.size()),
    					vd(vr.size()),
    					vq(vr.size()),
    					vm(vr.size());

    Clock	sClock;

    sClock.start("EigenIntegralSolverSerial");
    i = 0;
    getPreCondJacobi(&vm);
    vd = vm.cwiseProduct(vr);

    init_delta = newDelta = vd.dot(vr); //vd.dot(vr);
    dualOut << "gIteration starts... precision*init_delta =" << precision*init_delta <<endl;

    cout << minIter << " " << maxIter << endl;

    while (((newDelta > fabs(precision*precision*init_delta))||(i<minIter))&&(i<maxIter)) {
    	vq = SM*vd;
        alpha = newDelta/vd.dot(vq); //vd.dot(vq);

        vx += alpha*vd;
        vr -= alpha*vq;
        vs = vm.cwiseProduct(vr);

        oldDelta = newDelta;
        newDelta = vr.dot(vs); //vr.dot(vs);
        beta = newDelta / oldDelta;

        vd = vs + beta*vd;
        i++;
    }
    sClock.get();
    convergence = newDelta/precision/init_delta;
    dualOut << "Iteration cycle = " << i << " new_delta / (precision*init_delta) = " << convergence << endl;

    if (convergence > 50.0)
    {
        dualOut << "Sereve divergence, calculation stop!"<< endl;
        return false;
    }
    return true;
}

bool PCG_Eigen::integralSolveParallel (		const unsigned	minIter,
                                    		const unsigned	maxIter,
                                    		const double 	precision,
                                    		const bool		useCompSM)
{
    unsigned i;

    double	init_delta,newDelta,convergence;
    double 	alpha,beta;
    double	oldDelta;

    Eigen::VectorXd		vs(vr.size()),
    					vd(vr.size()),
    					vq(vr.size()),
    					vm(vr.size());

    char fileName[255];
    sprintf(fileName,"%s/vqCheck.vtk",OutputFolder);
    ofstream file(fileName, ios::out);

    Clock	sClock;

    sClock.start("EigenIntegralSolverParallel");
    i = 0;
    getPreCondJacobiParallel(&vm);
#pragma omp parallel for
   for (unsigned t=0; t<chunkNum; t++) {
    	if (t==chunkNum-1) {
    		vd.tail(chunkSizeEnd) = vm.tail(chunkSizeEnd).cwiseProduct(vr.tail(chunkSizeEnd));
    	} else {
    		vd.segment(t*chunkSize,chunkSize) = vm.segment(t*chunkSize,chunkSize).cwiseProduct(vr.segment(t*chunkSize,chunkSize));
    	}
    }

    init_delta = newDelta = parallelDot(&vd,&vr); //vd.dot(vr);
    dualOut << "gIteration starts... precision*init_delta =" << precision*init_delta <<endl;

    cout << minIter << " " << maxIter << endl;

    while (((newDelta > fabs(precision*precision*init_delta))||(i<minIter))&&(i<maxIter)) {
		#pragma omp parallel for
    	for (unsigned t=0; t<chunkNum; t++) {
    		if (t==chunkNum-1) {
    	        vq.tail(chunkSizeEnd) = bSM[t]*vd;
    	    } else {
    	        vq.segment(t*chunkSize,chunkSize) = bSM[t]*vd;
    	    }
    	}

        if (!useCompSM) {
            setRestrainZero	(&vq);
        }
        setConstrain(&vq);

        alpha = newDelta/parallelDot(&vd,&vq); //vd.dot(vq);

		#pragma omp parallel for
        for (unsigned t=0; t<chunkNum; t++) {
        	if (t==chunkNum-1) {
        		vx.tail(chunkSizeEnd) += vd.tail(chunkSizeEnd)*alpha;
        		vr.tail(chunkSizeEnd) -= vq.tail(chunkSizeEnd)*alpha;
        		vs.tail(chunkSizeEnd) = vm.tail(chunkSizeEnd).cwiseProduct(vr.tail(chunkSizeEnd));
        	} else {
        		vx.segment(t*chunkSize,chunkSize) += vd.segment(t*chunkSize,chunkSize)*alpha;
        		vr.segment(t*chunkSize,chunkSize) -= vq.segment(t*chunkSize,chunkSize)*alpha;
        		vs.segment(t*chunkSize,chunkSize) = vm.segment(t*chunkSize,chunkSize).cwiseProduct(vr.segment(t*chunkSize,chunkSize));
        	}
        }

        oldDelta = newDelta;
        newDelta = parallelDot(&vr,&vs); //vr.dot(vs);
        beta = newDelta / oldDelta;

	#pragma omp parallel for
        for (unsigned t=0; t<chunkNum; t++) {
        	if (t==chunkNum-1) {
        		vd.tail(chunkSizeEnd) = vs.tail(chunkSizeEnd) + beta*vd.tail(chunkSizeEnd);
            } else {
            	vd.segment(t*chunkSize,chunkSize) = vs.segment(t*chunkSize,chunkSize) + beta*vd.segment(t*chunkSize,chunkSize);
            }
        }
        i++;
    }
    sClock.get();
    convergence = newDelta/precision/init_delta;
    dualOut << "Iteration cycle = " << i << " new_delta / precision*init_delta = " << convergence << endl;

    if (convergence > 50.0) {
        dualOut << "Sereve divergence, calculation stop!"<< endl;
        return false;
    }
    return true;
}

bool PCG_Eigen::integralSolveParallel2 (    const unsigned  minIter,
                                            const unsigned  maxIter,
                                            const double    precision,
                                            const bool      useCompSM)
{
    unsigned i;

    double  init_delta,newDelta,convergence;
    double  alpha,beta;
    double  oldDelta;

    Eigen::VectorXd     vs(vr.size()),
                        vd(vr.size()),
                        vq(vr.size()),
                        vm(vr.size());

    char fileName[255];
    sprintf(fileName,"%s/vqCheck.vtk",OutputFolder);
    ofstream file(fileName, ios::out);

    Clock   sClock;

    sClock.start("EigenIntegralSolverParallel");
    i = 0;
    getPreCondJacobiParallel(&vm);
#pragma omp parallel for
   for (unsigned t=0; t<chunkNum; t++) {
       vd.segment(chunkPos[t],chunkNumList[t]) = vm.segment(chunkPos[t],chunkNumList[t]).cwiseProduct(vr.segment(chunkPos[t],chunkNumList[t]));
   }

    init_delta = newDelta = parallelDot2(&vd,&vr); //vd.dot(vr);
    dualOut << "gIteration starts... precision*init_delta =" << precision*init_delta <<endl;

    std::cout << minIter << " " << maxIter << endl;



    while (((newDelta > fabs(precision*precision*init_delta))||(i<minIter))&&(i<maxIter)) {
        #pragma omp parallel for
        for (unsigned t=0; t<chunkNum; t++) {
            vq.segment(chunkPos[t],chunkNumList[t]) = bSM[t]*vd;
        }


        if (!useCompSM) {
            setRestrainZero (&vq);
        }
        setConstrain(&vq);

        alpha = newDelta/parallelDot2(&vd,&vq); //vd.dot(vq);

        #pragma omp parallel for
        for (unsigned t=0; t<chunkNum; t++) {
            vx.segment(chunkPos[t],chunkNumList[t]) += vd.segment(chunkPos[t],chunkNumList[t])*alpha;
            vr.segment(chunkPos[t],chunkNumList[t]) -= vq.segment(chunkPos[t],chunkNumList[t])*alpha;
            vs.segment(chunkPos[t],chunkNumList[t]) = vm.segment(chunkPos[t],chunkNumList[t]).cwiseProduct(vr.segment(chunkPos[t],chunkNumList[t]));
        }

        oldDelta = newDelta;
        newDelta = parallelDot2(&vr,&vs); //vr.dot(vs);
        beta = newDelta / oldDelta;

    #pragma omp parallel for
        for (unsigned t=0; t<chunkNum; t++) {
            vd.segment(chunkPos[t],chunkNumList[t]) = vs.segment(chunkPos[t],chunkNumList[t]) + beta*vd.segment(chunkPos[t],chunkNumList[t]);
        }
        i++;
    }
    sClock.get();
    convergence = newDelta/precision/init_delta;
    dualOut << "Iteration cycle = " << i << " new_delta / precision*init_delta = " << convergence << endl;

    if (convergence > 50.0) {
        dualOut << "Sereve divergence, calculation stop!"<< endl;
        return false;
    }
    return true;
}

double PCG_Eigen::parallelDot			(	Eigen::VectorXd		*p_v1,
											Eigen::VectorXd		*p_v2)
{
	double sum = 0;
	#pragma omp parallel for reduction (+:sum)
    for (unsigned t=0; t<chunkNum; t++)
    {
    	if (t==chunkNum-1)
        {
    		sum += p_v1->tail(chunkSizeEnd).dot(p_v2->tail(chunkSizeEnd));
        }
        else
        {
        	sum += p_v1->segment(t*chunkSize,chunkSize).dot(p_v2->segment(t*chunkSize,chunkSize));
        }
    }

    return sum;
}

double PCG_Eigen::parallelDot2          (   Eigen::VectorXd     *p_v1,
                                            Eigen::VectorXd     *p_v2)
{
    double sum = 0;
    #pragma omp parallel for reduction (+:sum)
    for (unsigned t=0; t<chunkNum; t++) {
        sum += p_v1->segment(chunkPos[t],chunkNumList[t]).dot(p_v2->segment(chunkPos[t],chunkNumList[t]));
    }
    return sum;
}

bool PCG_Eigen::useEigenCGSolver	(	unsigned		maxIteration,
                                        double			precision)
{
    using namespace Eigen;
    Clock	sClock;
    sClock.start("EigenCGSolver");
    ConjugateGradient<SpMat>	ecg;
    ecg.setMaxIterations(maxIteration);
    ecg.setTolerance(precision);
    ecg.compute(SM);
    vx = ecg.solve(vr);
    sClock.get();

    dualOut << "Iteration No = " 		<< ecg.iterations() << endl;
    dualOut << "Estimated error = " 	<< ecg.error() << endl;

    if (ecg.error()>50*sqrt(precision))
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool PCG_Eigen::useBiCGSTAB_Solver	(	unsigned		maxIteration,
                                        double			precision)
{
	using namespace Eigen;
	Clock	sClock;
	sClock.start("EigenBiCGSTAB_Solver");
	BiCGSTAB<SpMat> solver;
	solver.setMaxIterations(maxIteration);
	solver.setTolerance(precision);
	solver.compute(SM);
	vx = solver.solve(vr);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
	return (solver.info()==Eigen::Success);
}

bool PCG_Eigen::sparseCholeskySolve ()
{
	Clock clock("sparseCholeskySolve");
	Eigen::SimplicialCholesky<SpMat>		sCholesky(SM);
	vx = sCholesky.solve(vr);
	clock.get();
	std::cout << "sCholesky.info() = " << sCholesky.info() << '\n';
	std::cout << "(sCholesky.info()==Eigen::Success) = " << (bool (sCholesky.info()==Eigen::Success)) << '\n';
	return (sCholesky.info()==Eigen::Success);
}

void PCG_Eigen::testSolvers	(	const unsigned	minIter,
                                const unsigned	maxIter,
                                const double 	precision,
                                const bool		useCompSM)
{
    Eigen::VectorXd tmp_vx,dx,tmp_vr;

    tmp_vr = vr;
    sparseCholeskySolve();;
    tmp_vx = vx;
    vx = Eigen::VectorXd::Zero(DOF_Num);
    vr = tmp_vr;
    integralSolveSerial(minIter,maxIter,precision);
    dx = tmp_vx - vx;

    cout << " Difference between two solvers = " << dx.dot(dx)/vx.dot(vx)/vx.cols() << endl;
}

void PCG_Eigen::initialize	(	const bool 					isRelax,
                                const bool					useCompSM)
{
    dualOut << "Ginitialize..." << endl;

    local_ttvNodeNum	= tNodeNum +tbNodeNum;
    local_tNodeNum		= tNodeNum +tbNodeNum;

    if (useCompSM)
    {
    	if (isUpdateDOF) {
    		updateDOFCompact();		//From Base
    		if ((!Use_OpenMP)||(ThreadNumUsed<2)) {
    			fillSMSimple_Full();
    		} else {
    			fillSMblock_full3();
    		}
    	}
    	fillConList();
        initDispVec_x_Compact(isRelax);
        initGradVec_r_Compact(isRelax);
    } else {
    	if (isUpdateDOF) {
    		updateDOFSimple();		//From Base
    		if ((!Use_OpenMP)||(ThreadNumUsed<2)) {
    		    fillSMSimple_Full();
    		} else {
    		    fillSMblock_full3();
    		}
    	}
    	fillConList();
        initDispVec_x_Simple(isRelax);
        initGradVec_r_Simple(isRelax);
    }
}

void PCG_Eigen::initGradVec_r_Simple	(	const bool				isRelax)
{
    unsigned i,j,NodeID,DOF;

    vr.resize(calNodeNum*DofPerNode);

    if (isRelax)
    {
        //No external force is applied

        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF) firstprivate (Adj_ResForce) shared (GNodeTab)
        for (i = 0; i < calNodeNum; i++)
        {
            NodeID = calNodeList[i];
            for (j=0; j<DofPerNode; j++)
            {
                DOF = i*DofPerNode+j;
                vr[DOF] = GNodeTab[NodeID].extF[j];
            }
        }
    }
    else
    {
        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF) firstprivate (Adj_ResForce) shared (GNodeTab)
        for (i = 0; i < freeNodeNum; i++)
        {
            NodeID = calNodeList[i];
            for (j=0; j<DofPerNode; j++)
            {
                DOF = i*DofPerNode+j;
                vr[DOF] = GNodeTab[NodeID].extdF[j] ; //- sum;
            }
        }

        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF) shared (GNodeTab)
        for (i = freeNodeNum; i < calNodeNum; i++)
        {
            NodeID = calNodeList[i];
            for (j=0; j<DofPerNode; j++)
            {
                DOF = i*DofPerNode+j;
                vr[DOF] = GNodeTab[NodeID].extdF[j]; //- sum;
            }
        }
    }
}

void PCG_Eigen::initGradVec_r_Compact	(	const bool				isRelax)
{
    unsigned i,j,d,NodeID,DOF;
    vr.resize(DOF_Num);
    if (isRelax)
    {
        //No external force is applied
        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF) firstprivate (Adj_ResForce) shared (GNodeTab)
        for (i = 0; i < freeNodeNum; i++) {
            NodeID = calNodeList[i];
            for (j=0; j<DofPerNode; j++) {
                DOF = i*DofPerNode+j;
                vr[DOF] = GNodeTab[NodeID].extF[j];
            }
        }

        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF,d) shared (GNodeTab)
        for (i = freeNodeNum; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            d=0;
            for (j=0; j<DofPerNode; j++) {
                if (!pDispNodeList[i-freeNodeNum].Dir[j]) {
                    vr[DOF+d] = GNodeTab[NodeID].extF[j]; //- sum;
                    d++;
                }
            }
        }
    } else {
        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF) firstprivate (Adj_ResForce) shared (GNodeTab)
        for (i = 0; i < freeNodeNum; i++)
        {
            NodeID = calNodeList[i];
            for (j=0; j<DofPerNode; j++)
            {
                DOF = i*DofPerNode+j;
                vr[DOF] = GNodeTab[NodeID].extdF[j]; //- sum;
            }
        }

        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF,d) shared (GNodeTab)
        for (i = freeNodeNum; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            d=0;
            for (j=0; j<DofPerNode; j++) {
                if (!pDispNodeList[i-freeNodeNum].Dir[j]) {
                    vr[DOF+d] = GNodeTab[NodeID].extdF[j]; //- sum;
                    d++;
                }
            }
        }
    }
    	if (Use_OpenMP) {
			#pragma omp parallel for
			for (unsigned t=0; t<chunkNum; t++) {
			    vr.segment(chunkPos[t],chunkNumList[t]) -= bSM[t]*vx;
			}
    	} else {
    		vr -= SM*vx;
    	}
    	std::cout << "initDispVec_r_Compact end" << std::endl;
}

void PCG_Eigen::initDispVec_x_Simple	(	const bool 				isRelax)
{
    unsigned i,NodeID,k;
    using Eigen::VectorXd;

    vx.resize(DofPerNode*calNodeNum);

    if (isRelax)
    {
        for (i = 0; i < calNodeNum; i++)
        {
            NodeID = calNodeList[i];
            for (k=0; k<DofPerNode; k++)
            {
                vx[i*DofPerNode+k] = 0.0; //(*p_nodeState)[NodeID].d[k]; //oldDelta_d[NodeID*DofPerNode+k];
            }
        }
    }
    else
    {
//        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,k) shared (GNodeTab)
        for (i = 0; i < calNodeNum; i++)
        {
            NodeID = calNodeList[i];
            for (k=0; k<DofPerNode; k++)
            {
                vx[i*DofPerNode+k] = 0.0; //(*p_nodeState)[NodeID].d[k];
            }
        }
    }
}

void PCG_Eigen::initDispVec_x_Compact	(	const bool 				isRelax)
{
    vx = Eigen::VectorXd::Zero(DOF_Num);
    if ((isUpdateDOF)||(!isRelax)) {
    }
    for (unsigned i=0; i<conList.size(); i++) {
        vx[conList[i].first] = conList[i].second;
    }
    if (false) {
        Eigen::VectorXd vx1 = Eigen::VectorXd::Zero(DOF_Num);
        std::cout << vx1 << std::endl;
        std::cout << "**********************" << std::endl;
        std::cout << conList.size() << std::endl;
        for (unsigned i=0; i<conList.size(); i++) {
            vx1[conList[i].first] = conList[i].second;
            std::cout << conList[i].first << ' ' << conList[i].second << std::endl;
        }
        std::cout << "++++++++++++++++++++++" << std::endl;
        Eigen::VectorXd vx2 = vx1 - vx;
        std::cout << vx2 << std::endl;
    }
    std::cout << "initDispVec_x_Compact end" << std::endl;
}

Eigen::Matrix3d    PCG_Eigen::getSpringLocalSM   (      const unsigned      NodeID,
                                                        const unsigned      dir)
{
    Eigen::Matrix3d            lSM;
    unsigned    LatticeID = GNodeTab[NodeID].nbLatticeID[dir];

    double C  = LatTab[LatticeID].factor*LatTab[LatticeID].k[0];
    Eigen::Vector3d        T(LatTab[LatticeID].axes[0][0],LatTab[LatticeID].axes[0][1],LatTab[LatticeID].axes[0][2]);
    lSM = C*T*T.transpose();
    return lSM;
}

Eigen::Matrix3d    PCG_Eigen::getSpringLocalSM   (      const unsigned      LatticeID)
{
    Eigen::Matrix3d            lSM;
    double C  = LatTab[LatticeID].factor*LatTab[LatticeID].k[0];
    Eigen::Vector3d        T(LatTab[LatticeID].axes[0][0],LatTab[LatticeID].axes[0][1],LatTab[LatticeID].axes[0][2]);
    lSM = C*T*T.transpose();
    return lSM;
}

Eigen::Matrix3d    PCG_Eigen::getShearLocalSM   (      const unsigned      NodeID,
                                                       const unsigned      dir)
{
    Eigen::Matrix3d            lSM;
    unsigned    LatticeID = GNodeTab[NodeID].nbLatticeID[dir];

    Eigen::Matrix3d             T;
    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];
    Eigen::Matrix3d             k;
    k << LatTab[LatticeID].k[0],                    0.0,                    0.0,
                            0.0, LatTab[LatticeID].k[1],                    0.0,
                            0.0,                    0.0, LatTab[LatticeID].k[1];
    lSM = LatTab[LatticeID].factor*T.transpose()*k*T;
    return lSM;
}

Eigen::Matrix3d    PCG_Eigen::getShearLocalSM   (      const unsigned      LatticeID)
{
    Eigen::Matrix3d            lSM;

    Eigen::Matrix3d             T;
    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];
    Eigen::Matrix3d             k;
    k << LatTab[LatticeID].k[0],                    0.0,                    0.0,
                            0.0, LatTab[LatticeID].k[1],                    0.0,
                            0.0,                    0.0, LatTab[LatticeID].k[1];
    lSM = LatTab[LatticeID].factor*T.transpose()*k*T;
    return lSM;
}

Eigen::Matrix3d    PCG_Eigen::getSpringLocalSM   (      const unsigned      LatticeID,
                                                        const double        factor)
{
    Eigen::Matrix3d            lSM;
    double C = (factor-LatTab[LatticeID].factor)*LatTab[LatticeID].k[0];
    Eigen::Vector3d        T(LatTab[LatticeID].axes[0][0],LatTab[LatticeID].axes[0][1],LatTab[LatticeID].axes[0][2]);
    lSM = C*T*T.transpose();
    return lSM;
}

Eigen::Matrix3d    PCG_Eigen::getShearLocalSM   (      const unsigned      LatticeID,
                                                       const double        factor)
{
    Eigen::Matrix3d            lSM;
    double C = (factor-LatTab[LatticeID].factor);
    Eigen::Matrix3d             T;
    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];
    Eigen::Matrix3d             k;
    k << LatTab[LatticeID].k[0],                    0.0,                    0.0,
                            0.0, LatTab[LatticeID].k[1],                    0.0,
                            0.0,                    0.0, LatTab[LatticeID].k[1];
    lSM = C*T.transpose()*k*T;
    return lSM;
}

Eigen::MatrixXd    PCG_Eigen::getFullLocalSM2   (       const unsigned      LatticeID,
                                                        const double        factor)
{
    double C = (factor-LatTab[LatticeID].factor);
    Eigen::MatrixXd T(2*DofPerNode,2*DofPerNode);

    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];

    Eigen::MatrixXd             k(DofPerNode,DofPerNode);
    k << LatTab[LatticeID].k[0],0.0,0.0,0.0,0.0,0.0,
         0.0,LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
         0.0,0.0,LatTab[LatticeID].k[1],0.0,0.0,0.0,
         0.0,0.0,0.0,(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]),0.0,0.0,
         0.0,0.0,0.0,0.0,LatTab[LatticeID].k[2],0.0,
         0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].k[3];

    Eigen::MatrixXd             B(DofPerNode,2*DofPerNode);
    double yc = LatTab[LatticeID].offset[0];
    double zc = LatTab[LatticeID].offset[1];
    double h2 = LatTab[LatticeID].length/2.0;
    double ip = 1.0; //std::sqrt((LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3])/LatTab[LatticeID].k[0]);
    double i1 = 1.0; //std::sqrt(LatTab[LatticeID].k[2]/LatTab[LatticeID].k[0]);
    double i2 = 1.0; //std::sqrt(LatTab[LatticeID].k[3]/LatTab[LatticeID].k[0]);
    B <<  -1.0,  0.0,  0.0,  0.0,  -zc,   yc,  1.0,  0.0,  0.0,  0.0,   zc,  -yc,
           0.0, -1.0,  0.0,   zc,  0.0,  -h2,  0.0,  1.0,  0.0,  -zc,  0.0,  -h2,
           0.0,  0.0, -1.0,  -yc,   h2,  0.0,  0.0,  0.0,  1.0,   yc,  +h2,  0.0,
           0.0,  0.0,  0.0,  -ip,  0.0,  0.0,  0.0,  0.0,  0.0,   ip,  0.0,  0.0,
           0.0,  0.0,  0.0,  0.0,  -i1,  0.0,  0.0,  0.0,  0.0,  0.0,   i1,  0.0,
           0.0,  0.0,  0.0,  0.0,  0.0,  -i2,  0.0,  0.0,  0.0,  0.0,  0.0,   i2;
    return C*T.transpose()*B.transpose()*k*B*T;
}

Eigen::MatrixXd    PCG_Eigen::getFullLocalSM2   (       const unsigned      LatticeID,
                                                        const double        factor,
                                                        const double        factorS)
{
    double C1 = factor;
    double S1 = factorS;
    double C = LatTab[LatticeID].factor;
    double S = LatTab[LatticeID].factorShear;
    Eigen::MatrixXd T(2*DofPerNode,2*DofPerNode);

    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];

    Eigen::MatrixXd             k(DofPerNode,DofPerNode);
    Eigen::MatrixXd             k1(DofPerNode,DofPerNode);
    k1<< C1*LatTab[LatticeID].k[0],0.0,0.0,0.0,0.0,0.0,
         0.0,S1*LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
         0.0,0.0,S1*LatTab[LatticeID].k[1],0.0,0.0,0.0,
         0.0,0.0,0.0,C1*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]),0.0,0.0,
         0.0,0.0,0.0,0.0,C1*LatTab[LatticeID].k[2],0.0,
         0.0,0.0,0.0,0.0,0.0,C1*LatTab[LatticeID].k[3];
    k << C*LatTab[LatticeID].k[0],0.0,0.0,0.0,0.0,0.0,
         0.0,S*LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
         0.0,0.0,S*LatTab[LatticeID].k[1],0.0,0.0,0.0,
         0.0,0.0,0.0,C*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]),0.0,0.0,
         0.0,0.0,0.0,0.0,C*LatTab[LatticeID].k[2],0.0,
         0.0,0.0,0.0,0.0,0.0,C*LatTab[LatticeID].k[3];

    Eigen::MatrixXd             B(DofPerNode,2*DofPerNode);
    double yc = LatTab[LatticeID].offset[0];
    double zc = LatTab[LatticeID].offset[1];
    double h2 = LatTab[LatticeID].length/2.0;
    double ip = 1.0;
    double i1 = 1.0;
    double i2 = 1.0;
    B <<  -1.0,  0.0,  0.0,  0.0,  -zc,   yc,  1.0,  0.0,  0.0,  0.0,   zc,  -yc,
           0.0, -1.0,  0.0,   zc,  0.0,  -h2,  0.0,  1.0,  0.0,  -zc,  0.0,  -h2,
           0.0,  0.0, -1.0,  -yc,   h2,  0.0,  0.0,  0.0,  1.0,   yc,  +h2,  0.0,
           0.0,  0.0,  0.0,  -ip,  0.0,  0.0,  0.0,  0.0,  0.0,   ip,  0.0,  0.0,
           0.0,  0.0,  0.0,  0.0,  -i1,  0.0,  0.0,  0.0,  0.0,  0.0,   i1,  0.0,
           0.0,  0.0,  0.0,  0.0,  0.0,  -i2,  0.0,  0.0,  0.0,  0.0,  0.0,   i2;
    return T.transpose()*B.transpose()*k1*B*T - T.transpose()*B.transpose()*k*B*T;
}

Eigen::MatrixXd    PCG_Eigen::getShearLocalSM2   (      const unsigned      LatticeID,
                                                        const double        factorS)
{
    double S = (factorS-LatTab[LatticeID].factorShear);
    Eigen::MatrixXd T(2*DofPerNode,2*DofPerNode);

    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];

    Eigen::MatrixXd             k(DofPerNode,DofPerNode);
    Eigen::MatrixXd             k1(DofPerNode,DofPerNode);
    k << 0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,S*LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
         0.0,0.0,S*LatTab[LatticeID].k[1],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0;

    Eigen::MatrixXd             B(DofPerNode,2*DofPerNode);
    double yc = LatTab[LatticeID].offset[0];
    double zc = LatTab[LatticeID].offset[1];
    double h2 = LatTab[LatticeID].length/2.0;
    double ip = 1.0; //std::sqrt((LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3])/LatTab[LatticeID].k[0]);
    double i1 = 1.0; //std::sqrt(LatTab[LatticeID].k[2]/LatTab[LatticeID].k[0]);
    double i2 = 1.0; //std::sqrt(LatTab[LatticeID].k[3]/LatTab[LatticeID].k[0]);
    B <<  -1.0,  0.0,  0.0,  0.0,  -zc,   yc,  1.0,  0.0,  0.0,  0.0,   zc,  -yc,
           0.0, -1.0,  0.0,   zc,  0.0,  -h2,  0.0,  1.0,  0.0,  -zc,  0.0,  -h2,
           0.0,  0.0, -1.0,  -yc,   h2,  0.0,  0.0,  0.0,  1.0,   yc,  +h2,  0.0,
           0.0,  0.0,  0.0,  -ip,  0.0,  0.0,  0.0,  0.0,  0.0,   ip,  0.0,  0.0,
           0.0,  0.0,  0.0,  0.0,  -i1,  0.0,  0.0,  0.0,  0.0,  0.0,   i1,  0.0,
           0.0,  0.0,  0.0,  0.0,  0.0,  -i2,  0.0,  0.0,  0.0,  0.0,  0.0,   i2;

    return T.transpose()*B.transpose()*k*B*T;
}

Eigen::MatrixXd    PCG_Eigen::getFullLocalSM_removeShear   (       const unsigned      LatticeID)
{
    Eigen::MatrixXd T(2*DofPerNode,2*DofPerNode);

    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];

    Eigen::MatrixXd             k(DofPerNode,DofPerNode);
    Eigen::MatrixXd             k1(DofPerNode,DofPerNode);
    double C = LatTab[LatticeID].factor;
    double S = LatTab[LatticeID].factorShear;
    double C1 = LatTab[LatticeID].factor;
    double S1 = Tiny;
    k1<< C1*LatTab[LatticeID].k[0],0.0,0.0,0.0,0.0,0.0,
         0.0,S1*LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
         0.0,0.0,S1*LatTab[LatticeID].k[1],0.0,0.0,0.0,
         0.0,0.0,0.0,C1*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]),0.0,0.0,
         0.0,0.0,0.0,0.0,C1*LatTab[LatticeID].k[2],0.0,
         0.0,0.0,0.0,0.0,0.0,C1*LatTab[LatticeID].k[3];
    k << C*LatTab[LatticeID].k[0],0.0,0.0,0.0,0.0,0.0,
         0.0,S*LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
         0.0,0.0,S*LatTab[LatticeID].k[1],0.0,0.0,0.0,
         0.0,0.0,0.0,C*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]),0.0,0.0,
         0.0,0.0,0.0,0.0,C*LatTab[LatticeID].k[2],0.0,
         0.0,0.0,0.0,0.0,0.0,C*LatTab[LatticeID].k[3];

    Eigen::MatrixXd             B(DofPerNode,2*DofPerNode);
    double yc = LatTab[LatticeID].offset[0];
    double zc = LatTab[LatticeID].offset[1];
    double h2 = LatTab[LatticeID].length/2.0;
    double ip = 1.0; //std::sqrt((LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3])/LatTab[LatticeID].k[0]);
    double i1 = 1.0; //std::sqrt(LatTab[LatticeID].k[2]/LatTab[LatticeID].k[0]);
    double i2 = 1.0; //std::sqrt(LatTab[LatticeID].k[3]/LatTab[LatticeID].k[0]);
    B <<  -1.0,  0.0,  0.0,  0.0,  -zc,   yc,  1.0,  0.0,  0.0,  0.0,   zc,  -yc,
           0.0, -1.0,  0.0,   zc,  0.0,  -h2,  0.0,  1.0,  0.0,  -zc,  0.0,  -h2,
           0.0,  0.0, -1.0,  -yc,   h2,  0.0,  0.0,  0.0,  1.0,   yc,  +h2,  0.0,
           0.0,  0.0,  0.0,  -ip,  0.0,  0.0,  0.0,  0.0,  0.0,   ip,  0.0,  0.0,
           0.0,  0.0,  0.0,  0.0,  -i1,  0.0,  0.0,  0.0,  0.0,  0.0,   i1,  0.0,
           0.0,  0.0,  0.0,  0.0,  0.0,  -i2,  0.0,  0.0,  0.0,  0.0,  0.0,   i2;
    return T.transpose()*B.transpose()*k1*B*T - T.transpose()*B.transpose()*k1*B*T;
}

Tensor    PCG_Eigen::getMemForce       (       const unsigned      LatticeID)
{
    unsigned NodeID = LatTab[LatticeID].nb[0];
    unsigned nbNodeID = LatTab[LatticeID].nb[1];
    Eigen::VectorXd vd(2*DofPerNode);
    for (unsigned d=0; d<DofPerNode; d++) {
        vd[d] = GNodeTab[NodeID].d[d];
        vd[d+DofPerNode] = GNodeTab[nbNodeID].d[d];
    }
    Eigen::MatrixXd T(2*DofPerNode,2*DofPerNode);
    Eigen::MatrixXd k(DofPerNode,DofPerNode);
    Eigen::VectorXd vf(2*DofPerNode);
    if (DofPerNode==6) {
		T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			 LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			 LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			 0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,
			 0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,
			 0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,
			 0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
			 0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
			 0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
			 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
			 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
			 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];
		Eigen::MatrixXd k(DofPerNode,DofPerNode);
		k << LatTab[LatticeID].factor*LatTab[LatticeID].k[0],0.0,0.0,0.0,0.0,0.0,
			 0.0,LatTab[LatticeID].factorShear*LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
			 0.0,0.0,LatTab[LatticeID].factorShear*LatTab[LatticeID].k[1],0.0,0.0,0.0,
			 0.0,0.0,0.0,LatTab[LatticeID].factor*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]),0.0,0.0,
			 0.0,0.0,0.0,0.0,LatTab[LatticeID].factor*LatTab[LatticeID].k[2],0.0,
			 0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].factor*LatTab[LatticeID].k[3];

		Eigen::MatrixXd             B(DofPerNode,2*DofPerNode);
		double yc = LatTab[LatticeID].offset[0];
		double zc = LatTab[LatticeID].offset[1];
		double h2 = LatTab[LatticeID].length/2.0;
		double ip = 1.0; // std::sqrt((LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3])/LatTab[LatticeID].k[0]);
		double i1 = 1.0; //std::sqrt(LatTab[LatticeID].k[2]/LatTab[LatticeID].k[0]);
		double i2 = 1.0; //std::sqrt(LatTab[LatticeID].k[3]/LatTab[LatticeID].k[0]);
		B <<  -1.0,  0.0,  0.0,  0.0,  -zc,   yc,  1.0,  0.0,  0.0,  0.0,   zc,  -yc,
			   0.0, -1.0,  0.0,   zc,  0.0,  -h2,  0.0,  1.0,  0.0,  -zc,  0.0,  -h2,
			   0.0,  0.0, -1.0,  -yc,   h2,  0.0,  0.0,  0.0,  1.0,   yc,  +h2,  0.0,
			   0.0,  0.0,  0.0,  -ip,  0.0,  0.0,  0.0,  0.0,  0.0,   ip,  0.0,  0.0,
			   0.0,  0.0,  0.0,  0.0,  -i1,  0.0,  0.0,  0.0,  0.0,  0.0,   i1,  0.0,
			   0.0,  0.0,  0.0,  0.0,  0.0,  -i2,  0.0,  0.0,  0.0,  0.0,  0.0,   i2;
		vf = B.transpose()*k*B*T*vd;
    } else {
    	T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
    		 LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
    		 LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
    		 0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
    		 0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
    		 0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];
    	Eigen::MatrixXd k(2*DofPerNode,2*DofPerNode);
    	k << -LatTab[LatticeID].k[0],0.0,0.0,LatTab[LatticeID].k[0],0.0,0.0,
    	     0.0,0.0,0.0,0.0,0.0,0.0,
    		 0.0,0.0,0.0,0.0,0.0,0.0,
    		 LatTab[LatticeID].k[0],0.0,0.0,-LatTab[LatticeID].k[0],0.0,0.0,
    		 0.0,0.0,0.0,0.0,0.0,0.0,
    		 0.0,0.0,0.0,0.0,0.0,0.0;
    	 vf = LatTab[LatticeID].factor*k*T*vd;
    }
    Tensor memForce;
    for (unsigned d=0; d<DofPerNode; d++) {
        memForce[d] = vf[d];
    }
    return memForce;
}

Eigen::MatrixXd    PCG_Eigen::getFullLocalSM2   (       const unsigned      LatticeID)
{
    Eigen::MatrixXd T(2*DofPerNode,2*DofPerNode);

    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];

    Eigen::MatrixXd             k(DofPerNode,DofPerNode);
    k << LatTab[LatticeID].factor*LatTab[LatticeID].k[0],0.0,0.0,0.0,0.0,0.0,
         0.0,LatTab[LatticeID].factorShear*LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
         0.0,0.0,LatTab[LatticeID].factorShear*LatTab[LatticeID].k[1],0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].factor*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]),0.0,0.0,
         0.0,0.0,0.0,0.0,LatTab[LatticeID].factor*LatTab[LatticeID].k[2],0.0,
         0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].factor*LatTab[LatticeID].k[3];

    Eigen::MatrixXd             B(DofPerNode,2*DofPerNode);
    double yc = LatTab[LatticeID].offset[0];
    double zc = LatTab[LatticeID].offset[1];
    double h2 = LatTab[LatticeID].length/2.0;
    double ip = 1.0; //std::sqrt((LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3])/LatTab[LatticeID].k[0]);
    double i1 = 1.0; //std::sqrt(LatTab[LatticeID].k[2]/LatTab[LatticeID].k[0]);
    double i2 = 1.0; //std::sqrt(LatTab[LatticeID].k[3]/LatTab[LatticeID].k[0]);
    B <<  -1.0,  0.0,  0.0,  0.0,  -zc,   yc,  1.0,  0.0,  0.0,  0.0,   zc,  -yc,
           0.0, -1.0,  0.0,   zc,  0.0,  -h2,  0.0,  1.0,  0.0,  -zc,  0.0,  -h2,
           0.0,  0.0, -1.0,  -yc,   h2,  0.0,  0.0,  0.0,  1.0,   yc,   h2,  0.0,
           0.0,  0.0,  0.0,  -ip,  0.0,  0.0,  0.0,  0.0,  0.0,   ip,  0.0,  0.0,
           0.0,  0.0,  0.0,  0.0,  -i1,  0.0,  0.0,  0.0,  0.0,  0.0,   i1,  0.0,
           0.0,  0.0,  0.0,  0.0,  0.0,  -i2,  0.0,  0.0,  0.0,  0.0,  0.0,   i2;
    return T.transpose()*B.transpose()*k*B*T;
}

Eigen::MatrixXd    PCG_Eigen::getFullLocalSM2   (       const unsigned      NodeID,
                                                        const unsigned      dir)
{
    Eigen::MatrixXd T(DofPerNode,DofPerNode);
    unsigned    LatticeID = GNodeTab[NodeID].nbLatticeID[dir];

    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];
    Eigen::MatrixXd             k(DofPerNode,DofPerNode);
    k << LatTab[LatticeID].factor*LatTab[LatticeID].k[0],0.0,0.0,0.0,0.0,0.0,
         0.0,LatTab[LatticeID].factorShear*LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
         0.0,0.0,LatTab[LatticeID].factorShear*LatTab[LatticeID].k[1],0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].factor*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]),0.0,0.0,
         0.0,0.0,0.0,0.0,LatTab[LatticeID].factor*LatTab[LatticeID].k[2],0.0,
         0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].factor*LatTab[LatticeID].k[3];

    Eigen::MatrixXd             B(DofPerNode,DofPerNode),B1(DofPerNode,DofPerNode);
    double yc = LatTab[LatticeID].offset[0];
    double zc = LatTab[LatticeID].offset[1];
    double h2 = LatTab[LatticeID].length/2.0;
    double ip = 1.0; //std::sqrt((LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3])/LatTab[LatticeID].k[0]);
    double i1 = 1.0; //std::sqrt(LatTab[LatticeID].k[2]/LatTab[LatticeID].k[0]);
    double i2 = 1.0; //std::sqrt(LatTab[LatticeID].k[3]/LatTab[LatticeID].k[0]);
    B <<  -1.0,  0.0,  0.0,  0.0,  -zc,   yc,
           0.0, -1.0,  0.0,   zc,  0.0,  -h2,
           0.0,  0.0, -1.0,  -yc,   h2,  0.0,
           0.0,  0.0,  0.0,  -ip,  0.0,  0.0,
           0.0,  0.0,  0.0,  0.0,  -i1,  0.0,
           0.0,  0.0,  0.0,  0.0,  0.0,  -i2;
    B1 <<  1.0,  0.0,  0.0,  0.0,   zc,  -yc,
           0.0,  1.0,  0.0,  -zc,  0.0,  -h2,
           0.0,  0.0,  1.0,   yc,  +h2,  0.0,
           0.0,  0.0,  0.0,   ip,  0.0,  0.0,
           0.0,  0.0,  0.0,  0.0,   i1,  0.0,
           0.0,  0.0,  0.0,  0.0,  0.0,   i2;

    Eigen::MatrixXd lSM(DofPerNode,2*DofPerNode);
    unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[dir];
    if (nbNodeID>NodeID) {
        lSM.block(0,0,DofPerNode,DofPerNode) = T.transpose()*B.transpose()*k*B*T;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = T.transpose()*B.transpose()*k*B1*T;
    } else {
        lSM.block(0,0,DofPerNode,DofPerNode) = T.transpose()*B1.transpose()*k*B1*T;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = T.transpose()*B1.transpose()*k*B*T;
    }
    return lSM;
}

Eigen::MatrixXd    PCG_Eigen::getFullLocalSM2   (       const unsigned      NodeID,
                                                        const unsigned      dir,
                                                        const double        factor)
{
    Eigen::MatrixXd T(DofPerNode,DofPerNode);
    unsigned    LatticeID = GNodeTab[NodeID].nbLatticeID[dir];
    double C = (factor-LatTab[LatticeID].factor);
    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];
    Eigen::MatrixXd             k(DofPerNode,DofPerNode);
    k << LatTab[LatticeID].factor*LatTab[LatticeID].k[0],0.0,0.0,0.0,0.0,0.0,
         0.0,LatTab[LatticeID].factorShear*LatTab[LatticeID].k[1],0.0,0.0,0.0,0.0,
         0.0,0.0,LatTab[LatticeID].factorShear*LatTab[LatticeID].k[1],0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].factor*(LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]),0.0,0.0,
         0.0,0.0,0.0,0.0,LatTab[LatticeID].factor*LatTab[LatticeID].k[2],0.0,
         0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].factor*LatTab[LatticeID].k[3];

    Eigen::MatrixXd             B(DofPerNode,DofPerNode),B1(DofPerNode,DofPerNode);
    double yc = LatTab[LatticeID].offset[0];
    double zc = LatTab[LatticeID].offset[1];
    double h2 = LatTab[LatticeID].length/2.0;
    double ip = 1.0; //std::sqrt((LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3])/LatTab[LatticeID].k[0]);
    double i1 = 1.0; //std::sqrt(LatTab[LatticeID].k[2]/LatTab[LatticeID].k[0]);
    double i2 = 1.0; //std::sqrt(LatTab[LatticeID].k[3]/LatTab[LatticeID].k[0]);
    B <<  -1.0,  0.0,  0.0,  0.0,  -zc,   yc,
           0.0, -1.0,  0.0,   zc,  0.0,  -h2,
           0.0,  0.0, -1.0,  -yc,   h2,  0.0,
           0.0,  0.0,  0.0,  -ip,  0.0,  0.0,
           0.0,  0.0,  0.0,  0.0,  -i1,  0.0,
           0.0,  0.0,  0.0,  0.0,  0.0,  -i2;
    B1 <<  1.0,  0.0,  0.0,  0.0,   zc,  -yc,
           0.0,  1.0,  0.0,  -zc,  0.0,  -h2,
           0.0,  0.0,  1.0,   yc,  +h2,  0.0,
           0.0,  0.0,  0.0,   ip,  0.0,  0.0,
           0.0,  0.0,  0.0,  0.0,   i1,  0.0,
           0.0,  0.0,  0.0,  0.0,  0.0,   i2;

    Eigen::MatrixXd lSM(DofPerNode,2*DofPerNode);
    unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[dir];
    if (nbNodeID>NodeID) {
        lSM.block(0,0,DofPerNode,DofPerNode) = C*T.transpose()*B.transpose()*k*B*T;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = C*T.transpose()*B.transpose()*k*B1*T;
    } else {
        lSM.block(0,0,DofPerNode,DofPerNode) = C*T.transpose()*B1.transpose()*k*B*T;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = C*T.transpose()*B1.transpose()*k*B1*T;
    }
    return lSM;
}

Eigen::MatrixXd    PCG_Eigen::getBeamSM   		(       const unsigned      NodeID,
                                                        const unsigned      dir)
{
    Eigen::MatrixXd T(DofPerNode,DofPerNode);
    unsigned    LatticeID = GNodeTab[NodeID].nbLatticeID[dir];

    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];

    double EA_L = LatTab[LatticeID].k[0];
    double EIz_L3 = 12.0*LatTab[LatticeID].k[3]/(LatTab[LatticeID].length*LatTab[LatticeID].length);
    double EIy_L3 = 12.0*LatTab[LatticeID].k[2]/(LatTab[LatticeID].length*LatTab[LatticeID].length);
    double EIz_L2 = 6.0*LatTab[LatticeID].k[3]/LatTab[LatticeID].length;
    double EIy_L2 = 6.0*LatTab[LatticeID].k[2]/LatTab[LatticeID].length;
    double EIz_L  = 4.0*LatTab[LatticeID].k[3];
    double EIy_L  = 4.0*LatTab[LatticeID].k[2];
    double GJ_L   = (LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]);
    Eigen::MatrixXd lSM(DofPerNode,2*DofPerNode);
    unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[dir];
    if (nbNodeID>NodeID) {
    	Eigen::MatrixXd   K1(DofPerNode,DofPerNode),K2(DofPerNode,DofPerNode);
        K1 <<  EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,
                0.0,  EIz_L3,     0.0,    0.0,     0.0,  EIz_L2,
                0.0,     0.0,  EIy_L3,    0.0, -EIy_L2,     0.0,
                0.0,     0.0,     0.0,   GJ_L,     0.0,     0.0,
                0.0,     0.0, -EIy_L2,    0.0,   EIy_L,     0.0,
                0.0,  EIz_L2,     0.0,    0.0,     0.0,   EIz_L;

        K2 << -EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,
                0.0, -EIz_L3,     0.0,    0.0,     0.0,  EIz_L2,
                0.0,     0.0, -EIy_L3,    0.0, -EIy_L2,     0.0,
                0.0,     0.0,     0.0,  -GJ_L,     0.0,     0.0,
                0.0,     0.0,  EIy_L2,    0.0, 2*EIy_L,     0.0,
                0.0, -EIz_L2,     0.0,    0.0,     0.0, 2*EIz_L;

        lSM.block(0,0,DofPerNode,DofPerNode) 			= LatTab[LatticeID].factor*T.transpose()*K1*T;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) 	= LatTab[LatticeID].factor*T.transpose()*K2*T;
    } else {
    	Eigen::MatrixXd   K3(DofPerNode,DofPerNode),K4(DofPerNode,DofPerNode);
    	K3 <<  -EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,
    	         0.0, -EIz_L3,     0.0,    0.0,     0.0, -EIz_L2,
    	         0.0,     0.0, -EIy_L3,    0.0,  EIy_L2,     0.0,
    	         0.0,     0.0,     0.0,  -GJ_L,     0.0,     0.0,
    	         0.0,     0.0, -EIy_L2,    0.0, 2*EIy_L,     0.0,
    	         0.0,  EIz_L2,     0.0,    0.0,     0.0, 2*EIz_L;
    	K4 <<   EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,
    	         0.0,  EIz_L3,     0.0,    0.0,     0.0, -EIz_L2,
    	         0.0,     0.0,  EIy_L3,    0.0,  EIy_L2,     0.0,
    	    	 0.0,     0.0,     0.0,   GJ_L,     0.0,     0.0,
    	    	 0.0,     0.0,  EIy_L2,    0.0,   EIy_L,     0.0,
    	    	 0.0, -EIz_L2,     0.0,    0.0,     0.0,   EIz_L;

        lSM.block(0,0,DofPerNode,DofPerNode)          = LatTab[LatticeID].factor*T.transpose()*K4*T;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = LatTab[LatticeID].factor*T.transpose()*K3*T;
    }

    return lSM;
}

Eigen::MatrixXd    PCG_Eigen::getBeamSM    (       const unsigned      LatticeID)
{
    Eigen::MatrixXd T(2*DofPerNode,2*DofPerNode);

    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];
    double EA_L = LatTab[LatticeID].k[0];
    double EIz_L3 = 12.0*LatTab[LatticeID].k[3]/(LatTab[LatticeID].length*LatTab[LatticeID].length);
    double EIy_L3 = 12.0*LatTab[LatticeID].k[2]/(LatTab[LatticeID].length*LatTab[LatticeID].length);
    double EIz_L2 = 6.0*LatTab[LatticeID].k[3]/LatTab[LatticeID].length;
    double EIy_L2 = 6.0*LatTab[LatticeID].k[2]/LatTab[LatticeID].length;
    double EIz_L  = 4.0*LatTab[LatticeID].k[3];
    double EIy_L  = 4.0*LatTab[LatticeID].k[3];
    double GJ_L   = (LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]);
    Eigen::MatrixXd K(2*DofPerNode,2*DofPerNode);
    K  <<  EA_L,     0.0,     0.0,    0.0,     0.0,     0.0, -EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,
            0.0,  EIz_L3,     0.0,    0.0,     0.0,  EIz_L2,   0.0, -EIz_L3,     0.0,    0.0,     0.0,  EIz_L2,
            0.0,     0.0,  EIy_L3,    0.0, -EIy_L2,     0.0,   0.0,     0.0, -EIy_L3,    0.0, -EIy_L2,     0.0,
            0.0,     0.0,     0.0,   GJ_L,     0.0,     0.0,   0.0,     0.0,     0.0,  -GJ_L,     0.0,     0.0,
            0.0,     0.0, -EIy_L2,    0.0,   EIy_L,     0.0,   0.0,     0.0,  EIy_L2,    0.0, 2*EIy_L,     0.0,
            0.0,  EIz_L2,     0.0,    0.0,     0.0,   EIz_L,   0.0, -EIz_L2,     0.0,    0.0,     0.0, 2*EIz_L,
          -EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,  EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,
            0.0, -EIz_L3,     0.0,    0.0,     0.0, -EIz_L2,   0.0,  EIz_L3,     0.0,    0.0,     0.0, -EIz_L2,
            0.0,     0.0, -EIy_L3,    0.0,  EIy_L2,     0.0,   0.0,     0.0,  EIy_L3,    0.0,  EIy_L2,     0.0,
            0.0,     0.0,     0.0,  -GJ_L,     0.0,     0.0,   0.0,     0.0,     0.0,   GJ_L,     0.0,     0.0,
            0.0,     0.0, -EIy_L2,    0.0, 2*EIy_L,     0.0,   0.0,     0.0,  EIy_L2,    0.0,   EIy_L,     0.0,
            0.0,  EIz_L2,     0.0,    0.0,     0.0, 2*EIz_L,   0.0, -EIz_L2,     0.0,    0.0,     0.0,   EIz_L;

    return LatTab[LatticeID].factor*T.transpose()*K*T;
}

Eigen::MatrixXd    PCG_Eigen::getBeamSM    (       const unsigned      LatticeID,
												   const double		   factor)
{
    Eigen::MatrixXd T(2*DofPerNode,2*DofPerNode);
    double C = (factor-LatTab[LatticeID].factor);

    T << LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2],0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[0][0], LatTab[LatticeID].axes[0][1], LatTab[LatticeID].axes[0][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[1][0], LatTab[LatticeID].axes[1][1], LatTab[LatticeID].axes[1][2],
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,LatTab[LatticeID].axes[2][0], LatTab[LatticeID].axes[2][1], LatTab[LatticeID].axes[2][2];
    double EA_L = LatTab[LatticeID].k[0];
    double EIz_L3 = 12.0*LatTab[LatticeID].k[3]/(LatTab[LatticeID].length*LatTab[LatticeID].length);
    double EIy_L3 = 12.0*LatTab[LatticeID].k[2]/(LatTab[LatticeID].length*LatTab[LatticeID].length);
    double EIz_L2 = 6.0*LatTab[LatticeID].k[3]/LatTab[LatticeID].length;
    double EIy_L2 = 6.0*LatTab[LatticeID].k[2]/LatTab[LatticeID].length;
    double EIz_L  = 4.0*LatTab[LatticeID].k[3];
    double EIy_L  = 4.0*LatTab[LatticeID].k[3];
    double GJ_L   = (LatTab[LatticeID].k[2]+LatTab[LatticeID].k[3]);
    Eigen::MatrixXd K(2*DofPerNode,2*DofPerNode);
    K  <<  EA_L,     0.0,     0.0,    0.0,     0.0,     0.0, -EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,
            0.0,  EIz_L3,     0.0,    0.0,     0.0,  EIz_L2,   0.0, -EIz_L3,     0.0,    0.0,     0.0,  EIz_L2,
            0.0,     0.0,  EIy_L3,    0.0, -EIy_L2,     0.0,   0.0,     0.0, -EIy_L3,    0.0, -EIy_L2,     0.0,
            0.0,     0.0,     0.0,   GJ_L,     0.0,     0.0,   0.0,     0.0,     0.0,  -GJ_L,     0.0,     0.0,
            0.0,     0.0, -EIy_L2,    0.0,   EIy_L,     0.0,   0.0,     0.0,  EIy_L2,    0.0, 2*EIy_L,     0.0,
            0.0,  EIz_L2,     0.0,    0.0,     0.0,   EIz_L,   0.0, -EIz_L2,     0.0,    0.0,     0.0, 2*EIz_L,
          -EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,  EA_L,     0.0,     0.0,    0.0,     0.0,     0.0,
            0.0, -EIz_L3,     0.0,    0.0,     0.0, -EIz_L2,   0.0,  EIz_L3,     0.0,    0.0,     0.0, -EIz_L2,
            0.0,     0.0, -EIy_L3,    0.0,  EIy_L2,     0.0,   0.0,     0.0,  EIy_L3,    0.0,  EIy_L2,     0.0,
            0.0,     0.0,     0.0,  -GJ_L,     0.0,     0.0,   0.0,     0.0,     0.0,   GJ_L,     0.0,     0.0,
            0.0,     0.0, -EIy_L2,    0.0, 2*EIy_L,     0.0,   0.0,     0.0,  EIy_L2,    0.0,   EIy_L,     0.0,
            0.0,  EIz_L2,     0.0,    0.0,     0.0, 2*EIz_L,   0.0, -EIz_L2,     0.0,    0.0,     0.0,   EIz_L;

    return C*T.transpose()*K*T;
}



void PCG_Eigen::fillSMSimple_Full	()
{
	unsigned DOF,i_pDisp;
	unsigned NodeID,nbNum,nbDOF;
	unsigned isFree[DofPerNode],isNbFree[DofPerNode];
	std::array<unsigned,DofPerNode>    pos,nbPos;
	std::cout << "Starting fillSMSimple_Full..." << std::endl;

	using Eigen::VectorXd;
	SM.resize(DOF_Num,DOF_Num);
	SM.reserve(VectorXd::Constant(DOF_Num,maxCoordNo*DofPerNode));

	for (DOF=0; DOF<freeNodeNum*DofPerNode; DOF+=DofPerNode) {
		NodeID 	= calNodeList[DOF/DofPerNode];
		nbNum 	= GNodeTab[NodeID].nbLatticeID.size();
		for (unsigned k=0; k<nbNum; k++) {
			nbDOF = NodeIDtoDOF[GNodeTab[NodeID].nbNodeID[k]];
			Eigen::MatrixXd     lSM(DofPerNode, 2*DofPerNode);
			if (KNum==1) {
				Eigen::MatrixXd     lSM2 = getSpringLocalSM(NodeID,k);
				lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
				lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
			} else if (KNum==2) {
				Eigen::MatrixXd     lSM2 = getShearLocalSM(NodeID,k);
				lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
				lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
			} else {
				lSM = getFullLocalSM2(NodeID,k);
			}
			for (unsigned row=0; row<DofPerNode; row++) {
				for (unsigned col=0; col<DofPerNode; col++) {
					SM.coeffRef(DOF+row,DOF+col) += lSM(row,col);
			    }
			}
			if (nbDOF<freeNodeNum*DofPerNode) {
				for (unsigned row=0; row<DofPerNode; row++) {
					for (unsigned col=0; col<DofPerNode; col++) {
						SM.coeffRef(DOF+row,nbDOF+col) += lSM(row,DofPerNode+col);
					}
				}
			}
			else if (nbDOF<DOF_Num) {
				i_pDisp = (lower_bound(pDispDOF.begin(),pDispDOF.end(),nbDOF)-pDispDOF.begin());
				if (!isFullRestrain(pDispNodeList[i_pDisp].Dir)) {
					for (unsigned d=0; d<DofPerNode; d++) {
						isNbFree[d]=!pDispNodeList[i_pDisp].Dir[d];
					}
					unsigned cnt = 0;
					for (unsigned i=0; i<DofPerNode; i++) {
						nbPos[i] = cnt*isNbFree[i];
					    cnt += isNbFree[i];
					}
					for (unsigned row=0; row<DofPerNode; row++) {
						for (unsigned col=0; col<DofPerNode; col++) {
					    	SM.coeffRef(DOF+row,nbDOF+nbPos[col]) += lSM(row,DofPerNode+col)*isNbFree[col];
					    }
					}
				}
			}
		}
	}
	DOF = freeNodeNum*DofPerNode;
	for (unsigned iDOF=freeNodeNum; iDOF<calNodeNum; iDOF++) {
		NodeID 	= calNodeList[iDOF];
		if (!isFullRestrain(pDispNodeList[iDOF-freeNodeNum].Dir)) {
			for (unsigned d=0; d<DofPerNode; d++) {
				isFree[d] = !pDispNodeList[iDOF-freeNodeNum].Dir[d];
			}
			unsigned cnt = 0;
			for (unsigned i=0; i<DofPerNode; i++) {
				pos[i] = cnt*isFree[i];
			    cnt += isFree[i];
			}
			nbNum 	= GNodeTab[NodeID].nbLatticeID.size();
			for (unsigned k=0; k<nbNum; k++) {
				nbDOF = NodeIDtoDOF[GNodeTab[NodeID].nbNodeID[k]];
				Eigen::MatrixXd     lSM(DofPerNode, 2*DofPerNode);
				if (KNum==1) {
					Eigen::MatrixXd     lSM2 = getSpringLocalSM(NodeID,k);
					lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
					lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
				} else if (KNum==2) {
					Eigen::MatrixXd     lSM2 = getShearLocalSM(NodeID,k);
					lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
					lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
				} else {
					lSM = getFullLocalSM2(NodeID,k);
				}
				for (unsigned row=0; row<DofPerNode; row++) {
					for (unsigned col=0; col<DofPerNode; col++) {
						SM.coeffRef(DOF+pos[row],DOF+pos[col]) += lSM(row,col)*isFree[row]*isFree[col];
				    }
				}
				if (nbDOF<freeNodeNum*DofPerNode) {
					for (unsigned row=0; row<DofPerNode; row++) {
						for (unsigned col=0; col<DofPerNode; col++) {
					    	SM.coeffRef(DOF+pos[row],nbDOF+col) += lSM(row,DofPerNode+col)*isFree[row];
					    }
					}

				}
				else if (nbDOF<DOF_Num)
				{
					i_pDisp = (lower_bound(pDispDOF.begin(),pDispDOF.end(),nbDOF)-pDispDOF.begin());
					if (!isFullRestrain(pDispNodeList[i_pDisp].Dir)) {
						for (unsigned d=0; d<DofPerNode; d++) {
							isNbFree[d] = !pDispNodeList[i_pDisp].Dir[d];
						}
						cnt = 0;
						for (unsigned i=0; i<DofPerNode; i++) {
							nbPos[i] = cnt*isNbFree[i];
						    cnt += isNbFree[i];
						}
						for (unsigned row=0; row<DofPerNode; row++) {
							for (unsigned col=0; col<DofPerNode; col++) {
								SM.coeffRef(DOF+pos[row],nbDOF+nbPos[col]) += lSM(row,DofPerNode+col)*isFree[row]*isNbFree[col];
						    }
						}
					}
				}
			}
		}
		DOF +=getFreeDOFNum(pDispNodeList[iDOF-freeNodeNum].Dir);
	}
	SM.makeCompressed();
	std::cout << "Starting fillSMSimple_Full finished" << std::endl;
}


Tensor PCG_Eigen::getMemForce			(	const unsigned			NodeID,
											const unsigned			k)
{
	Eigen::MatrixXd     lSM(DofPerNode, 2*DofPerNode);
	if (KNum==1) {
		Eigen::MatrixXd     lSM2 = getSpringLocalSM(NodeID,k);
		lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
		lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
	} else if (KNum==2) {
		Eigen::MatrixXd     lSM2 = getShearLocalSM(NodeID,k);
		lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
		lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
	} else {
		lSM = getFullLocalSM2(NodeID,k);
	}
	Eigen::VectorXd		vd(2*DofPerNode);
	unsigned nbNodeID = GNodeTab[NodeID].nbNodeID[k];
	for (unsigned d=0; d<DofPerNode; d++) {
		vd[d] = GNodeTab[NodeID].d[d]+GNodeTab[NodeID].d0[d];
		vd[d+DofPerNode] = GNodeTab[nbNodeID].d[d]+GNodeTab[NodeID].d0[d];
	}
	Eigen::VectorXd vf = lSM*vd;
	std::array<double,DofPerNode> f;
	for (unsigned d=0; d<DofPerNode; d++) {
		f[d] = vf[d];
	}
	return f;
}

Tensor PCG_Eigen::getResidual			(	const unsigned			NodeID)
{

	std::array<double,DofPerNode> residual;
	double maxLatForce = 0.0, maxLatMoment = 0.0;
	Geometry geo;
	for (unsigned d=0; d<DofPerNode; d++) {
		residual[d] = 0.0;
	}
	for (unsigned k=0; k<GNodeTab[NodeID].nbNodeID.size(); k++) {
		std::array<double,DofPerNode> memForce = getMemForce(NodeID,k);
		for (unsigned d=0; d<DofPerNode; d++) {
			residual[d] += memForce[d];
		}
		double forceMag=0.0, momentMag=0.0;
		for (unsigned d=0; d<Dim; d++) {
			forceMag += memForce[d]*memForce[d];
			momentMag += memForce[d+Dim]*memForce[d+Dim];
		}
		forceMag = std::sqrt (forceMag);
		momentMag = std::sqrt (momentMag);
		if (forceMag > maxLatForce) {
			maxLatForce = forceMag;
		}
		if (momentMag > maxLatMoment) {
			maxLatMoment = momentMag;
		}
	}
	for (unsigned d=0; d<DofPerNode; d++) {
		residual[d] -= GNodeTab[NodeID].extF[d];
	}
	return residual;
}

unsigned PCG_Eigen::getMaxCoordNum		( 		const unsigned 		shift,
												const unsigned		endDOF)
{
	unsigned NodeID;
	unsigned max = 0;
	unsigned start = shift/DofPerNode;
	unsigned end = start+endDOF/DofPerNode;
	for (unsigned i=start; i<end; i++) {
		NodeID = calNodeList[i];
		if (GNodeTab[NodeID].n>max) {
			max = GNodeTab[NodeID].n;
		}
	}
	cout << max << " ";
	return max;
}

void PCG_Eigen::fillSMblock_full2     ()
{
    using Eigen::VectorXd;

    Clock clock;
    clock.start("fillSMblock");
    chunkNum = ThreadNumUsed*((unsigned) (calNodeNum/nodeNumLimit+1));
    if (DOF_Num/(max(DOF_Num-freeNodeNum*DofPerNode,(unsigned) 1))>=chunkNum) {
        chunkSize = (DOF_Num/chunkNum)/DofPerNode*DofPerNode;
        chunkSizeEnd = DOF_Num - chunkSize*(chunkNum-1);
    } else {
        chunkSize = (freeNodeNum/(chunkNum-1))*DofPerNode;
        chunkSizeEnd = DOF_Num - chunkSize*(chunkNum-1);
    }
    cout << DOF_Num << " " << chunkNum << " " << chunkSize << " " << DOF_Num-freeNodeNum*DofPerNode << " " << chunkSizeEnd << endl;
    bSM.clear();
    bSM.resize(chunkNum);
    for (unsigned t=0; t<chunkNum; t++) {
        unsigned    isFree[DofPerNode],isNbFree[DofPerNode];
        std::array<unsigned,DofPerNode>    pos,nbPos;
        unsigned NodeID,nbNum,nbDOF;
        unsigned shift = t*chunkSize;
        unsigned endDOF;

        if (t==chunkNum-1) {
            bSM[t].resize(chunkSizeEnd,DOF_Num);
            endDOF = freeNodeNum*DofPerNode-chunkSize*(chunkNum-1);
            bSM[t].reserve(VectorXd::Constant(DOF_Num,(getMaxCoordNum(shift,endDOF))*DofPerNode));
        } else {
            bSM[t].resize(chunkSize,DOF_Num);
            endDOF = chunkSize;
            bSM[t].reserve(VectorXd::Constant(DOF_Num,maxCoordNo*DofPerNode));
        }
        for (unsigned DOF=0; DOF<endDOF; DOF+=DofPerNode) {
            NodeID  = calNodeList[(DOF+shift)/DofPerNode];
            nbNum   = GNodeTab[NodeID].nbLatticeID.size();
            for (unsigned k=0; k<nbNum; k++) {
            	Eigen::MatrixXd     lSM(DofPerNode, 2*DofPerNode);
            	if (KNum==1) {
            		Eigen::MatrixXd     lSM2 = getSpringLocalSM(NodeID,k);
            		lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
            		lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
            	} else if (KNum==2) {
            		Eigen::MatrixXd     lSM2 = getShearLocalSM(NodeID,k);
            		lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
            		lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
            	} else {
            		lSM = getFullLocalSM2(NodeID,k);
            	}
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t].coeffRef(DOF+row,DOF+shift+col) += lSM(row,col);
                    }
                }

                nbDOF = NodeIDtoDOF[GNodeTab[NodeID].nbNodeID[k]];
                if (nbDOF<freeNodeNum*DofPerNode) {
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t].coeffRef(DOF+row,nbDOF+col) += lSM(row,DofPerNode+col);
                        }
                    }
                } else if (nbDOF<DOF_Num) {
                    unsigned i_pDisp = (lower_bound(pDispDOF.begin(),pDispDOF.end(),nbDOF)-pDispDOF.begin());
                    if (!isFullRestrain(pDispNodeList[i_pDisp].Dir)) {
                        for (unsigned d=0; d<DofPerNode; d++) {
                            isNbFree[d]=!pDispNodeList[i_pDisp].Dir[d];
                        }
                        unsigned cnt = 0;
                        for (unsigned i=0; i<DofPerNode; i++) {
                            nbPos[i] = cnt*isNbFree[i];
                            cnt += isNbFree[i];
                        }
                        for (unsigned row=0; row<DofPerNode; row++) {
                            for (unsigned col=0; col<DofPerNode; col++) {
                                bSM[t].coeffRef(DOF+row,nbDOF+nbPos[col]) += lSM(row,DofPerNode+col)*isNbFree[col];
                            }
                        }
                    }
                }
            }
        }

        if (t == chunkNum-1) {
            unsigned DOF = freeNodeNum*DofPerNode-chunkSize*(chunkNum-1);
            for (unsigned iDOF=freeNodeNum; iDOF<calNodeNum; iDOF++) {
                NodeID  = calNodeList[iDOF];
                if (!isFullRestrain(pDispNodeList[iDOF-freeNodeNum].Dir)) {
                    for (unsigned d=0; d<DofPerNode; d++) {
                        isFree[d] = !pDispNodeList[iDOF-freeNodeNum].Dir[d];
                    }
                    unsigned cnt = 0;
                    for (unsigned i=0; i<DofPerNode; i++) {
                        pos[i] = cnt*isFree[i];
                        cnt += isFree[i];
                    }
                    nbNum   = GNodeTab[NodeID].nbLatticeID.size();
                    for (unsigned k=0; k<nbNum; k++) {
                        nbDOF = NodeIDtoDOF[GNodeTab[NodeID].nbNodeID[k]];
                        Eigen::MatrixXd     lSM(DofPerNode, 2*DofPerNode);
                        if (KNum==1) {
                        	Eigen::MatrixXd     lSM2 = getSpringLocalSM(NodeID,k);
                            lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
                            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
                        } else if (KNum==2) {
                            Eigen::MatrixXd     lSM2 = getShearLocalSM(NodeID,k);
                            lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
                            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
                        } else {
                            lSM = getFullLocalSM2(NodeID,k);
                        }
                        for (unsigned row=0; row<DofPerNode; row++) {
                            for (unsigned col=0; col<DofPerNode; col++) {
                                bSM[t].coeffRef(DOF+pos[row],DOF+shift+pos[col]) += lSM(row,col)*isFree[row]*isFree[col];
                            }
                        }
                        if (nbDOF<freeNodeNum*DofPerNode) {
                            for (unsigned row=0; row<DofPerNode; row++) {
                                for (unsigned col=0; col<DofPerNode; col++) {
                                    bSM[t].coeffRef(DOF+pos[row],nbDOF+col) += lSM(row,DofPerNode+col)*isFree[row];
                                }
                            }
                        } else if (nbDOF<DOF_Num) {
                            unsigned i_pDisp = (lower_bound(pDispDOF.begin(),pDispDOF.end(),nbDOF)-pDispDOF.begin());
                            if (!isFullRestrain(pDispNodeList[i_pDisp].Dir)) {
                                for (unsigned d=0; d<DofPerNode; d++) {
                                    isNbFree[d] = !pDispNodeList[i_pDisp].Dir[d];
                                }
                                cnt = 0;
                                for (unsigned i=0; i<DofPerNode; i++) {
                                    nbPos[i] = cnt*isNbFree[i];
                                    cnt += isNbFree[i];
                                }
                                for (unsigned row=0; row<DofPerNode; row++) {
                                    for (unsigned col=0; col<DofPerNode; col++) {
                                        bSM[t].coeffRef(DOF+pos[row],nbDOF+nbPos[col]) += lSM(row,DofPerNode+col)*isFree[row]*isNbFree[col];
                                    }
                                }
                            }
                        }
                    }
                }
                DOF +=getFreeDOFNum(pDispNodeList[iDOF-freeNodeNum].Dir);
            }
        }
        bSM[t].makeCompressed();
    }
    omp_set_num_threads(ThreadNumUsed);
    clock.get();
}

void PCG_Eigen::fillChunkNumList ()
{
    chunkNum = ThreadNumUsed*((unsigned) (calNodeNum/nodeNumLimit+1));
    unsigned chunkSizeAvg = (DOF_Num/chunkNum)/DofPerNode*DofPerNode;
    freeLocation = std::min(freeNodeNum*DofPerNode / chunkSizeAvg,chunkNum-1);
    chunkNumList.clear();
    chunkPos.clear();
    for (unsigned t=0; t<freeLocation; t++) {
        chunkNumList.push_back(chunkSizeAvg);
    }
    unsigned chunkCnt = freeNodeNum*DofPerNode - chunkSizeAvg*freeLocation;
    unsigned iDOF=0;
    for (unsigned t=freeLocation; t<chunkNum; t++) {
        if (t==chunkNum-1) {
            chunkSizeAvg = UINT_MAX;
        }
        while ((chunkCnt<chunkSizeAvg)&&(iDOF<pDispNodeList.size())) {
            chunkCnt+=getFreeDOFNum(pDispNodeList[iDOF].Dir);
            iDOF++;
        }
        chunkNumList.push_back(chunkCnt);
        chunkCnt = 0;
    }
    unsigned cnt = 0;
    dualOut << "ChunkNumList: ";
    for (unsigned t=0; t<chunkNum; t++) {
        chunkPos.push_back(cnt);
        cnt += chunkNumList[t];
    }
    dualOut << std::endl;
    dualOut << "freeLocation = " << freeLocation << std::endl;
    dualOut << "total = " << cnt << " , DOF_Num = " << DOF_Num << std::endl;
    std::cout << DOF_Num << " " << chunkNum  << " " << DOF_Num-freeNodeNum*DofPerNode << std::endl;
}

std::array<unsigned,2> PCG_Eigen::getPosDOF (     const unsigned      DOF)
{
    std::array<unsigned,2> pos;
    pos[0] = 0;
    pos[1] = DOF;
    while (pos[1]>=chunkNumList[pos[0]]) {
        pos[1] -= chunkNumList[pos[0]];
        pos[0]++;
    }
    return pos;
}

bool PCG_Eigen::updateSM_full              (      const unsigned          LatticeID,
                                                  const double            factor)
{
    static bool isRepeat = false;
    if (isUpdateDOF) {
        if (!isRepeat) {
            dualOut << "The stiffness matrix will be reassembled, no update is required!" << std::endl;
            isRepeat = true;
        }
        return true;
    }
    unsigned NodeID1 = LatTab[LatticeID].nb[0];
    unsigned NodeID2 = LatTab[LatticeID].nb[1];
    if ((GNodeTab[NodeID1].type== Type_Unstable)||(GNodeTab[NodeID2].type== Type_Unstable)) {
        tripOut << "[PCG_Eigen::updateSM]: One of the nodes is unstable, please check!" << std::endl;
        return false;
    }
    unsigned DOF1 = NodeIDtoDOF [NodeID1];
    unsigned DOF2 = NodeIDtoDOF [NodeID2];
    std::array<unsigned,2> pos = getPosDOF (DOF1);
    unsigned t1 = pos[0];
    unsigned tDOF1 = pos[1];
    pos = getPosDOF (DOF2);
    unsigned t2 = pos[0];
    unsigned tDOF2 = pos[1];

    Eigen::MatrixXd lSM(2*DofPerNode,2*DofPerNode);
    if (KNum==1) {
        Eigen::MatrixXd partLSM = getSpringLocalSM(LatticeID,factor);
        lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
        lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
        lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
    } else if (KNum==2) {
        Eigen::MatrixXd partLSM = getShearLocalSM(LatticeID,factor);
        lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
        lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
        lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
    } else if (DofPerNode==6) {
        lSM = getFullLocalSM2(LatticeID,factor);
    }
    if ((DOF1<freeNodeNum*DofPerNode)&&(DOF2<freeNodeNum*DofPerNode)) {
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t1].coeffRef(tDOF1+row,DOF1+col) += lSM(row,col);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t1].coeffRef(tDOF1+row,DOF2+col) += lSM(row,col+DofPerNode);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t2].coeffRef(tDOF2+row,DOF1+col) += lSM(row+DofPerNode,col);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t2].coeffRef(tDOF2+row,DOF2+col) += lSM(row+DofPerNode,col+DofPerNode);
            }
        }
    } else if ((DOF1<DOF_Num)&&(DOF2<DOF_Num)){
        std::array<unsigned,DofPerNode>    isFree1,isFree2;
        std::array<unsigned,DofPerNode>    pos1,pos2;
        for (unsigned i=0; i<DofPerNode; i++) {
            isFree1[i] = isFree2[i] = true;
            pos1[i] = pos2[i] = i;
        }
        bool isFullRestrain1 = false;
        bool isFullRestrain2 = false;
        if (DOF1>=freeNodeNum*DofPerNode) {
            unsigned iDisp1=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF1)-pDispDOF.begin());
            isFullRestrain1 = isFullRestrain(pDispNodeList[iDisp1].Dir);
            if (!isFullRestrain1) {
                for (unsigned d=0; d<DofPerNode; d++) {
                    isFree1[d]=!pDispNodeList[iDisp1].Dir[d];
                }
                unsigned cnt = 0;
                for (unsigned i=0; i<DofPerNode; i++) {
                    pos1[i] = cnt*isFree1[i];
                    cnt += isFree1[i];
                }
            }
        }
        if (DOF2>=freeNodeNum*DofPerNode) {
            unsigned iDisp2=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF2)-pDispDOF.begin());
            isFullRestrain2 = isFullRestrain(pDispNodeList[iDisp2].Dir);
            if (!isFullRestrain2) {
                for (unsigned d=0; d<DofPerNode; d++) {
                    isFree2[d]=!pDispNodeList[iDisp2].Dir[d];
                }
                unsigned cnt = 0;
                for (unsigned i=0; i<DofPerNode; i++) {
                    pos2[i] = cnt*isFree2[i];
                    cnt += isFree2[i];
                }
            }
        }
        if (!isFullRestrain1) {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+pos1[row],DOF1+pos1[col]) += lSM(row,col)*isFree1[row]*isFree1[col];
                }
            }
            if (!isFullRestrain2) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t1].coeffRef(tDOF1+pos1[row],DOF2+pos2[col]) += lSM(row,col+DofPerNode)*isFree1[row]*isFree2[col];
                    }
                }
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF1+pos1[col]) += lSM(row+DofPerNode,col)*isFree2[row]*isFree1[col];
                    }
                }
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                    }
                }

            }
        } else if (!isFullRestrain2) {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                }
            }
        }
    }
    isRepeat = false;
    return true;
}

void PCG_Eigen::fillSMblock_full3     ()
{
    using Eigen::VectorXd;

    Clock clock;
    clock.start("fillSMblock3");
    fillChunkNumList ();
    bSM.clear();
    bSM.resize(chunkNum);

    unsigned shift = 0;
    unsigned iDOF = 0;
    for (unsigned t=0; t<chunkNum; t++) {
        unsigned    isFree[DofPerNode],isNbFree[DofPerNode];
        std::array<unsigned,DofPerNode>    pos,nbPos;
        unsigned NodeID,nbNum,nbDOF;
        bSM[t].resize(chunkNumList[t],DOF_Num);
        bSM[t].reserve(VectorXd::Constant(DOF_Num,maxCoordNo*DofPerNode));
        if (t<=freeLocation) {
            unsigned endDOF = std::min (chunkNumList[t],freeNodeNum*DofPerNode-shift);
            for (unsigned DOF=0; DOF<endDOF; DOF+=DofPerNode) {
                NodeID  = calNodeList[(DOF+shift)/DofPerNode];
                nbNum   = GNodeTab[NodeID].nbLatticeID.size();
                for (unsigned k=0; k<nbNum; k++) {
                    nbDOF = NodeIDtoDOF[GNodeTab[NodeID].nbNodeID[k]];
                    Eigen::MatrixXd     lSM(DofPerNode, 2*DofPerNode);
                    if (KNum==1) {
                        Eigen::MatrixXd     lSM2 = getSpringLocalSM(NodeID,k);
                        lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
                        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
                    } else if (KNum==2) {
                        Eigen::MatrixXd     lSM2 = getShearLocalSM(NodeID,k);
                        lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
                        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
                    } else {
                        lSM = getFullLocalSM2(NodeID,k);
                    }
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t].coeffRef(DOF+row,DOF+shift+col) += lSM(row,col);
                        }
                    }

                    if (nbDOF<freeNodeNum*DofPerNode) {
                        for (unsigned row=0; row<DofPerNode; row++) {
                            for (unsigned col=0; col<DofPerNode; col++) {
                                bSM[t].coeffRef(DOF+row,nbDOF+col) += lSM(row,DofPerNode+col);
                            }
                        }
                    } else if (nbDOF<DOF_Num) {
                        unsigned i_pDisp = (lower_bound(pDispDOF.begin(),pDispDOF.end(),nbDOF)-pDispDOF.begin());
                        if (!isFullRestrain(pDispNodeList[i_pDisp].Dir)) {
                            for (unsigned d=0; d<DofPerNode; d++) {
                                isNbFree[d]=!pDispNodeList[i_pDisp].Dir[d];
                            }
                            unsigned cnt = 0;
                            for (unsigned i=0; i<DofPerNode; i++) {
                                nbPos[i] = cnt*isNbFree[i];
                                cnt += isNbFree[i];
                            }
                            for (unsigned row=0; row<DofPerNode; row++) {
                                for (unsigned col=0; col<DofPerNode; col++) {
                                    bSM[t].coeffRef(DOF+row,nbDOF+nbPos[col]) += lSM(row,DofPerNode+col)*isNbFree[col];
                                }
                            }
                        }
                    }
                }
            }
        }
        unsigned DOF = 0;
        if (t == freeLocation) {
            DOF = freeNodeNum*DofPerNode-shift;
        }
        if (t>=freeLocation) {
            while ((DOF<chunkNumList[t])&&(iDOF<calNodeNum-freeNodeNum)) {
                NodeID  = calNodeList[iDOF+freeNodeNum];
                if (!isFullRestrain(pDispNodeList[iDOF].Dir)) {
                    for (unsigned d=0; d<DofPerNode; d++) {
                        isFree[d] = !pDispNodeList[iDOF].Dir[d];
                    }
                    unsigned cnt = 0;
                    for (unsigned i=0; i<DofPerNode; i++) {
                        pos[i] = cnt*isFree[i];
                        cnt += isFree[i];
                    }
                    nbNum   = GNodeTab[NodeID].nbLatticeID.size();
                    for (unsigned k=0; k<nbNum; k++) {
                        nbDOF = NodeIDtoDOF[GNodeTab[NodeID].nbNodeID[k]];
                        Eigen::MatrixXd     lSM(DofPerNode, 2*DofPerNode);
                        if (KNum==1) {
                            Eigen::MatrixXd     lSM2 = getSpringLocalSM(NodeID,k);
                            lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
                            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
                        } else if (KNum==2) {
                            Eigen::MatrixXd     lSM2 = getShearLocalSM(NodeID,k);
                            lSM.block(0,0,DofPerNode,DofPerNode) = lSM2;
                            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -lSM2;
                        } else {
                            lSM = getFullLocalSM2(NodeID,k);
                        }
                        for (unsigned row=0; row<DofPerNode; row++) {
                            for (unsigned col=0; col<DofPerNode; col++) {
                                bSM[t].coeffRef(DOF+pos[row],DOF+shift+pos[col]) += lSM(row,col)*isFree[row]*isFree[col];
                            }
                        }
                        if (nbDOF<freeNodeNum*DofPerNode) {
                            for (unsigned row=0; row<DofPerNode; row++) {
                                for (unsigned col=0; col<DofPerNode; col++) {
                                    bSM[t].coeffRef(DOF+pos[row],nbDOF+col) += lSM(row,DofPerNode+col)*isFree[row];
                                }
                            }
                        } else if (nbDOF<DOF_Num) {
                            unsigned i_pDisp = (lower_bound(pDispDOF.begin(),pDispDOF.end(),nbDOF)-pDispDOF.begin());
                            if (!isFullRestrain(pDispNodeList[i_pDisp].Dir)) {
                                for (unsigned d=0; d<DofPerNode; d++) {
                                    isNbFree[d] = !pDispNodeList[i_pDisp].Dir[d];
                                }
                                cnt = 0;
                                for (unsigned i=0; i<DofPerNode; i++) {
                                    nbPos[i] = cnt*isNbFree[i];
                                    cnt += isNbFree[i];
                                }
                                for (unsigned row=0; row<DofPerNode; row++) {
                                    for (unsigned col=0; col<DofPerNode; col++) {
                                        bSM[t].coeffRef(DOF+pos[row],nbDOF+nbPos[col]) += lSM(row,DofPerNode+col)*isFree[row]*isNbFree[col];
                                    }
                                }
                            }
                        }
                    }
                }
                DOF += getFreeDOFNum(pDispNodeList[iDOF].Dir);
                iDOF++;
            }
        }
        bSM[t].makeCompressed();
        shift += chunkNumList[t];
    }
    omp_set_num_threads(ThreadNumUsed);
    clock.get();
}

void PCG_Eigen::fillSMblock_full4              ()
{
    Clock clock;
    clock.start("fillSMblock4");
    fillChunkNumList ();
    bSM.clear();
    bSM.resize(chunkNum);
    for (unsigned t=0; t<chunkNum; t++) {
        bSM[t].resize(chunkNumList[t],DOF_Num);
        bSM[t].reserve(Eigen::VectorXd::Constant(DOF_Num,maxCoordNo*DofPerNode));
    }

    for (unsigned LatticeID = 0; LatticeID < tLatticeNum; LatticeID++) {
        unsigned NodeID1 = LatTab[LatticeID].nb[0];
        unsigned NodeID2 = LatTab[LatticeID].nb[1];
        if ((GNodeTab[NodeID1].type== Type_Unstable)||(GNodeTab[NodeID2].type== Type_Unstable)) {
            tripOut << "[PCG_Eigen::updateSM]: One of the nodes is unstable, please check!" << std::endl;
            continue;
        }
        unsigned DOF1 = NodeIDtoDOF [NodeID1];
        unsigned DOF2 = NodeIDtoDOF [NodeID2];

        std::array<unsigned,2> pos = getPosDOF (DOF1);
        unsigned t1 = pos[0];
        unsigned tDOF1 = pos[1];
        pos = getPosDOF (DOF2);
        unsigned t2 = pos[0];
        unsigned tDOF2 = pos[1];
        Eigen::MatrixXd lSM(2*DofPerNode,2*DofPerNode);
        if (KNum==1) {
            Eigen::MatrixXd partLSM = getSpringLocalSM(LatticeID);
            lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
        } else if (KNum==2) {
            Eigen::MatrixXd partLSM = getShearLocalSM(LatticeID);
            lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
        } else if (DofPerNode==6) {
            lSM = getFullLocalSM2(LatticeID);
        }
        if ((DOF1<freeNodeNum*DofPerNode)&&(DOF2<freeNodeNum*DofPerNode)) {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+row,DOF1+col) += lSM(row,col);
                }
            }
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+row,DOF2+col) += lSM(row,col+DofPerNode);
                }
            }
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+row,DOF1+col) += lSM(row+DofPerNode,col);
                }
            }
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+row,DOF2+col) += lSM(row+DofPerNode,col+DofPerNode);
                }
            }
        } else if ((DOF1<DOF_Num)&&(DOF2<DOF_Num)){
            std::array<unsigned,DofPerNode>    isFree1,isFree2;
            std::array<unsigned,DofPerNode>    pos1,pos2;
            for (unsigned i=0; i<DofPerNode; i++) {
                isFree1[i] = isFree2[i] = true;
                pos1[i] = pos2[i] = i;
            }
            bool isFullRestrain1 = false;
            bool isFullRestrain2 = false;
            if (DOF1>=freeNodeNum*DofPerNode) {
                unsigned iDisp1=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF1)-pDispDOF.begin());
                isFullRestrain1 = isFullRestrain(pDispNodeList[iDisp1].Dir);
                if (!isFullRestrain1) {
                    for (unsigned d=0; d<DofPerNode; d++) {
                        isFree1[d]=!pDispNodeList[iDisp1].Dir[d];
                    }
                    unsigned cnt = 0;
                    for (unsigned i=0; i<DofPerNode; i++) {
                        pos1[i] = cnt*isFree1[i];
                        cnt += isFree1[i];
                    }
                }
            }
            if (DOF2>=freeNodeNum*DofPerNode) {
                unsigned iDisp2=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF2)-pDispDOF.begin());
                isFullRestrain2 = isFullRestrain(pDispNodeList[iDisp2].Dir);
                if (!isFullRestrain2) {
                    for (unsigned d=0; d<DofPerNode; d++) {
                        isFree2[d]=!pDispNodeList[iDisp2].Dir[d];
                    }
                    unsigned cnt = 0;
                    for (unsigned i=0; i<DofPerNode; i++) {
                        pos2[i] = cnt*isFree2[i];
                        cnt += isFree2[i];
                    }
                }
            }
            if (!isFullRestrain1) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t1].coeffRef(tDOF1+pos1[row],DOF1+pos1[col]) += lSM(row,col)*isFree1[row]*isFree1[col];
                    }
                }
                if (!isFullRestrain2) {
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t1].coeffRef(tDOF1+pos1[row],DOF2+pos2[col]) += lSM(row,col+DofPerNode)*isFree1[row]*isFree2[col];
                        }
                    }
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t2].coeffRef(tDOF2+pos2[row],DOF1+pos1[col]) += lSM(row+DofPerNode,col)*isFree2[row]*isFree1[col];
                        }
                    }
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                        }
                    }

                }
            } else if (!isFullRestrain2) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                    }
                }
            }
        }
    }
    clock.get();
    return;
}

bool PCG_Eigen::updateSM_full              (      const unsigned          LatticeID,
                                                  const double            factor,
                                                  const double            factorS)
{
    std::vector<SpMat>          bSM_temp = bSM;
    static bool isRepeat = false;
    if (isUpdateDOF) {
        if (!isRepeat) {
            dualOut << "The stiffness matrix will be reassembled, no update is required!" << std::endl;
            isRepeat = true;
        }
        return true;
    }
    unsigned NodeID1 = LatTab[LatticeID].nb[0];
    unsigned NodeID2 = LatTab[LatticeID].nb[1];
    if ((GNodeTab[NodeID1].type== Type_Unstable)||(GNodeTab[NodeID2].type== Type_Unstable)) {
        tripOut << "[PCG_Eigen::updateSM]: One of the nodes is unstable, please check!" << std::endl;
        return false;
    }
    unsigned DOF1 = NodeIDtoDOF [NodeID1];
    unsigned DOF2 = NodeIDtoDOF [NodeID2];
    std::array<unsigned,2> pos = getPosDOF (DOF1);
    unsigned t1 = pos[0];
    unsigned tDOF1 = pos[1];
    pos = getPosDOF (DOF2);
    unsigned t2 = pos[0];
    unsigned tDOF2 = pos[1];

    Eigen::MatrixXd lSM(2*DofPerNode,2*DofPerNode);
    if (KNum==1) {
        Eigen::MatrixXd partLSM = getSpringLocalSM(LatticeID,factor);
        lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
        lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
        lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
    } else if (KNum==2) {
        Eigen::MatrixXd partLSM = getShearLocalSM(LatticeID,factor);
        lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
        lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
        lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
        lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
    } else if (DofPerNode==6) {
        lSM = getFullLocalSM2(LatticeID,factor,factorS);
    }
    if ((DOF1<freeNodeNum*DofPerNode)&&(DOF2<freeNodeNum*DofPerNode)) {
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t1].coeffRef(tDOF1+row,DOF1+col) += lSM(row,col);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t1].coeffRef(tDOF1+row,DOF2+col) += lSM(row,col+DofPerNode);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t2].coeffRef(tDOF2+row,DOF1+col) += lSM(row+DofPerNode,col);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t2].coeffRef(tDOF2+row,DOF2+col) += lSM(row+DofPerNode,col+DofPerNode);
            }
        }
    } else if ((DOF1<DOF_Num)&&(DOF2<DOF_Num)){
        std::array<unsigned,DofPerNode>    isFree1,isFree2;
        std::array<unsigned,DofPerNode>    pos1,pos2;
        for (unsigned i=0; i<DofPerNode; i++) {
            isFree1[i] = isFree2[i] = true;
            pos1[i] = pos2[i] = i;
        }
        bool isFullRestrain1 = false;
        bool isFullRestrain2 = false;
        if (DOF1>=freeNodeNum*DofPerNode) {
            unsigned iDisp1=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF1)-pDispDOF.begin());
            isFullRestrain1 = isFullRestrain(pDispNodeList[iDisp1].Dir);
            if (!isFullRestrain1) {
                for (unsigned d=0; d<DofPerNode; d++) {
                    isFree1[d]=!pDispNodeList[iDisp1].Dir[d];
                }
                unsigned cnt = 0;
                for (unsigned i=0; i<DofPerNode; i++) {
                    pos1[i] = cnt*isFree1[i];
                    cnt += isFree1[i];
                }
            }
        }
        if (DOF2>=freeNodeNum*DofPerNode) {
            unsigned iDisp2=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF2)-pDispDOF.begin());
            isFullRestrain2 = isFullRestrain(pDispNodeList[iDisp2].Dir);
            if (!isFullRestrain2) {
                for (unsigned d=0; d<DofPerNode; d++) {
                    isFree2[d]=!pDispNodeList[iDisp2].Dir[d];
                }
                unsigned cnt = 0;
                for (unsigned i=0; i<DofPerNode; i++) {
                    pos2[i] = cnt*isFree2[i];
                    cnt += isFree2[i];
                }
            }
        }
        if (!isFullRestrain1) {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+pos1[row],DOF1+pos1[col]) += lSM(row,col)*isFree1[row]*isFree1[col];
                }
            }
            if (!isFullRestrain2) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t1].coeffRef(tDOF1+pos1[row],DOF2+pos2[col]) += lSM(row,col+DofPerNode)*isFree1[row]*isFree2[col];
                    }
                }
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF1+pos1[col]) += lSM(row+DofPerNode,col)*isFree2[row]*isFree1[col];
                    }
                }
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                    }
                }

            }
        } else if (!isFullRestrain2) {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                }
            }
        }
    }
    isRepeat = false;
    return true;
}

bool PCG_Eigen::updateSM_full              (    std::vector<unsigned>&     latList,
                                                const double             _reConRatio,
                                                const double             gamma)
{
    static bool isRepeat = false;
    if (isUpdateDOF) {
        if (!isRepeat) {
            dualOut << "The stiffness matrix will be reassembled, no update is required!" << std::endl;
            isRepeat = true;
        }
        return true;
    }
    bool isTiny = true;
    if ((_reConRatio<10*Tiny)&&(gamma<10*Tiny)) {
        isTiny = false;
    }
    #pragma omp parallel for if (Use_OpenMP) shared (LatTab,GNodeTab,latList) firstprivate (ReConRatio,isTiny)
    for (unsigned i=0; i<latList.size(); i++) {
        unsigned LatticeID = latList[i];
        double factorS = gamma;
        if (isTiny) {
            Tensor memForce = getMemForce(LatticeID);
            double sigma_s = std::sqrt(memForce[1]*memForce[1]+memForce[2]*memForce[2])/LatTab[LatticeID].area;
            factorS *= std::max(std::min(1.0,gamma*LatTab[LatticeID].fShear/sigma_s),0.0);
        }
        unsigned NodeID1 = LatTab[LatticeID].nb[0];
        unsigned NodeID2 = LatTab[LatticeID].nb[1];
        if ((GNodeTab[NodeID1].type== Type_Unstable)||(GNodeTab[NodeID2].type== Type_Unstable)) {
            tripOut << "[PCG_Eigen::updateSM]: One of the nodes is unstable, please check!" << std::endl;
            continue;
        }
        unsigned DOF1 = NodeIDtoDOF [NodeID1];
        unsigned DOF2 = NodeIDtoDOF [NodeID2];

        std::array<unsigned,2> pos = getPosDOF (DOF1);
        unsigned t1 = pos[0];
        unsigned tDOF1 = pos[1];
        pos = getPosDOF (DOF2);
        unsigned t2 = pos[0];
        unsigned tDOF2 = pos[1];

        Eigen::MatrixXd lSM(2*DofPerNode,2*DofPerNode);
        if (KNum==1) {
            Eigen::MatrixXd partLSM = getSpringLocalSM(LatticeID,_reConRatio);
            lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
        } else if (KNum==2) {
            Eigen::MatrixXd partLSM = getShearLocalSM(LatticeID,_reConRatio);
            lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
        } else if (DofPerNode==6) {
            lSM = getFullLocalSM2(LatticeID,_reConRatio,factorS);
        }

        if ((DOF1<freeNodeNum*DofPerNode)&&(DOF2<freeNodeNum*DofPerNode)) {
        #pragma omp critical
            {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+row,DOF1+col) += lSM(row,col);
                }
            }
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+row,DOF2+col) += lSM(row,col+DofPerNode);
                }
            }
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+row,DOF1+col) += lSM(row+DofPerNode,col);
                }
            }
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+row,DOF2+col) += lSM(row+DofPerNode,col+DofPerNode);
                }
            }
            }
        } else if ((DOF1<DOF_Num)&&(DOF2<DOF_Num)){
            std::array<unsigned,DofPerNode>    isFree1,isFree2;
            std::array<unsigned,DofPerNode>    pos1,pos2;
            for (unsigned i=0; i<DofPerNode; i++) {
                isFree1[i] = isFree2[i] = true;
                pos1[i] = pos2[i] = i;
            }
            bool isFullRestrain1 = false;
            bool isFullRestrain2 = false;
            if (DOF1>=freeNodeNum*DofPerNode) {
                unsigned iDisp1=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF1)-pDispDOF.begin());
                isFullRestrain1 = isFullRestrain(pDispNodeList[iDisp1].Dir);
                if (!isFullRestrain1) {
                    for (unsigned d=0; d<DofPerNode; d++) {
                        isFree1[d]=!pDispNodeList[iDisp1].Dir[d];
                    }
                    unsigned cnt = 0;
                    for (unsigned i=0; i<DofPerNode; i++) {
                        pos1[i] = cnt*isFree1[i];
                        cnt += isFree1[i];
                    }
                }
            }
            if (DOF2>=freeNodeNum*DofPerNode) {
                unsigned iDisp2=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF2)-pDispDOF.begin());
                isFullRestrain2 = isFullRestrain(pDispNodeList[iDisp2].Dir);
                if (!isFullRestrain2) {
                    for (unsigned d=0; d<DofPerNode; d++) {
                        isFree2[d]=!pDispNodeList[iDisp2].Dir[d];
                    }
                    unsigned cnt = 0;
                    for (unsigned i=0; i<DofPerNode; i++) {
                        pos2[i] = cnt*isFree2[i];
                        cnt += isFree2[i];
                    }
                }
            }
            #pragma omp critical
            {
            if (!isFullRestrain1) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t1].coeffRef(tDOF1+pos1[row],DOF1+pos1[col]) += lSM(row,col)*isFree1[row]*isFree1[col];
                    }
                }
                if (!isFullRestrain2) {
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t1].coeffRef(tDOF1+pos1[row],DOF2+pos2[col]) += lSM(row,col+DofPerNode)*isFree1[row]*isFree2[col];
                        }
                    }
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t2].coeffRef(tDOF2+pos2[row],DOF1+pos1[col]) += lSM(row+DofPerNode,col)*isFree2[row]*isFree1[col];
                        }
                    }
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                        }
                    }

                }
            } else if (!isFullRestrain2) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                    }
                }
            }
            }
        }
        LatTab[LatticeID].factor = _reConRatio;
        LatTab[LatticeID].factorShear = factorS;
    }
    isRepeat = false;
    return true;
}

bool PCG_Eigen::updateSM_full              (    std::vector<unsigned>&       latList,
                                                const double                 _reConRatio,
                                                std::vector<double>&         gammaList)
{
    static bool isRepeat = false;
    if (isUpdateDOF) {
        if (!isRepeat) {
            dualOut << "The stiffness matrix will be reassembled, no update is required!" << std::endl;
            isRepeat = true;
        }
        return true;
    }

    #pragma omp parallel for if (Use_OpenMP) shared (LatTab,GNodeTab,gammaList,latList)
    for (unsigned i=0; i<latList.size(); i++) {
        unsigned LatticeID = latList[i];
        double gamma = gammaList[i];

        unsigned NodeID1 = LatTab[LatticeID].nb[0];
        unsigned NodeID2 = LatTab[LatticeID].nb[1];
        if ((GNodeTab[NodeID1].type== Type_Unstable)||(GNodeTab[NodeID2].type== Type_Unstable)) {
            tripOut << "[PCG_Eigen::updateSM]: One of the nodes is unstable, please check!" << std::endl;
            continue;
        }
        unsigned DOF1 = NodeIDtoDOF [NodeID1];
        unsigned DOF2 = NodeIDtoDOF [NodeID2];

        std::array<unsigned,2> pos = getPosDOF (DOF1);
        unsigned t1 = pos[0];
        unsigned tDOF1 = pos[1];
        pos = getPosDOF (DOF2);
        unsigned t2 = pos[0];
        unsigned tDOF2 = pos[1];

        Eigen::MatrixXd lSM(2*DofPerNode,2*DofPerNode);
        if (KNum==1) {
            Eigen::MatrixXd partLSM = getSpringLocalSM(LatticeID,_reConRatio);
            lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
        } else if (KNum==2) {
            Eigen::MatrixXd partLSM = getShearLocalSM(LatticeID,_reConRatio);
            lSM.block(0,0,DofPerNode,DofPerNode) = partLSM;
            lSM.block(0,DofPerNode,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,0,DofPerNode,DofPerNode) = -partLSM;
            lSM.block(DofPerNode,DofPerNode,DofPerNode,DofPerNode) = partLSM;
        } else if (DofPerNode==6) {
            lSM = getFullLocalSM2(LatticeID,_reConRatio,gamma);
        }

        if ((DOF1<freeNodeNum*DofPerNode)&&(DOF2<freeNodeNum*DofPerNode)) {
        #pragma omp critical
            {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+row,DOF1+col) += lSM(row,col);
                }
            }
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+row,DOF2+col) += lSM(row,col+DofPerNode);
                }
            }
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+row,DOF1+col) += lSM(row+DofPerNode,col);
                }
            }
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+row,DOF2+col) += lSM(row+DofPerNode,col+DofPerNode);
                }
            }
            }
        } else if ((DOF1<DOF_Num)&&(DOF2<DOF_Num)){
            std::array<unsigned,DofPerNode>    isFree1,isFree2;
            std::array<unsigned,DofPerNode>    pos1,pos2;
            for (unsigned i=0; i<DofPerNode; i++) {
                isFree1[i] = isFree2[i] = true;
                pos1[i] = pos2[i] = i;
            }
            bool isFullRestrain1 = false;
            bool isFullRestrain2 = false;
            if (DOF1>=freeNodeNum*DofPerNode) {
                unsigned iDisp1=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF1)-pDispDOF.begin());
                isFullRestrain1 = isFullRestrain(pDispNodeList[iDisp1].Dir);
                if (!isFullRestrain1) {
                    for (unsigned d=0; d<DofPerNode; d++) {
                        isFree1[d]=!pDispNodeList[iDisp1].Dir[d];
                    }
                    unsigned cnt = 0;
                    for (unsigned i=0; i<DofPerNode; i++) {
                        pos1[i] = cnt*isFree1[i];
                        cnt += isFree1[i];
                    }
                }
            }
            if (DOF2>=freeNodeNum*DofPerNode) {
                unsigned iDisp2=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF2)-pDispDOF.begin());
                isFullRestrain2 = isFullRestrain(pDispNodeList[iDisp2].Dir);
                if (!isFullRestrain2) {
                    for (unsigned d=0; d<DofPerNode; d++) {
                        isFree2[d]=!pDispNodeList[iDisp2].Dir[d];
                    }
                    unsigned cnt = 0;
                    for (unsigned i=0; i<DofPerNode; i++) {
                        pos2[i] = cnt*isFree2[i];
                        cnt += isFree2[i];
                    }
                }
            }
            #pragma omp critical
            {
            if (!isFullRestrain1) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t1].coeffRef(tDOF1+pos1[row],DOF1+pos1[col]) += lSM(row,col)*isFree1[row]*isFree1[col];
                    }
                }
                if (!isFullRestrain2) {
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t1].coeffRef(tDOF1+pos1[row],DOF2+pos2[col]) += lSM(row,col+DofPerNode)*isFree1[row]*isFree2[col];
                        }
                    }
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t2].coeffRef(tDOF2+pos2[row],DOF1+pos1[col]) += lSM(row+DofPerNode,col)*isFree2[row]*isFree1[col];
                        }
                    }
                    for (unsigned row=0; row<DofPerNode; row++) {
                        for (unsigned col=0; col<DofPerNode; col++) {
                            bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                        }
                    }

                }
            } else if (!isFullRestrain2) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                    }
                }
            }
            }
        }
        LatTab[LatticeID].factor = _reConRatio;
        LatTab[LatticeID].factorShear = gamma;
    }
    isRepeat = false;
    return true;
}


bool PCG_Eigen::updateSM_shear              (     const unsigned          LatticeID,
                                                  const double            factorS)
{
    std::vector<SpMat>          bSM_temp = bSM;
    static bool isRepeat = false;
    if (isUpdateDOF) {
        if (!isRepeat) {
            dualOut << "The stiffness matrix will be reassembled, no update is required!" << std::endl;
            isRepeat = true;
        }
        return true;
    }
    unsigned NodeID1 = LatTab[LatticeID].nb[0];
    unsigned NodeID2 = LatTab[LatticeID].nb[1];
    if ((GNodeTab[NodeID1].type== Type_Unstable)||(GNodeTab[NodeID2].type== Type_Unstable)) {
        tripOut << "[PCG_Eigen::updateSM]: One of the nodes is unstable, please check!" << std::endl;
        return false;
    }
    unsigned DOF1 = NodeIDtoDOF [NodeID1];
    unsigned DOF2 = NodeIDtoDOF [NodeID2];
    std::array<unsigned,2> pos = getPosDOF (DOF1);
    unsigned t1 = pos[0];
    unsigned tDOF1 = pos[1];
    pos = getPosDOF (DOF2);
    unsigned t2 = pos[0];
    unsigned tDOF2 = pos[1];

    Eigen::MatrixXd lSM(2*DofPerNode,2*DofPerNode);
    if (DofPerNode==6) {
        lSM = getShearLocalSM2(LatticeID,factorS);
    }
    if ((DOF1<freeNodeNum*DofPerNode)&&(DOF2<freeNodeNum*DofPerNode)) {
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t1].coeffRef(tDOF1+row,DOF1+col) += lSM(row,col);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t1].coeffRef(tDOF1+row,DOF2+col) += lSM(row,col+DofPerNode);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t2].coeffRef(tDOF2+row,DOF1+col) += lSM(row+DofPerNode,col);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t2].coeffRef(tDOF2+row,DOF2+col) += lSM(row+DofPerNode,col+DofPerNode);
            }
        }
    } else if ((DOF1<DOF_Num)&&(DOF2<DOF_Num)){
        std::array<unsigned,DofPerNode>    isFree1,isFree2;
        std::array<unsigned,DofPerNode>    pos1,pos2;
        for (unsigned i=0; i<DofPerNode; i++) {
            isFree1[i] = isFree2[i] = true;
            pos1[i] = pos2[i] = i;
        }
        bool isFullRestrain1 = false;
        bool isFullRestrain2 = false;
        if (DOF1>=freeNodeNum*DofPerNode) {
            unsigned iDisp1=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF1)-pDispDOF.begin());
            isFullRestrain1 = isFullRestrain(pDispNodeList[iDisp1].Dir);
            if (!isFullRestrain1) {
                for (unsigned d=0; d<DofPerNode; d++) {
                    isFree1[d]=!pDispNodeList[iDisp1].Dir[d];
                }
                unsigned cnt = 0;
                for (unsigned i=0; i<DofPerNode; i++) {
                    pos1[i] = cnt*isFree1[i];
                    cnt += isFree1[i];
                }
            }
        }
        if (DOF2>=freeNodeNum*DofPerNode) {
            unsigned iDisp2=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF2)-pDispDOF.begin());
            isFullRestrain2 = isFullRestrain(pDispNodeList[iDisp2].Dir);
            if (!isFullRestrain2) {
                for (unsigned d=0; d<DofPerNode; d++) {
                    isFree2[d]=!pDispNodeList[iDisp2].Dir[d];
                }
                unsigned cnt = 0;
                for (unsigned i=0; i<DofPerNode; i++) {
                    pos2[i] = cnt*isFree2[i];
                    cnt += isFree2[i];
                }
            }
        }
        if (!isFullRestrain1) {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+pos1[row],DOF1+pos1[col]) += lSM(row,col)*isFree1[row]*isFree1[col];
                }
            }
            if (!isFullRestrain2) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t1].coeffRef(tDOF1+pos1[row],DOF2+pos2[col]) += lSM(row,col+DofPerNode)*isFree1[row]*isFree2[col];
                    }
                }
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF1+pos1[col]) += lSM(row+DofPerNode,col)*isFree2[row]*isFree1[col];
                    }
                }
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                    }
                }

            }
        } else if (!isFullRestrain2) {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                }
            }
        }
    }
    isRepeat = false;
    return true;
}

bool PCG_Eigen::updateSM_shear              (      const unsigned          LatticeID)
{
    static bool isRepeat = false;
    if (isUpdateDOF) {
        if (!isRepeat) {
            dualOut << "The stiffness matrix will be reassembled, no update is required!" << std::endl;
            isRepeat = true;
        }
        return true;
    }
    unsigned NodeID1 = LatTab[LatticeID].nb[0];
    unsigned NodeID2 = LatTab[LatticeID].nb[1];
    if ((GNodeTab[NodeID1].type== Type_Unstable)||(GNodeTab[NodeID2].type== Type_Unstable)) {
        tripOut << "[PCG_Eigen::updateSM_shear]: One of the nodes is unstable, please check!" << std::endl;
        return false;
    }
    unsigned DOF1 = NodeIDtoDOF [NodeID1];
    unsigned DOF2 = NodeIDtoDOF [NodeID2];
    std::array<unsigned,2> pos = getPosDOF (DOF1);
    unsigned t1 = pos[0];
    unsigned tDOF1 = pos[1];
    pos = getPosDOF (DOF2);
    unsigned t2 = pos[0];
    unsigned tDOF2 = pos[1];

    Eigen::MatrixXd lSM(2*DofPerNode,2*DofPerNode);
    if (KNum==1) {
        tripOut << "[PCG_Eigen::updateSM_shear]: KNum = 1, no shear spring to be removed" << std::endl;
        return false;
    } else if (KNum==2) {
        tripOut << "[PCG_Eigen::updateSM_shear]: KNum = 2, no shear spring to be removed" << std::endl;
        return false;
    } else if (DofPerNode==6) {
        lSM=getFullLocalSM_removeShear(LatticeID);
    }
    if ((DOF1<freeNodeNum*DofPerNode)&&(DOF2<freeNodeNum*DofPerNode)) {
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t1].coeffRef(tDOF1+row,DOF1+col) += lSM(row,col);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t1].coeffRef(tDOF1+row,DOF2+col) += lSM(row,col+DofPerNode);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t2].coeffRef(tDOF2+row,DOF1+col) += lSM(row+DofPerNode,col);
            }
        }
        for (unsigned row=0; row<DofPerNode; row++) {
            for (unsigned col=0; col<DofPerNode; col++) {
                bSM[t2].coeffRef(tDOF2+row,DOF2+col) += lSM(row+DofPerNode,col+DofPerNode);
            }
        }
    } else if ((DOF1<DOF_Num)&&(DOF2<DOF_Num)) {
        std::array<unsigned,DofPerNode>    isFree1,isFree2;
        std::array<unsigned,DofPerNode>    pos1,pos2;
        for (unsigned i=0; i<DofPerNode; i++) {
            isFree1[i] = isFree2[i] = true;
            pos1[i] = pos2[i] = i;
        }
        bool isFullRestrain1 = false;
        bool isFullRestrain2 = false;
        if (DOF1>=freeNodeNum*DofPerNode) {
            unsigned iDisp1=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF1)-pDispDOF.begin());
            isFullRestrain1 = isFullRestrain(pDispNodeList[iDisp1].Dir);
            if (!isFullRestrain1) {
                for (unsigned d=0; d<DofPerNode; d++) {
                    isFree1[d]=!pDispNodeList[iDisp1].Dir[d];
                }
                unsigned cnt = 0;
                for (unsigned i=0; i<DofPerNode; i++) {
                    pos1[i] = cnt*isFree1[i];
                    cnt += isFree1[i];
                }
            }
        }
        if (DOF2>=freeNodeNum*DofPerNode) {
            unsigned iDisp2=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF2)-pDispDOF.begin());
            isFullRestrain2 = isFullRestrain(pDispNodeList[iDisp2].Dir);
            if (!isFullRestrain2) {
                unsigned iDisp2=(std::lower_bound(pDispDOF.begin(),pDispDOF.end(),DOF2)-pDispDOF.begin());
                for (unsigned d=0; d<DofPerNode; d++) {
                    isFree2[d]=!pDispNodeList[iDisp2].Dir[d];
                }
                unsigned cnt = 0;
                for (unsigned i=0; i<DofPerNode; i++) {
                    pos2[i] = cnt*isFree2[i];
                    cnt += isFree2[i];
                }
            }
        }
        if (!isFullRestrain1) {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t1].coeffRef(tDOF1+pos1[row],DOF1+pos1[col]) += lSM(row,col)*isFree1[row]*isFree1[col];
                }
            }
            if (!isFullRestrain2) {
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t1].coeffRef(tDOF1+pos1[row],DOF2+pos2[col]) += lSM(row,col+DofPerNode)*isFree1[row]*isFree2[col];
                    }
                }
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF1+pos1[col]) += lSM(row+DofPerNode,col)*isFree2[row]*isFree1[col];
                    }
                }
                for (unsigned row=0; row<DofPerNode; row++) {
                    for (unsigned col=0; col<DofPerNode; col++) {
                        bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                    }
                }

            }
        } else if (!isFullRestrain2) {
            for (unsigned row=0; row<DofPerNode; row++) {
                for (unsigned col=0; col<DofPerNode; col++) {
                    bSM[t2].coeffRef(tDOF2+pos2[row],DOF2+pos2[col]) += lSM(row+DofPerNode,col+DofPerNode)*isFree2[row]*isFree2[col];
                }
            }
        }
    }
    isRepeat = false;
    return true;
}

void PCG_Eigen::getPreCondJacobi			(	Eigen::VectorXd*	p_vm)
{
	if (DOF_Num != SM.rows()) {
        tripOut << " Stiffness Matrix size is not correct!" << endl;
        tripOut << "DOF_Num = " << DOF_Num << " SM.rows() = " << SM.rows() << std::endl;
    }
    p_vm -> resize(DOF_Num);
    for (unsigned i=0; i<DOF_Num; i++) {
    	if (std::fabs(SM.coeff(i,i))<Tiny) {
    	   std::cout << "ill-conditioned stiffness matrix! " << i << "\n";
    	}
        (*p_vm)[i] = 1.0/SM.coeff(i,i);
        if ((*p_vm)[i] !=(*p_vm)[i]) {
        	std::cout << "vector m contain -nan!" << std::endl;
        }
    }
}

void PCG_Eigen::getPreCondJacobiParallel	(	Eigen::VectorXd*	p_vm)
{
    unsigned sumDOF = 0;
    for (unsigned t=0; t<chunkNum; t++) {
    	sumDOF += bSM[t].rows();
    }
    if ((DOF_Num != sumDOF)||(DOF_Num != bSM[0].cols())) {
        dualOut << " Stiffness Matrix size is not correct!" << std::endl;
        tripOut << "DOF_Num = " << DOF_Num << " sumDOF = " << sumDOF << " bSM[0].cols() = " << bSM[0].cols() << std::endl;
    }
    p_vm -> resize(DOF_Num);
    unsigned shift = 0;
    unsigned cnt = 0, cnt1 = 0;
    for (unsigned t=0; t<chunkNum; t++) {
    	for (unsigned i=0; i<bSM[t].rows(); i++) {
    	    if (std::fabs(bSM[t].coeff(i,i+shift))<Tiny) {
    	        std::cout << "ill-conditioned stiffness matrix! " << t << ' ' << i+shift << "\n";
    	        (*p_vm)[i+shift] = 1.0;
    	        cnt++;
    	    }
    		(*p_vm)[i+shift] = 1.0/(bSM[t].coeff(i,i+shift));
    		if ((*p_vm)[i+shift] !=(*p_vm)[i+shift]) {
    		    tripOut << "vector m contain -nan! bSM[t].coeff(i,i+shift) = " << bSM[t].coeff(i,i+shift) << std::endl;
    		    unsigned NodeID = DOFtoNodeID[i+shift];
    		    if (NodeID!=tNodeNum) {
    		        cnt1++;
    		    }
    		}
    	}
    	shift += bSM[t].rows();
    }
    if (cnt1) {
        tripOut << "Problematic node number = " << cnt1 << std::endl;
    }

    for (unsigned i=0; i<conList.size(); i++) {
        (*p_vm)[conList[i].first] = 0.0;
    }
}

void PCG_Eigen::getPreCondIdentity			(	Eigen::VectorXd*	p_vm)
{
    *p_vm = Eigen::VectorXd::Ones(DOF_Num);
}

void PCG_Eigen::setRestrainZero				(	Eigen::VectorXd		*p_vq)
{
    unsigned i,iDOF;
    #pragma omp parallel for private (iDOF,i)
    for (iDOF = freeNodeNum; iDOF < calNodeNum; iDOF++)
    {
        for (i=0; i<DofPerNode; i++)
        {
            if (pDispNodeList[iDOF-freeNodeNum].Dir[i])
            {
                (*p_vq)[iDOF*DofPerNode+i] = 0.0;
            }
        }
    }
}

void PCG_Eigen::setConstrain             (   Eigen::VectorXd     *p_vq)
{
    for (unsigned i=0; i < conList.size(); i++) {
        (*p_vq)[conList[i].first] = 0.0;
    }
}

void PCG_Eigen::updateNodeStateSimple	(	const unsigned			step,
                                            const bool 				isRelax)
{
    unsigned i,k,NodeID,DOF;

    if (isRelax)
    {
        for (i = 0; i < calNodeNum; i++)
        {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (k = 0; k< DofPerNode ; k++)
            {
                GNodeTab[NodeID].d[k] = vx[DOF+k];
            }

        }
    }
    else
    {
        for (i = 0; i < calNodeNum; i++)
        {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (k = 0; k< DofPerNode ; k++)
            {
                GNodeTab[NodeID].d[k] += vx[DOF+k];			//For spring model
            }
        }
    }

    updateOldDeltaSimple();

}

void PCG_Eigen::updateNodeStateCompact	(	const unsigned			step,
                                            const bool 				isRelax)
{
    unsigned d,NodeID,DOF;

    if (isRelax){
        for (unsigned i = 0; i < freeNodeNum; i++){
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (unsigned k = 0; k< DofPerNode ; k++){
                GNodeTab[NodeID].d[k] = vx[DOF+k];
            }

        }

        for (unsigned i = freeNodeNum; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            d = 0;
            for (unsigned k = 0; k< DofPerNode ; k++) {
                if (!pDispNodeList[i-freeNodeNum].Dir[k]) {
                    GNodeTab[NodeID].d[k] = vx[DOF+d];
                    d++;
                } else {
                    GNodeTab[NodeID].d[k] = 0.0;
                }
            }
        }

    } else {
        for (unsigned i = 0; i < freeNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (unsigned k = 0; k< DofPerNode ; k++) {
                GNodeTab[NodeID].d[k] += vx[DOF+k];			//For spring model
            }
        }

        for (unsigned i = freeNodeNum; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            d = 0;
            for (unsigned k = 0; k< DofPerNode ; k++) {
                if (!pDispNodeList[i-freeNodeNum].Dir[k]) {
                    GNodeTab[NodeID].d[k] += vx[DOF+d];			//For spring model
                    d++;
                } else {
                    GNodeTab[NodeID].d[k] = 0.0;
                }
            }
        }
    }

    updateOldDeltaCompact();

}

void PCG_Eigen::updateConstrainForce    ()
{
    std::cout << "updateConstrainForce... " << std::endl;
    if (nodeLists.constrain.size()!=0) {
        #pragma omp parallel for
        for (unsigned t=0; t<chunkNum; t++) {
            vr.segment(chunkPos[t],chunkNumList[t]) = bSM[t]*vx;
        }

        for (unsigned i=0; i<nodeLists.constrain.size(); i++) {
            unsigned NodeID = nodeLists.constrain[i].NodeID;
            std::vector<Restrain>::iterator it = std::find(pDispNodeList.begin(),pDispNodeList.end(),NodeID);
            if (it==pDispNodeList.end()) {
                for (unsigned k=0; k<nodeLists.constrain[i].conDisp.size(); k++) {
                    GNodeTab[NodeID].extF[nodeLists.constrain[i].conDir[k]] = vr[NodeIDtoDOF[NodeID]+nodeLists.constrain[i].conDir[k]];
                }
            } else {
                for (unsigned k=0; k<nodeLists.constrain[i].conDisp.size(); k++) {
                    unsigned pos = 0;
                    for (unsigned l=0; l<nodeLists.constrain[i].conDir[k]; l++) {
                        pos += !pDispNodeList[it-pDispNodeList.begin()].Dir[l];
                    }
                    GNodeTab[NodeID].extF[nodeLists.constrain[i].conDir[k]] = vr[NodeIDtoDOF[NodeID]+pos];
                }
            }
        }
    }
}


std::vector<unsigned> PCG_Eigen::testSM  ()
{
    static unsigned cnt = 0;
    std::vector<unsigned> latList;
    std::vector<SpMat>          bSM_temp = bSM;
    fillSMblock_full3();
    dualOut << "[PCG_Eigen::testSM]..." << std::endl;
    bool isPass = true;
    for (unsigned t=0; t<chunkNum; t++) {
        if (bSM[t].isApprox(bSM_temp[t],1e-5)) {
            dualOut << "***********Pass " << t << std::endl;
        } else {
            char fileName[255];
            sprintf(fileName,"%s/bSM_Check_a%04d.vtk",OutputSubFolder,cnt);
            ofstream file1(fileName, ios::out);
            tripOut << "***********Fail " << t << ' ' << (bSM[t]-bSM_temp[t]).norm() << std::endl;
            SpMat          bSM_diff = bSM_temp[t] - bSM[t];
            file1 << bSM_diff << std::endl;
            for (int k=0; k<bSM_diff.outerSize(); ++k) {
                for (SpMat::InnerIterator it(bSM_diff,k); it; ++it) {
                    if (std::fabs(it.value())>1e-5) {
                        unsigned DOF1 =  it.row()+chunkPos[t];
                        unsigned DOF2 =  it.col();
                        if ((DOFtoNodeID[DOF1]<tNodeNum)&&(DOFtoNodeID[DOF2]<tNodeNum)
                                &&(DOFtoNodeID[DOF1]!=DOFtoNodeID[DOF2])) {
                            unsigned NodeID1 = DOFtoNodeID[DOF1];
                            tripOut << NodeID1 << ',';
//                            GNodeTab[NodeID1].type = 999;
                            unsigned NodeID2 = DOFtoNodeID[DOF2];
                            tripOut << NodeID2 << ' ';
//                            GNodeTab[NodeID2].type = 999;
                            for (unsigned n=0; n<GNodeTab[NodeID1].nbLatticeID.size(); n++) {
                                if (GNodeTab[NodeID1].nbNodeID[n]==NodeID2) {
                                    latList.push_back(GNodeTab[NodeID1].nbLatticeID[n]);
                                    dualOut << GNodeTab[NodeID1].nbLatticeID[n] << ' ' << std::endl;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            isPass = false;
        }
    }
    if (!isPass) {
        cnt++;
    }
    return latList;
}

std::vector<unsigned> PCG_Eigen::checkAllLatticeSM()
{
    dualOut << "[checkAllLatticeSM]..." << std::endl;
    static unsigned cnt = 0;
    std::vector<unsigned> latList;
    testSM();
    for (unsigned LatticeID = 0; LatticeID<tLatticeNum; LatticeID++) {
        double factor = LatTab[LatticeID].factor/2.0;
        double factorS = LatTab[LatticeID].factorShear/2.0;
        updateSM_full (LatticeID,factor,factorS);
        std::vector<SpMat>          bSM_temp = bSM;
        Eigen::MatrixXd lSM = getFullLocalSM2(LatticeID,factor,factorS);
        Eigen::MatrixXd lSM1 = getFullLocalSM2(LatticeID);
        LatTab[LatticeID].factor = factor;
        LatTab[LatticeID].factorShear = factorS;
        fillSMblock_full4();
        for (unsigned t=0; t<chunkNum; t++) {
            if (!bSM[t].isApprox(bSM_temp[t],1e-5)) {
                dualOut << LatticeID << ' ' << LatTab[LatticeID].factor << ' ' << LatTab[LatticeID].factorShear << std::endl;
                dualOut << LatTab[LatticeID].nb[0] << ' ' << LatTab[LatticeID].nb[1] << std::endl;
                SpMat          bSM_diff = bSM[t]-bSM_temp[t];
                dualOut << "***********Fail " << t << ' ' << bSM_diff.norm() << std::endl;
                char fileName[255];
                sprintf(fileName,"%s/bSM_Check_a%04d.vtk",OutputSubFolder,cnt);
                ofstream file1(fileName, ios::out |ios::app);
                std::vector<SpMat>          bSM_temp = bSM;
                Eigen::MatrixXd lSM2 = getFullLocalSM2(LatticeID);
                file1 << t << ' ' << LatticeID << ' ' << LatTab[LatticeID].factor << ' ' << LatTab[LatticeID].factorShear << std::endl;
                file1 << LatTab[LatticeID].nb[0] << ' ' << LatTab[LatticeID].nb[1] << ' ' <<  bSM_diff.norm() << std::endl;
                file1 << bSM_diff << std::endl;
                file1 << std::endl <<"-------------------------" << std::endl;
                file1 << lSM << std::endl;
                file1 << std::endl <<"-------------------------" << std::endl;
                file1 << lSM1 << std::endl;
                file1 << std::endl <<"-------------------------" << std::endl;
                file1 << lSM2 << std::endl;
                file1 << std::endl <<"-------------------------" << std::endl;
                file1 << lSM1 +lSM - lSM2 << std::endl;
                file1 << std::endl <<"*************************" << std::endl;
                cnt++;
                latList.push_back(LatticeID);
            }
        }
    }
    for (unsigned LatticeID = 0; LatticeID<tLatticeNum; LatticeID++) {
        LatTab[LatticeID].factor *= 2.0;
        LatTab[LatticeID].factorShear *= 2.0;
    }
    fillSMblock_full3();
    dualOut << "Abnormal lattice no = " << latList.size() << std::endl;
    return latList;
}

void PCG_Eigen::updateOldDeltaSimple()
{
    unsigned NodeID,DOF;

    oldDelta_d.clear();
    oldDelta_d.resize(local_tNodeNum*DofPerNode,0.0);

    for (unsigned i=0; i < calNodeNum; i++) {
        NodeID = calNodeList[i];
        DOF = NodeIDtoDOF[NodeID];
        for (unsigned k = 0; k< DofPerNode ; k++) {
            oldDelta_d[NodeID*DofPerNode+k] = vx[DOF+k];
        }
    }
}

void PCG_Eigen::updateOldDeltaCompact()
{
    unsigned i,NodeID,DOF;

    oldDelta_d.clear();
    oldDelta_d.resize(local_tNodeNum*DofPerNode,0.0);

    for (i=0; i < freeNodeNum; i++) {
        NodeID = calNodeList[i];
        DOF = NodeIDtoDOF[NodeID];
        for (unsigned k = 0; k< DofPerNode ; k++) {
            oldDelta_d[NodeID*DofPerNode+k] = vx[DOF+k];
        }
    }

    for (i=freeNodeNum; i < calNodeNum; i++) {
        NodeID = calNodeList[i];
        DOF = NodeIDtoDOF[NodeID];
        unsigned d = 0;
        for (unsigned k = 0; k< DofPerNode ; k++) {
            if (!pDispNodeList[i-freeNodeNum].Dir[k]) {
                oldDelta_d[NodeID*DofPerNode+k] = vx[DOF+d];
                d++;
            } else {
                oldDelta_d[NodeID*DofPerNode+k] = 0.0;
            }
        }
    }
}

void PCG_Eigen::clean ()
{
    vr.resize(0);
    vx.resize(0);

    SM.resize(0,0);
    bSM.clear();

    freeNodeList.clear();
    pDispNodeList.clear();
    calNodeList.clear();
    NodeIDtoDOF.clear();
    pDispDOF.clear();
}
