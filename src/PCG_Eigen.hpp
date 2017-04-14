/*
 * PCG_Eigen.hpp
 *
 *  Created on: Feb 15, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef PCG_EIGEN_HPP_
#define PCG_EIGEN_HPP_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "MiscTools.hpp"
//#include "Output.hpp"
//#include "Testing.hpp"
#include "GConjGrad.hpp"

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

typedef Eigen::SparseMatrix<double> 	SpMat;
typedef Eigen::MatrixXf                 DenseMat;
typedef Eigen::Triplet<double> 			Triplet;

class PCG_Eigen : public PreCondConjGrad {
	friend class Stress;
	friend class EnergyCal;
	friend class GLatRemoval;
public:
	PCG_Eigen ();
    PCG_Eigen(	unsigned 	_maxCoordNo,
    			unsigned	_nodeNumLimit)			{ 	maxCoordNo = _maxCoordNo;
    													nodeNumLimit = _nodeNumLimit;
    													isUpdateDOF = true;
    													freeLocation = 0;
    													cnt_PCG = 0;};
    virtual ~PCG_Eigen();
    void updateDOF					(	const bool 						_isUpdateDOF) {
    	isUpdateDOF = _isUpdateDOF;
    };

    bool solve											(	const unsigned 					step,
                                        					const bool 						isRelax,
                                        					const bool						useCompSM);

    bool updateSM_full              					(   const unsigned                  LatticeID,
                                        					const double                    factor);

    bool updateSM_full                                  (   const unsigned                  LatticeID,
                                                            const double                    factor,
                                                            const double                    factorS);

    bool updateSM_full                                  (   std::vector<unsigned>&          latList,
                                                            const double                    _reConRatio,
                                                            const double                    gamma);

    bool updateSM_full                                  (   std::vector<unsigned>&          latList,
                                                            const double                    _reConRatio,
                                                            std::vector<double>&            gammaList);

    bool updateSM_shear                                 (   const unsigned                  LatticeID);

    bool updateSM_shear                                 (   const unsigned                  LatticeID,
                                                            const double                    factorS);

    Tensor getMemForce                     				(   const unsigned                  LatticeID);

    Tensor getResidual			                        (	const unsigned					NodeID);

    std::vector<unsigned>   testSM                      ();

    std::vector<unsigned>   checkAllLatticeSM           ();

private:
    bool						isUpdateDOF;
    unsigned                    cnt_PCG;
    unsigned					maxCoordNo;
    unsigned					nodeNumLimit;
    unsigned					chunkNum,chunkSize,chunkSizeEnd;
    std::vector<unsigned>       chunkNumList;
    std::vector<unsigned>       chunkPos;
    unsigned                    freeLocation;
    double						SMentryTolerance;
    std::vector<unsigned>       conDOFList;
    Eigen::VectorXd 			vr;
    Eigen::VectorXd 			vx;

    SpMat						SM;
    std::vector<SpMat>			bSM;


    bool integralSolve2				(	const unsigned				minIter,
                                        const unsigned				maxIter,
                                        const double 				precision,
                                        const bool					useCompSM);

    bool integralSolveSerial 		(	const unsigned				minIter,
                                        const unsigned				maxIter,
                                        const double 				precision);

    bool integralSolveParallel 		(	const unsigned				minIter,
                                        const unsigned				maxIter,
                                        const double 				precision,
                                        const bool					useCompSM);

    bool integralSolveParallel2     (   const unsigned              minIter,
                                        const unsigned              maxIter,
                                        const double                precision,
                                        const bool                  useCompSM);

    bool useEigenCGSolver			(	unsigned					maxIteration,
                                        double						precision);

    bool sparseCholeskySolve 		();

    bool useBiCGSTAB_Solver			(	unsigned					maxIteration,
                                        double						precision);

    void testSolvers				(	const unsigned				minIter,
                                        const unsigned				maxIter,
                                        const double 				precision,
                                        const bool					useCompSM);

    void initialize					(	const bool 					isRelax,
                                        const bool					useCompSM);

    unsigned getMaxCoordNum			( 	const unsigned 				shift,
    									const unsigned				endDOF);

    void fillChunkNumList           ();

    std::array<unsigned,2> getPosDOF(  const unsigned               DOF);

    void fillSMSimple_Full			();

    void fillSMblock_full2          ();

    void fillSMblock_full3          ();

    void fillSMblock_full4          ();

    Eigen::Matrix3d getSpringLocalSM    (   const unsigned              NodeID,
                                            const unsigned              dir);

    Eigen::Matrix3d getSpringLocalSM    (   const unsigned              LatticeID,
                                            const double                factor);

    Eigen::Matrix3d getSpringLocalSM    (   const unsigned      		LatticeID);

    Eigen::Matrix3d getShearLocalSM     (   const unsigned              NodeID,
                                            const unsigned              dir);

    Eigen::Matrix3d getShearLocalSM     (   const unsigned              LatticeID,
                                            const double                factor);

    Eigen::Matrix3d getShearLocalSM   	(   const unsigned      		LatticeID);

    Eigen::MatrixXd getFullLocalSM2     (   const unsigned          NodeID,
                                            const unsigned          dir);

    Eigen::MatrixXd getFullLocalSM2     (   const unsigned          LatticeID,
                                            const double            factor);

    Eigen::MatrixXd getFullLocalSM2     (   const unsigned          LatticeID,
                                            const double            factor,
                                            const double            factorS);

    Eigen::MatrixXd getFullLocalSM2     (   const unsigned          LatticeID);

    Eigen::MatrixXd getFullLocalSM2     (   const unsigned          NodeID,
                                            const unsigned          dir,
                                            const double            factor);

    Eigen::MatrixXd getShearLocalSM2   (    const unsigned          LatticeID,
                                            const double            factorS);

    Eigen::MatrixXd getFullLocalSM_removeShear   (     const unsigned      LatticeID);

    Eigen::MatrixXd getBeamSM   		(   const unsigned      	NodeID,
                                            const unsigned      	dir);

    Eigen::MatrixXd getBeamSM    		(   const unsigned      	LatticeID);

    Eigen::MatrixXd getBeamSM    		(   const unsigned      	LatticeID,
    										const double			factor);

    Tensor getMemForce					(	const unsigned			NodeID,
    										const unsigned			k);

    void getPreCondJacobi			(	Eigen::VectorXd*			p_vm);

    void getPreCondJacobiParallel	(	Eigen::VectorXd*			p_vm);

    void getPreCondIdentity			(	Eigen::VectorXd*			p_vm);

    double parallelDot				(	Eigen::VectorXd				*p_v1,
    									Eigen::VectorXd				*p_v2);

    double parallelDot2             (   Eigen::VectorXd             *p_v1,
                                        Eigen::VectorXd             *p_v2);

    void initGradVec_r_Simple		(   const bool					isRelax);

    void initGradVec_r_Compact		(	const bool					isRelax);

    void initDispVec_x_Simple		(	const bool 					isRelax);

    void initDispVec_x_Compact		(	const bool 					isRelax);

    void setRestrainZero			(	Eigen::VectorXd*			p_vq);

    void setConstrain               (   Eigen::VectorXd*            p_vq);

    void updateNodeStateSimple		(	const unsigned				step,
                                        const bool 					isRelax);

    void updateNodeStateCompact		(	const unsigned				step,
                                        const bool 					isRelax);
    void updateConstrainForce       ();

    void updateOldDeltaSimple       ();

    void updateOldDeltaCompact      ();

    void clean                      ();
};

#endif /* PCG_EIGEN_HPP_ */
