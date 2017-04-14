/*
 * GConjGrad.hpp
 *
 *  Created on: Nov 30, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef GCONJGRAD_HPP_
#define GCONJGRAD_HPP_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "MiscTools.hpp"
#include "ProProcess.hpp"
//#include "Output.hpp"
//#include "Testing.hpp"
#include "SparseMatrix.hpp"

#include <eigen3/Eigen/Sparse>




class PreCondConjGrad {
public:
    PreCondConjGrad();
    virtual ~PreCondConjGrad();

    virtual bool solve		(	const unsigned 					step,
                                const bool 						isRelax,
                                const bool						useCompSM);

protected:
    std::vector<double>				oldDelta_d;				//Index arranged by NodeID
    std::vector<unsigned>			freeNodeList;
    std::vector<Restrain>			pDispNodeList;
    std::vector<unsigned>			calNodeList;
    std::vector<unsigned>			NodeIDtoDOF;
    std::vector<unsigned>           DOFtoNodeID;
    std::vector<unsigned>			pDispDOF;
    std::vector<IDAndInfo>          conList;

    unsigned						freeNodeNum;
    unsigned						pDispNodeNum;
    unsigned						calNodeNum;
    unsigned						dumDOF;
    unsigned						DOF_Num;

    unsigned						local_tNodeNum, local_ttvNodeNum;

    void adjSMCal			( 	const unsigned					NodeID,
                                const unsigned 					dir,
                                double*							p_C,
                                double*							p_Cx,
                                double*							p_Cy,
                                double*							p_Cz);
    void adjSMCal			( 	const unsigned					NodeID,
                                const unsigned 					nbNodeID,
                                const double					keChange,
                                double*							p_C,
                                double*							p_Cx,
                                double*							p_Cy,
                                double*							p_Cz);

    virtual void    updateDOFSimple();

    virtual void    updateDOFCompact();

    virtual void    fillConList();

    virtual void    simpleUpdateConList         (   const double            disp);


    virtual void 	initDispVec_x_Simple		(	const bool 				isRelax);

    virtual void 	initGradVec_r_Simple		(	const bool						isRelax);

    virtual void 	fillCompSM					(	StiffnessMatrix*			p_SM);

    virtual void 	fillSMSimple				(	StiffnessMatrix*			p_SM);

    virtual void 	updateNodeStateSimple		(	const unsigned			step,
            										const bool 				isRelax);

    virtual void 	updateNodeStateCompact		(	const unsigned			step,
            										const bool 				isRelax);

    virtual void 	updateOldDeltaSimple();

    virtual void 	updateOldDeltaCompact();

    unsigned getFreeDOFNum	(	const std::array<bool,DofPerNode>	Restrain) {
        unsigned cnt = 0;
        for (unsigned i=0; i<DofPerNode; i++) {
            if (!Restrain[i]) {
                cnt++;
            }
        }
        return cnt;
    }

    bool isFullRestrain	(	const std::array<bool,DofPerNode>	Restrain) {
        bool flag = true;
        for (unsigned i=0; i<DofPerNode; i++) {
            flag = flag && Restrain[i];
        }
        return flag;
    }

private:

    std::vector<double>				m;	//Preconditional vector
    std::vector<double>				s;
    std::vector<double>				r;
    std::vector<double>				x;
    std::vector<double>				q;
    std::vector<double>				d;


    unsigned 						i_preCalMat[Dim][Dim];

    void fillGPreCalMat		();

    void genPreCondJacobi 	(	StiffnessMatrix*	p_SM);

    void genPreCondIC		(	StiffnessMatrix*	p_SM,
                                Preconditioner*		p_PC);



    void simpleSMCal		( 	const unsigned					NodeID,
                                const unsigned 					dir,
                                double*							p_C,
                                double*							p_Cx,
                                double*							p_Cy,
                                double*							p_Cz);

    void 	ginitialize		(	StiffnessMatrix*				p_SM,
                                const bool 						isRelax);

    void initGradVec_d();

    void integralSolve 		(	double*							p_newDelta,
                                StiffnessMatrix*				p_SM);

    void moduleSolve	 	(	double*							p_newDelta);

    void gupdateVec_q       ();

    void newUpdateVec_q 	(	StiffnessMatrix*				p_SM);

    void gupdateVec_s       ();

    void gupdateVec_d 		(	double						beta);

    void gupdateVecs_xrs	(	double						alpha);

    void tmpUpdateNodeState		(	const bool 				isRelax);

    void undoTmpUpdateNodeState	(	const bool 					isRelax);

    unsigned	getSizeOfData	();

    void 	QuickMatMulti	(	StiffnessMatrix*			p_SM);

    void PCGclean ();

    double dot (std::vector<double> *p_v1,std::vector<double> *p_v2);

};


#endif /* GCONJGRAD_HPP_ */
