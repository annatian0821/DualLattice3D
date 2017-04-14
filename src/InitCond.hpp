/*
 * InitCond.hpp
 *
 *  Created on: Oct 12, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */
#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "Loading.hpp"
#include "GLatRemoval.hpp"
#include "MiscTools.hpp"
#include "ConnectLat.hpp"
#include "PCG_Eigen.hpp"
#include "ProProcess.hpp"
#include "Output.hpp"
#include <random>

#ifndef INITCOND_HPP_
#define INITCOND_HPP_

class InitCond
{
public:
	InitCond				() {};

	~InitCond				() {};

	InitCond				(	GLatRemoval*					p_glatRemoval);

    void setInitCond		(	ConnectLat*						p_conLat,
            					Vertices*						p_verti,
            					PCG_Eigen*						p_EigenPCG,
            					std::vector<unsigned>*			p_pxCrackLatList,
            					const std::array<unsigned,2>	pxCrackPerimeterNum,
            					const unsigned					maxCoordNo);

    void setInitCond		(	ConnectLat*						p_conLat,
                                GLatRemoval*                    p_glatRemoval,
            					Vertices*						p_verti,
            					FluidCal*						p_fCal,
            					PCG_Eigen*						p_EigenPCG,
            					std::vector<unsigned>*			p_pxCrackLatList,
								const std::array<unsigned,2>	pxCrackPerimeterNum,
            					const double					initP,
            					const unsigned					maxCoordNo);

    unsigned stableLat 		(	ConnectLat*						p_conLat,
        						Vertices*						p_verti,
        						PCG_Eigen*						p_EigenPCG,
        						bool&							isStable,
        						unsigned&						relaxStep,
        						unsigned&						totalStep,
        						const unsigned					maxIter,
        						const unsigned					minLatTolerance);

    void    setBoreHoleNodeList (   std::vector<unsigned>*      p_latList);

    void    setInitCrackAperture (  std::vector<unsigned>*      p_pxCrackLatList);

private:
    GLatRemoval*				p_glatRemoval;
    std::vector<unsigned>       boreHoleNodeList;
    void initNodeAndLatState	();

    void applyBoundaryLoading	();

    void applyBoundaryStress	( 	const Tensor6		initStress);

    void applyinitialCondition	(	const Tensor6					initStress,
    								Vertices*						p_verti,
									ConnectLat*						p_conLat,
									PCG_Eigen*						p_EigenPCG,
    								const unsigned					maxCoordNo);

    void applyBoundaryStressFree( 	const Tensor6					initStress);

    void applyMinRestrain ();

    void setBoreHoleConstrain   ();

    bool stablizeInitModel 		(	Vertices*						p_verti,
    								ConnectLat*						p_conLat,
    								PCG_Eigen*						p_EigenPCG,
									const unsigned					maxCoordNo,
									const unsigned					maxBreakNum,
    								const unsigned					maxStep);

    void genGPreExRegCrack		(	const double 				xc,
                                    const double 				yc,
                                    const double 				zc,
                                    const double 				a,
                                    const double 				b,
                                    const double 				c,
                                    ConnectLat*					p_conLat,
                                    Vertices*					p_verti);

    void genPriInjectWell		(	Point							centre,
                                    const double					r,
                                    const double					injectRate,
                                    const double					pressure,
                                    ConnectLat*						p_conLat,
                                    Vertices*						p_verti,
                                    FluidCal*						p_fCal);

    unsigned genPriInjectWell	(	const std::vector<unsigned>*	p_pxCrackLatList,
                                    const double					injectRate,
                                    const double					pressure,
                                    ConnectLat*						p_conLat,
                                    GLatRemoval*                    p_glatRemoval,
                                    Vertices*						p_verti,
                                    FluidCal*						p_fCal);

    void genPriInjectWellSingle	(	const unsigned					LatticeID,
                                    const double					injectRate,
                                    const double					pressure,
                                    ConnectLat*						p_conLat,
                                    GLatRemoval*                    p_glatRemoval,
                                    Vertices*						p_verti,
                                    FluidCal*						p_fCal);

    void applyinitialConstrain  (   const unsigned                  face,
                                    const unsigned                  dir,
                                    const double                    disp);

    void applyUnbreakableNb		(	const std::vector<unsigned>*	p_pxCrackLatList);

    void applyUnbreakableNbAll	(	const std::vector<unsigned>*	p_pxCrackLatList,
									const std::array<unsigned,2>	pxCrackPerimeterNum);

    void modifyTensileStrain	(	std::vector<unsigned>*			p_latList,
    								double							value);

    void applyUnbreakableOuterLat	();

    void applySmallDistortPxCrack(	const std::vector<unsigned>*   p_pxCrackLatList);

    void applyBoundaryRestrain	( 	const RestrainMat			   restrainMat);

    RestrainMat  getRestrainMat	( 	const Tensor6				   initStress);

    void rmPxCrack              (   ConnectLat*                    p_conLat,
                                    Vertices*                      p_verti,
                                    PCG_Eigen*                     p_EigenPCG,
                                    std::vector<unsigned>*         p_pxCrackLatList);

    void modifyLatStiffnessAndThreshold(const unsigned LatticeID, const double ratio);

    void modifyLatStiffness(const unsigned LatticeID, const double ratio);

    void modifyLatThreshold(const unsigned LatticeID, const double ratio);

    void setZero				();

    void outputBoundaryRestrain    ();

    bool isConstrain			(	const std::array<bool,DofPerNode>	bConstrain	)
    {
        bool flag = false;
        for (unsigned i=0; i<DofPerNode; i++) {
            flag = flag || bConstrain[i];
            if (flag) return flag;
        }
    	return flag;
    }

};

#endif /* INITCOND_HPP_ */
