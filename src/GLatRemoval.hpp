/*
 * GLatRemoval.h
 *
 *  Created on: Nov 29, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef GLATREMOVAL_H_
#define GLATREMOVAL_H_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "Statistics.hpp"
#include "ConnectLat.hpp"
#include "ProProcess.hpp"
#include "MiscTools.hpp"
#include "PCG_Eigen.hpp"

/*
bool operator== (const iDouble &n1, const iDouble &n2)
{
	return n1.data == n2.data;
}
*/

inline bool reverseSortDir (double ii, double jj) { return (ii>jj); }

class GLatRemoval {
	friend class ConnectLat;
public:
    GLatRemoval					();

    GLatRemoval					(	const double	_chkCoplanarTol);
    virtual ~GLatRemoval();

    void setRecord				(	bool inBool ) {
        isRecord = inBool;
    };

    unsigned remove				(	ConnectLat*					p_conLat,
                                    Vertices*					p_verti,
                                    FailStats*					p_fStats,
                                    PCG_Eigen*					p_EigenPCG,
                                    EnergyCal*                  p_eCal,
                                    const unsigned				step,
                                    const double				time,
                                    double* 					maxLatNormForce,
                                    const unsigned 				maxBreak);

    unsigned remove				(	ConnectLat*					p_conLat,
                                    Vertices*					p_verti,
                                    PCG_Eigen*					p_EigenPCG,
                                    const unsigned				maxBreak,
                                    const unsigned 				step,
                                    const double				time);

    void removeUnstableNode		(	const unsigned 				NodeID,
    								const unsigned				step,
    								const double				time,
                                    ConnectLat*					p_conLat,
                                    Vertices*					p_verti,
                                    PCG_Eigen*					p_EigenPCG,
                                    const std::vector<unsigned> rmDirList);

    void removeUnstableNode		(	const unsigned 				NodeID,
            						NodeLists*					p_nodeLists,
                                    std::vector<unsigned>*		p_latList);

    void removeLattice			(	const unsigned 				LatticeID);

    void removeLattice			(	const unsigned 				LatticeID,
                                    ConnectLat*					p_conLat,
                                    Vertices*					p_verti,
                                    FailStats*					p_fStats,
                                    PCG_Eigen*					p_EigenPCG,
                                    EnergyCal*                  p_eCal,
                                    const unsigned				step,
                                    const double				time);

    void removeLattice			(	const unsigned 				LatticeID,
                                    ConnectLat*					p_conLat,
                                    Vertices*					p_verti,
                                    PCG_Eigen*					p_EigenPCG,
                                    const unsigned				step,
                                    const double				time);

    void removeLatShear         (   const unsigned              LatticeID,
                                    const double                factorS,
                                    const double                fShear,
                                    ConnectLat*                 p_conLat,
                                    Vertices*                   p_verti,
                                    FailStats*                  p_fStats,
                                    PCG_Eigen*                  p_EigenPCG,
                                    EnergyCal*                  p_eCal,
                                    const unsigned              step,
                                    const double                time);

    double getMaxLatNormForce 	(	unsigned 					maxBreak);

    unsigned reConnectLat		(	ConnectLat*					p_conLat,
    								PCG_Eigen*					p_EigenPCG,
    								const unsigned              maxNo);

    unsigned disConnectLat		(	Vertices*					p_verti,
        							ConnectLat*					p_conLat,
        							PCG_Eigen*					p_EigenPCG,
        							const unsigned				step,
        							const double				time,
                                    const unsigned              maxNo);

    void disConnectLat			(	Vertices*					p_verti,
    								ConnectLat*					p_conLat,
    								PCG_Eigen*					p_EigenPCG,
    								const unsigned				step,
    								const double				time,
    								const std::vector<unsigned>*p_latList);

    bool updateReConLatShearStiffness  (    ConnectLat*         p_conLat,
                                            PCG_Eigen*          p_EigenPCG,
                                            Vertices*           p_verti,
                                            const unsigned      step,
                                            const double        time);

    bool updateReConLatShearStiffness  (    ConnectLat*         p_conLat,
                                            PCG_Eigen*          p_EigenPCG,
                                            Vertices*           p_verti,
                                            const double        gamma,
                                            const unsigned      step,
                                            const double        time);
    unsigned getLatRmNum		      () {  return cntLatRm;        }
    void setNonReConLatList           (   std::vector<unsigned>*          p_latList);

    void setNonReConLat					(	unsigned 			LatticeID);
    unsigned checkAllLatticeSM  (   PCG_Eigen*                  p_EigenPCG,
                                    ConnectLat*                 p_conLat);

private:

    unsigned 			    cntNodeRm;
    unsigned 			    cntLatRm;
    unsigned                cntLatShearRm;
    bool 				    isRecord;
    double 				    accSE;
    double				    chkCoplanarTol;
    std::vector<unsigned>   nonReConList;

    void findFailLatSimple		(	std::vector<iDouble>*		p_exceed,
                                    double* 					maxLatNormForce);

    void findFailLatSimple		(	std::vector<iDouble>*		p_exceed,
                                    std::vector<double>* 		p_latNormForceList);

    void findFailLatTensile     (   std::vector<iDouble>*       p_exceed,
                                    std::vector<double>*        p_latNormForceList);

    void findFailLatMC          (   std::vector<iTriple>*       p_exceed,
                                    std::vector<double>*        p_latNormForceList);

    void findFailLatTensionAndShear  (  std::vector<iTriple>*   p_exceed,
                                        std::vector<double>*    p_latNormForceList);

    void findFailLatResultant   (   std::vector<iDouble>*   	p_exceed,
                                    std::vector<double>*    	p_latNormForceList);

    void checkNbNode			(	const unsigned 				LatticeID,
    								const unsigned				step,
    								const double				time,
    								ConnectLat*					p_conLat,
    								Vertices*					p_verti,
    								PCG_Eigen*					p_EigenPCG);

    double reConnectLatSingle	(	unsigned					LatticeID,
    								PCG_Eigen*					p_EigenPCG,
    								ConnectLat*                 p_conLat);

    void disConnectLatSingle	(	const unsigned				LatticeID,
        							const unsigned				step,
        							const double				time,
        							Vertices*					p_verti,
        							ConnectLat*					p_conLat,
        							PCG_Eigen*					p_EigenPCG);

    bool isLatAtBoundary 		(	const std::vector<iDouble> 		rmList);

    void recordFailLatStressState  (   const unsigned      LatticeID);
};
#endif /* GLATREMOVAL_H_ */
