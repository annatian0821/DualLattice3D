/*
 * Stress.hpp
 *
 *  Created on: 31 Aug 2015
 *      Author: john
 */

#ifndef STRESS_HPP_
#define STRESS_HPP_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "PCG_Eigen.hpp"
#include "MiscTools.hpp"

class Stress
{
	friend class PCG_Eigen;
public:
    Stress() {};
    ~Stress() {};

    void 	calBasic			();

    void 	calPriStress		();

    double	getStress			(	const unsigned 				NodeID,
                                    const unsigned 				k);

    double 	getPriStress		(	const unsigned 				NodeID,
                                    const unsigned 				k);

    double  getStress_p         (   const unsigned              NodeID);

    double  getStress_q         (   const unsigned              NodeID);

    double	getPriDir			(	const unsigned 				NodeID,
                                    const unsigned 				d,
                                    const unsigned 				k);

    Vec     getPriDir           (   const unsigned              NodeID,
                                    const unsigned              d);

private:
    unsigned								nNum;
    std::vector<std::array<double,6> > 		stress;
    std::vector<std::array<double,3> > 		pStress;
    std::vector<std::array<double,3*3> > 	pDir;

    void calPDir				(	const unsigned				i,
                                    std::vector<iDouble>*		p_tmpStress);
};


#endif /* STRESS_HPP_ */
