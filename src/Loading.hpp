/*
 * Loading.hpp
 *
 *  Created on: Oct 25, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef LOADING_HPP_
#define LOADING_HPP_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "ConnectLat.hpp"
#include "GLatRemoval.hpp"
#include "FluidNetwork.hpp"
#include "MiscTools.hpp"
#include "FluidCal.hpp"

class Loading
{
public:
	Loading () {};
	Loading							(	double					_incrRatio,
										double					_decrRatio)
	{
		incrRatio = _incrRatio;
		decrRatio = _decrRatio;
		sumF = 0.0;
	}

	double	getSumF					()	{ return sumF; };

	double	getDeltaF				()	{ return deltaF; };

	double  getFacePressure         (   const unsigned          face,
	                                    const unsigned          dir);

	void applyFluidPressure			(	FluidCal*				p_fCal,
                                        ConnectLat*             p_conLat);

	void addLoad					(	const double			maxNorm,
            							const unsigned 			face,
            							const unsigned 			dir);

	void addLoad					(	const double			maxNorm,
		            					ConnectLat*				p_conLat);

	void addLoad                    (   const double            maxNorm);

	void reduceLoad					(	const double            maxNorm,
	                                    const unsigned 			face,
	            						const unsigned 			dir);

	void reduceLoad					(	const double			maxNorm,
	            						ConnectLat*				p_conLat);

	void reduceLoad                 (   const double            maxNorm);

	void initLoad					(	const double			deltaF,
	            						const unsigned 			face,
	            						const unsigned 			dir);

	void initialDisp                (   const unsigned           face,
	                                    const unsigned           dir,
	                                    const double             disp);

	void initLoad					(	const double			_deltaF,
	            						ConnectLat*				p_conLat);

	void setOldMaxNorm				(	const unsigned			maxBreak);

    void assignPresDisp				(	const unsigned 			NodeID,
                                        const double 			dx,
                                        const double 			dy,
                                        const double 			dz,
                                        const bool* 			flag);

    void setForceZero 				();
    void setConForceZero            ();

    void applyGUniPressureOnBoundary (  const double            pressure,
                                        const unsigned          face,
                                        const unsigned          dir);

    void applyConstFluidPressure	(	std::vector<unsigned> 	latList,
    									const double			sumP);

    void logSumF					(	const unsigned 			step );

    void outputForceNodes           (   const unsigned          face);

private:
    double			incrRatio,decrRatio;
    double			sumF,deltaF,oldMaxNorm;
    void applySingleConstFluidPressure (	const unsigned		LatticeID,
        									const double		sumP);

    void applySingleConstFluidPressure (	const unsigned		LatticeID,
        									const double		deltaP,
        									const double		sumP);

    void assignNodeExtForce				(	const unsigned 			NodeID,
                                            const double 			fx,
                                            const double 			fy,
                                            const double 			fz);

    void clearAllLoading				();

};

#endif /* LOADING_HPP_ */
