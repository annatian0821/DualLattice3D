/*
 * Testing.hpp
 *
 *  Created on: Dec 31, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef TESTING_HPP_
#define TESTING_HPP_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "MiscTools.hpp"
#include "ProProcess.hpp"
#include "Output.hpp"

class Testing
{
public:
    Testing() {};
    Testing				(	const unsigned 			face,
                            char* 					inPath,
                            char* 					inPrefix)
    {
        unsigned b,d;
        unsigned f[6] = {Face1, Face2, Face3, Face4, Face5, Face6};
        old_bDisp.resize(nodeLists.boundary[face].size(),0.0);

        for (b=0; b<6; b++)
        {
            for (d=0; d<Dim; d++)
            {
                old_AllbDisp[b*Dim+d].resize(nodeLists.boundary[f[b]].size(),0.0);
                old_AllbForce[b*Dim+d].resize(nodeLists.boundary[f[b]].size(),0.0);
            }
        }
        accWD = 0.0;
        accAllWD = 0.0;

        path = inPath;
        prefix = inPrefix;

        if (inPrefix[0]!='\0')
        {
            strcat(prefix,"-");
        }
    }

    ~Testing() {};

    void 	check();
    void	check					(	const double				accSE,
                                        const unsigned				totalStep,
                                        const bool					isRelax);

    bool 	testMinCoordNum			(	const unsigned				MinCoordNum);

    double	calAllResidual			();

    double 	insideReaction			(   const unsigned 					step);

private:

    std::vector<double>			old_bDisp;
    std::vector<double>			old_AllbDisp[6*Dim];
    std::vector<double>			old_AllbForce[6*Dim];
    double						accWD, accAllWD;
    char*						path;
    char*						prefix;

    double	checkEnergy				(	const double				accSE,
                                        const unsigned				totalStep,
                                        const bool					isRelax);

    double	calSingleResidual		(	const unsigned				NodeID,
                                        std::vector<unsigned>		calDir);

    double	getTotalStrainEnergy	();

    void 	totalExtWorkDone		();

    double 	getTotalExtWorkDoneSimple	(	const unsigned 				face,
                                            const unsigned 				dir,
                                            const bool					isRelax);

    double	getTotalExtWorkDoneComplete (	const unsigned				step);

    void 	writeBoundaryReaction	(	const unsigned 				step);

    void 	getBoundaryReaction		(	std::vector<std::vector<double> >*	p_bForce,
                                        const unsigned						step);



    void 	checkForce				(	std::vector<std::vector<double> >*	p_bForce,
                                        const unsigned 						face,
                                        const unsigned						step);
};

#endif /* TESTING_HPP_ */
