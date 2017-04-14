
#ifndef _OUTPUT_HPP
#define _OUTPUT_HPP

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "ProProcess.hpp"
#include "ConnectLat.hpp"
#include "Vertices.hpp"
#include "MiscTools.hpp"
#include "FluidNetwork.hpp"
#include "FluidCal.hpp"
#include "PCG_Eigen.hpp"
#include "Stress.hpp"

class Output
{
public:
    Output 		() {};

    ~Output 	() {};

    Output 				            (	char* 					inPath,
                                        char* 					inPrefix);

    void writeAll					(	const unsigned 					step,
                                        const double                    extForce,
    	                               	FluidCal*						p_fCal,
    	                               	ConnectLat*						p_conLat,
    	                               	Vertices*						p_verti,
    	                                PCG_Eigen*                      p_EigenPCG);

    void writeLatAndNetwork			(	const unsigned 					step,
                                        const double                    extForce,
                                    	ConnectLat*						p_conLat,
                                    	Vertices*						p_verti,
                                        PCG_Eigen*                      p_EigenPCG);

    void writeNodeAndLattice		(	const unsigned 					step,
                                        PCG_Eigen*                      p_EigenPCG);

    void writeForceChain            (   const unsigned                  step,
                                        const double                    extForce,
                                        PCG_Eigen*                      p_EigenPCG);

    void writeFluidNetwork			(	const unsigned 					step,
                                    	FluidCal*						p_fCal);

    void writeMSplot				(	const unsigned 					step,
                                        ConnectLat*						p_conLat,
                                        Vertices*						p_verti);

    void writePriDir 				(	const unsigned 					step,
    		                            ConnectLat*						p_conLat,
    		                            Vertices*						p_verti);

    void writeNode 		(	const unsigned 					step);

    void writeGNode		(	const unsigned 					step);

    void writeCentroid 	();

    void writeGLattice	(	const unsigned 					step);

    void writeFailSurface 	(	const unsigned 					step,
                                const std::vector<unsigned>*	p_failLatList,
                                Vertices*						p_verti);

    void writeFailNetwork	(	const unsigned 					step,
                                ConnectLat*						p_conLat,
                                Vertices*						p_verti,
                                PCG_Eigen*                      p_EigenPCG);

    void writetFailSurface 	();

    void writeAllFailSurface   ( Vertices*                   p_verti);

    void writettFailSurface ();

    void writeStress		(	const unsigned 				step);

    void writeLatAxes       (   const unsigned              step,
                                ConnectLat*                 p_conLat);
    void writeLatAxesAll    ();

    void writeIntersect     (   const unsigned              step,
                                ConnectLat*                 p_conLat);

    void writeOuterBox			();
    void writeInnerBox			();

    void test			    (   const char* fileName);

    void writettFailSurface 	(Vertices*					p_verti);

    void writeRestrain     ();

    void printSingleCell	();

    void writeRawNode         (     const unsigned          totalStep);

    void writeRawLattice      (     PCG_Eigen*              p_EigenPCG,
                                    const unsigned          totalStep);

private:

    char*				path;
    char*				prefix;

    std::ofstream nodeFile;
    std::ofstream latticeFile;

    unsigned tDispNode;
    unsigned tDispLattice;
    unsigned tDispfNode;
    unsigned tDispPipe;

    void Node_d 			(	std::ostream*					p_file);

    void Node_ID 			(	std::ostream*                   p_file);

    void Node_extF 			(	std::ostream*					p_file);

    void Node_extF_Mag 		(	std::ostream*					p_file);

    void Node_stress		(	std::ostream*					p_file);

    void Node_PriStress 	(	std::ostream*					p_file);

    void GNode_ne 			(	std::ostream*					p_file);

    void GNode_v 			(	const char* 					fileName);

    void Node_residual 		(	std::ostream*					p_file);

    void Node_isFree        (   std::ostream*                   p_file);

    void Node_type          (   std::ostream*                   p_file);

    void Lattice_ke			(	std::ostream*					p_file,
                                const char* 					varName);

    void Lattice_F 			(	std::ostream*					p_file);

    void Lattice_fLocal     (   std::ostream*                   p_file,
                                PCG_Eigen*                      p_EigenPCG);

    void Lattice_fLocal_simp(   std::ostream*                   p_file,
                                PCG_Eigen*                      p_EigenPCG);

    void Lattice_fStress    (   std::ostream*                   p_file,
                                PCG_Eigen*                      p_EigenPCG);

    void Lattice_fPerimeter (   ConnectLat*                     p_conLat,
                                const unsigned                  step);

    void Lattice_offset 	(	std::ostream*			        p_file);

    void Lattice_isBreak 	(	const char* 					fileName,
                                const char* 					varName);

    void GLattice_normF		(	std::ostream*					p_file,
                                const char* 					varName);

    void GLattice_ke		(	std::ostream*					p_file,
                                const char* 					varName);

    void Lattice_ShearToAxialRatio 		(	std::ostream*		p_file);

    void GLattice_e			(	std::ostream*					p_file,
    							const char* 					varName);

    void GLattice_area		(	const char* 					fileName,
                                const char* 					varName);

    void GLattice_Ft		(	std::ostream*					p_file,
    							const char* 					varName);

    void Fluid_width 		(	std::ostream*					p_file,
                                FluidCal*	 					p_fCal);

    void Fluid_pressure 	(	std::ostream*					p_file,
                                FluidCal*	 					p_fluid);

    void Fluid_flowRate 	(	std::ostream*					p_file,
                                FluidCal* 						p_fluid);

    void Fluid_leakOff      (   std::ostream*                   p_file,
                                FluidCal*                       p_fCal);

//	void Lattice_ang(const char* fileName, const std::vector<LatState>* latState);
};

#endif
