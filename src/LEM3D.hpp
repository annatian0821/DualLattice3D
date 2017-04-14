#pragma once
#ifndef _LEM3D_HPP
#define _LEM3D_HPP

#include <iostream>
#include <fstream>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

#include <climits>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <list>
#include <vector>
#include <array>
#include <set>
#include <algorithm>
#include <time.h>
#include <stddef.h>
#include <sys/time.h>
#include <string>

//#include "dataclass.hpp"
//Geometry
const unsigned 			Dim		 				= 3;						//2D case - 2, 3D case - 3
extern unsigned         ProblemID;
extern unsigned         GeometryID;
extern int	 			Nx;	 					//= 175;                 		//Number of points in x direction 175
extern int	 			Ny;	 					//= 101;               		//Number of points in y direction 101
extern int	 			Nz;						//= 101;						//Number of points in z direction
extern double			ScaleBoundary;
extern double			ScaleOuterLength;
const unsigned 			NeighbourNum		 	= 6+6*(Dim-2);				//
const unsigned 			BoundaryNum		 		= 4;
const unsigned			DofPerNode				= 6;						//DOF per node, 3 for spring model
const unsigned          KNum                    = 4;

const double 			Pi								= 3.14159265358979;

typedef std::array<double,Dim>								Point;
typedef std::array<double,Dim>								UniVec;
typedef std::array<double,Dim>								Vec;
typedef std::array<UniVec,Dim>								BasisVecs;
typedef	std::array<Point,2>									Line;
typedef std::vector<Point>									Surface;
typedef std::array<float,Dim>								PointF;
typedef std::vector<PointF>									SurfaceF;
typedef	std::array<double,DofPerNode>						Tensor;
typedef	std::array<double,6>								Tensor6;

typedef std::array<std::array<bool,DofPerNode>,6> 			RestrainMat;
typedef	std::pair<unsigned,double>							IDAndInfo;	//Store fNodeID and Injection rate
typedef std::pair<unsigned,std::array<double,3> >           IDAndInfo3;
typedef std::pair<unsigned,unsigned>						IDpair;

typedef boost::iostreams::tee_device<std::ostream,std::ofstream> TeeDispFile;
typedef boost::iostreams::tee_device<boost::iostreams::stream<TeeDispFile>,std::ofstream> TeeDispFileFile;

//Lattice Properties
extern double 					UnitLength;
extern double					CrackRadius;
extern double 					NormLatStiffness; 			       		//ks, Material stiffness
extern double                   ShearToAxialStiffnessRatio;              // alpha = ks/kn, default = 1.0
extern double					MomentToAxialStiffnessRatio;			// default = 1.0
extern std::array<double,2> 	MicroAxialStrength; 						//Force need to break a lattice
extern double                   MicroShearStrength;
extern double                   StrengthFactor;                         //Modification factor for strength
extern double                   Phi;
extern double 					Theta;							// Max. principal stress to break a lattice

extern double					PackingDensity;							// Percent from most dense HCF packing 0.7405
extern double					MinArea;
extern double					MinApecture;							// Minimum aperture to define closure of crack

extern double					LoadIncrFactor;							//Determine the min. increment of loading
extern double					LoadDecrFactor;							//Reduction of loading during relaxation, value < 1.0
extern double                   FluidFactor;                            //Factor to modify fluid pressure
extern double					InitForce;								//Specify the initial force
extern double					ReConRatio;						        //Critical Surface Energy
extern double					FractureToughness;
extern double					MaxPressure;							//Max fluid pressure allowed
extern double					Viscosity;
extern double                   LeakOffCoeff;
extern double					InjectRate;								//Input constant discharge rate
extern double					TimePerStep;								//TimeStep
extern double					InitTime;
extern double                   avgD;
//Control parameters
extern bool                     Use_UnevenMesh;
extern bool 			        Use_Eigen;
extern bool 			        Use_CompactStiffnessMatrix;
extern bool				        Use_LEFM;
extern bool				        Use_UnbreakableOuterMesh;
extern bool                     Use_CubicFlow;
extern bool                     Use_ReconnectLattice;
extern bool                     Use_FullPlasticModel;
extern bool                     Use_FrictionModel;
extern bool                     Use_RandomNode;


extern bool 			Use_RanLatStiffness;
extern bool             DispControl;

extern bool 			Use_RanLatStrength;
extern bool 			isTestCase;

extern bool				Use_OpenMP;

extern bool				Adj_ResForce;
extern bool             Use_const_NodeNum;

extern double			PCG_Precision;							//Parameter determining the precision for Preconditioned Conjugate Gradient Solver
extern double           Sigma;
extern double 			Disp_Scale;								//The scale factor for displacement displayed in output

extern unsigned         LatFailureCriterion;
extern unsigned         StrengthDistModel;
extern unsigned         LoadDir;
extern unsigned         LoadFace;
extern unsigned 		MaxPerStep_BrokenLattice;				//Control no of broken lattice each step
extern unsigned 		MaxPerStep_SplitNode;

extern unsigned 		BoundaryLoadCase;						//Loading exerted on boundary, 0 : no boundary disp applied
extern unsigned 		MaxLoadStep;
extern unsigned			MaxRelaxStep;
extern unsigned			MaxTotalStep;

extern unsigned			SampleRate;								// The frequency of printing full node and lattice information (1 means to output every time)

//Parallel processing
extern unsigned 		ThreadNumUsed;

const unsigned 		Type_0 						= 0;								//real node
const unsigned 		Type_Virt 					= 1;								//virtual node
const unsigned 		Type_Outer					= 2;    							//spilt node
const unsigned 		Type_Boundary 				= 3;								//edge node
const unsigned 		Type_Unstable 				= 4;								//removed node
const unsigned 		Type_Dummy 					= 5;								//dummy


const unsigned		Face1						= 0;
const unsigned		Face2						= 1;
const unsigned		Face3						= 2;
const unsigned		Face4						= 3;
const unsigned		Face5						= 4;
const unsigned		Face6						= 5;
//Boundary Fixity

extern RestrainMat					BoundaryMat;
extern Tensor6						InitStress;
//Shift to reserve space for virtual nodes

extern unsigned		xShift;
extern unsigned		yShift;
extern unsigned 	zShift;

//Constants
const double 		Tiny 							= 1.e-10;
const double		tiny							= 1.e-20;
const double 		Huge 							= 1.e+10;

//Dummy values
extern unsigned 			DumNodeID;
extern unsigned 			DumLatticeID;
extern unsigned				OffsetDumLatticeID;
const unsigned 				DumDir 					= NeighbourNum;

//Temp variable
extern unsigned				tInNodeNum;									//Total number of node in inner fine mesh
extern unsigned 			tNodeNum;									//Total number of real node
extern unsigned             crackNodeNum;                               //Total number of node representing pre-existing fracture
extern unsigned             outerCrackNodeNum;                           //Total number of node representing weak plane
extern unsigned				tbNodeNum;									//Total number of boundary nodes
extern unsigned 			tvNodeNum;									//Total number of virtual node
extern unsigned 			tLatticeNum;								//Total number of real lattice
extern unsigned				tbLatticeNum;								//Total number of boundary lattice
extern unsigned 			tvLatticeNum;								//Total number of virtual lattice
extern unsigned             PriInjLattice;                               //LatticeID of primany injection point

extern unsigned				initRmNode;									//Number of node removed to form initial crack
extern unsigned				initRmLat;									//Number of lattice removed to form initial crack

extern char*				OutputFolder;								//Folder that store output files
extern char*				OutputSubFolder;							//Sub Folder that store statistics files
extern char*				OutputFilePrefix;							//Prefix of all output files
extern char*				RunTimeCommendFile;							//File that
extern std::string          OutputRawFolder;

//For printing log file and cout
extern std::ofstream 								ofs;
extern TeeDispFile 									tee1;
extern boost::iostreams::stream<TeeDispFile> 		dualOut;

extern std::ofstream 								errFile;
extern TeeDispFile 									teeErr;
extern boost::iostreams::stream<TeeDispFile> 		errOut;
extern TeeDispFileFile 								tee2;
extern boost::iostreams::stream<TeeDispFileFile> 	tripOut;


#endif
