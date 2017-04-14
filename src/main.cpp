/*
 * main.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: john
 */
#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "Input.hpp"
#include "InitCond.hpp"
#include "Statistics.hpp"
#include "Loading.hpp"
#include "GConjGrad.hpp"
#include "ProProcess.hpp"
#include "Testing.hpp"
#include "GLatRemoval.hpp"
#include "Output.hpp"
#include "MiscTools.hpp"
#include "PCG_Eigen.hpp"
#include "Tesellation.hpp"
#include "Vertices.hpp"
#include "ConnectLat.hpp"
#include "FluidNetwork.hpp"
#include "FluidCal.hpp"
#include "voro++.hh"

NodeLists nodeLists;


std::ofstream ofs;
std::ofstream errFile;

TeeDispFile tee1(std::cout,ofs);
boost::iostreams::stream<TeeDispFile> dualOut(tee1);
TeeDispFile teeErr(std::cerr,ofs);
boost::iostreams::stream<TeeDispFile> errOut(teeErr);
TeeDispFileFile tee2(errOut,errFile);
boost::iostreams::stream<TeeDispFileFile> tripOut(tee2);

int hydraulicFracturing		(	Clock* 							p_gClock);

int ucs						(	Clock* 							p_gClock);

std::array<bool,2>  solidFluidCoupling     (	ConnectLat*						p_conLat,
                                                Vertices*						p_verti,
                                                FluidCal*						p_fCal,
                                                PCG_Eigen*						p_EigenPCG,
                                                Loading*						p_loading,
                                                RTCommand*						p_rtc,
                                                Output*							p_output,
                                                GLatRemoval*					p_glatRemoval,
                                                bool& 							isRunTimeStop,
                                                unsigned&						relaxStep);

int main (int argc, char** argv) {
    if(argc != 2) {
        std::cerr << "Incorrect number of arguments" << std::endl;
        std::cerr << "[USAGE:] ./lem /path/to/inputFile/inputFile.txt" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    Clock* 					p_gClock		= new Clock();
    Input*					p_input			= new Input();
    p_input         -> initLogFile ();
    p_input         -> configFile(argv);

    dualOut << "**************************************************************************" << std::endl;
    dualOut << "*****************************START OF PROGRAM*****************************" << std::endl;
    dualOut << "**************************************************************************" << std::endl;
    dualOut << "Program Staring time = " << p_gClock -> getDateAndTime() << std::endl;
    p_gClock 		-> start("Global");
    unsigned status = 0;
    if (p_input -> isRead) {
        delete	p_input;
        if (ProblemID==0) {
            ucs(p_gClock);
        }
        if ((ProblemID==1)||(ProblemID==2)) {
            hydraulicFracturing (p_gClock);
        }
        if (ProblemID==999) {
            status = ucs(p_gClock);
        }
        if (ProblemID==9999) {
            status = hydraulicFracturing (p_gClock);
        }
        if (status == 0) {
            dualOut << "Program terminated successfully..." << std::endl;
        } else if (status == 1) {
            dualOut << "Program terminated because of unstable model..." << std::endl;
        } else if (status ==3) {
            dualOut << "Program terminated because of specified totalStep is reached..." << std::endl;
        } else if (status == 4) {
            dualOut << "Program terminated because of specified relaxStep is reached..." << std::endl;
        }
        dualOut << "Program ending time = " << p_gClock -> getDateAndTime() << std::endl;
        dualOut << "**************************************************************************" << std::endl;
        dualOut << "******************************END OF PROGRAM******************************" << std::endl;
        dualOut << "**************************************************************************" << std::endl;
        p_gClock			-> get();
    }
    else {
        dualOut << "Cannot read input file, program stop..." << std::endl;
    }
    delete		p_gClock;
    return 0;
}


int ucs						(	Clock* 					p_gClock)
{
    unsigned runStatus = 0;
	RTCommand*					p_rtc 			= new RTCommand (RunTimeCommendFile);
	Vertices*					p_verti			= new Vertices();
	ConnectLat*					p_conLat		= new ConnectLat(0.0);
	Output*						p_output		= new Output(OutputFolder,OutputFilePrefix);
	unsigned					maxCoordNo = 0;
	VoroTesell*					p_voro 			= new VoroTesell();
	p_voro                      -> generate(1.0*UnitLength,ScaleBoundary,
	                                ScaleOuterLength,&maxCoordNo,p_verti);
	std::vector<unsigned>*      p_pxCrackLatList= new std::vector<unsigned> (p_voro -> getPreExCrackLatID ());
	std::vector<unsigned>       weakPlaneLatList= p_voro -> getWeakPlaneLatID ();
	std::array<unsigned,2> 		pxCrackPerimeterNum = p_voro->getPxCrackPerimeterNum();
	Statistics*					p_stat 			= new Statistics(OutputSubFolder,OutputFilePrefix,maxCoordNo,
	                                                    2*p_pxCrackLatList->size(),p_voro->getRefineNode());
	if (GeometryID==2) {
	    p_stat      -> setCrackMonthLatList (0.05*Nx*UnitLength,p_pxCrackLatList);
	}
	p_stat		-> writeInLatAll(2*p_pxCrackLatList->size(),0,100,p_conLat);

	delete 		p_voro;
	p_output	-> writeOuterBox();
	p_output	-> writeInnerBox();

	GLatRemoval*		p_glatRemoval	= new GLatRemoval(1e-4);
	InitCond*			p_initCond 		= new InitCond(p_glatRemoval);
	PCG_Eigen*			p_EigenPCG		= new PCG_Eigen(maxCoordNo,150000);
	p_initCond 	-> setInitCond	(p_conLat,p_verti,p_EigenPCG,p_pxCrackLatList,pxCrackPerimeterNum,maxCoordNo);
	p_initCond  -> setInitCrackAperture(p_pxCrackLatList);
	delete		p_pxCrackLatList;
	if (ProblemID==999) {
	        p_output                    -> writeLatAndNetwork (0,InitForce,p_conLat,p_verti,p_EigenPCG);
	        p_stat                      -> writeWeakPlaneStress (p_conLat,p_EigenPCG,&weakPlaneLatList,InitStress[2],0);
	        p_stat                      -> writeStressStrain  (InitStress[2],0);
	        p_stat                      -> fillLatPriDirList (60,0);
	        p_stat                      -> writeStress(&weakPlaneLatList,-InitStress[2],0,0);
	        p_stat                      -> writeLatStress (InitStress[2],0,p_conLat,p_EigenPCG);
	        p_stat                      -> writelatPolar (180,0,InitStress[2],p_conLat,p_EigenPCG);
	        p_stat                      -> writelatPolar (180,9,0,InitStress[2],p_conLat,p_EigenPCG);
	        return 0;
	}

	GLatForce*			p_glatForce 	= new GLatForce();
	FailStats*			p_fStats 		= new FailStats(OutputSubFolder,OutputFilePrefix,maxCoordNo);
	omp_set_num_threads(ThreadNumUsed);
	PreCondConjGrad*	p_customPCG 	= new PreCondConjGrad();
	p_glatRemoval		-> setRecord(true);
	Loading*			p_loading		= new Loading(LoadIncrFactor,LoadDecrFactor);

	bool 				isRunTimeStop 	= false;

	unsigned				relaxStep = 0,totalStep = 0;
	bool					relaxCheck 	= true;
	bool					isStable	= true;
	double					maxNorm 	= Huge;
	p_output            -> writeRestrain();
	const unsigned 			loadDir = LoadDir;
	if (DispControl) {
	    p_loading                   -> initialDisp    (LoadFace,loadDir,InitForce);
	} else {
	    p_loading					-> initLoad(InitForce,LoadFace,loadDir);
	}
	isStable = p_EigenPCG 		-> solve(totalStep,true,Use_CompactStiffnessMatrix);
	if (!isStable) {
		tripOut << "***Model unstable, quit calculation!***" << std::endl;
		runStatus = 1;
		p_output		-> writeLatAndNetwork (totalStep,InitForce,p_conLat,p_verti,p_EigenPCG);
		p_fStats 		-> write (totalStep+1,30);
		return 1;
	}
	p_conLat					-> setLogOffset(true);
	EnergyCal*       			p_eCal = new EnergyCal();
	p_eCal                      -> fillForceNodeList (LoadFace,loadDir);
	unsigned rmLatNum = 0;
	unsigned loadCnt = 0;
	for (unsigned loadStep = 0; loadStep < MaxLoadStep; loadStep++) {
		relaxStep = 0;
		relaxCheck = true;
		tripOut << "-------------------LoadStep " << loadStep << " -------------------" << std::endl;
		double sumF;
		do {
			tripOut << "--------RelaxStep " << relaxStep << "--------TotalStep " << totalStep <<" --------" << std::endl;
			p_rtc->update ();
			isRunTimeStop = p_rtc->isStop();
			if (isRunTimeStop) {
				break;
			}
			if (Use_ReconnectLattice) {
			    p_initCond		-> stableLat (p_conLat,p_verti,p_EigenPCG,isStable,relaxStep,totalStep,10,20);
			} else {
			    isStable = p_EigenPCG       -> solve(totalStep,true,Use_CompactStiffnessMatrix);
			    totalStep++;
			}
			if (!isStable) {
			    break;
			}
			if (DispControl) {
			    sumF = p_loading-> getFacePressure(LoadFace,loadDir);
			    p_loading                   -> setConForceZero();
			} else {
			    sumF = p_loading->getSumF();
			}
			p_stat                      -> writeStressStrain    (sumF,loadCnt);
			p_stat						-> writeStressStrain	(2,2,sumF,loadCnt);
			p_stat						-> writeStressStrain	(1,1,sumF,loadCnt);
			p_stat						-> writeStressStrain	(0,0,sumF,loadCnt);
			p_gClock->get();
			if ((CrackRadius>Tiny)&&(ProblemID==0)) {
			    p_conLat                    -> logFractureVol (sumF,loadCnt);
			}
			if (GeometryID==2) {
			    p_stat                      -> writeCrackMouthDispHist(sumF);
			}
			p_eCal                      -> record();
			rmLatNum = p_glatRemoval 	-> remove	(p_conLat,p_verti,p_fStats,p_EigenPCG,p_eCal,
			                                            loadCnt,0.0,&maxNorm,MaxPerStep_BrokenLattice);
			if (rmLatNum) {
			    isStable = p_EigenPCG       -> solve(loadCnt,true,Use_CompactStiffnessMatrix);
			}
			p_stat                      -> writeLatClusterAgg (p_conLat,loadCnt);
            p_eCal -> accumulate();

			std::cout << "loadCnt = " << loadCnt << std::endl;
			if (loadCnt%SampleRate==0) {
			    p_output    -> writeLatAndNetwork (loadCnt,sumF,p_conLat,p_verti,p_EigenPCG);
			    p_fStats    -> write (loadCnt,30);
			    p_stat      -> writeWeakPlaneStress (p_conLat,p_EigenPCG,&weakPlaneLatList,sumF,loadCnt);
			    p_stat      -> writeLatCluster (p_conLat,100,loadCnt);
			    p_stat      -> writelatPolar (180,loadCnt,sumF,p_conLat,p_EigenPCG);
			    p_stat      -> writelatPolar (180,5,loadCnt,sumF,p_conLat,p_EigenPCG);
			    p_stat      -> writeFailLatPolar (180,p_conLat->getAllGroupVec(),loadCnt,sumF,"all",p_EigenPCG);
			    p_stat      -> writeFailLatPolar (180,p_conLat->getAllLat(0),loadCnt,sumF,"tension",p_EigenPCG);
			    p_stat      -> writeFailLatPolar (60,p_conLat->getAllLat(1),loadCnt,sumF,"shear",p_EigenPCG);
			    p_stat      -> writeFailLatPolar (60,p_conLat->getReConLat(),loadCnt,sumF,"reCon",p_EigenPCG);
			    p_stat      -> writeFailLatPolar (60,p_conLat->getDisConLat(),loadCnt,sumF,"disCon",p_EigenPCG);
			    p_stat      -> fillLatPriDirList (60,loadCnt);
			    p_stat      -> writeStress(&weakPlaneLatList,-sumF,loadCnt,0);
			    p_stat      -> writeLatStress (sumF,loadCnt,p_conLat,p_EigenPCG);
			}
			if (rmLatNum == 0) {
			    if (DispControl) {
			        p_loading   -> addLoad  (maxNorm);
			    } else {
			        p_loading   -> addLoad  (maxNorm,LoadFace,loadDir);
			    }
			    p_eCal -> compute (loadCnt,sumF);
			    loadCnt++;
			} else {
			    if (DispControl) {
			        p_loading   -> reduceLoad   (maxNorm);
			    } else {
				    p_loading	-> reduceLoad	(maxNorm,LoadFace,loadDir);
			    }
				p_conLat 	-> logConOffset (loadCnt,p_verti);
				p_conLat 	-> logOffset (loadCnt);
				loadCnt++;
			}
			relaxCheck = ((rmLatNum>0) && (relaxStep < MaxRelaxStep) && (totalStep < MaxTotalStep));
		} while (relaxCheck);
		if (!isStable) {
			tripOut << "***LEM Model unstable, quit calculation!***" << std::endl;
			p_output	-> writeLatAndNetwork (loadCnt,sumF,p_conLat,p_verti,p_EigenPCG);
			p_fStats 	-> write (loadCnt,30);
			runStatus = 1;
			break;
		}
		if (isRunTimeStop) {
		    p_output    -> writeLatAndNetwork (loadCnt,sumF,p_conLat,p_verti,p_EigenPCG);
		    p_fStats    -> write (loadCnt,30);
		    runStatus = 2;
		    break;
		}
		if (totalStep >= MaxTotalStep) {
			dualOut << "No of step reaches the specified maximum, calculation stop!" << std::endl;
			runStatus = 3;
			break;
		}
		if (relaxStep >= MaxRelaxStep) {
		    tripOut << "No of relax step reach the specified maximum = " << MaxRelaxStep << " calculation stop!" << std::endl;
		    runStatus = 4;
		}
	}

	delete					p_initCond;
	delete					p_rtc;
	delete					p_output;
	delete					p_fStats;
	delete					p_glatForce;
	delete					p_stat;
	delete					p_loading;
	delete					p_glatRemoval;
	delete					p_EigenPCG;
	delete					p_customPCG;
	delete					p_verti;
	delete					p_conLat;
	delete                  p_eCal;
	delete[] 				OutputFolder;
	delete[]				OutputFilePrefix;

	return runStatus;
}

int hydraulicFracturing		(	Clock* 					p_gClock)
{
    unsigned runStatus = 0;
	RTCommand*					p_rtc 			= new RTCommand (RunTimeCommendFile);
	Vertices*					p_verti			= new Vertices();
	ConnectLat*					p_conLat		= new ConnectLat(MinApecture);
	Output*						p_output		= new Output(OutputFolder,OutputFilePrefix);
	unsigned	maxCoordNo = 0;
	VoroTesell*	p_voro = new VoroTesell();
	p_voro -> generate(1.0*UnitLength,ScaleBoundary,ScaleOuterLength,&maxCoordNo,p_verti);
	std::vector<unsigned>*      p_pxCrackLatList= new std::vector<unsigned> (p_voro -> getPreExCrackLatID ());
	std::vector<unsigned>*      p_boreHoleNodeList = new std::vector<unsigned> (p_voro -> getBoreHoleNodeList ());
	std::vector<unsigned>       weakPlaneLatList= p_voro -> getWeakPlaneLatID ();
	std::array<unsigned,2> pxCrackPerimeterNum = p_voro->getPxCrackPerimeterNum();
	Statistics*	p_stat = new Statistics(OutputSubFolder,OutputFilePrefix,maxCoordNo,
	                                    2*pxCrackPerimeterNum.size(),p_voro->getRefineNode());
	p_stat		-> writeInLatAll(2*p_pxCrackLatList->size(),0,100,p_conLat);
	delete 		p_voro;
	p_output	-> writeOuterBox();
	p_output	-> writeInnerBox();
	p_stat 		-> compute(30);

	const double 		minP = 0; //std::max(1.15*InitStress[2],0.0);//1.1*std::min(std::min(InitStress[0],InitStress[1]),InitStress[2]);
    FluidCal*			p_fCal			= new FluidCal(Viscosity,InitTime,TimePerStep,1e-3,
                                                1e-4,MinApecture,LoadIncrFactor,LoadDecrFactor,FluidFactor,minP,MaxPressure,25000);
    GLatRemoval*		p_glatRemoval	= new GLatRemoval(MinApecture);
	InitCond*			p_initCond 		= new InitCond(p_glatRemoval);
	PCG_Eigen*			p_EigenPCG		= new PCG_Eigen(maxCoordNo,20000);

	p_initCond      -> setBoreHoleNodeList (p_boreHoleNodeList);
	p_initCond 	    -> setInitCond	(p_conLat,p_glatRemoval,p_verti,p_fCal,p_EigenPCG,
									p_pxCrackLatList,pxCrackPerimeterNum,minP,maxCoordNo);
	delete      p_boreHoleNodeList;

	GLatForce*			p_glatForce = new GLatForce();
	FailStats*			p_fStats 	= new FailStats(OutputSubFolder,OutputFilePrefix,maxCoordNo);
	Testing*			p_test		= new Testing(Face6,OutputSubFolder,OutputFilePrefix);

	omp_set_num_threads(ThreadNumUsed);

	PreCondConjGrad*	p_customPCG 	= new PreCondConjGrad();
	p_glatRemoval		-> setRecord(true);
	Loading*			p_loading		= new Loading(LoadIncrFactor,LoadDecrFactor);

	bool isRunTimeStop = false;
	unsigned				relaxStep = 0,totalStep = 0;
	bool					relaxCheck 			= true;
	bool					isStable			= true;
	double					maxNorm 			= Huge;

	if (ProblemID==2) {
	    p_fCal                  ->  setInitialList(p_conLat);
	}
	p_initCond                  -> setInitCrackAperture(p_pxCrackLatList);
	delete      p_pxCrackLatList;
	p_conLat					-> setLogOffset(true);
	p_loading					-> initLoad(InitForce,p_conLat);
//	p_fCal                      -> applyConstFluidPressure(InitForce);
	isStable = p_EigenPCG 		-> solve(totalStep,true,Use_CompactStiffnessMatrix);
	p_stat						-> writeStress(&weakPlaneLatList,p_fCal->getPriInjP(),0,0.0);

	EnergyCal*       p_eCal = new EnergyCal();
	unsigned rmLatAcc = 0;
	double time = 0.0;
	unsigned rmLatNum = 0;

	p_fCal      -> updateFluidNetwork (p_verti,p_conLat,0.0);

	if (ProblemID==2) {
	    std::vector<unsigned>   fNetwork = p_conLat->getGroup(0);
	    unsigned initFNodeNum = fNetwork.size();
        do {
            dualOut << "Search for breakdown pressure..." << std::endl;
            rmLatNum = p_glatRemoval    -> remove   (p_conLat,p_verti,p_fStats,p_EigenPCG,p_eCal,
                                                     totalStep,time,&maxNorm,MaxPerStep_BrokenLattice);
            if (rmLatNum!=0) {
                isStable = p_EigenPCG       -> solve(totalStep,true,Use_CompactStiffnessMatrix);
                if (Use_FullPlasticModel) {
                    p_glatRemoval   -> updateReConLatShearStiffness  (p_conLat,p_EigenPCG,p_verti,ReConRatio,totalStep,time);
                } else if (Use_FrictionModel) {
                    p_glatRemoval   -> updateReConLatShearStiffness  (p_conLat,p_EigenPCG,p_verti,totalStep,time);
                }
                totalStep++;
            }
            if (false) {
                p_output    -> writeNodeAndLattice(totalStep,p_EigenPCG);
                p_output    -> writeFailNetwork(totalStep,p_conLat,p_verti,p_EigenPCG);
            }
        } while (rmLatNum!=0);
        dualOut << "Finish searching for breakdown pressure" << std::endl;
        p_fCal      ->  setBHInitStorage();
        if (GeometryID==4) {
            p_fCal      ->  setNotchList(&weakPlaneLatList);
        }
	}
	p_fCal      -> setInitialArea(p_conLat);
	p_loading   -> initLoad(InitForce,p_conLat);
	p_EigenPCG  -> solve(0,true,Use_CompactStiffnessMatrix);
    if (true) {
        p_output    -> writeAll (0,1.0,p_fCal,p_conLat,p_verti,p_EigenPCG);
        p_stat      -> writeLatCluster (p_conLat,100,0);
        double minR = avgD/(std::min(std::min(UnitLength*Nx,UnitLength*Ny),UnitLength*Nz));
        if (((GeometryID==0)||(GeometryID==3))&&(CrackRadius>minR)) {
        }

        p_stat      -> writeStress(&weakPlaneLatList,p_conLat,1.0,totalStep,time);
        p_stat      -> writeWeakPlaneStress (p_conLat,p_EigenPCG,&weakPlaneLatList,InitForce,0);
        if (StrengthFactor<2.0) {
            p_stat      -> writelatPolar(180,0,1.0,p_conLat,p_EigenPCG);
            p_stat      -> writelatPolar(180,5,0,1.0,p_conLat,p_EigenPCG);
            p_stat      -> writeFailLatPolar (180,p_conLat->getAllGroupVec(),totalStep,1.0,"all",p_EigenPCG);
            p_stat      -> writeFailLatPolar (180,p_conLat->getAllLat(0),totalStep,1.0,"tension",p_EigenPCG);
            p_stat      -> writeFailLatPolar (60,p_conLat->getAllLat(1),totalStep,1.0,"shear",p_EigenPCG);
            p_stat      -> writeFailLatPolar (180,p_conLat->getReConLat(),totalStep,1.0,"reCon",p_EigenPCG);
            p_stat      -> writeFailLatPolar (180,p_conLat->getDisConLat(),totalStep,1.0,"disCon",p_EigenPCG);
            p_stat      -> fillLatPriDirList (60,totalStep);
        }
    }
	for (unsigned timeStep = 0; timeStep <MaxLoadStep; timeStep++) {
		relaxStep = 0;
		relaxCheck = true;
		time = p_fCal->getTime();
		tripOut << "-------------------TimeStep " << timeStep
				<< ", time = " << time << " -------------------" << std::endl;
		std::array<bool,2> SF_coupling_status;
		if ((SF_coupling_status[1])&&(timeStep!=0)) {
		    //Fracture propagation stage
			do {
			    p_eCal->record(p_fCal);
				rmLatNum = p_glatRemoval 	-> remove	(p_conLat,p_verti,p_fStats,p_EigenPCG,p_eCal,
														totalStep,time,&maxNorm,MaxPerStep_BrokenLattice);
				if (rmLatNum!=0) {
				    isStable = p_EigenPCG   -> solve(totalStep,true,Use_CompactStiffnessMatrix);
				    p_eCal->accumulate();
				}
				if (Use_ReconnectLattice) {
				    p_initCond				-> stableLat (p_conLat,p_verti,p_EigenPCG,isStable,relaxStep,totalStep,10,20);
				}
				isRunTimeStop = p_rtc		-> isStop();
				p_glatForce                 -> calculate();
				p_rtc->update ();
				rmLatAcc += rmLatNum;
				relaxStep++;
			} while ((rmLatNum>0) && isStable && !isRunTimeStop);
			p_stat                      -> writeLatClusterAgg (p_conLat,timeStep);
			if ((isStable)&&(!isRunTimeStop)&&(rmLatAcc>0)) {
			    p_eCal-> compute (time);
			}
			dualOut << "Total lattice removed in this time step = " << rmLatAcc << std::endl;
			relaxStep = 0;
		}
		//Fluid flow stage
		SF_coupling_status[0] = true;
		SF_coupling_status[1] = false;
		while (relaxCheck) {
			tripOut << "-------------------TotalStep " << totalStep << "-------------------" << std::endl;
			maxNorm = 0.;
			SF_coupling_status =   solidFluidCoupling (p_conLat,p_verti,p_fCal,p_EigenPCG,p_loading,p_rtc,p_output,
			                                            p_glatRemoval,isRunTimeStop,relaxStep);
			dualOut << "Step = " << timeStep << " relaxStep = " << relaxStep << " totalStep = " << totalStep << " maxNorm = " << maxNorm
					<< " deltaF = " << p_loading->getDeltaF() << " sumF = " << p_loading->getSumF() << std::endl;
			p_gClock->get();
			isRunTimeStop = p_rtc->isStop();
			p_rtc->update ();
			relaxCheck = (SF_coupling_status[0] && (relaxStep < MaxRelaxStep) && (totalStep < MaxTotalStep)
							&&(!isRunTimeStop)&&(SF_coupling_status[1]));
		}
		if (!SF_coupling_status[0]) {
			p_output	-> writeAll (totalStep,1.0,p_fCal,p_conLat,p_verti,p_EigenPCG);
			p_fStats 	-> write (timeStep+1,30);
			runStatus = 1;
			break;
		}
		if (totalStep >= MaxTotalStep) {
			dualOut << "No of step reaches the specified maximum, calculation stop!" << std::endl;
//			p_fCal  -> updateCalInfo(qStatus);
			runStatus = 3;
			break;
		}
		if (relaxStep >= MaxRelaxStep) {
		    dualOut << "The max no of relaxStep has been reach, proceed to next step" << std::endl;
		}
		p_conLat 	-> logFLat (time,rmLatAcc);
		p_conLat 	-> logConOffset (time,p_verti);
		p_conLat 	-> logOffset (time);
		rmLatAcc =0;
		if (isRunTimeStop) {
			p_glatForce	-> calculate ();
			p_output	-> writeAll (totalStep,1.0,p_fCal,p_conLat,p_verti,p_EigenPCG);
			p_fStats 	-> write (totalStep,30);
			runStatus = 2;
			break;
		}
		p_stat          -> Fluid_pressure (p_fCal,timeStep);
		if ((timeStep+1)%SampleRate==0) {
			p_output	-> writeAll (timeStep+1,1.0,p_fCal,p_conLat,p_verti,p_EigenPCG);
			double minR = avgD/(std::min(std::min(UnitLength*Nx,UnitLength*Ny),UnitLength*Nz));
			if (((GeometryID==0)||(GeometryID==3))&&(CrackRadius>minR)) {
			}
			p_stat      -> writeStress(&weakPlaneLatList,p_conLat,1.0,timeStep,time);
			p_stat      -> writeWeakPlaneStress (p_conLat,p_EigenPCG,&weakPlaneLatList,p_fCal->getPriInjP(),timeStep);
			if (StrengthFactor<2.0) {
			    p_fCal      -> printConnectivity(timeStep);
			    p_stat      -> writeLatCluster   (p_conLat,100,timeStep);
			    p_stat      -> writelatPolar     (180,timeStep,1.0,p_conLat,p_EigenPCG);
			    p_stat      -> writelatPolar     (180,3,timeStep,1.0,p_conLat,p_EigenPCG);
			    p_stat      -> writeFailLatPolar (180,p_conLat->getAllGroupVec(),timeStep,1.0,"all",p_EigenPCG);
			    p_stat      -> writeFailLatPolar (180,p_conLat->getAllLat(0),timeStep,1.0,"tension",p_EigenPCG);
			    p_stat      -> writeFailLatPolar (60,p_conLat->getAllLat(1),timeStep,1.0,"shear",p_EigenPCG);
			    p_stat      -> writeFailLatPolar (60,p_conLat->getReConLat(),timeStep,1.0,"reCon",p_EigenPCG);
			    p_stat      -> writeFailLatPolar (60,p_conLat->getDisConLat(),timeStep,1.0,"disCon",p_EigenPCG);
			    p_stat      -> fillLatPriDirList (60,timeStep);
			}
		}
		if ((timeStep==10)&&(SampleRate==1)) {
		    SampleRate *= 10;
		}
	}
	delete					p_initCond;
	delete					p_rtc;
	delete					p_output;
	delete					p_fStats;
	delete					p_fCal;
	delete					p_glatForce;
	delete					p_stat;
	delete					p_loading;
	delete					p_glatRemoval;
	delete					p_EigenPCG;
	delete					p_customPCG;
	delete					p_test;
	delete					p_verti;
	delete					p_conLat;
	delete                  p_eCal;
	delete[] 				OutputFolder;
	delete[]				OutputFilePrefix;
	return runStatus;
}

std::array<bool,2> solidFluidCoupling			 (	ConnectLat*						p_conLat,
    											    Vertices*						p_verti,
    											    FluidCal*						p_fCal,
    											    PCG_Eigen*						p_EigenPCG,
    											    Loading*						p_loading,
    											    RTCommand*						p_rtc,
    											    Output*							p_output,
    											    GLatRemoval*					p_glatRemoval,
    											    bool& 							isRunTimeStop,
    											    unsigned&						relaxStep)
{
	static bool isFirst = true;
	unsigned i = 0;
//	static std::pair<bool,bool> 		qStatus(false,false);
	std::array<bool,2>  out;
	static std::array<bool,2>  qStatus = {true,false};
	if (p_rtc->isStop()) {
	    out[0] = false;
	    out[1] = false;
	    return out;
	}
	bool isFluidStable = true;
	if (isFirst) {
	    p_fCal                      -> updateFluidNetwork(p_verti,p_conLat,InitTime);
	    p_fCal                      -> setP0(InitForce);
	    p_fCal                      -> setTime (InitTime);
//	    p_loading                   -> applyFluidPressure(p_fCal,p_conLat);
	    isFirst = false;
	}

	tripOut << "-------------------Start Solid-Fluid Dual Lattice Modelling-------------------" << std::endl;

	if (Use_FullPlasticModel) {
	    p_glatRemoval           -> updateReConLatShearStiffness  (p_conLat,p_EigenPCG,p_verti,ReConRatio,
	                                p_fCal->getTimeStep(),p_fCal->getTime());
	} else if (Use_FrictionModel) {
	    p_glatRemoval   -> updateReConLatShearStiffness  (p_conLat,p_EigenPCG,p_verti,p_fCal->getTimeStep(),p_fCal->getTime());
	}

	if ((qStatus[1]==true)||(!Use_CubicFlow)) {
	    std::vector<unsigned>   latList = p_fCal->getLatID();
	    p_loading                   -> applyConstFluidPressure(latList,p_fCal->getPriInjP());
	} else {
	    p_loading                   -> applyFluidPressure(p_fCal,p_conLat);
	}

/*	if ((p_fCal->isPChange())&&(Use_CubicFlow)) {
	    p_loading                   -> applyFluidPressure(p_fCal,p_conLat);
	}*/
	double isSolidStable   = p_EigenPCG         -> solve(p_fCal->getTimeStep(),true,Use_CompactStiffnessMatrix);
	qStatus                = p_fCal             -> updateQ();
	dualOut << "Q status : " << qStatus[0] << ' ' <<  qStatus[1] << std::endl;
	if (qStatus[0]) {
	    qStatus[0] = qStatus[1] = false;
	    dualOut << " Mass balance achieve, proceed to next time step" << qStatus[0] << " convergence of Q = " << qStatus[1];
	    p_fCal->initizeNextStep ();
	    out[0]= true, out[1] = true;
	    return out;
	} else {
	    out[1]= false;
	}

/*	if (qStatus[0]) {
	    p_fCal->initizeNextStep ();
	    out[0]= true, out[1] = true;
	    return true;
	}*/

	bool is_p_converge = false;
	if (Use_CubicFlow) {
	    isFluidStable     = p_fCal     -> solve();
	    is_p_converge      = p_fCal     -> interpolate_pi();
	    if (is_p_converge) {
	        p_fCal      -> interpolatingP0 ();
	    }
	} else {
	    p_fCal          ->  interpolatingP0 ();
	    p_fCal          ->  updatingPressureProfile((!Use_CubicFlow)||(is_p_converge));
	    qStatus[1] = true;
	    out[0]= true;
	}

//	p_fCal						    -> updateFluidNetwork(p_verti,p_conLat,time);
/*	if (Use_CubicFlow) {
	    fluidStatus.second = p_fCal -> calculate(qStatus,timeStep);
	} else {
	    fluidStatus.second = p_fCal	-> calculate_static(qStatus,timeStep);
	}

	if (!fluidStatus.second) {
		tripOut << "[stableLatFluid] - Fluid model unstable, quit calculation" << std::endl;
		return false;
	}*/
//	p_loading 					-> applyFluidPressure(p_fCal,p_conLat);

//    fluidStatus.first  = p_fCal -> updateCalInfo(qStatus);

    p_rtc->update ();
//    isStable = p_EigenPCG 		-> solve(timeStep,true,Use_CompactStiffnessMatrix);
//    p_fCal -> updateRes();
    if (!isSolidStable) {
        tripOut << "Solid Lattice model unstable, quit calculation" << std::endl;
        out[0]= false;
        return out;
    }
    if (!isFluidStable) {
        tripOut << "Fluid Lattice model unstable, quit calculation" << std::endl;
        out[0]= false;
        return out;
    }
		relaxStep++;
		i++;

	tripOut << "-------------------Finished Solid-Fluid Dual Lattice Modelling-------------------" << std::endl;
	return out;
}
