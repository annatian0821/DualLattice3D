/*
 * Input.cpp
 *
 *  Created on: Oct 9, 2013
 *      Author: john
 */

#include "Input.hpp"

using namespace std;

//Global Variable to be declared

unsigned tInNodeNum,tNodeNum,tbNodeNum,tvNodeNum,tLatticeNum,tbLatticeNum,tvLatticeNum;
unsigned crackNodeNum,outerCrackNodeNum;
//Global Variable to be inputted
unsigned        ProblemID;
unsigned        GeometryID;
int 			Nx,Ny,Nz;
double			ScaleBoundary;
double			ScaleOuterLength;
double			CrackRadius;

//Lattice Properties
double 					UnitLength;
double 					NormLatStiffness;
double                  ShearToAxialStiffnessRatio;              // alpha = ks/kn, default = 1.0
double					MomentToAxialStiffnessRatio;
std::array<double,2>	MicroAxialStrength;
double                  MicroShearStrength;
double                  Phi;
double                  StrengthFactor;                         //Modification factor for strength
double					Theta;
double					TimePerStep;
double					InitTime;

double				Disorder;

double				PackingDensity;
double				MinArea;
double				MinApecture;

double				LoadIncrFactor;							//Determine the min. increment of loading
double				LoadDecrFactor;
double              FluidFactor;
double				InitForce;								//Specify the initial force
double				ReConRatio;
double				Viscosity;
double              LeakOffCoeff;
double				MaxPressure;
double				InjectRate;								//Input constant discharge rate
double				PCG_Precision;							//Parameter determining the precision for Preconditioned Conjugate Gradient Solver
double              Sigma;
double              Disp_Scale;

//unsigned			dofPerNode;
unsigned            LatFailureCriterion;
unsigned            StrengthDistModel;
unsigned            LoadDir;
unsigned            LoadFace;
unsigned 			MaxPerStep_BrokenLattice;				//Control no of broken lattice each step
unsigned			MaxRelaxStep;
unsigned			MaxTotalStep;

unsigned 			BoundaryLoadCase;				//Loading exerted on boundary, 0 : no boundary disp applied
unsigned 			MaxLoadStep;

unsigned			SampleRate;						// The frequency of printing full node and lattice information (1 means to output every time)
unsigned 			ThreadNumUsed;

bool                Use_UnevenMesh;
bool 				Use_Eigen;
bool 				Use_CompactStiffnessMatrix;
bool				Use_LEFM;
bool				Use_UnbreakableOuterMesh;
bool                Use_CubicFlow;
bool                Use_ReconnectLattice;
bool                Use_RandomNode;

bool 				Use_RanLatStiffness;
bool                DispControl;

bool				Use_FrictionModel;
bool                Use_FullPlasticModel;
bool 				Use_RanLatStrength;
bool 				isTestCase;
bool                isStable;

bool				Use_OpenMP;

bool				Adj_ResForce;
bool                Use_const_NodeNum;



RestrainMat			BoundaryMat;
Tensor6				InitStress;

char*						OutputFolder;								//Folder that store output files
char*						OutputSubFolder;
char*						OutputFilePrefix;							//Prefix of all output files
char*						RunTimeCommendFile;
std::string                 OutputRawFolder;

bool Input::initLogFile () {
    time_t t = time(NULL);
    tm* timePtr = localtime(&t);
    if (isTestCase) {
        LogPathName = "/home/john/output/log/log-test.log";
        ErrLogPathName = "/home/john/output/log/errLog-test.log";
    } else {
        char logPathName[1000];
        char errLogPathName[1000];
        std::sprintf(logPathName,"/home/john/output/log/log-%02d%02d%02d-%02d%02d%02d.log",
                timePtr->tm_year-100,timePtr->tm_mon+1,timePtr->tm_mday,timePtr->tm_hour,timePtr->tm_min,timePtr->tm_sec);
        std::sprintf(errLogPathName,"/home/john/output/log/errLog-%02d%02d%02d-%02d%02d%02d.log",
                timePtr->tm_year-100,timePtr->tm_mon+1,timePtr->tm_mday,timePtr->tm_hour,timePtr->tm_min,timePtr->tm_sec);
        LogPathName = logPathName;
        ErrLogPathName = errLogPathName;
    }
    ofs.open (LogPathName, std::ios::out | std::ios::app);
    errFile.open (ErrLogPathName, std::ios::out | std::ios::app);
}

int Input::mkpath(std::string s,mode_t mode) {
    size_t pre=0,pos;
    std::string dir;
    int mdret;

    if(s[s.size()-1]!='/') {
        // force trailing / so we can handle everything in loop
        s+='/';
    }

    while((pos=s.find_first_of('/',pre))!=std::string::npos) {
        dir=s.substr(0,pos++);
        pre=pos;
        if(dir.size()==0) continue; // if leading / first time is 0 length
        if((mdret=mkdir(dir.c_str(),mode)) && errno!=EEXIST) {
            return mdret;
        }
    }
    return mdret;
}

RTCommand::RTCommand	(char*	_fileName) {
	fileName = _fileName;
}

bool RTCommand::isStop ()
{
	infile.open(fileName, ios::in);
	if (infile) {
		infile.ignore (10000,'#');
		infile.ignore (255, '=');
		unsigned	inBool;
		infile >> inBool;
		if (inBool) {
			dualOut << "A commend from run time interface file " << fileName << " requesting termination of calculation."
					<< " Program will terminate after completion of all calculations in current step" << endl;
		}
		infile.close();
		return inBool;
	} else {
		cout << "Cannot read file, return false and continue running the program" << endl;
		infile.close();
		return false;
	}
	return false;
}

bool RTCommand::update ()
{
	static unsigned cnt = 0;
	infile.open(fileName, ios::in);
	if (infile) {
		infile.ignore (10000,'#');
		infile.ignore (255, '=');
		unsigned	inBool;
		infile >> inBool;
		unsigned			dataUnsigned;
		infile.ignore (255, '=');
		infile >> dataUnsigned;
		if (dataUnsigned<=cnt) {
			return false;
		} else {
			Clock clock;
			dualOut << "-----------Runtime commend No. "<< dataUnsigned << " received on " << clock.getDateAndTime() << endl;
			vector<unsigned>	inUnsigned;
			vector<unsigned>	isChange;
			unsigned			changeFlag;
			for (unsigned i=0; i<3; i++)
			{
				infile.ignore (255, '=');
				infile >> changeFlag;
				isChange.push_back(changeFlag);
				infile >> dataUnsigned;
				inUnsigned.push_back(dataUnsigned);
			}
			double 			dataDouble;
			vector<double>	inDouble;
			for (unsigned i=0; i<3; i++)
			{
				infile.ignore (255, '=');
				infile >> changeFlag;
				isChange.push_back(changeFlag);
				infile >> dataDouble;
				inDouble.push_back(dataDouble);
			}
			if (!isChange.empty()) {
				if (isChange.front()) {
					SampleRate = inUnsigned[0];
					dualOut << "SampleRate has changed to " << SampleRate << endl;
				}
				isChange.erase(isChange.begin());
			}
			if (!isChange.empty()) {
				if (isChange.front()) {
					MaxTotalStep = inUnsigned[1];
					dualOut << "MaxTotalStep has changed to " << MaxTotalStep << endl;
				}
				isChange.erase(isChange.begin());
			}
			if (!isChange.empty()) {
				if (isChange.front()) {
					MaxPerStep_BrokenLattice = inUnsigned[2];
					dualOut << "MaxPerStep_BrokenLattice has changed to " << MaxPerStep_BrokenLattice << endl;
				}
				isChange.erase(isChange.begin());
			}
			if (!isChange.empty()) {
				if (isChange.front()) {
					LoadIncrFactor = inDouble[0];
					dualOut << "LoadIncrFactor has changed to " << LoadIncrFactor << endl;
				}
				isChange.erase(isChange.begin());
			}
			if (!isChange.empty()) {
				if (isChange.front()) {
					LoadDecrFactor = inDouble[1];
					dualOut << "LoadDecrFactor has changed to " << LoadDecrFactor << endl;
				}
				isChange.erase(isChange.begin());
			}
			if (!isChange.empty()) {
			    if (isChange.front()) {
			        TimePerStep = inDouble[2];
			        dualOut << "TimePerStep has changed to " << TimePerStep << endl;
			    }
			    isChange.erase(isChange.begin());
			}
			if (!isChange.empty()) {
				dualOut << "Extra input is detected and is discarded!" << endl;
			}
			dualOut <<"----------------End of runtime commend No. "<< ++cnt << " --------------------" << endl;
			infile.close();
			return true;
		}
	} else {
		cout << "Cannot read file, return false and continue running the program" << endl;
		infile.close();
		return false;
	}
		return false;
}

void Input::configFile (char** argv)
{
    unsigned 		dataUnsigned;
    int 			dataInt;
    double			dataDouble;
    bool			dataBool;
    unsigned        ignoreLen = 1000;
    string			tString;

    RestrainMat		tempMat;

    std::string fileName = argv[1];
    std::cout << "File read as: " << fileName << std::endl;

    infile.open(fileName, ios::in);
    if (infile) {
        infile.ignore(10000,'#');
        {
            for (unsigned i=0; i<5; i++) {
                infile.ignore (ignoreLen, '=');
                infile >> dataInt;
                inInt.push_back(dataInt);
            }

            for (unsigned i=0; i<30; i++) {
                infile.ignore (ignoreLen, '=');
                infile >> dataDouble;
                inDouble.push_back(dataDouble);
            }

            for (unsigned i=0; i<11; i++) {
                infile.ignore (ignoreLen, '=');
                infile >> dataUnsigned;
                inUnsigned.push_back(dataUnsigned);
            }

            for (unsigned i=0; i<17; i++) {
                infile.ignore (ignoreLen, '=');
                infile >> dataBool;
                inBool.push_back(dataBool);
            }

            for (unsigned i=0; i<6; i++) {
                for (unsigned j=0; j<6; j++) {
                    infile.ignore (ignoreLen, '=');
                    if (j<DofPerNode) {
                        infile >> tempMat[i][j];
                    }
                }
            }

            for (unsigned i=0; i<6; i++) {
                infile.ignore (ignoreLen, '=');
                infile >> InitStress[i];
            }

            infile.ignore (ignoreLen, '=');
            getline(infile,tString,'#');
            OutputFolder = new char[tString.size()+1];
            copy(tString.begin(), tString.end(), OutputFolder);
            OutputFolder[tString.size()] = '\0';
            cout << OutputFolder;


            int status = this->mkpath(OutputFolder,0755);
            std::cout << "Creating directory tree: " << OutputFolder << std::endl;

            infile.ignore (ignoreLen, '=');
            getline(infile,tString,'#');
            OutputSubFolder = new char[tString.size()+1];
            copy(tString.begin(), tString.end(), OutputSubFolder);
            OutputSubFolder[tString.size()] = '\0';

            status = this->mkpath(OutputSubFolder,0755);
            std::cout << "Creating directory tree: " << OutputSubFolder << std::endl;

            infile.ignore (ignoreLen, '=');
            getline(infile,tString,'#');
            OutputFilePrefix = new char[tString.size()+1];
            copy(tString.begin(), tString.end(), OutputFilePrefix);
            OutputFilePrefix[tString.size()] = '\0';

            cout << OutputFilePrefix;

            infile.ignore (ignoreLen, '=');
            getline(infile,tString,'#');
            RunTimeCommendFile = new char[tString.size()+1];
            copy(tString.begin(), tString.end(), RunTimeCommendFile);
            RunTimeCommendFile[tString.size()] = '\0';

            OutputRawFolder = OutputFolder;
            if (OutputRawFolder.back()=='/') {
                OutputRawFolder += "Raw";
            } else {
                OutputRawFolder += "/Raw";
            }
            this->mkpath(OutputRawFolder,0755);
            std::cout << "Creating directory tree: " << OutputRawFolder << std::endl;

            cout << "Read Config.txt completed" <<endl;
        }
        ProblemID                   = inInt[0];
        GeometryID                  = inInt[1];
        Nx 							= inInt[2];
        Ny 							= inInt[3];
        Nz 							= inInt[4];

        UnitLength					= inDouble[0];
        ScaleBoundary				= inDouble[1];
        ScaleOuterLength			= inDouble[2];
        CrackRadius					= inDouble[3];
        NormLatStiffness 			= inDouble[4];              		//Material stiffness
        ShearToAxialStiffnessRatio  = inDouble[5];
        MomentToAxialStiffnessRatio = inDouble[6];
        MicroAxialStrength[0] 		= inDouble[7];						//Force need to break a lattice (Tension)
        MicroAxialStrength[1] 		= inDouble[8];						//Force need to break a lattice (Compression)
        MicroShearStrength          = inDouble[9];
        Phi                         = inDouble[10];
        StrengthFactor              = inDouble[11];
        Theta					    = inDouble[12];						//Max. principal stress to break a lattice
        MaxPressure					= inDouble[13];
        InjectRate					= inDouble[14];
        Viscosity					= inDouble[15];
        LeakOffCoeff                = inDouble[16];
        TimePerStep					= inDouble[17];
        InitTime					= inDouble[18];
        PackingDensity				= inDouble[19];
        MinArea						= inDouble[20];
        MinApecture					= inDouble[21];
        LoadIncrFactor				= inDouble[22];						//Determine the min. increment of loading
        LoadDecrFactor				= inDouble[23];
        FluidFactor                 = inDouble[24];
        InitForce					= inDouble[25];						//Specify the initial force
        ReConRatio			        = inDouble[26];
        PCG_Precision				= inDouble[27];
        Sigma                       = inDouble[28];
        Disp_Scale                  = inDouble[29];

        LatFailureCriterion         = inUnsigned[0];
        StrengthDistModel           = inUnsigned[1];
        LoadDir                     = inUnsigned[2];
        LoadFace                    = inUnsigned[3];
        MaxPerStep_BrokenLattice	= inUnsigned[4];					//Control no of broken lattice each step
        BoundaryLoadCase 			= inUnsigned[5];					//Loading exerted on boundary, 0 : no boundary disp applied
        MaxLoadStep 				= inUnsigned[6];
        MaxRelaxStep				= inUnsigned[7];
        MaxTotalStep				= inUnsigned[8];
        SampleRate					= inUnsigned[9];
        ThreadNumUsed				= inUnsigned[10];

        Use_UnevenMesh              = inBool[0];
        Use_Eigen 					= inBool[1];						//  1 =  Use Matrix Library, 0 use custom library
        Use_CompactStiffnessMatrix	= inBool[2];						//  1 =  Reduce Stiffness Matrix size by dropping the restrain DOF
        Use_LEFM					= inBool[3];
        Use_UnbreakableOuterMesh	= inBool[4];
        Use_CubicFlow               = inBool[5];
        Use_ReconnectLattice        = inBool[6];
        Use_RandomNode              = inBool[7];

        Use_RanLatStiffness 	    = inBool[8];
        DispControl                 = inBool[9];

        Use_FrictionModel			= inBool[10];
        Use_FullPlasticModel        = inBool[11];

        Use_RanLatStrength			= inBool[12];
        isTestCase			        = inBool[13];

        Use_OpenMP					= inBool[14];

        Adj_ResForce				= inBool[15];
        Use_const_NodeNum           = inBool[16];


        BoundaryMat	 				= tempMat;

        logConfig ();
        isRead = true;

    }
    else
    {
        dualOut << "Cannot open file: Config.txt";
        isRead = false;
    }
}

void Input::logConfig ()
{
    Clock gClock;

    dualOut << " ----------Configurations----------" << std::endl;
    dualOut << " ProblemID = " << ProblemID << std::endl;
    dualOut << " GeometryID = " << GeometryID << std::endl;
    dualOut << " Dimensions of model :" << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " (Nx,Ny,Nz) = " << Nx << " " << Ny << " " << Nz << std::endl;
    dualOut << " Use_RandomNode = " << Use_RandomNode << std::endl;
    dualOut << " Use_UnevenMesh = " << Use_UnevenMesh << std::endl;
    dualOut << " ScaleBoundary = " << ScaleBoundary << std::endl;
    dualOut << " ScaleOuterLength = " << ScaleOuterLength << std::endl;
    dualOut << " CrackRadius = " << CrackRadius	 <<std::endl;
    dualOut << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Lattice Parameters :" << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " DofPerNode = " << DofPerNode << std::endl;
    dualOut << " KNum = " << KNum << std::endl;
    dualOut << " UnitLength = " << UnitLength << std::endl;
    dualOut << " CrackRadius = " << CrackRadius << std::endl;
    dualOut << " Theta = " << Theta << std::endl;
    dualOut << " NormLatStiffness = " << NormLatStiffness << std::endl;
    dualOut << " ShearToAxialStiffnessRatio = " << ShearToAxialStiffnessRatio << std::endl;
    dualOut << " MomentToAxialStiffnessRatio = " << MomentToAxialStiffnessRatio << std::endl;
    dualOut << " Use_LEFM = " << Use_LEFM << std::endl;
    dualOut << " MicroAxialStrength (max,min) = " << MicroAxialStrength[0] << " " << MicroAxialStrength[1] << std::endl;
    dualOut << " MicroShearStrength = " << MicroShearStrength << std::endl;
    dualOut << " Phi = " << Phi << std::endl;
    dualOut << " StrengthFactor = " << StrengthFactor << std::endl;
    dualOut << " StrengthDistModel = " << StrengthDistModel << std::endl;
    dualOut << " Use_RanLatStiffness = " << Use_RanLatStiffness << std::endl;
    dualOut << " Use_RanLatStrength = " << Use_RanLatStrength << std::endl;
    dualOut << " Sigma = " << Sigma << std::endl;
    dualOut << " LatFailureCriterion = " << LatFailureCriterion << std::endl;
    dualOut << " ReConRatio = " << ReConRatio << std::endl;
    dualOut << " Disorder = " << Disorder << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Fluid Parameters :" << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Use_CubicFlow = "  << Use_CubicFlow << std::endl;
    dualOut << " Max fluid pressure = " << MaxPressure << std::endl;
    dualOut << " Fluid viscosity = " << Viscosity << std::endl;
    dualOut << " LeakOffCoeff = " << LeakOffCoeff << std::endl;
    dualOut << " Injection rate = " << InjectRate << std::endl;
    dualOut << " FluidFactor = "    << FluidFactor << std::endl;
    dualOut << " MinApecture = "    << MinApecture << std::endl;
    dualOut << " TimeStep = "		<< TimePerStep << std::endl;
    dualOut << " InitTime = "		<< InitTime << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Voronoi Tesellation :" << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " PackingDensity = " << PackingDensity << std::endl;
    dualOut << " MinArea = " << MinArea << std::endl;
    dualOut << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Control parameters :" << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " InitForce = "      << InitForce << std::endl;
    dualOut << " MaxLoadstep = " 	<< MaxLoadStep << std::endl;
    dualOut << " MaxRelaxStep = " 	<< MaxRelaxStep << std::endl;
    dualOut << " MaxTotalStep = " 	<< MaxTotalStep << std::endl;
    dualOut << " MaxPerStep_BrokenLattice = " << MaxPerStep_BrokenLattice << std::endl;
    dualOut << " LoadIncrFactor = " << LoadIncrFactor << std::endl;
    dualOut << " LoadDecrFactor = " << LoadDecrFactor << std::endl;
    dualOut << " Sample Rate = "	<< SampleRate << std::endl;
    dualOut << " InitForce = " 		<< InitForce << std::endl;
    dualOut << " LoadDir = "        << LoadDir << std::endl;
    dualOut << " LoadFace = "       << LoadFace << std::endl;
    dualOut << " DispControl = "    << DispControl << std::endl;
    dualOut << " Use_FrictionModel = " << Use_FrictionModel << std::endl;
    dualOut << " Use_FullPlasticModel = " << Use_FullPlasticModel << std::endl;
    dualOut << " Adj_ResForce = " 	<< Adj_ResForce << std::endl;
    dualOut << " MinApecture = "	<< MinApecture << std::endl;
    dualOut << " Use_ReconnectLattice = " << Use_ReconnectLattice << std::endl;
    dualOut << " Use_const_NodeNum = " << Use_const_NodeNum << std::endl;
    dualOut << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Solver parameters :" << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " PCG_Precision = " << PCG_Precision << std::endl;
    dualOut << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Output parameters :" << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Disp_Scale = " << Disp_Scale << std::endl;
    dualOut << " isTestCase = " << isTestCase << std::endl;
    dualOut << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Boundary Conditions : " << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Boundary fixity (x,y,z) :" << std::endl;
    dualOut << " Face1 = " << BoundaryMat[0][0] << " " << BoundaryMat[0][1] << " " << BoundaryMat[0][2];
    if (DofPerNode==6) {
        dualOut << ' ' << BoundaryMat[0][3] << " " << BoundaryMat[0][4] << " " << BoundaryMat[0][5];
    }
    dualOut << std::endl;
    dualOut << " Face2 = " << BoundaryMat[1][0] << " " << BoundaryMat[1][1] << " " << BoundaryMat[1][2];
    if (DofPerNode==6) {
        dualOut << ' ' << BoundaryMat[1][3] << " " << BoundaryMat[2][4] << " " << BoundaryMat[3][5];
    }
    dualOut << std::endl;
    dualOut << " Face3 = " << BoundaryMat[2][0] << " " << BoundaryMat[2][1] << " " << BoundaryMat[2][2];
    if (DofPerNode==6) {
        dualOut << ' ' << BoundaryMat[2][3] << " " << BoundaryMat[2][4] << " " << BoundaryMat[2][5];
    }
    dualOut << std::endl;
    dualOut << " Face4 = " << BoundaryMat[3][0] << " " << BoundaryMat[3][1] << " " << BoundaryMat[3][2];
    if (DofPerNode==6) {
        dualOut << ' ' << BoundaryMat[3][3] << " " << BoundaryMat[3][4] << " " << BoundaryMat[3][5];
    }
    dualOut << std::endl;
    dualOut << " Face5 = " << BoundaryMat[4][0] << " " << BoundaryMat[4][1] << " " << BoundaryMat[4][2];
    if (DofPerNode==6) {
        dualOut << ' ' << BoundaryMat[4][3] << " " << BoundaryMat[4][4] << " " << BoundaryMat[4][5];
    }
    dualOut << std::endl;
    dualOut << " Face6 = " << BoundaryMat[5][0] << " " << BoundaryMat[5][1] << " " << BoundaryMat[5][2];
    if (DofPerNode==6) {
        dualOut << ' ' << BoundaryMat[5][3] << " " << BoundaryMat[5][4] << " " << BoundaryMat[5][5];
    }
    dualOut << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " Initial Stress :" << std::endl;
    dualOut << " (Sigma_xx,Sigma_yy,Sigma_zz)  = " << InitStress[0] << " " << InitStress[1] << " " << InitStress[2] << std::endl;
    dualOut << " (Sigma_xy,Sigma_xz,Sigma_yz)  = " << InitStress[3] << " " << InitStress[4] << " " << InitStress[5] << std::endl;
    dualOut << std::endl;


    dualOut << " ------------------------" << std::endl;
    dualOut << " Misc. parameters : " << std::endl;
    dualOut << " ------------------------" << std::endl;
    dualOut << " ThreadNumUsed = " 		<< ThreadNumUsed << std::endl;
    dualOut << " Output path = "   		<< OutputFolder  << std::endl;
    dualOut << " Output SubFolder = " 	<< OutputSubFolder << std::endl;
    dualOut << " Prefix for output = "  << OutputFilePrefix  << std::endl;
    dualOut << " RunTimeCommendFile = " << RunTimeCommendFile  << std::endl;
    dualOut << " LogPath = " << LogPathName << std::endl;
    dualOut << " ErrLogPath = " << ErrLogPathName << std::endl;
    dualOut << " Time = " << gClock.getDateAndTime() << std::endl;
    dualOut << std::endl;

    dualOut << " ----------End of Configurations----------" << std::endl;

    char fileName[255]; //filename
    sprintf(fileName,"%s/ConfigLog.vtk",OutputFolder);
    ofstream disk(fileName, ios::out);

    sprintf(fileName,"%s/ConfigLog.vtk",OutputSubFolder);
    ofstream dropbox(fileName, ios::out);

    TeeDispFile tee1(disk,dropbox);
    boost::iostreams::stream<TeeDispFile> file(tee1);


    file << " ----------Configurations----------" << std::endl;
    file << " ProblemID = " << ProblemID << std::endl;
    file << " GeometryID = " << GeometryID << std::endl;
    file << " Dimensions of model :" << std::endl;
    file << " ------------------------" << std::endl;
    file << " (Nx,Ny,Nz) = " << Nx << " " << Ny << " " << Nz << std::endl;
    file << " Use_RandomNode = " << Use_RandomNode << std::endl;
    file << " Use_UnevenMesh = " << Use_UnevenMesh << std::endl;
    file << " ScaleBoundary = " << ScaleBoundary << std::endl;
    file << " ScaleOuterLength = " << ScaleOuterLength << std::endl;
    file << " CrackRadius = " << CrackRadius  <<std::endl;
    file << std::endl;
    file << " ------------------------" << std::endl;
    file << " Lattice Parameters :" << std::endl;
    file << " ------------------------" << std::endl;
    file << " DofPerNode = " << DofPerNode << std::endl;
    file << " KNum = " << KNum << std::endl;
    file << " UnitLength = " << UnitLength << std::endl;
    file << " CrackRadius = " << CrackRadius << std::endl;
    file << " Theta = " << Theta << std::endl;
    file << " NormLatStiffness = " << NormLatStiffness << std::endl;
    file << " ShearToAxialStiffnessRatio = " << ShearToAxialStiffnessRatio << std::endl;
    file << " MomentToAxialStiffnessRatio = " << MomentToAxialStiffnessRatio << std::endl;
    file << " Use_LEFM = " << Use_LEFM << std::endl;
    file << " MicroAxialStrength (max,min) = " << MicroAxialStrength[0] << " " << MicroAxialStrength[1] << std::endl;
    file << " MicroShearStrength = " << MicroShearStrength << std::endl;
    file << " Phi = " << Phi << std::endl;
    file << " StrengthFactor = " << StrengthFactor << std::endl;
    file << " StrengthDistModel = " << StrengthDistModel << std::endl;
    file << " Use_RanLatStiffness = " << Use_RanLatStiffness << std::endl;
    file << " Use_RanLatStrength = " << Use_RanLatStrength << std::endl;
    file << " Sigma = " << Sigma << std::endl;
    file << " LatFailureCriterion = " << LatFailureCriterion << std::endl;
    file << " ReConRatio = " << ReConRatio << std::endl;
    file << " Disorder = " << Disorder << std::endl;
    file << " ------------------------" << std::endl;
    file << " Fluid Parameters :" << std::endl;
    file << " ------------------------" << std::endl;
    file << " Use_CubicFlow = "  << Use_CubicFlow << std::endl;
    file << " Max fluid pressure = " << MaxPressure << std::endl;
    file << " Fluid viscosity = " << Viscosity << std::endl;
    file << " LeakOffCoeff = " << LeakOffCoeff << std::endl;
    file << " Injection rate = " << InjectRate << std::endl;
    file << " MinApecture = "    << MinApecture << std::endl;
    file << " TimeStep = "       << TimePerStep << std::endl;
    file << " InitTime = "       << InitTime << std::endl;
    file << " ------------------------" << std::endl;
    file << " Voronoi Tesellation :" << std::endl;
    file << " ------------------------" << std::endl;
    file << " PackingDensity = " << PackingDensity << std::endl;
    file << " MinArea = " << MinArea << std::endl;
    file << std::endl;
    file << " ------------------------" << std::endl;
    file << " Control parameters :" << std::endl;
    file << " ------------------------" << std::endl;
    file << " MaxLoadstep = "    << MaxLoadStep << std::endl;
    file << " MaxRelaxStep = "   << MaxRelaxStep << std::endl;
    file << " MaxTotalStep = "   << MaxTotalStep << std::endl;
    file << " MaxPerStep_BrokenLattice = " << MaxPerStep_BrokenLattice << std::endl;
    file << " LoadIncrFactor = " << LoadIncrFactor << std::endl;
    file << " LoadDecrFactor = " << LoadDecrFactor << std::endl;
    file << " Sample Rate = "    << SampleRate << std::endl;
    file << " InitForce = "      << InitForce << std::endl;
    file << " LoadDir = "        << LoadDir << std::endl;
    file << " LoadFace = "       << LoadFace << std::endl;
    file << " DispControl = "    << DispControl << std::endl;
    file << " Use_FrictionModel = " << Use_FrictionModel << std::endl;
    file << " Use_FullPlasticModel = " << Use_FullPlasticModel << std::endl;
    file << " Adj_ResForce = "   << Adj_ResForce << std::endl;
    file << " MinApecture = "    << MinApecture << std::endl;
    file << " Use_ReconnectLattice = " << Use_ReconnectLattice << std::endl;
    file << std::endl;
    file << " ------------------------" << std::endl;
    file << " Solver parameters :" << std::endl;
    file << " ------------------------" << std::endl;
    file << " PCG_Precision = " << PCG_Precision << std::endl;
    file << std::endl;
    file << " ------------------------" << std::endl;
    file << " Output parameters :" << std::endl;
    file << " ------------------------" << std::endl;
    file << " Disp_Scale = " << Disp_Scale << std::endl;
    file << " isTestCase = " << isTestCase << std::endl;
    file << std::endl;
    file << " ------------------------" << std::endl;
    file << " Boundary Conditions : " << std::endl;
    file << " ------------------------" << std::endl;
    file << " Boundary fixity (x,y,z) :" << std::endl;
    file << " Face1 = " << BoundaryMat[0][0] << " " << BoundaryMat[0][1] << " " << BoundaryMat[0][2];
    if (DofPerNode==6) {
        file << ' ' << BoundaryMat[0][3] << " " << BoundaryMat[0][4] << " " << BoundaryMat[0][5];
    }
    file << std::endl;
    file << " Face2 = " << BoundaryMat[1][0] << " " << BoundaryMat[1][1] << " " << BoundaryMat[1][2];
    if (DofPerNode==6) {
        file << ' ' << BoundaryMat[1][3] << " " << BoundaryMat[2][4] << " " << BoundaryMat[3][5];
    }
    file << std::endl;
    file << " Face3 = " << BoundaryMat[2][0] << " " << BoundaryMat[2][1] << " " << BoundaryMat[2][2];
    if (DofPerNode==6) {
        file << ' ' << BoundaryMat[2][3] << " " << BoundaryMat[2][4] << " " << BoundaryMat[2][5];
    }
    file << std::endl;
    file << " Face4 = " << BoundaryMat[3][0] << " " << BoundaryMat[3][1] << " " << BoundaryMat[3][2];
    if (DofPerNode==6) {
        file << ' ' << BoundaryMat[3][3] << " " << BoundaryMat[3][4] << " " << BoundaryMat[3][5];
    }
    file << std::endl;
    file << " Face5 = " << BoundaryMat[4][0] << " " << BoundaryMat[4][1] << " " << BoundaryMat[4][2];
    if (DofPerNode==6) {
        file << ' ' << BoundaryMat[4][3] << " " << BoundaryMat[4][4] << " " << BoundaryMat[4][5];
    }
    file << std::endl;
    file << " Face6 = " << BoundaryMat[5][0] << " " << BoundaryMat[5][1] << " " << BoundaryMat[5][2];
    if (DofPerNode==6) {
        file << ' ' << BoundaryMat[5][3] << " " << BoundaryMat[5][4] << " " << BoundaryMat[5][5];
    }
    file << std::endl;
    file << " ------------------------" << std::endl;
    file << " Initial Stress :" << std::endl;
    file << " (Sigma_xx,Sigma_yy,Sigma_zz)  = " << InitStress[0] << " " << InitStress[1] << " " << InitStress[2] << std::endl;
    file << " (Sigma_xy,Sigma_xz,Sigma_yz)  = " << InitStress[3] << " " << InitStress[4] << " " << InitStress[5] << std::endl;
    file << std::endl;

    file << " ------------------------" << std::endl;
    file << " Misc. parameters : " << std::endl;
    file << " ------------------------" << std::endl;
    file << " ThreadNumUsed = "      << ThreadNumUsed << std::endl;
    file << " Output path = "        << OutputFolder  << std::endl;
    file << " Output SubFolder = "   << OutputSubFolder << std::endl;
    file << " Prefix for output = "  << OutputFilePrefix  << std::endl;
    file << " RunTimeCommendFile = " << RunTimeCommendFile  << std::endl;
    file << " LogPath = " << LogPathName << std::endl;
    file << " ErrLogPath = " << ErrLogPathName << std::endl;
    file << " Time = " << gClock.getDateAndTime() << std::endl;
    file << std::endl;

    file << " ----------End of Configurations----------" << std::endl;


}
