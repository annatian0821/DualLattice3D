/*
 * GConjGrad.cpp
 *
 *  Created on: Nov 30, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "GConjGrad.hpp"

using namespace std;

PreCondConjGrad::PreCondConjGrad() {
    // TODO Auto-generated constructor stub

}

PreCondConjGrad::~PreCondConjGrad() {
    // TODO Auto-generated destructor stub
}


bool PreCondConjGrad::solve		(	const unsigned 					step,
                                    const bool 						isRelax,
                                    const bool						useCompSM)
{
    unsigned i;
    unsigned i_min = (unsigned) (1.5*max(max(Nx,Ny),Nz));
    unsigned i_max = 15*i_min;
    double precision = PCG_Precision; ///sqrt(tNodeNum);

    double							init_delta;
    double							new_delta;
    double 							convergence;
    Clock							lClock;
    Clock							sClock;
    GLatForce						lforce;
    StiffnessMatrix 				SM;

    lClock.start("CG");
    ginitialize(&SM,isRelax);

    lClock.get();
    i = 0;
    init_delta = new_delta = dot (&d,&r);

    cout << "Size of SM = " << ((double) SM.getSizeOf())/1024/1024 << "Mb" << endl;
    cout << "Size of data = " << ((double) getSizeOfData())/1024/1024 << "Mb" << endl;

    dualOut << "gIteration starts... precision*init_delta =" << precision*init_delta <<endl;

    sClock.start("Solver");
    while (((new_delta > fabs(precision*init_delta))||(i<i_min))&&(i<i_max))
    {
        integralSolve (&new_delta,&SM);
        i++;
    }
    sClock.get();
    updateNodeStateSimple (step,isRelax);
    convergence = new_delta / precision/init_delta;
    dualOut << "Iteration cycle = " << i << " new_delta / precision*init_delta = " << convergence << endl;

    PCGclean();

    if (convergence > 50.0)
    {
        dualOut << "Sereve divergence, calculation stop!"<< endl;
        lClock.get();
        return false;
    }

    lClock.get();
    return true;
}

void PreCondConjGrad::integralSolve (	double*						p_newDelta,
                                        StiffnessMatrix*			p_SM)
{
    double alpha,beta;
    unsigned 	i,j,k,i_max,dir;
    unsigned	sDOF,nb_DOF,kStart;
    double		dot,oldDelta;
    Clock							lClock;
    i_max = (calNodeNum)*DofPerNode;
    p_SM->MatVecProduct(&d,&q);
    dot = 0.0;
    #pragma omp parallel for firstprivate (i_max) reduction (+:dot)
    for (unsigned i=0; i<i_max; i++) {
        dot += d[i]*q[i];
    }
    alpha = *p_newDelta/dot;
    dot = 0;
    #pragma omp parallel for firstprivate (i_max,alpha) reduction (+:dot)
    for (unsigned i= 0; i< i_max; i++) {
        x[i] = x[i] + alpha*d[i];
        r[i] = r[i] - alpha*q[i];
        s[i] = m[i]*r[i];
        dot += r[i]*s[i];
    }
    oldDelta = *p_newDelta;
    *p_newDelta = dot;
    beta = *p_newDelta / oldDelta;
    #pragma omp parallel for firstprivate (i_max,beta)
    for (unsigned i= 0; i< i_max; i++) {
        d[i] = s[i] + beta*d[i];
    }
}

void	PreCondConjGrad::QuickMatMulti	(	StiffnessMatrix*	p_SM)
{
    unsigned 	iDOF,i,j,k;
    unsigned	sDOF,nb_DOF;
    unsigned	CoordNum;		//coordination no
    double 		sum;

    #pragma omp parallel for default (none) private (iDOF,i,j,k,sum,sDOF,nb_DOF,CoordNum) shared (p_SM)
    for (iDOF = 0; iDOF < freeNodeNum; iDOF++)
    {
        CoordNum = p_SM->getCoordNum(iDOF);
        sDOF = iDOF*DofPerNode;

        for (i=0; i<DofPerNode; i++)
        {
            sum =0.;
            for (k = 0; k < CoordNum; k++)
            {
                nb_DOF = p_SM->getNbDOF(iDOF,k);
                for (j=0; j<DofPerNode; j++)
                {
                    sum += (d[sDOF+j]-d[nb_DOF+j])*p_SM->getEntry(iDOF,k,i,j);
                }
            }
            q[sDOF+i] = sum;
        }
    }

    #pragma omp parallel for default (none) private (iDOF,i,j,k,sum,sDOF,nb_DOF,CoordNum) firstprivate(nodeLists) shared (p_SM)
    for (iDOF = freeNodeNum; iDOF < calNodeNum; iDOF++)
    {
        CoordNum = p_SM->getCoordNum(iDOF);
        sDOF = iDOF*DofPerNode;

        for (i=0; i<DofPerNode; i++)
        {
            sum = 0.;
            if (pDispNodeList[iDOF-freeNodeNum].Dir[i])
            {
                q[sDOF+i] = 0.;
            }
            else
            {
                for (k = 0; k < CoordNum; k++)
                {
                    nb_DOF = p_SM->getNbDOF(iDOF,k);

                    for (j=0; j<DofPerNode; j++)
                    {
                        sum += (d[sDOF+j]-d[nb_DOF+j])*p_SM->getEntry(iDOF,k,i,j);
                    }
                }

                q[sDOF+i] = sum;
            }
        }
    }
}

void PreCondConjGrad::moduleSolve (	double*							p_newDelta)
{
    double							alpha;
    double							beta;
    double							old_delta;


    alpha = *p_newDelta/(dot(&d,&q));					// alpha = delta_new/(d dot q)

    gupdateVecs_xrs	(alpha);						// x = x + alpha*d; r = r - alpha*q; s = M-1 * r
    old_delta = *p_newDelta;
    *p_newDelta = dot(&r,&s);							// new_delta = g dot g
    beta = *p_newDelta / old_delta;
    gupdateVec_d (beta);							// d = r + beta*d
}

void PreCondConjGrad::ginitialize	(	StiffnessMatrix*			p_SM,
                                        const bool 					isRelax)
{
    dualOut << "Ginitialize..." << endl;

    local_ttvNodeNum	= tNodeNum +tbNodeNum;
    local_tNodeNum		= tNodeNum +tbNodeNum;
    updateDOFSimple();		//From Base
    d.resize(DofPerNode*(calNodeNum+1),0.0);
    q.resize(DofPerNode*(calNodeNum+1),0.0);
    r.resize(DofPerNode*(calNodeNum+1),0.0);
    x.resize(DofPerNode*(calNodeNum+1),0.0);
    s.resize(DofPerNode*(calNodeNum+1),0.0);
    fillCompSM(p_SM);
    genPreCondJacobi(p_SM);
    // +1 for dummy DOF
    initDispVec_x_Simple(isRelax);
    initGradVec_r_Simple(isRelax);
    initGradVec_d();
}

void PreCondConjGrad::updateDOFSimple()
{
    calNodeList.clear();
    for (std::vector<unsigned>::iterator it = nodeLists.free.begin(); it!=nodeLists.free.end(); ++it) {
        calNodeList.push_back(*it);
    }

    freeNodeNum = nodeLists.free.size();
    dualOut << " calNodeList.size() = " << calNodeList.size() << endl;

    pDispNodeNum = nodeLists.restrain.size();
    dualOut << " pDispNodeNum = " << pDispNodeNum << endl;

    for (std::list<Restrain>::iterator rit=nodeLists.restrain.begin(); rit!=nodeLists.restrain.end(); ++rit) {
        calNodeList.push_back((*rit).NodeID);
        pDispNodeList.push_back(*rit);
    }
    dualOut << " calNodeList.size() becomes = " << calNodeList.size() << endl;

    calNodeNum = calNodeList.size();
    dumDOF = calNodeList.size()*DofPerNode;
    NodeIDtoDOF.resize(local_ttvNodeNum+1,dumDOF);

    for (unsigned i=0; i < calNodeList.size(); i++) {
        NodeIDtoDOF[calNodeList[i]] = DofPerNode*i;
    }
}

void PreCondConjGrad::updateDOFCompact()
{
	std::cout << "updateDOFCompact starts...." << std::endl;
    calNodeList.clear();
    for (std::vector<unsigned>::iterator it = nodeLists.free.begin(); it!=nodeLists.free.end(); ++it) {
        calNodeList.push_back(*it);
    }
    freeNodeNum = nodeLists.free.size();
    dualOut << " calNodeList.size() = " << calNodeList.size() << endl;
    pDispNodeNum = nodeLists.restrain.size();
    dualOut << " pDispNodeNum = " << pDispNodeNum << endl;
    pDispNodeList.clear();
    for (std::list<Restrain>::iterator rit=nodeLists.restrain.begin(); rit!=nodeLists.restrain.end(); ++rit) {
//    	cout << (*rit).NodeID << '\n';
    	calNodeList.push_back(rit->NodeID);
        pDispNodeList.push_back(*rit);
    }
    dualOut << " calNodeList.size() becomes = " << calNodeList.size() << endl;
    calNodeNum = calNodeList.size();
    dumDOF = calNodeList.size()*DofPerNode;
    NodeIDtoDOF.clear();
    NodeIDtoDOF.resize(local_ttvNodeNum+1,dumDOF);
    for (unsigned i=0; i < freeNodeNum; i++) {
        NodeIDtoDOF[calNodeList[i]] = DofPerNode*i;
    }
    unsigned DOF = DofPerNode*freeNodeNum;
    pDispDOF.clear();
    for (unsigned i=freeNodeNum; i < calNodeNum; i++) {
        NodeIDtoDOF[calNodeList[i]] = DOF;
        pDispDOF.push_back(DOF);
        DOF +=getFreeDOFNum(pDispNodeList[i-freeNodeNum].Dir);
    }
    DOF_Num = DOF;
    DOFtoNodeID.clear();
    DOFtoNodeID.resize(DOF_Num,tNodeNum);
    for (unsigned NodeID=0; NodeID<local_ttvNodeNum; NodeID++) {
        if (NodeIDtoDOF[NodeID]<DOFtoNodeID.size()){
            DOFtoNodeID[NodeIDtoDOF[NodeID]] = NodeID;
        }
    }
    std::cout << "updateDOFCompact Finished." << std::endl;
}

void PreCondConjGrad::fillConList()
{
    std::cout << "fillConList starts...." << std::endl;
    conList.clear();
    for (unsigned i=0; i<nodeLists.constrain.size(); i++) {
        unsigned NodeID = nodeLists.constrain[i].NodeID;
        IDAndInfo       con;
        std::vector<Restrain>::iterator it = std::find(pDispNodeList.begin(),pDispNodeList.end(),NodeID);
        if (it==pDispNodeList.end()) {
            for (unsigned k=0; k<nodeLists.constrain[i].conDisp.size(); k++) {
                con.first = NodeIDtoDOF[NodeID] + nodeLists.constrain[i].conDir[k];
                con.second = nodeLists.constrain[i].conDisp[k];
                conList.push_back(con);
            }
        } else {
            for (unsigned k=0; k<nodeLists.constrain[i].conDisp.size(); k++) {
                if (!isFullRestrain(it->Dir)) {
                    unsigned pos =0;
                    for (unsigned l=0; l<nodeLists.constrain[i].conDir[k]; l++) {
                        pos += !pDispNodeList[it-pDispNodeList.begin()].Dir[l];
                    }
                    con.first = NodeIDtoDOF[NodeID] + pos;
                    con.second = nodeLists.constrain[i].conDisp[k];
                    conList.push_back(con);
                }
            }
        }
    }
    std::cout << "fillConList Finished." << std::endl;
}

void PreCondConjGrad::simpleUpdateConList   (   const double        disp)
{
    for (unsigned i=0; i<conList.size(); i++) {
        conList[i].second = disp;
    }
}

void PreCondConjGrad::updateNodeStateSimple	(	const unsigned			step,
        										const bool 				isRelax)
{
    unsigned i,k,NodeID,DOF;
    double sum[3]= {0.0};
    double sumN[3]= {0.0};

    if (isRelax) {
        for (i = 0; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (k = 0; k< DofPerNode ; k++) {
                sumN[k] += (x[DOF+k] - GNodeTab[NodeID].d[k])/(GNodeTab[NodeID].d[k]+tiny);
                sum[k] += fabs(x[DOF+k] - GNodeTab[NodeID].d[k]);
                GNodeTab[NodeID].d[k] = x[DOF+k];
            }

        }
        cout << sum[0]/calNodeNum << ' '<< sum[1]/calNodeNum  << ' '<< sum[2]/calNodeNum << endl;
        cout << sumN[0]/calNodeNum << ' '<< sumN[1]/calNodeNum  << ' '<< sumN[2]/calNodeNum << endl;
    } else {
        for (i = 0; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (k = 0; k< DofPerNode ; k++) {
                GNodeTab[NodeID].d[k] += x[DOF+k];			//For spring model
            }
        }
    }
    updateOldDeltaSimple();
}

void PreCondConjGrad::updateNodeStateCompact	(	const unsigned			step,
        											const bool 				isRelax)
{
    unsigned i,d,NodeID,DOF;
    double sum[3]= {0.0};
    double sumN[3]= {0.0};

    if (isRelax) {
        for (i = 0; i < freeNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (unsigned k = 0; k< DofPerNode ; k++) {
                GNodeTab[NodeID].d[k] = x[DOF+k];
            }

        }

        for (i = freeNodeNum; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            d = 0;
            for (unsigned k = 0; k< DofPerNode ; k++) {
                if (!pDispNodeList[i-freeNodeNum].Dir[k]) {
                    GNodeTab[NodeID].d[k] = x[DOF+d];
                    d++;
                } else {
                    GNodeTab[NodeID].d[k] = 0.0;
                }
            }
        }

    } else {
        for (i = 0; i < freeNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (unsigned k = 0; k< DofPerNode ; k++) {
                GNodeTab[NodeID].d[k] += x[DOF+k];			//For spring model
            }
        }

        for (i = freeNodeNum; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            d = 0;
            for (unsigned k = 0; k< DofPerNode ; k++) {
                if (!pDispNodeList[i-freeNodeNum].Dir[k]) {
                    GNodeTab[NodeID].d[k] += x[DOF+d];			//For spring model
                    d++;
                }
            }
        }
    }

    updateOldDeltaCompact();
}

void PreCondConjGrad::updateOldDeltaSimple()
{
    unsigned i,k,NodeID,DOF;

    oldDelta_d.clear();
    oldDelta_d.resize(local_tNodeNum*DofPerNode,0.0);

    for (i=0; i < calNodeNum; i++)
    {
        NodeID = calNodeList[i];
        DOF = NodeIDtoDOF[NodeID];
        for (k = 0; k< DofPerNode ; k++)
        {
            oldDelta_d[NodeID*DofPerNode+k] = x[DOF+k];
        }
    }
}

void PreCondConjGrad::updateOldDeltaCompact()
{
    unsigned i,k,d,NodeID,DOF;

    oldDelta_d.clear();
    oldDelta_d.resize(local_tNodeNum*DofPerNode,0.0);

    for (i=0; i < freeNodeNum; i++)
    {
        NodeID = calNodeList[i];
        DOF = NodeIDtoDOF[NodeID];
        for (k = 0; k< DofPerNode ; k++)
        {
            oldDelta_d[NodeID*DofPerNode+k] = x[DOF+k];
        }
    }

    for (i=freeNodeNum; i < calNodeNum; i++)
    {
        NodeID = calNodeList[i];
        DOF = NodeIDtoDOF[NodeID];
        d = 0;
        for (k = 0; k< DofPerNode ; k++)
        {
            if (!pDispNodeList[i-freeNodeNum].Dir[k])
            {
                oldDelta_d[NodeID*DofPerNode+k] = x[DOF+d];
                d++;
            }
            else
            {
                oldDelta_d[NodeID*DofPerNode+k] = 0.0;
            }
        }
    }
}

void PreCondConjGrad::initDispVec_x_Simple	(	const bool 				isRelax)
{
    unsigned NodeID;
    if (isRelax)
    {
        #pragma omp parallel for if (Use_OpenMP) default (none) private (NodeID)
        for (unsigned i = 0; i < calNodeNum; i++)
        {
            NodeID = calNodeList[i];
            for (unsigned k=0; k<DofPerNode; k++) {
                x[i*DofPerNode+k] = 0.0;
            }
        }
    } else {
        #pragma omp parallel for if (Use_OpenMP) default (none) private (NodeID)
        for (unsigned i = 0; i < calNodeNum; i++)
        {
            NodeID = calNodeList[i];
            for (unsigned k=0; k<DofPerNode; k++)
            {
                x[i*DofPerNode+k] = 0.0;
            }
        }
    }
}

void PreCondConjGrad::initGradVec_d()			// b - Ax
{
    unsigned i;
    unsigned i_max = (calNodeNum)*DofPerNode;

    #pragma omp parallel for if (Use_OpenMP) default (none) private (i) firstprivate(i_max)
    for (i = 0; i < i_max; i++)
    {
        d[i] = m[i]*r[i];
    }
}

void PreCondConjGrad::initGradVec_r_Simple(	const bool				isRelax)
{
    unsigned i,j,NodeID,DOF;
    if (isRelax) {
        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF) firstprivate (Adj_ResForce) shared (GNodeTab)
        for (i = 0; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            for (j=0; j<DofPerNode; j++) {
                DOF = i*DofPerNode+j;
                r[DOF] = GNodeTab[NodeID].extF[j];
            }
        }
    } else {
        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF) firstprivate (Adj_ResForce) shared (GNodeTab)
        for (i = 0; i < freeNodeNum; i++) {
            NodeID = calNodeList[i];
            for (j=0; j<DofPerNode; j++) {
                DOF = i*DofPerNode+j;
                r[DOF] = GNodeTab[NodeID].extdF[j];
            }
        }

        #pragma omp parallel for if (Use_OpenMP) default (none) private (i,NodeID,j,DOF) shared (GNodeTab)
        for (i = freeNodeNum; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            for (j=0; j<DofPerNode; j++) {
                DOF = i*DofPerNode+j;
                r[DOF] = GNodeTab[NodeID].extdF[j];
            }
        }
    }
}

void PreCondConjGrad::genPreCondJacobi 	(	StiffnessMatrix*	p_SM)
{
    unsigned kStart,nbNum,kPos,iNbDOF;

    m.resize(DofPerNode*(calNodeNum+1),0.0);

    for (unsigned iDOF=0; iDOF<calNodeNum; iDOF++){
        nbNum = p_SM->getCoordNum(iDOF);
        kStart = p_SM->getStart(iDOF);

        for (unsigned k=0; k<kStart; k++) {
            kPos = p_SM->getNbPos(iDOF,k);
            iNbDOF = p_SM->getNbDOF(iDOF,k)/DofPerNode;

            for (unsigned i=0; i<Dim; i++) {
                m[iDOF*DofPerNode+i] += p_SM->getEntry(iNbDOF,kPos,i,i);
            }
        }

        for (unsigned k=0; k<nbNum-kStart; k++) {
            for (unsigned i=0; i<Dim; i++) {
                m[iDOF*DofPerNode+i] += p_SM->getEntry(iDOF,k,i,i);
            }
        }

        for (unsigned i=0; i<Dim; i++) {
            m[iDOF*DofPerNode+i] = 1.0/m[iDOF*DofPerNode+i];
        }
    }
}

void PreCondConjGrad::fillSMSimple	(	StiffnessMatrix*			p_SM)
{
    unsigned NodeID;
    unsigned nbNum;

    double	C,Cx,Cy,Cz;
    Mat<double> mat;

    p_SM->initialize(calNodeNum,freeNodeNum,&pDispNodeList);

    for (unsigned i=0; i<calNodeNum; i++)			//For spring model only
    {
        NodeID 	= calNodeList[i];
        nbNum 	= GNodeTab[NodeID].nbLatticeID.size();
        mat.resize(nbNum,2*DofPerNode);

        for (unsigned k=0; k<nbNum; k++)
        {
            adjSMCal(NodeID,k,&C,&Cx,&Cy,&Cz);
            p_SM->setNb(i,NodeIDtoDOF[GNodeTab[NodeID].nbNodeID[k]]);

            mat.set(k,0,C*Cx*Cx);
            mat.set(k,1,C*Cy*Cy);
            mat.set(k,2,C*Cz*Cz);
            mat.set(k,3,C*Cx*Cy);
            mat.set(k,4,C*Cx*Cz);
            mat.set(k,5,C*Cy*Cz);
        }
        p_SM->putMat(&mat);
    }
}

void PreCondConjGrad::fillCompSM	(		StiffnessMatrix*			p_SM)
{
    unsigned NodeID,iNbDOF;
    unsigned nbNum,cnt;
    bool isFind;

    double	C,Cx,Cy,Cz;
    Mat<double> mat;

    p_SM->initialize(calNodeNum,freeNodeNum,&pDispNodeList);

    for (unsigned i=0; i<calNodeNum; i++)			//For spring model only
    {
        NodeID 	= calNodeList[i];
        nbNum 	= GNodeTab[NodeID].nbNodeID.size();

        unsigned kStart	= GNodeTab[NodeID].start;

        p_SM->setStart(i,kStart);
        mat.resize(nbNum-kStart,2*DofPerNode);

        for (unsigned k=0; k<nbNum; k++) {
            p_SM->setNb(i,NodeIDtoDOF[GNodeTab[NodeID].nbNodeID[k]]);
        }

        for (unsigned k=kStart; k<nbNum; k++) {
            adjSMCal(NodeID,k,&C,&Cx,&Cy,&Cz);

            mat.set(k-kStart,0,C*Cx*Cx);
            mat.set(k-kStart,1,C*Cy*Cy);
            mat.set(k-kStart,2,C*Cz*Cz);
            mat.set(k-kStart,3,C*Cx*Cy);
            mat.set(k-kStart,4,C*Cx*Cz);
            mat.set(k-kStart,5,C*Cy*Cz);
        }
        p_SM->putMat(&mat);
    }

    for (unsigned i=0; i<calNodeNum; i++) {			//For spring model only
        unsigned kStart	= p_SM->getStart(i);
        for (unsigned k=0; k<kStart; k++) {
            iNbDOF = p_SM->getNbDOF(i,k)/DofPerNode;
            isFind = false;

            for (unsigned m=p_SM->getStart(iNbDOF); m<p_SM->getCoordNum(iNbDOF) ; m++) {
                if (i*DofPerNode==p_SM->getNbDOF(iNbDOF,m)) {
                    p_SM->setNbPos(i,m-p_SM->getStart(iNbDOF));
                    isFind = true;
                    break;
                }
            }

            if (!isFind) {
                p_SM->eraseNbDOF(i,k);
                p_SM->setStart(i,p_SM->getStart(i)-1);
                k--;
                kStart--;
            }
        }
    }
}

void PreCondConjGrad::adjSMCal		( 	const unsigned				NodeID,
                                        const unsigned 				dir,
                                        double*						p_C,
                                        double*						p_Cx,
                                        double*						p_Cy,
                                        double*						p_Cz)
{
    unsigned 	nbNodeID = GNodeTab[NodeID].nbNodeID[dir];
    unsigned 	LatticeID = GNodeTab[NodeID].nbLatticeID[dir];
    double		sqLen = 0.0;
    for (unsigned d = 0; d < Dim; d++) {
        sqLen += std::pow(( GNodeTab[NodeID].coord[d] - GNodeTab[nbNodeID].coord[d]),2);
    }
    *p_C  = LatTab[LatticeID].k[0]/sqLen;
    *p_Cx = (GNodeTab[nbNodeID].coord[0]-GNodeTab[NodeID].coord[0]);
    *p_Cy = (GNodeTab[nbNodeID].coord[1]-GNodeTab[NodeID].coord[1]);
    *p_Cz = (GNodeTab[nbNodeID].coord[2]-GNodeTab[NodeID].coord[2]);
}

void PreCondConjGrad::adjSMCal		( 	const unsigned				NodeID,
                                        const unsigned 				nbNodeID,
                                        const double				keChange,
                                        double*						p_C,
                                        double*						p_Cx,
                                        double*						p_Cy,
                                        double*						p_Cz)
{
    double		sqLen = 0.0;
    for (unsigned d = 0; d < Dim; d++) {
        sqLen += std::pow(( GNodeTab[NodeID].coord[d] - GNodeTab[nbNodeID].coord[d]),2);
    }
    *p_C  = keChange/sqLen;

    *p_Cx = (GNodeTab[nbNodeID].coord[0]-GNodeTab[NodeID].coord[0]);
    *p_Cy = (GNodeTab[nbNodeID].coord[1]-GNodeTab[NodeID].coord[1]);
    *p_Cz = (GNodeTab[nbNodeID].coord[2]-GNodeTab[NodeID].coord[2]);
}

void PreCondConjGrad::simpleSMCal		( 	const unsigned				NodeID,
        									const unsigned 				dir,
        									double*						p_C,
        									double*						p_Cx,
        									double*						p_Cy,
        									double*						p_Cz)
{
    unsigned nbNodeID,LatticeID;
    nbNodeID = GNodeTab[NodeID].nbNodeID[dir];
    LatticeID = GNodeTab[NodeID].nbLatticeID[dir];
    *p_C = LatTab[LatticeID].k[0]/LatTab[LatticeID].length/LatTab[LatticeID].length;
    *p_Cx = (GNodeTab[nbNodeID].coord[0]-GNodeTab[NodeID].coord[0]);
    *p_Cy = (GNodeTab[nbNodeID].coord[1]-GNodeTab[NodeID].coord[1]);
    *p_Cz = (GNodeTab[nbNodeID].coord[2]-GNodeTab[NodeID].coord[2]);
}

void PreCondConjGrad::gupdateVec_s()
{
    unsigned i;
    unsigned i_max = (calNodeNum)*DofPerNode;

    #pragma omp parallel for if (Use_OpenMP) default (none) private (i) firstprivate(i_max)
    for (i= 0; i< i_max; i++) {
        s[i] = m[i]*r[i];
    }
}

void PreCondConjGrad::gupdateVecs_xrs	(	double	alpha)
{
    unsigned i;
    unsigned i_max = (calNodeNum)*DofPerNode;


    #pragma omp parallel for if (Use_OpenMP) default (none) private (i) firstprivate(alpha,i_max)
    for (i= 0; i< i_max; i++)
    {
        x[i] = x[i] + alpha*d[i];
        r[i] = r[i] - alpha*q[i];
        s[i] = m[i]*r[i];
    }
}

void PreCondConjGrad::gupdateVec_d	(	double 				beta)
{
    unsigned i;
    unsigned i_max = (calNodeNum)*DofPerNode;

    #pragma omp parallel for if (Use_OpenMP) default (none) private (i) firstprivate(beta,i_max)
    for (i= 0; i< i_max; i++)
    {
        d[i] = s[i] + beta*d[i];
    }
}

void PreCondConjGrad::newUpdateVec_q 	(	StiffnessMatrix*	p_SM)
{
    unsigned 	i,j,k,iDOF;
    unsigned	sDOF,nb_DOF;
    unsigned	CoordNum;		//coordination no
    double		sum;

    #pragma omp parallel for if (Use_OpenMP) default (none) private (iDOF,i,j,k,sum,sDOF,nb_DOF,CoordNum) shared (p_SM)
    for (iDOF = 0; iDOF < freeNodeNum; iDOF++)
    {
        CoordNum = p_SM->getCoordNum(iDOF);
        sDOF = iDOF*DofPerNode;

        for (i=0; i<DofPerNode; i++)
        {
            sum =0.;
            for (k = 0; k < CoordNum; k++)
            {
                nb_DOF = p_SM->getNbDOF(iDOF,k);
                for (j=0; j<DofPerNode; j++)
                {
                    sum += (d[sDOF+j]-d[nb_DOF+j])*p_SM->getEntry(iDOF,k,i,j);		// Error
                }
            }
            q[sDOF+i] = sum;
        }
    }

    #pragma omp parallel for if (Use_OpenMP) default (none) private (iDOF,i,j,k,sum,sDOF,nb_DOF,CoordNum) firstprivate(nodeLists) shared (p_SM)
    for (iDOF = freeNodeNum; iDOF < calNodeNum; iDOF++)
    {
        CoordNum = p_SM->getCoordNum(iDOF);
        sDOF = iDOF*DofPerNode;

        for (i=0; i<DofPerNode; i++)
        {
            sum = 0.;
            if (pDispNodeList[iDOF-freeNodeNum].Dir[i])
            {
                q[sDOF+i] = 0.;
            }
            else
            {
                for (k=0; k<CoordNum; k++)
                {
                    nb_DOF = p_SM->getNbDOF(iDOF,k);
                    for (j=0; j<DofPerNode; j++)
                    {
                        sum += (d[sDOF+j]-d[nb_DOF+j])*p_SM->getEntry(iDOF,k,i,j);
                    }
                }
                q[sDOF+i] = sum;
            }
        }
    }
}

void PreCondConjGrad::tmpUpdateNodeState	(	const bool 				isRelax)
{
    unsigned NodeID,DOF;

    if (isRelax) {
        for (unsigned i = 0; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (unsigned k = 0; k< DofPerNode ; k++) {
                GNodeTab[NodeID].d[k] = x[DOF+k];
            }
        }
    } else {
        for (unsigned i = 0; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (unsigned k = 0; k< DofPerNode ; k++) {
                GNodeTab[NodeID].d[k] = x[DOF+k];
            }
        }
    }
}

void PreCondConjGrad::undoTmpUpdateNodeState	(	const bool 				isRelax)
{
    unsigned NodeID,DOF;

    if (isRelax)
    {
        for (unsigned i = 0; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (unsigned k = 0; k< DofPerNode ; k++) {
                GNodeTab[NodeID].d[k] -= x[DOF+k] - oldDelta_d[NodeID*DofPerNode+k];			//For spring model
            }
        }
    } else {
        for (unsigned i = 0; i < calNodeNum; i++) {
            NodeID = calNodeList[i];
            DOF = NodeIDtoDOF[NodeID];
            for (unsigned k = 0; k< DofPerNode ; k++) {
                GNodeTab[NodeID].d[k] -= x[DOF+k];			//For spring model
            }
        }
    }
}

unsigned	PreCondConjGrad::getSizeOfData	()
{
    return sizeof(double)*(d.size() + q.size() + x.size() + r.size() + s.size() + m.size());
}

void PreCondConjGrad::PCGclean ()
{
    r.clear();
    x.clear();
    q.clear();
    d.clear();
    m.clear();	//Preconditional vector
    s.clear();


    freeNodeList.clear();
    pDispNodeList.clear();
    calNodeList.clear();
    NodeIDtoDOF.clear();
    pDispDOF.clear();
}

double PreCondConjGrad::dot (vector<double> *p_v1,vector<double> *p_v2)
{
    unsigned i,i_max;
    double out = 0.;

    i_max = min(p_v1->size(),p_v2->size());

    #pragma omp parallel for if (Use_OpenMP) default (none) reduction (+:out) firstprivate(i_max) private (i) shared (p_v1,p_v2)
    for (i=0; i< i_max; i++)
    {
        out += (*p_v1)[i]*(*p_v2)[i];
    }
    return out;
}
