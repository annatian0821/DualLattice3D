/*
 * SparseMatrix.cpp
 *
 *  Created on: Feb 8, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "SparseMatrix.hpp"

using namespace std;

SparseMatrix::SparseMatrix() {
    // TODO Auto-generated constructor stub

}

SparseMatrix::~SparseMatrix() {
    // TODO Auto-generated destructor stub
}


void	SparseMatrix::insertRow	(	const unsigned					iDOF,
                                    const vector<double>*			p_inVec)
{
    sSM[iDOF].insert(p_inVec);
}

void 	StiffnessMatrix::initialize		(	const unsigned				NodeSize,
        const unsigned				freeNodeNum,
        const vector<Restrain>*		p_pDispNodeList)
{
    sizeTotal = NodeSize;
    sizeFree = freeNodeNum;

    unsigned size = p_pDispNodeList->size();

    pDispNodeList.clear();

    for (unsigned i=0; i<size; i++)
    {
        pDispNodeList.push_back((*p_pDispNodeList)[i]);
    }

    nbDOF.resize(NodeSize+1);
    start.resize(NodeSize+1);
    nbPos.resize(NodeSize+1);

    index[0*DofPerNode+0] = 0;
    index[0*DofPerNode+1] = 3;
    index[0*DofPerNode+2] = 4;
    index[1*DofPerNode+0] = 3;
    index[1*DofPerNode+1] = 1;
    index[1*DofPerNode+2] = 5;
    index[2*DofPerNode+0] = 4;
    index[2*DofPerNode+1] = 5;
    index[2*DofPerNode+2] = 2;
}

void 	StiffnessMatrix::MatVecProduct	(	const vector<double>*		p_d,
        vector<double>*				p_q)
{
    unsigned 	iDOF,i,j,k;
    unsigned	sDOF,nb_DOF,kStart,kPos;
    unsigned	CoordNum;		//coordination no
    double		sum;

    Clock							lClock;

    #pragma omp parallel for default (none) private (iDOF,i,j,k,sum,sDOF,nb_DOF,CoordNum,kStart,kPos) shared (p_d,p_q)
    for (iDOF = 0; iDOF < sizeFree; iDOF++)
    {
        CoordNum = nbDOF[iDOF].size();  		//getCoordNum(iDOF);
        sDOF = iDOF*DofPerNode;

        kStart = start[iDOF]; 					//p_SM->getStart(iDOF);
        for (i=0; i<DofPerNode; i++)
        {
            sum =0.;
            for (k = 0; k < kStart; k++)
            {
                nb_DOF = nbDOF[iDOF][k]; 		//p_SM->getNbDOF(iDOF,k);
                kPos = nbPos[iDOF][k];			//p_SM->getNbPos(iDOF,k);
                for (j=0; j<DofPerNode; j++)
                {
                    sum += ((*p_d)[sDOF+j]-(*p_d)[nb_DOF+j])*sSM[nb_DOF/DofPerNode].get(kPos,index[i*DofPerNode+j]);  		//p_SM->getEntry(nb_DOF/DofPerNode,kPos,i,j);
                }
            }

            for (k = 0; k < CoordNum-kStart; k++)
            {
                nb_DOF = nbDOF[iDOF][k+kStart]; //p_SM->getNbDOF(iDOF,k+kStart);
                for (j=0; j<DofPerNode; j++)
                {
                    sum += ((*p_d)[sDOF+j]-(*p_d)[nb_DOF+j])*sSM[iDOF].get(k,index[i*DofPerNode+j]);	//sum += (d[sDOF+j]-d[nb_DOF+j])*p_SM->getEntry(iDOF,k,i,j);
                }
            }
            (*p_q)[sDOF+i] = sum;
        }
    }

    #pragma omp parallel for default (none) private (iDOF,i,j,k,sum,sDOF,nb_DOF,CoordNum,kStart,kPos) shared (p_d,p_q) //firstprivate (freeNodeNum,calNodeNum)
    for (iDOF = sizeFree; iDOF < sizeTotal; iDOF++)
    {
        CoordNum = nbDOF[iDOF].size();
        sDOF = iDOF*DofPerNode;
        kStart = start[iDOF];

        for (i=0; i<DofPerNode; i++)
        {
            sum = 0.;
            if (pDispNodeList[iDOF-sizeFree].Dir[i])
            {
                (*p_q)[sDOF+i] = 0.;
            }
            else
            {
                for (k = 0; k < kStart; k++)
                {
                    nb_DOF = nbDOF[iDOF][k]; 				//p_SM->getNbDOF(iDOF,k);
                    kPos = nbPos[iDOF][k];					//p_SM->getNbPos(iDOF,k);
                    for (j=0; j<DofPerNode; j++)
                    {
                        sum += ((*p_d)[sDOF+j]-(*p_d)[nb_DOF+j])*sSM[nb_DOF/DofPerNode].get(kPos,index[i*DofPerNode+j]);		//p_SM->getEntry(nb_DOF/DofPerNode,kPos,i,j);
                    }
                }

                for (k=0; k<CoordNum-kStart; k++)
                {
                    nb_DOF = nbDOF[iDOF][k+kStart]; 		//p_SM->getNbDOF(iDOF,k+kStart); //p_sSM->getNbDOF(iDOF,k);
                    for (j=0; j<DofPerNode; j++)
                    {
                        sum += ((*p_d)[sDOF+j]-(*p_d)[nb_DOF+j])*sSM[iDOF].get(k,index[i*DofPerNode+j]);			//p_SM->getEntry(iDOF,k,i,j);
                    }
                }
                (*p_q)[sDOF+i] = sum;
            }
        }
    }
}

unsigned	StiffnessMatrix::getSizeOf	()
{
    unsigned i;
    unsigned cntUnsigned = 0, cntDouble = 0;

    cntUnsigned = index.size();
    for (i=0; i<nbDOF.size(); i++)
    {
        cntUnsigned += nbDOF[i].size();
    }

    for (i=0; i<sSM.size(); i++)
    {
        cntDouble += sSM[i].getSize();
    }

    return cntUnsigned*sizeof(unsigned)+cntDouble*sizeof(double);
}

void	StiffnessMatrix::printAll	(	unsigned outSize)
{
    unsigned i,j;
    char	fileName[255];

    sprintf(fileName,"%s/SM.txt",OutputSubFolder);

    ofstream SMFile(fileName); //, ios::out | ios::app);

    for (i=0; i<outSize; i++)
    {
        SMFile << "--------iDOF = " << i << " --------------"<< endl;
        SMFile << "start =  " << start[i] << endl;
        SMFile << "inbDOF : ";
        for (j=0; j<nbDOF[i].size(); j++)
        {
            SMFile << nbDOF[i][j]/3 << " ";
        }

        SMFile << endl;
        SMFile << "nbPos : ";
        for (j=0; j<nbPos[i].size(); j++)
        {
            SMFile << nbPos[i][j] << " ";
        }

        SMFile << endl;
        SMFile << "sSM row size = " << sSM[i].getRowNum() << endl;

        if (nbPos[i].size()!=start[i])
        {
            cout << "nbPos.size!=start, i = " << i << endl;
        }

        if (nbDOF[i].size()!=nbPos[i].size()+sSM[i].getRowNum())
        {
            cout << "nbDOF.size!=nbPos.size()+sSM.getRow, i = " << i << endl;
        }
    }
}

void 	Preconditioner::initialize		()
{
    index[0*3+0] = 0;
    index[0*3+1] = 1;
    index[0*3+2] = 2;
    index[1*3+0] = 3;
    index[1*3+1] = 4;
    index[1*3+2] = 5;
    index[2*3+0] = 6;
    index[2*3+1] = 7;
    index[2*3+2] = 8;
};

unsigned	Preconditioner::getSizeOf	()
{
    unsigned i;
    unsigned cntDouble = 0;

    for (i=0; i<sSM.size(); i++)
    {
        cntDouble += sSM[i].getSize();
    }

    return cntDouble*sizeof(double);
}

void		Preconditioner::copyAll			(	SparseMat*		p_SMat)
{
    unsigned i;
    for (i=0; i<p_SMat->size(); i++)
    {
        sSM.push_back((*p_SMat)[i]);
    }
}

void		Preconditioner::convertFromSM	(	StiffnessMatrix*		p_SM)
{
    Mat<double>	mat;
    double var;
    unsigned iDOF,i,j,k,nNum;
    for (iDOF=0; iDOF<p_SM->size(); iDOF++)
    {
        nNum = p_SM->getCoordNum(iDOF);
        mat.resize(nNum+1,Dim*Dim);		//All set to zero also

        for (k=1; k<nNum+1; k++)
        {
            for (i=0; i<Dim; i++)
            {
                for (j=0; j<Dim; j++)
                {
                    var = p_SM->getEntry(iDOF,k-1,i,j);
                    mat.set(k,i*Dim+j,var);
                    mat.add(0,i*Dim+j,var);
                }
            }
        }
        sSM.push_back(mat);
    }
}

void		Preconditioner::convertFromCompSM	(	StiffnessMatrix*		p_SM)
{
    double var;
    unsigned iDOF,i,j,k,nNum,kPos,kStart,iNbDOF;

    Mat<double>	mat;
    for (iDOF=0; iDOF<p_SM->size(); iDOF++)
    {
        nNum = p_SM->getCoordNum(iDOF);
        kStart = p_SM->getStart(iDOF);
        mat.resize(nNum-kStart+1,Dim*Dim);		//All set to zero also

        for (k=1; k<kStart+1; k++)
        {
            kPos = p_SM->getNbPos(iDOF,k-1);
            iNbDOF = p_SM->getNbDOF(iDOF,k-1)/DofPerNode;
            for (i=0; i<Dim; i++)
            {
                for (j=0; j<Dim; j++)
                {
                    mat.add(0,i*Dim+j,p_SM->getEntry(iNbDOF,kPos,i,j));
                }
            }
        }

        for (k=1; k<nNum-kStart+1; k++)
        {
            for (i=0; i<Dim; i++)
            {
                for (j=0; j<Dim; j++)
                {
                    var = p_SM->getEntry(iDOF,k-1,i,j);
                    mat.set(k,i*Dim+j,var);
                    mat.add(0,i*Dim+j,var);
                }
            }
        }
        sSM.push_back(mat);
    }
}

void 	Preconditioner::genPreCondIC					(	StiffnessMatrix*	p_SM)
{
    unsigned i,j,k,ii,jj,kk,nNum;
    vector<double>		vec;
    double				Akk,Ajk;
    Clock				clock;

    unsigned			calNodeNum = p_SM->size();
    clock.start("Incomplete Cholesky");

    initialize();
    convertFromCompSM(p_SM);

    #pragma omp parallel for firstprivate (calNodeNum) shared (p_SM) // schedule (static chunkSize)
    for (k=0; k<calNodeNum; k++)
    {
        for (kk=0; kk<Dim; kk++)
        {
            Akk = sqrt(sSM[k].get(0,kk*Dim+kk)); 	//(p_PC->getEntry(k,0,kk,kk));
            sSM[k].set(0,kk*Dim+kk,Akk);			//p_PC->setEntry(k,0,kk,kk,Akk);
            for (ii=kk+1; ii<Dim; ii++)
            {
                sSM[k].div(0,ii*Dim+kk,Akk);		//p_PC->divEntry(k,0,ii,kk,Akk);
            }

            nNum = sSM[k].getRowNum();				//p_PC->getRowNum(k);
            for (i=1; i<nNum; i++)
            {
                for (ii=0; ii<Dim; ii++)
                {
                    sSM[k].div(i,ii*Dim+kk,Akk);	//p_PC->divEntry(k,i,ii,kk,Akk);
                }
            }

            for (jj=kk+1; jj<Dim; jj++)
            {
                Ajk = sSM[k].get(0,jj*Dim+kk);		// Ajk = p_PC->getEntry(k,0,jj,kk);
                for (ii=jj; ii<Dim; ii++)
                {
                    sSM[k].add(0,ii*Dim+jj,-Ajk*sSM[k].get(0,ii*Dim+kk));		//	p_PC->addEntry(k,0,ii,jj,-Ajk*p_PC->getEntry(k,0,ii,kk));
                }
                nNum = sSM[k].getRowNum();  		//p_PC->getRowNum(k);
                for (i=1; i<nNum; i++)
                {
                    for (ii=0; ii<Dim; ii++)
                    {
                        sSM[k].add(i,ii*Dim+jj,-Ajk*sSM[k].get(i,ii*Dim+kk));				//p_PC->addEntry(k,i,ii,jj,-Ajk*p_PC->getEntry(k,i,ii,kk));
                    }
                }
            }

            for (j=k+1; j<calNodeNum; j++)
            {
                for (jj=0; jj<Dim; jj++)
                {
                    Ajk = sSM[j].get(0,jj*Dim+kk);															//Ajk = p_PC->getEntry(j,0,jj,kk);
                    nNum = sSM[j].getRowNum();																	// nNum = p_PC->getRowNum(j);
                    for (i=0; i<nNum; i++)			//should start from nbNodeID > NodeID
                    {
                        for (ii=0; ii<Dim; ii++)
                        {
                            sSM[j].add(i,ii*Dim+jj,-Ajk*sSM[j].get(i,ii*Dim+kk));								//p_PC->addEntry(j,i,ii,jj,-Ajk*p_PC->getEntry(j,i,ii,kk));
                        }
                    }
                }
            }
        }
    }
    clock.get("Finished");
}
