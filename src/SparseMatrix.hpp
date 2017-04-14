/*
 * SparseMatrix.hpp
 *
 *  Created on: Feb 8, 2014
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef SPARSEMATRIX_HPP_
#define SPARSEMATRIX_HPP_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "MiscTools.hpp"
#include <eigen3/Eigen/Sparse>

typedef	std::vector<Mat<double> > SparseMat;

class SparseMatrix {
public:
    SparseMatrix();
    virtual ~SparseMatrix();

    void	setEntry		(	const unsigned	iDOF,
                                const unsigned	k,
                                const unsigned	i,
                                const unsigned	j,
                                const double	var)
    {
        sSM[iDOF].set(k,index[i*DofPerNode+j],var);
    }

    void	addEntry		(	const unsigned	iDOF,
                                const unsigned	k,
                                const unsigned	i,
                                const unsigned	j,
                                const double	var)
    {
        sSM[iDOF].add(k,index[i*DofPerNode+j],var);
    }

    void	multiEntry		(	const unsigned	iDOF,
                                const unsigned	k,
                                const unsigned	i,
                                const unsigned	j,
                                const double	var)
    {
        sSM[iDOF].multi(k,index[i*DofPerNode+j],var);
    }

    void	divEntry		(	const unsigned	iDOF,
                                const unsigned	k,
                                const unsigned	i,
                                const unsigned	j,
                                const double	var)
    {
        sSM[iDOF].div(k,index[i*DofPerNode+j],var);
    }


    void 			putMat			(	Mat<double>*	p_mat)
    {
        sSM.push_back(*p_mat);
    }

    double 			getEntry		(	const unsigned	iDOF,
                                        const unsigned	k,
                                        const unsigned	i,
                                        const unsigned	j)
    {
        return sSM[iDOF].get(k,index[i*DofPerNode+j]);
    }

    unsigned		size			()
    {
        return sSM.size();
    };

    unsigned	getRowNum	(	unsigned	iDOF)
    {
        return sSM[iDOF].getRowNum();
    }

    SparseMat*	getSMptr	()	{
        return &sSM;
    }

    void	insertRow		(	const unsigned					iDOF,
                                const std::vector<double>*		p_inVec);


protected:
    SparseMat									sSM;
    std::array<unsigned,DofPerNode*DofPerNode >	index;
private:
};

class StiffnessMatrix : public SparseMatrix
{
public:
    StiffnessMatrix			() {};

    ~StiffnessMatrix		() {};

    inline StiffnessMatrix			(	const unsigned	NodeSize,
                                        const unsigned	LatSize);

    void 		initialize		(	const unsigned					NodeSize,
                                    const unsigned					freeNodeNum,
                                    const std::vector<Restrain>*	p_pDispNodeList);

    void 		MatVecProduct	(	const std::vector<double>*		p_d,
                                    std::vector<double>*			p_q);

    void 		setNb		(	unsigned		iDOF,
                                unsigned		nb_DOF) {
        nbDOF[iDOF].push_back(nb_DOF);
    }

    unsigned	getNbDOF	(	unsigned		iDOF,
                                unsigned		k)		{
        return nbDOF[iDOF][k];
    }

    void		rmNbDOF		(	unsigned		iDOF,
                                unsigned		rmDOF)	{
        nbDOF[iDOF].erase(std::remove(nbDOF[iDOF].begin(),nbDOF[iDOF].end(),rmDOF),nbDOF[iDOF].end());
    }

    void		eraseNbDOF	(	unsigned		iDOF,
                                unsigned		pos)	{
        nbDOF[iDOF].erase(nbDOF[iDOF].begin()+pos);
    }

    void		setStart 	(	unsigned		iDOF,
                                unsigned		_start)	{
        start[iDOF] = _start;
    }

    unsigned	getStart 	(	unsigned		iDOF)	{
        return start[iDOF];
    }

    void		setNbPos	(	unsigned		iDOF,
                                unsigned		k)		{
        nbPos[iDOF].push_back(k);
    }

    unsigned	getNbPos	(	unsigned		iDOF,
                                unsigned		k)		{
        return nbPos[iDOF][k];
    }

    unsigned	getCoordNum	(	unsigned		iDOF)	{
        return nbDOF[iDOF].size();
    }

    unsigned	getSizeOf	();

    double		outInitDelta	();

    void		printAll	(	unsigned		outSize);

private:

    std::vector<std::vector<unsigned> >			nbDOF;
    std::vector<unsigned>						start;
    std::vector<std::vector<unsigned> >			nbPos;

    std::vector<Restrain>						pDispNodeList;

    unsigned									sizeFree;
    unsigned									sizeTotal;

};


class Preconditioner : public SparseMatrix
{
public:

    Preconditioner			() {};
    ~Preconditioner			() {};

    void		initialize		();

    void		copyAll				(	SparseMat*			p_SMat);
    void		convertFromSM		(	StiffnessMatrix*	p_SM);
    void		convertFromCompSM	(	StiffnessMatrix*	p_SM);
    void 		genPreCondIC		(	StiffnessMatrix*	p_SM);

    unsigned	getSizeOf		();

private:
};

#endif /* SPARSEMATRIX_HPP_ */
