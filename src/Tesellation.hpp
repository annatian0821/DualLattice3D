/*
 * Tesellation.h
 *
 *  Created on: Nov 18, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef TESELLATION_H_
#define TESELLATION_H_

#include "LEM3D.hpp"
#include "Dataclass.hpp"
#include "MiscTools.hpp"
#include "Vertices.hpp"
#include "Statistics.hpp"
//#include "truncated_normal.hpp"
//#include "log_normal_truncated_ab.hpp"

#include <random>
#include "voro++.hh"

typedef CustomPair<unsigned,unsigned>	NodeIDPair;
class VoroTesell {
public:
    VoroTesell();
    virtual ~VoroTesell();

    void                        generate	(	const double			minDist,
    									        const double			scaleBoundary,
    									        double					outerScaleLength,
    									        unsigned*				p_maxCoordNo,
    									        Vertices*				p_verti);

    std::vector<unsigned> 		getPreExCrackLatID ();

    std::vector<unsigned>       getWeakPlaneLatID ();

    std::array<unsigned,2>		getPxCrackPerimeterNum ();	//return how many node around perimeter of pre-existing crack ;

    std::vector<unsigned>       getBoreHoleNodeList ();

    std::array<unsigned,2>      getRefineNode() {
        return refineNode;
    }

private:
    unsigned					Nnum;
    double                      coarseness;
    std::array<unsigned,2>      refineNode;
    std::array<bool,6>           isRegFace;
    unsigned                    GeoID;
    unsigned					Npx;
    unsigned					Npy;
    unsigned					Npz;
    unsigned                    nx,ny,nz;
    unsigned                    ex,ey,ez;
    Point                       bottomLocation;
    double                      BHdepth;
    std::vector<unsigned>		preExistingCrack;
    std::vector<unsigned>       weakPlaneLatID;
    std::vector<unsigned>       boreHoleNodeList;
    std::vector<Point>          tmpNodeList;
    unsigned					bStartEnd[6][2];
    std::array<unsigned,2>		pxPerimeterNum;
    std::vector<unsigned>       unBreakableNodeList;
    unsigned                    BHNodeNum;
    double                      scaleX,scaleY,scaleZ;
    double						minD;
//    double                      avgD;
    double                      notchLen;
    double                      strengthAdj;
    double                      anisotropy;
    struct			VoroInfo	{
    public:
        VoroInfo				();

        VoroInfo				(	const unsigned			nNum)
        {
            nbList.resize(nNum);
            area.resize(nNum);
//			vectices.resize(nNum);
            vCoord.resize(nNum);
            vList.resize(nNum);
            vListIndex.resize(nNum);
        };
        ~VoroInfo				() {};
        std::vector<unsigned>								IDList;
        std::vector<double>									vol;
        std::vector<std::array<double,3> >					centroid;
        std::vector<std::vector<int> >						nbList;
        std::vector<std::vector<double> >					area;
        std::vector<std::vector<double> >					vCoord;
        std::vector<std::vector<std::vector<int> > >		vList;
        std::vector<std::vector<int> >						vListIndex;
    };

    VoroInfo*											p_voroInfo;
    std::vector<NodeIDPair>								weakPlane;
    std::vector<NodeIDPair>                             unBreakablePair;

    void initBoundary                       ();

    std::vector<NodeIDPair> genRanNode		(	voro::container_poly* 	p_conp,
                                    			unsigned* 				p_NodeID);

    void genRanNode							(	unsigned* 				p_NodeID,
    											voro::container_poly* 	p_conp,
    											Mat3Dvec*				p_NodeInPartition,
    											const Point 			centre,
    											const Vec     			normVec,
    											const double			lt,
    											const double			ld,
    											const double			weakPlaneRadius,
    											const bool				isConjNodeListEmpty);

    void genRanNodeReduced                   (   unsigned*               p_NodeID,
                                                voro::container_poly*   p_conp,
                                                Mat3Dvec*               p_NodeInPartition,
                                                const Point             centre,
                                                const Vec               normVec,
                                                const double            lt,
                                                const double            ld,
                                                const double            weakPlaneRadius,
                                                const unsigned          dir,
                                                double                  rRatio,
                                                const bool              isConjNodeListEmpty);

    void genRanNodeReduced                   (   unsigned*           p_NodeID,
                                                 voro::container_poly*     p_conp,
                                                 Mat3Dvec*           p_NodeInPartition,
                                                 const Point         centre,
                                                 const Vec           normVec,
                                                 const double        lt,
                                                 const double        ld,
                                                 const double        weakPlaneRadius,
                                                 const unsigned      dir,
                                                 double              rRatio,
                                                 double              anisotropy,
                                                 const bool          isConjNodeListEmpty);

    void genRanNodeReduced                  (   unsigned*               p_NodeID,
                                                voro::container_poly*   p_conp,
                                                Mat3Dvec*               p_NodeInPartition,
                                                const double            lmin,
                                                const unsigned          dir,
                                                double                  rRatio);

    void genRanNodeAnisotropic              (	unsigned* 				p_NodeID,
    											voro::container_poly* 	p_conp,
    											double					aniZone);

    std::vector<NodeIDPair> genRegNodes     (   voro::container_poly*      p_conp,
                                                unsigned*                  p_NodeID);

    void genSimpleCubicNodes				(	voro::container_poly*      p_conp,
    											unsigned* 					p_NodeID);

    void genBCC_Nodes						(	voro::container_poly* 		p_conp,
    											unsigned* 					p_NodeID);

    void genFCC_Nodes						(	voro::container_poly* 		p_conp,
    											unsigned* 					p_NodeID);

    void genHCP_Nodes                       (   voro::container_poly*       p_conp,
                                                unsigned*                   p_NodeID);

    void genRegBoundaryNodes				(	voro::container_poly* 		p_conp,
    											unsigned* 					p_NodeID);

    std::vector<NodeIDPair> genPennyCrack	(	const Point				centre,
    											const double			radius,
    											unsigned* 				p_NodeID,
    											voro::container_poly* 	p_conp,
    											Mat3Dvec*				p_NodeInPartition);

    std::vector<NodeIDPair> genPennyCrackReg	(	const Point				centre,
    												const double			radius,
    												const double            relaxD,
    												const double            meanZ,
    												unsigned* 				p_NodeID,
    												voro::container_poly* 	p_conp,
    												Mat3Dvec*				p_NodeInPartition);

    std::vector<NodeIDPair> genPennyCrackRan    (   const Point             centre,
                                                    const double            radius,
                                                    const double            ld,
                                                    const double            lt,
                                                    unsigned*               p_NodeID,
                                                    voro::container_poly*   p_conp,
                                                    Mat3Dvec*               p_NodeInPartition);

   bool genInclinedPennyCrackRan                (   const Point                 centre,
                                                    const double                radius,
                                                    Vec                         normVec,
                                                    const double                ld,
                                                    const double                lt,
                                                    unsigned*                   p_NodeID,
                                                    std::vector<NodeIDPair>*    p_conjNodeList,
                                                    voro::container_poly*       p_conp,
                                                    Mat3Dvec*                   p_NodeInPartition);

    void genIncOuterPennyCrackRan               (   const Point             centre,
                                                    const double            inRadius,
                                                    const double            outRadius,
                                                    Vec                     normVec,
                                                    const double            ld,
                                                    const double            lt,
                                                    unsigned*               p_NodeID,
                                                    voro::container_poly*   p_conp,
                                                    Mat3Dvec*               p_NodeInPartition);

    void  genOuterPennyCrack					(	const Point				centre,
    												const double			crackRadius,
    												const double			outRadius,
    												const double            relaxD,
    												const double            meanZ,
    												unsigned* 				p_NodeID,
													voro::container_poly* 	p_conp,
    												Mat3Dvec*				p_NodeInPartition);

    void genRoughPennyCrack						(	const Point				centre,
    												const double			crackRadius,
    												const double            lt);

    void genNotches                             (   const double                width,
                                                    const double                ld,
                                                    const double                lt,
                                                    unsigned*                   p_NodeID,
                                                    std::vector<NodeIDPair>*    p_conjNodeList,
                                                    voro::container_poly*       p_conp,
                                                    Mat3Dvec*                   p_NodeInPartition);

    void genNotchWeakPlane                      (   const Point             centre,
                                                    const double            notchLen,
                                                    const double            lt,
                                                    unsigned*               p_NodeID,
                                                    voro::container_poly*   p_conp,
                                                    Mat3Dvec*               p_NodeInPartition);

    void genBorehole                            (   const Point                 centre,
                                                    const double                radius,
                                                    const double                t,
                                                    const unsigned              resolution,
                                                    Vec                         normVec,
                                                    unsigned*                   p_NodeID,
                                                    voro::container_poly*       p_conp,
                                                    Mat3Dvec*                   p_NodeInPartition,
                                                    std::vector<NodeIDPair>*    p_conjNodeList);

    void genBorehole2                           (    const Point                 centre,
                                                     const double                radius,
                                                     const double                t,
                                                     const unsigned              resolution,
                                                     Vec                         normVec,
                                                     unsigned*                   p_NodeID,
                                                     voro::container_poly*       p_conp,
                                                     Mat3Dvec*                   p_NodeInPartition,
                                                     std::vector<NodeIDPair>*    p_conjNodeList);

    void genOuterNode							(	voro::container_poly* 	p_conp,
    												unsigned* 				p_NodeID,
													double					scaleLength,
													const double			scaleBoundary,
													const unsigned			layer);

    void genOuterNode                           (   voro::container_poly*   p_conp,
                                                    unsigned*               p_NodeID,
                                                    const double            Lmin,
                                                    double                  scaleLength,
                                                    const double            scaleBoundary,
                                                    unsigned                layer);

    void genOuterNodeAni						(	voro::container_poly* 	p_conp,
    												unsigned* 				p_NodeID,
    												double					scaleLength,
    												const double			scaleBoundary,
    												const unsigned			layer);

    void inputOuterBoundary		                (	voro::container_poly* 	p_conp,
                                                    unsigned* 				p_NodeID,
                                                    double					scaleLength,
                                                    const double			scaleBoundary,
                                                    std::vector<double>*	p_Di);

    Mat3Dvec  getNodeInPartition                (   const double            lmin);

    bool testInside (const Point                        rTest)
    {
        return (((rTest[0]>0.0)&&(rTest[0]<Nx*UnitLength)
                ||(rTest[1]>0.0)&&(rTest[1]<Ny*UnitLength)
                ||(rTest[2]>0.0)&&(rTest[2]<Nz*UnitLength)));
    }

    unsigned testOverlap		(	Mat3Dvec* 						p_3Dvec,
                                    const double                    lmin,
                                    const Point						rTest,
                                    const unsigned					NodeID,
                                    const std::array<unsigned,Dim> 	p);

   unsigned testOverlap        (    Mat2Dvec*                           p_2Dvec,
                                    const Point                         rTest,
                                    const unsigned                      NodeID,
                                    const std::array<unsigned,2>        p);

    unsigned testOverlapAni    (    Mat3Dvec*                           p_3Dvec,
                                    const double                        lmin,
                                    const double                        anistropy,
                                    const unsigned                      dir,
                                    const Point                         rTest,
                                    const unsigned                      NodeID,
                                    const std::array<unsigned,Dim>      p);

    bool isRestrictedZone       (   const Point         testPoint,
                                    const Point         centre,
                                    const Vec           normVec,
                                    const double        offsetTol,
                                    const double        radius);

    bool isRestrictedZone       (   const Point         testPoint,
                                    const double        width,
                                    const double        ld,
                                    const double        t);
    bool isRestrictedZoneBH     (   const Point         testPoint,
                                    const Point         bottom,
                                    const double        radius,
                                    const double        depth);
    bool isInRange              (   const Point         point);

    void inputBoundary 			(	voro::container_poly* 	p_conp,
                                    unsigned* 				p_NodeID,
                                    std::vector<double>*	p_Di);

    void fillInBoundary				();
    std::vector<unsigned> getInBoundary	(	const unsigned		layerNum);
    void applyUnbreakableLattice	(	std::vector<unsigned>*	p_nodeList);
    void applyUnbreakableInBoundary	();
    void applyUnbreakableCoarseMesh	();
    void applyUnbreakableNode ();

    void updateVoroInfo			(	voro::container_poly* 				p_con);

    void getNbListAndArea		(	voro::container_poly* 				p_con,
                                    std::vector<std::vector<int> >* 	p_nbList,
                                    std::vector<std::vector<double> >*	p_area);

    void getNbListAndArea		(	voro::container_poly* 				p_con,
                                    std::vector<std::vector<int> >* 	p_nbList,
                                    std::vector<std::vector<double> >*	p_area,
                                    std::vector<unsigned>* 				p_IDList );

    void findNeighbour	        (	const double						scaleLength,
    								std::vector<NodeIDPair>*			p_conjNodeList);

    void findPreExistingCrack   (   std::vector<NodeIDPair>*            p_conjNodeList);
    void findWeakPlane          ();
    void setUnbreakablePair     ();

    bool testSameBoundary		(	const unsigned					Node1,
                                    const unsigned					Node2);

    bool testMinCoordNum		(	const unsigned					MinCoordNum);

    void updateBoundaryList		(	const unsigned 					NodeID,
                                    const int 						vFace,
                                    const double                    area);

    void updateNodeVoroVol		();

    void putNode				(	voro::container_poly* 			p_conp,
                                    unsigned* 						p_NodeID,
                                    const double 					rx,
                                    const double 					ry,
                                    const double 					rz,
                                    const unsigned 					type);

    void putNode				(	voro::container_poly* 			p_conp,
                                    unsigned*						p_NodeID,
                                    const Point						p,
                                    const unsigned 					type);

    bool insertNode             (   const Point                         inPoint,
                                    unsigned*                           p_NodeID,
                                    voro::container_poly*               p_conp,
                                    Mat3Dvec*                           p_NodeInPartition);

    void putLattice				(	const unsigned 					NodeID,
                                    const unsigned					nbNodeID,
                                    const double					latStiffness,
                                    const std::array<double,2>		latBreakStress,
    								const double					factor,
                                    const double					latArea,
                                    unsigned*						p_LatticeID,
                                    const unsigned					vIndex,
                                    const bool						useLEFM);

    void putBeam                (   const unsigned              NodeID,
                                    const unsigned              nbNodeID,
                                    const double                latStiffness,
                                    const std::array<double,2>  latBreakStress,
                                    const double                factor,
                                    const double                alpha,
                                    const double                latArea,
                                    unsigned*                   p_LatticeID,
                                    const unsigned              vIndex,
                                    const bool                  useLEFM);
    void assignRandomStrength   (   const unsigned              type,
                                    const double                sd,
                                    const double                min);

    void 	updateVertices		(	Vertices*						p_verti);

    void    updateLatInfo       (   Vertices*                       p_verti);

    void	adjCellArea			(	std::vector<unsigned>*			p_NodeList);

    double 	getMeanCellArea		(	const unsigned 					start,
    								const unsigned 					end);

    double 	getCellAreaSD		(	const unsigned 					start,
    								const unsigned 					end,
    								const double					mean);

    void    writeCrackNodeIDPair(   std::vector<NodeIDPair>*        p_conjNodeList);

    void    writeWeakPlaneNodeIDPair ();

    void    setBreakableLatticeAtBoundary ();

    void    setRegFace ();

    unsigned 	updateNodeCoordinationNo();
    void 		calLatParameters();
    void 		reset 			(	voro::container_poly* 			p_conp);
};

class RanPoints {
public:
    RanPoints() {};
    ~RanPoints() {};
private:
    double x,y,z;
};

#endif /* TESELLATION_H_ */
