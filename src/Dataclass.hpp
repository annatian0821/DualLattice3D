
#include "LEM3D.hpp"

#ifndef _DATACLASS_HPP
#define _DATACLASS_HPP

class Restrain
{
public:
    Restrain(const unsigned ID, const std::array<bool,DofPerNode> rArray) {
        NodeID = ID;
        for (unsigned i=0; i<DofPerNode; i++) {
            Dir[i] = rArray[i];
        }
    };
    Restrain(const unsigned ID) {
        NodeID = ID;
        for (unsigned i=0; i<DofPerNode; i++) {
            Dir[i] = false;
        }
    };

    unsigned	NodeID;									//Store restrained nodes (boundaries, preDisp, etc)
    std::array<bool,DofPerNode>		Dir;

    void		remove(const unsigned ID);

    void        addRestrain (const std::array<bool,DofPerNode> rArray) {
        for (unsigned i=0; i<DofPerNode; i++) {
            Dir[i] = Dir[i] || rArray[i];
        }
    }

    bool operator== (const Restrain &rhs) const {
        bool flag = false;
        for (unsigned i=0; i<DofPerNode; i++) {
            flag = flag || Dir[i];
            if (flag) {
                return NodeID == rhs.NodeID;
            }
        }
        return false;
    }

    bool operator== (const unsigned &inNodeID) const {
        return NodeID == inNodeID;
    }

    bool operator< (const Restrain &rhs) const { return NodeID < rhs.NodeID; }

private:

};

struct Constrain
{
public:
    Constrain() {
        NodeID = 0;
    };
    Constrain(const unsigned ID, const std::vector<unsigned> dirVec, const std::vector<double>    dispVec) {
        if (dirVec.size()!=dispVec.size()) {
            tripOut << "[Constrain] Sizes of dirVec and dispVec inputted do not match!" << std::endl;
        }
        if (dirVec.size()>DofPerNode) {
            tripOut << "[Constrain] Sizes of dirVec exceed DofPerNode!" << std::endl;
        }
        if (dirVec.size()>DofPerNode) {
            tripOut << "[Constrain] Sizes of dirVec exceed DofPerNode!" << std::endl;
        }
        NodeID = ID;
        conDir = dirVec;
        conDisp = dispVec;
    };

    Constrain(const unsigned ID, const unsigned dir, const double    disp) {
        NodeID = ID;
        conDir.push_back(dir);
        conDisp.push_back(disp);
    };

    void assign (const unsigned dir, const double disp) {
        std::vector<unsigned>::iterator it = std::find(conDir.begin(),conDir.end(),dir);
        if (it==conDir.end()) {
            tripOut << "[Constrain::assign] Input dir cannot be found, nothing is done!" << std::endl;
        } else {
            conDisp[it - conDir.begin()] = disp;
        }
    }

    void add (const unsigned dir, const double    disp) {
        std::vector<unsigned>::iterator it = std::find(conDir.begin(), conDir.end(), dir);
        if (it==conDir.end()) {
            conDir.push_back(dir);
            conDisp.push_back(disp);
        } else {
            conDisp[it-conDir.begin()] += disp;
        }
    };

    unsigned    NodeID;                                 //Store restrained nodes (boundaries, preDisp, etc)
    std::vector<unsigned>       conDir;
    std::vector<double>       conDisp;
private:

};

class GNode
{
public:
    //Dafault Constructor
    GNode() {};
    //Overloaded Constructor - general info
    GNode	(	const unsigned 		type_input,
                const double 		coord_x,
                const double 		coord_y,
                const double 		coord_z)
    {
        type     = type_input;
        coord[0] = coord_x;
        coord[1] = coord_y;
        coord[2] = coord_z;
        n=0;
        ne=0;
        v=0.0;
        start = 0;
        for (unsigned i=0; i<DofPerNode; i++) {
            d[i] = 0.0;
            d0[i] = 0.0;
            extF[i] = 0.0;
            extdF[i] = 0.0;
        }
    };

    GNode	(	const unsigned 		type_input,
                const Point			_coord)
    {
        type    = type_input;
        coord[0]   = _coord[0];
        coord[1]   = _coord[1];
        coord[2]   = _coord[2];
        n       = 0;
        ne      = 0;
        v       = 0.0;
        start = 0;
        for (unsigned i=0; i<DofPerNode; i++) {
            d[i] = 0.0;
            d0[i] = 0.0;
            extF[i] = 0.0;
            extdF[i] = 0.0;
        }
    };
    // Destructor
    ~GNode() {};
    //Variable
    unsigned 				type;                              				//node type
    Point 					coord;                          				//node coordinates x, y
    unsigned				n;												//Initial coordination number, no of neighborhood
    unsigned				ne;												//Effective coordination number, excluding broken lattice
    unsigned                start;
    double					v;												//Volume of Voronoi cell
    Point					centroid;							            //Coordinates of centroid of voronoi cell
    Tensor                  d;
    Tensor                  d0;
    Tensor                  extdF;
    Tensor                  extF;
    //neighbouring nodes
    std::vector<unsigned> 				nbNodeID;                      			//ID of neighboring nodes
    std::vector<unsigned> 				nbLatticeID;							//ID of neighboring lattices
    inline void 			setNeighbour 	(	const unsigned 		nb_NodeID,
            									const unsigned 		LatticeID);
private:

};

void        GNode::setNeighbour (   const unsigned nb_NodeID,
                                    const unsigned LatticeID)
{
    nbNodeID.push_back(nb_NodeID);
    nbLatticeID.push_back(LatticeID);
    n++;
}

class Lattice
{
public:
    Lattice()	{};
    inline Lattice	(	const unsigned 				NodeID,
                        const unsigned 				nbNodeID,
                        const double 				l,
                        const double                _area,
                        const double 				kn,
                        const std::array<double,2> 	breakDisp);

    inline Lattice	(   const unsigned              NodeID,
                        const unsigned              nbNodeID,
                        const double                l,
                        const double                _area,
                        const double                kn,
                        const double                ks,
                        const std::array<double,2>  breakDisp);

    inline Lattice  (   const unsigned              NodeID,
                        const unsigned              nbNodeID,
                        const double                l,
                        const double                _area,
                        const double                kn,
                        const double                ks,
                        const double                k_theta,
                        const std::array<double,2>  breakDisp,
                        const double                _shearStr);

    ~Lattice()	{};

    std::array<unsigned,2> 				nb;						//NodeID of Nodes on both ends of lattice
    double 								length;					//length of lattice
    double                              factor;
    double                              factorShear;
    Point                               centroid;               //Centroid of facet
    std::array<UniVec,Dim>              axes;                   //Local axes (e0 = axial, e1 major principal, e2 minor principal
    std::array<double,KNum> 			k;						//Stiffness
    std::array<double,2>				et;						//Breaking elongation (0) or contraction (1)
    double                              shearStr;               //Breaking shear stress
    double                              fShear;                 // shear stress at failure
    Point                               intersect;               //intersect point of lattice and facet
    std::array<double,Dim-1>            offset;                 //offset of intersect point to centroid in local coordinates (y,z)
    double								area;					//Area of lattice
    double                              initOpening;            //Initial opening of fracture
    bool                                isBreak;
private:
};

Lattice::Lattice                (   const unsigned              NodeID,
                                    const unsigned              nbNodeID,
                                    const double                l,
                                    const double                _area,
                                    const double                kn,
                                    const std::array<double,2>  breakDisp)
{
    nb[0] = NodeID;
    nb[1] = nbNodeID;
    length = l;
    k[0] = kn;
    et = breakDisp;
    area = _area;
    isBreak = false;
    factor = 1.0;
    factorShear = 1.0;
    shearStr = Huge;
    fShear = Huge;
    initOpening = 0.0;
}

Lattice::Lattice                (   const unsigned              NodeID,
                                    const unsigned              nbNodeID,
                                    const double                l,
                                    const double                _area,
                                    const double                kn,
                                    const double                ks,
                                    const std::array<double,2>  breakDisp)
{
    nb[0] = NodeID;
    nb[1] = nbNodeID;
    length = l;
    k[0] = kn;
    k[1] = ks;
    et = breakDisp;
    area = _area;
    isBreak = false;
    factor = 1.0;
    factorShear = 1.0;
    shearStr = Huge;
    initOpening = 0.0;
    fShear = Huge;
}

Lattice::Lattice                (   const unsigned              NodeID,
                                    const unsigned              nbNodeID,
                                    const double                l,
                                    const double                _area,
                                    const double                kn,
                                    const double                ks,
                                    const double                k_theta,
                                    const std::array<double,2>  breakDisp,
                                    const double                _shearStr)
{
    nb[0] = NodeID;
    nb[1] = nbNodeID;
    length = l;
    k[0] = kn;
    k[1] = ks;
    k[2] = k_theta;
    k[3] = k_theta;
    et = breakDisp;
    area = _area;
    isBreak = false;
    factor = 1.0;
    factorShear = 1.0;
    shearStr = _shearStr;
    initOpening = 0.0;
    fShear = Huge;
}
// Using Born spring as Lattice element. Both normal and shear forces between cells are considered but no moment between them.
class Beam {
public:
	Beam()	{};

    inline Beam     	(	const unsigned              NodeID,
                            const unsigned              nbNodeID,
                            const double                l,
                            const double                _area,
                            const double                kn,
                            const double                ks,
                            const double                i11,
                            const double                i22,
                            const double                i33,
                            const std::array<double,2>  _et);

    inline Beam	        (	const unsigned 				NodeID,
                            const unsigned 				nbNodeID,
                            const double 				l,
                            const double 				_area);

    ~Beam()	{};
    std::array<unsigned,2>				nb;													//NodeID of Nodes on both ends of lattice
    double 								length;												//length of lattice
    std::array<double,6>				k;													//Stiffness: 0:axial 1: shear, 2: major bending, 3:minor bending, 4:torsion
    std::array<double,2>				et;												    //Breaking elongation (0) or contraction (1)
    double								area;
    bool                                isBreak;
private:
};

Beam::Beam              (   const unsigned              NodeID,
                            const unsigned              nbNodeID,
                            const double                l,
                            const double                _area,
                            const double                kn,
                            const double                ks,
                            const double                i11,
                            const double                i22,
                            const double                i33,
                            const std::array<double,2>  _et)
{
    nb[0] = NodeID;
    nb[1] = nbNodeID;
    length = l;
    area = _area;
    k[0] = kn;
    k[1] = k[2] = ks;
    k[3] = i11;
    k[4] = i22;
    k[5] = i33;
    isBreak = false;
}

Beam::Beam              (   const unsigned              NodeID,
                            const unsigned              nbNodeID,
                            const double                l,
                            const double                _area)
{
    nb[0] = NodeID;
    nb[1] = nbNodeID;
    length = l;
    area = _area;
    for (unsigned i=0; i<DofPerNode; i++) {
        k[i] = 0.0;
    }
    isBreak = false;

}


class				NodeID_Dim
{
public:
    NodeID_Dim();
    ~NodeID_Dim();

    void 			resize	(unsigned 	nx, unsigned ny,	unsigned nz, unsigned nxx, unsigned nyy, unsigned nzz, unsigned DumID);
    void			assign	(unsigned 	x, 	unsigned y, 	unsigned z, unsigned NodeID);
    unsigned		get		(unsigned 	x, 	unsigned y, 	unsigned z);

private:
    std::vector<unsigned>	arrNodeID;
    int						size_x, size_y, size_z;
};

class NodeLists
{
public:
    std::array<std::vector<unsigned>,6> 	boundary;								//Store boundary nodes
    std::array<std::vector<double>,6>       boundaryArea;                            //Store area of boundary
    std::array<std::vector<unsigned>,6>	    inBoundary;
    std::vector<unsigned>	free;										//Store real nodes that do not restrained and stable
    std::vector<unsigned>	virt;										//Store virtual nodes
    std::list<Restrain>		restrain;									//Store nodes that are restrainted
    std::vector<Constrain>  constrain;                                  //Store nodes that are constrained
    std::vector<unsigned>	unstable;									//Store nodes that are unstable
    bool removeRestrain (const unsigned ID);
private:

};

extern std::vector			<GNode>				GNodeTab;
extern std::vector			<Lattice>		    LatTab;
extern std::vector			<Beam>		        LatBeamTab;

extern NodeID_Dim								nodeID_Dim;
extern NodeLists								nodeLists;





#endif
