/*
 * MiscTools.hpp
 *
 *  Created on: Oct 24, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#ifndef MISCTOOLS_HPP_
#define MISCTOOLS_HPP_

#include "LEM3D.hpp"
#include <unistd.h>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
typedef std::array<double,2>								Point2D;

std::array<bool,6> getFaceFixity ();
class Clock
{
public:
    Clock()
    {

    };

    Clock(std::string clockName);

    ~Clock()
    {

    };

    void start(std::string clockName);
    void get();
    double get(std::string taskName);
    std::string getDateAndTime();
    double seed ();
private:
    std::string		name;
    struct 			timeval tv;
    double 			wall,cpu,wallStart,cpuStart;
};

template <class T>
class Mat
{
public:
    Mat() {};
    Mat	( const unsigned size_x, const unsigned size_y)
    {
        nx = size_x;
        ny = size_y;
        vec.resize(nx*ny,0);
    };
    ~Mat() {};

    double						nx,ny;

    void	resize (const unsigned size_x, const unsigned size_y)		{
        nx = size_x;
        ny = size_y;
        vec.resize(nx*ny,0);
    };

    void 	set 	(const unsigned i,const unsigned j,const double var) 	{
        vec[i*ny+j] = var;
    };
    void	add		(const unsigned i,const unsigned j,const double var) 	{
        vec[i*ny+j] += var;
    };
    void	multi	(const unsigned i,const unsigned j,const double var) 	{
        vec[i*ny+j] *= var;
    };
    void	div		(const unsigned i,const unsigned j,const double var) 	{
        vec[i*ny+j] /= var;
    };

    void	insert	(const std::vector<T>* p_inVec)						{
        for (unsigned i=0; i<p_inVec->size(); i++)
        {
            vec.insert(vec.begin(),p_inVec->begin(),p_inVec->end());
        }
        nx +=1;
    };

    unsigned	getRowNum	()	{   return nx; };

    double	get	(const unsigned i,const unsigned j) 					{
        return vec[i*ny+j];
    };
    unsigned getSize ()													{
        return nx*ny;
    }

private:
    std::vector<T>			vec;
};

class				Mat3Dvec
{
public:
    Mat3Dvec() {};

    ~Mat3Dvec() {};

    void 				resize		(	const unsigned 			nx,
                                        const unsigned 			ny,
                                        const unsigned 			nz);

    inline void			push_back	(	const unsigned 			x,
                                        const unsigned 			y,
                                        const unsigned 			z,
                                        const unsigned 			NodeID)
    {
        vecNodeID[x*(size_y*size_z)+ y*size_z + z].push_back(NodeID);
    }

    inline void         push_back   (   const std::array<unsigned,3>          p,
                                        const unsigned                        NodeID)
    {
        vecNodeID[p[0]*(size_y*size_z)+ p[1]*size_z + p[2]].push_back(NodeID);
    }

    inline int 		get			        (	unsigned 			    x,
                                            unsigned 			    y,
                                            unsigned 			    z,
                                            std::vector<unsigned>*	p_IDList)
    {
        if (x>=size_x) {
            tripOut << "[Mat3Dvec.get] - x = " << x << " >= " << size_x << " return -1 \n";
            return -1;
        }
        if (y>=size_y) {
            tripOut << "[Mat3Dvec.get] - y = " << y << " >= " << size_y << " return -1 \n";
            return -1;
        }
        if (z>=size_z) {
            tripOut << "[Mat3Dvec.get] - z = " << z << " >= " << size_z << " return -1 \n";
            return -1;
        }
        if (vecNodeID[x*(size_y*size_z)+ y*size_z + z].empty())
        {
            return 0;
        }
        else
        {
            //insert vector at the end
            p_IDList -> clear();
            *p_IDList = vecNodeID[x*(size_y*size_z)+ y*size_z + z];
            return 1;
        }
    }

    inline bool			append		(	const unsigned 			x,
                                        const unsigned 			y,
                                        const unsigned 			z,
                                        std::vector<unsigned>*	p_IDList)
    {
        if (vecNodeID[x*(size_y*size_z)+ y*size_z + z].empty()) {
            return false;
        } else {
            //insert vector at the end
            p_IDList -> insert (p_IDList -> end(), vecNodeID[x*(size_y*size_z)+ y*size_z + z].begin(), vecNodeID[x*(size_y*size_z)+ y*size_z + z].end());
            return true;
        }
    };
    std::array<unsigned,3>      getDims ()
    {
        std::array<unsigned,3> dims;
        dims[0] = size_x;
        dims[1] = size_y;
        dims[2] = size_z;
        return dims;
    }

private:
    std::vector<std::vector<unsigned> >			vecNodeID;
    unsigned									size_x, size_y, size_z;
};

class               Mat2Dvec
{
public:
    Mat2Dvec() {
        size_x = size_y = 0;
    };
    ~Mat2Dvec() {};

    void                resize      (   const unsigned          nx,
                                        const unsigned          ny);

    inline void         push_back   (   const unsigned          x,
                                        const unsigned          y,
                                        const unsigned          NodeID)
    {
        vecNodeID[x*size_y+ y].push_back(NodeID);
    }

    inline int          get         (   const unsigned          x,
                                        const unsigned          y,
                                        std::vector<unsigned>*  p_IDList)
    {
        if (x>=size_x) {
            tripOut << "[Mat2Dvec.get] - x = " << x << " >= " << size_x << " return -1 \n";
            return -1;
        }
        if (y>=size_y) {
            tripOut << "[Mat2Dvec.get] - y = " << y << " >= " << size_y << " return -1 \n";
            return -1;
        }
        if (vecNodeID[x*size_y+ y].empty())
        {
            return 0;
        }
        else
        {
            //insert vector at the end
            p_IDList -> clear();
            *p_IDList = vecNodeID[x*size_y+ y];
            return true;
        }
    }

    std::array<unsigned,2>      getDims ()
    {
        std::array<unsigned,2> dims;
        dims[0] = size_x;
        dims[1] = size_y;
        return dims;
    }

    inline bool         append      (   const unsigned          x,
                                        const unsigned          y,
                                        std::vector<unsigned>*  p_IDList)
    {
        if (vecNodeID[x*size_y+ y].empty()) {
            return false;
        } else {
            //insert vector at the end
            p_IDList -> insert (p_IDList -> end(), vecNodeID[x*size_y+ y].begin(), vecNodeID[x*size_y+ y].end());
            return true;
        }
    };

private:
    std::vector<std::vector<unsigned> >         vecNodeID;
    unsigned                                    size_x, size_y;
};

class iDouble
{
public:
    iDouble (const double input, const unsigned i)
    {
        data = input;
        index = i;
    }

    ~iDouble () {};

    //For sorting from greatest to smallest
    inline bool operator< (const iDouble& rhs) const
    {
        return data > rhs.data;
    }

//	friend bool operator== (const iDouble &n1, const iDouble &n2);

    double		data;
    unsigned	index;
private:

};

class iTriple
{
public:
    iTriple (const double input1, const double input2, const unsigned i)
    {
        data1 = input1;
        data2 = input2;
        index = i;
    }

    ~iTriple () {};

    //For sorting from greatest to smallest
    inline bool operator< (const iTriple& rhs) const
    {
        return data1 > rhs.data1;
    }

//  friend bool operator== (const iDouble &n1, const iDouble &n2);

    double      data1;
    unsigned    data2;
    unsigned    index;
private:

};

class Geometry
{
public:
    Geometry() {};
    ~Geometry() {};
//    template <class T1,class T2>
    double dist(Point N1, Point N2) {
        return sqrt((N1[0]-N2[0])*(N1[0]-N2[0])+
                	(N1[1]-N2[1])*(N1[1]-N2[1])+
                	(N1[2]-N2[2])*(N1[2]-N2[2]));
    };
    double dist(Point N1, PointF N2) {
        return sqrt((N1[0]-N2[0])*(N1[0]-N2[0])+
                    (N1[1]-N2[1])*(N1[1]-N2[1])+
                    (N1[2]-N2[2])*(N1[2]-N2[2]));
    };

    double mag (Point 		N) {
    	return sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
    }

    double mag2 (std::array<double,DofPerNode>       N) {
        return sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
    }

    Vec     rotateVec   (   const double        x1,
                            const double        y1,
                            const UniVec        n);

    Vec	crossProd		(	const Vec			v1,
    						const Vec			v2)
    {
    	Vec out;
    	out[0] = v1[1]*v2[2]-v1[2]*v2[1];
    	out[1] = v1[2]*v2[0]-v1[0]*v2[2];
    	out[2] = v1[0]*v2[1]-v1[1]*v2[0];
    	return out;
    }

    inline double linePointDist	(	Point	point,
    								Line	line);

    double dot				(	const Point				d0,
    							const Point				d1)
    {
    	double var = 0.0;
    	for (unsigned i=0; i<Dim; i++) {
    		var += d0[i]*d1[i];
    	}
    	return var;
    };
    double dot              (   const std::array<float,Dim>             d0,
                                const Point                             d1)
   {
        double var = 0.0;
        for (unsigned i=0; i<Dim; i++) {
            var += d0[i]*d1[i];
        }
        return var;
   };
   std::array<double,2> distOffsetProj        (   Point N, Point N0, Vec     uniVec)
   {
       std::array<double,2> distOffset;
       Point intersect,intersect_local;
       Point d;
       for (unsigned i=0; i<Dim; i++) {
           d[i] = N0[i]-N[i];
       }
       double k = dot(d,uniVec);
       for (unsigned i=0; i<Dim; i++) {
           intersect[i] =N[i] + k*uniVec[i];
           intersect_local[i] = N0[i] - intersect[i];
       }
       distOffset[0] = dist(intersect,N0);
       distOffset[1] = k;
       return distOffset;
    }
    Point subtract			(	const Point				d0,
								const Point				d1)
    {
    	Point result;
    	for (unsigned i=0; i<Dim; i++) {
    		result[i] = d0[i] - d1[i];
    	}
    	return result;
    }
    Point			getCentroid						(	const double                    latArea,
                                                        const Surface					&v);			//vertices coordinates

    double 			getArea3D_Polygon				(	Surface 						V,
    													Point 							N);

    std::array<double,3>     getSecMomentOfArea2D   (   std::vector<std::array<double,2> >      &v2D);

    inline double 	getMaxCentroidVertiDist			(	const Point 					centroid,
    													const Surface					&v);
    inline bool		isCoplanar						(	const std::vector<Point>*		p_ptList,
    													const double					tolerance);

    inline bool		isCoplanar						(	const Point						p1,
    													const Point						p2,
    													const Point						p3,
    													const Point						p4,
    													const double					tolerance);

    BasisVecs		getBasisVec					    (	const Point						nbCoord0,
														const Point						nbCoord1,
    													const Surface					&v,				//vertices coordinates
    													const Point						centroid);

    Point 			projectPoint					(	const Point						x,
    													const Point						p,
    													const Point						n);

private:
    bool			getCentroid2D					(	const Surface					&v,
    													Point							&centroid);

};

struct		Plane
{
    Plane () {};
    ~Plane () {};
    Plane		(	Point		_pt,
          			Vec			_norm) {
    	pt = _pt;
    	//Normalize the normal vector
        double	normMag = std::sqrt(_norm[0]*_norm[0]+_norm[1]*_norm[1]+_norm[2]*_norm[2]);
        if (normMag>Tiny) {
			for (unsigned i=0; i<Dim; i++) {
				norm[i] = _norm[i]/normMag;
			}
        } else {
        	tripOut << "Plane constructor - Norm input is too small!" << std::endl;
        }
    };
    Point		pt;
    Vec			norm;
    double		torNorm;

    inline double 		getPtPlaneDist			(	const Point			point);

    inline Point		getProjPoint			(	const Point			point);

    inline Point		getProjPoint			(	const Point			point,
    												double*				p_dist);

    inline Point		getintersectPoint		(	const Point			p1,
													const Point			p2);

    inline bool			isNormClose				(	const Point			p1,
													const Point			p2,
													const double		tol);
};

double 		Plane::getPtPlaneDist			(	const Point			point) {
    Vec	ptSubtract;
    for (unsigned i=0; i<Dim; i++) {
    	ptSubtract[i] = point[i] - pt[i];
    }
    return norm[0]*ptSubtract[0]+norm[1]*ptSubtract[1]+norm[2]*ptSubtract[2];
}

Point		Plane::getProjPoint			(	const Point			point,
											double*				p_dist)
{
	*p_dist = 0.0;
    for (unsigned i=0; i<Dim; i++) {
    	*p_dist += (point[i]-pt[i])*norm[i];
    }
    Point	result;
   	for (unsigned i=0; i<Dim; i++) {
    	result[i] = point[i] - *p_dist*norm[i];
    }
    return result;
}

Point		Plane::getProjPoint			(	const Point			point) {
    double 	dist = 0.0;
    for (unsigned i=0; i<Dim; i++) {
    	dist += (point[i]-pt[i])*norm[i];
    }
    Point	result;
   	for (unsigned i=0; i<Dim; i++) {
    	result[i] = point[i] - dist*norm[i];
    }
    return result;
}

Point		Plane::getintersectPoint	(	const Point			p1,
											const Point			p2)
{
	Point result;
	Geometry geo;
	Point ca = geo.subtract(pt,p1);
	Point ba = geo.subtract(p2,p1);
	double d=geo.dot(ca,norm)/geo.dot(ba,norm);
	for (unsigned i=0; i<Dim; i++) {
		result[i] = p1[i]+d*(p2[i]-p1[i]);
	}
	return result;
}

bool		Plane::isNormClose			(	const Point			p1,
											const Point			p2,
											const double		tol)
{
	Geometry geo;
	Vec		pNorm = geo.subtract(p1,p2);
	double	magNorm = geo.mag(pNorm);
	for (unsigned i=0; i<Dim; i++) {
		pNorm[i] /= magNorm;
	}
	double npDot = geo.dot(norm,pNorm);
	return (std::fabs(npDot)>tol);
}

double Geometry::linePointDist		(	const Point		point,
										const Line		line)
{
	Point	d0,d1;
	for (unsigned i=0; i<Dim; i++) {
		d0[i] = line[0][i]-point[i];
		d1[i] = line[1][i]-line[0][i];
	}
	double d0d0 = dot(d0,d0);
	double d1d1 = dot(d1,d1);
	double d0d1 = dot(d0,d1);
	return (d0d0*d1d1-d0d1*d0d1)/d1d1;
}

double Geometry::getMaxCentroidVertiDist		(	const Point 	centroid,
													const Surface	&v)
{
	double maxDist = 0.0;
	unsigned size = v.size();
	Geometry geo;
	double dist;
	for (unsigned i=0; i<size; i++) {
		dist = geo.dist(centroid,v[i]);
		if (dist>maxDist) {
			maxDist = dist;
		}
	}
	return maxDist;
}

inline bool		Geometry::isCoplanar						(	const std::vector<Point>*		p_ptList,
																const double					tolerance)
{
	if (p_ptList->size()<4) {
		return false;
	} else {
		for (unsigned i=3; i<p_ptList->size(); i++) {
			if (!(isCoplanar((*p_ptList)[0],(*p_ptList)[1],(*p_ptList)[2],(*p_ptList)[i],tolerance))) {
				return false;
			}
		}
		return true;
	}
}
inline bool		Geometry::isCoplanar						(	const Point						p1,
    															const Point						p2,
    															const Point						p3,
    															const Point						p4,
    															const double					tolerance)
{
	Eigen::Vector3d		v1,v2,v3,v4;
	v1 << p1[0] , p1[1] , p1[2];
	v2 << p2[0] , p2[1] , p2[2];
	v3 << p3[0] , p3[1] , p3[2];
	v4 << p4[0] , p4[1] , p4[2];
	return (std::fabs((v3-v1).dot((v2-v1).cross(v4-v3)))<tolerance);
}

template <typename T1, typename T2>
struct CustomPair {
	CustomPair () {};
	~CustomPair () {};
	CustomPair (	const T1 &data1, const T2 &data2) {
		i1 = data1;
		i2 = data2;
	}
	T1 i1;
	T2 i2;
	bool operator< 	(	const CustomPair& rhs) const { return i1 < rhs.i1; }
	bool operator== (	const CustomPair& rhs) const { return i1 == rhs.i1; }
	bool operator== (	const T1 rhs) const { return i1 == rhs; }
	bool operator< (   const T1 rhs) const { return i1 < rhs; }
};

template <class T1, class T2>
void sortPair	(	std::vector<T1>*	p_vec1,
                    std::vector<T2>*	p_vec2)
{
    unsigned i;
    std::vector<std::pair<T1,T2> >	vPair;
    std::pair<T1,T2>	newPair;
    if ( p_vec1->size()!=p_vec2->size())
    {
        dualOut << "Two input vectors are of different size, nothing is done!" << std::endl;
        return;
    }

    for (i=0; i<p_vec1->size(); i++)
    {
        newPair = std::make_pair((*p_vec1)[i],(*p_vec2)[i]);
        vPair.push_back(newPair);
    }

    std::sort(vPair.begin(),vPair.end());

    for (i=0; i<p_vec1->size(); i++)
    {
        (*p_vec1)[i] = vPair[i].first;
        (*p_vec2)[i] = vPair[i].second;
    }
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <class T>
void uniqueVector (std::vector<T>*	p_vec)
{
    std::sort(p_vec->begin(), p_vec->end());
    p_vec->erase(std::unique(p_vec->begin(),p_vec->end()), p_vec->end());
}

template <class T>
void uniqueList (std::list<T>*	p_list)
{
	p_list->sort();
	p_list->unique();
}

template<typename T>
inline void remove(std::vector<T> & v, const T & item)
{
    v.erase(std::remove(v.begin(), v.end(), item), v.end());
}

template<class T1,class T2>
inline bool comparePairFirst (const std::pair<T1,T2> &p1,const std::pair<T1,T2> &p2)
{
	return (p1.first==p2.first);
}

template<class T1,class T2>
inline bool comparePairSecond (const std::pair<T1,T2> &p1,const std::pair<T1,T2> &p2)
{
	return (p1.second==p2.second);
}

template<typename T>
inline void setDiffVector (	std::vector<T> &in1, std::vector<T> &in2, std::vector<T> &out)
{
	std::sort(in1.begin(),in1.end());
	std::sort(in2.begin(),in2.end());
	std::set_symmetric_difference (in1.begin(),in1.end(),in2.begin(),in2.end(),std::back_inserter(out));
}

template<typename T>
inline void setDiffVector ( std::set<T> &in1, std::vector<T> &in2, std::vector<T> &out)
{
    std::sort(in1.begin(),in1.end());
    std::sort(in2.begin(),in2.end());
    std::set_symmetric_difference (in1.begin(),in1.end(),in2.begin(),in2.end(),std::back_inserter(out));
}

template<typename T>
inline std::vector<T> setDiffVector (	std::vector<T> &in1, std::vector<T> &in2)
{
	std::vector<T> out;
	std::sort(in1.begin(),in1.end());
	std::sort(in2.begin(),in2.end());
	std::set_symmetric_difference (in1.begin(),in1.end(),in2.begin(),in2.end(),std::back_inserter(out));
	return out;
}

template<typename T>
void eraseValue(std::vector<T> &vec, T value)
{
    auto pr = std::equal_range(std::begin(vec), std::end(vec), value);
    vec.erase(pr.first, pr.second);
}

template<typename T>
void eraseValue(std::list<T> &vec, T value)
{
    auto pr = std::equal_range(std::begin(vec), std::end(vec), value);
    vec.erase(pr.first, pr.second);
}

template<typename T>
inline void printData (const T  data) {
	std::cout << ' ' << data;
}

template<typename T>
inline void dualOutData (const T  data) {
	dualOut << ' ' << data;
}

template<typename T>
struct WriteData {
	WriteData () {};
	~WriteData () {};
	WriteData (std::ofstream*	_p_file) {
		p_file = _p_file;
	}
	std::ofstream*	p_file;
	void write (	const T  data) {
		*p_file << ' ' << data;
	}
};

#endif /* MISCTOOLS_HPP_ */
