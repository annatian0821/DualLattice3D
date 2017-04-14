/*
 * MiscTools.cpp
 *
 *  Created on: Oct 24, 2013
 *      Author: John K. W. Wong, kwjw2@cam.ac.uk
 */

#include "MiscTools.hpp"

using namespace std;

Point	Geometry::getCentroid						(	const double        latArea,
                                                        const Surface       &v)			//vertices coordinates
{
	Point	centroid;

	unsigned vNum = v.size();
	double xy,xz,yz;
	double Axy = 0.0;
	double Axz = 0.0;
	double Ayz = 0.0;

	const double rTol = 1e-4;
	const double tolerance = rTol*latArea;

	std::array<double,2> centroid_xy;
	std::array<double,2> centroid_xz;
	std::array<double,2> centroid_yz;

	centroid_xy[0] = centroid_xy[1] = 0.0;
	centroid_xz[0] = centroid_xz[1] = 0.0;
	centroid_yz[0] = centroid_yz[1] = 0.0;

	centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;
	// For all vertices except last
	for (unsigned i=0; i<vNum-1; i++) {
		xy = v[i][0]*v[i+1][1]-v[i+1][0]*v[i][1];
		xz = v[i][0]*v[i+1][2]-v[i+1][0]*v[i][2];
		yz = v[i][1]*v[i+1][2]-v[i+1][1]*v[i][2];

		Axy +=xy;
		Axz +=xz;
		Ayz +=yz;

		centroid_xy[0] += (v[i][0]+v[i+1][0])*xy;
		centroid_xy[1] += (v[i][1]+v[i+1][1])*xy;
		centroid_xz[0] += (v[i][0]+v[i+1][0])*xz;
		centroid_xz[1] += (v[i][2]+v[i+1][2])*xz;
        centroid_yz[0] += (v[i][1]+v[i+1][1])*yz;
        centroid_yz[1] += (v[i][2]+v[i+1][2])*yz;

	}
	// Do last vertex
	xy = v[vNum-1][0]*v[0][1]-v[0][0]*v[vNum-1][1];
	xz = v[vNum-1][0]*v[0][2]-v[0][0]*v[vNum-1][2];
	yz = v[vNum-1][1]*v[0][2]-v[0][1]*v[vNum-1][2];

	Axy +=xy;
	Axz +=xz;
	Ayz +=yz;

	centroid_xy[0] += (v[vNum-1][0]+v[0][0])*xy;
	centroid_xy[1] += (v[vNum-1][1]+v[0][1])*xy;
	centroid_xz[0] += (v[vNum-1][0]+v[0][0])*xz;
	centroid_xz[1] += (v[vNum-1][2]+v[0][2])*xz;
	centroid_yz[0] += (v[vNum-1][1]+v[0][1])*yz;
	centroid_yz[1] += (v[vNum-1][2]+v[0][2])*yz;
	//To handle geometry that are closed to a line after projection on xy or xz plane

	if (std::fabs(Axy)<tolerance) {
	    if (std::fabs(Axz)<tolerance){
	        centroid[1] = centroid_yz[0]/(3.0*Ayz);
	        centroid[2] = centroid_yz[1]/(3.0*Ayz);
	        unsigned d = 0;
	        for (unsigned i=0; i<vNum; i++) {
	            centroid[d] += v[i][d];
	        }
	        centroid[d] /= vNum;

	    } else if (std::fabs(Ayz)<tolerance){
	        centroid[0] = centroid_xz[0]/(3.0*Axz);
	        centroid[2] = centroid_xz[1]/(3.0*Axz);
	        unsigned d = 1;
	        for (unsigned i=0; i<vNum; i++) {
	            centroid[d] += v[i][d];
	        }
	        centroid[d] /= vNum;
	    } else {
	        centroid[0] = centroid_xz[0]/(3.0*Axz);
	        centroid[1] = centroid_yz[0]/(3.0*Ayz);
	        centroid[2] = centroid_yz[1]/(3.0*Ayz);
	    }
	} else if (std::fabs(Axz)<tolerance) {
	    if (std::fabs(Ayz)<tolerance) {
	       centroid[0] = centroid_xy[0]/(3.0*Axy);
	       centroid[1] = centroid_xy[1]/(3.0*Axy);
	       unsigned d = 2;
	       for (unsigned i=0; i<vNum; i++) {
	           centroid[d] += v[i][d];
	       }
	       centroid[d] /= vNum;
	    } else {
	        centroid[0] = centroid_xy[0]/(3.0*Axy);
	        centroid[1] = centroid_xy[1]/(3.0*Axy);
	        centroid[2] = centroid_yz[1]/(3.0*Ayz);
	    }
	} else if (std::fabs(Ayz)<tolerance) {
	    centroid[0] = centroid_xy[0]/(3.0*Axy);
	    centroid[1] = centroid_xy[1]/(3.0*Axy);
	    centroid[2] = centroid_xz[1]/(3.0*Axz);
	} else {
	    centroid[0] = (centroid_xy[0]/Axy+centroid_xz[0]/Axz)/6.0;
	    centroid[1] = (centroid_xy[1]/Axy+centroid_yz[0]/Ayz)/6.0;
	    centroid[2] = (centroid_xz[1]/Axz+centroid_yz[1]/Ayz)/6.0;
	}
	return centroid;
}

//===================================================================


// area3D_Polygon(): compute the area of a 3D planar polygon
//  Input:  int n = the number of vertices in the polygon
//          Point* V = an array of n+1 points in a 2D plane with V[n]=V[0]
//          Point N = a normal vector of the polygon's plane
//  Return: the (float) area of the polygon

double Geometry::getArea3D_Polygon		(	std::vector<Point> 			V,
											Point 						N)
{
    double area = 0;
    double an, ax, ay, az; // abs value of normal and its coords
    int  coord;           // coord to ignore: 1=x, 2=y, 3=z
    int  i, j, k;         // loop indices
    unsigned n = V.size();
    if (n < 3) {
    	tripOut << "[Geometry::getArea3D_Polygon] polygon inputted has less than 3 vertices" << std::endl;
    	return 0;  // a degenerate polygon
    }
    V.push_back(V[0]);

    // select largest abs coordinate to ignore for projection
    ax = (N[0]>0 ? N[0] : -N[0]);    // abs x-coord
    ay = (N[1]>0 ? N[1] : -N[1]);    // abs y-coord
    az = (N[2]>0 ? N[2] : -N[2]);    // abs z-coord

    if (std::fabs(N[0])>Tiny) {
    	coord = 3;
    }
    coord = 3;                    // ignore z-coord
    if (ax > ay) {
        if (ax > az) coord = 1;   // ignore x-coord
    }
    else if (ay > az) coord = 2;  // ignore y-coord

    // compute area of the 2D projection
    switch (coord) {
      case 1:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i][1] * (V[j][2] - V[k][2]));
        break;
      case 2:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i][2] * (V[j][0] - V[k][0]));
        break;
      case 3:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i][0] * (V[j][1] - V[k][1]));
        break;
    }
    switch (coord) {    // wrap-around term
      case 1:
        area += (V[n][1] * (V[1][2] - V[n-1][2]));
        break;
      case 2:
        area += (V[n][2] * (V[1][0] - V[n-1][0]));
        break;
      case 3:
        area += (V[n][0] * (V[1][1] - V[n-1][1]));
        break;
    }

    // scale to get area before projection
    an = std::sqrt( ax*ax + ay*ay + az*az); // length of normal vector
    switch (coord) {
      case 1:
        area *= (an / (2 * N[0]));
        break;
      case 2:
        area *= (an / (2 * N[1]));
        break;
      case 3:
        area *= (an / (2 * N[2]));
        break;
    }
    return area;
}

Point Geometry::projectPoint							(	const Point			x,
															const Point			p,
															const Point			n)
{
	double t = (n[0]*(p[0]-x[0]) + n[1]*(p[1]-x[1]) + n[2]*(p[2]-x[2])) / (n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	Point out;
	out[0] = x[0] + t*n[0];
	out[1] = x[1] + t*n[1];
	out[2] = x[2] + t*n[2];
	return out;
}

BasisVecs	Geometry::getBasisVec						(	const Point			nbCoord0,
															const Point			nbCoord1,
															const Surface		&v,		//vertices coordinates
															const Point			centroid)
{
	UniVec		e1,e2,e3;
	double		magnitude = 0.0;
	Vec			diff;
	for (unsigned j=0; j<Dim; j++) {
		diff[j]= nbCoord0[j]-nbCoord1[j];
	}
	magnitude = mag(diff);
	for (unsigned j=0; j<Dim; j++) {
		e3[j] = diff[j]/ magnitude;
	}
	unsigned vNum = v.size();
	std::vector<Point>		nVerti;
	for (unsigned i=0; i<vNum; i++) {
		Point pt;
		for (unsigned j=0; j<Dim; j++) {
			pt[j] = v[i][j] - centroid[j];
		}
		nVerti.push_back(pt);
	}
	//Calculate 1st axis representing one of the transverse directions
	magnitude = mag(nVerti[0]);
	for (unsigned j=0; j<Dim; j++) {
		e1[j] = nVerti[0][j]/magnitude;
	}
	e2 = crossProd(e1,e3);
	cout << e1[0] << ' ' << e1[1] << ' ' << e1[2] << '\n';
	cout << e2[0] << ' ' << e2[1] << ' ' << e2[2] << '\n';
	cout << e3[0] << ' ' << e3[1] << ' ' << e3[2] << '\n';
	//Calcaulate 2D coordinates with reference to surface defined by e1 and e2
	std::vector<Point2D>	verti2D;
	Point2D	pt2D;
	for (unsigned i=0; i<vNum; i++) {
		pt2D[0] = dot(nVerti[i],e1);
		pt2D[1] = dot(nVerti[i],e2);
		verti2D.push_back(pt2D);
	}
	std::array<UniVec,3>		priDir;
	double Ixx = 0.0, Iyy = 0.0, Ixy = 0.0;
	double	var;
	for (unsigned i=0; i<vNum-1; i++) {
		var = (verti2D[i][0]*verti2D[i+1][1]-verti2D[i+1][0]*verti2D[i][1]);
		Ixx += (verti2D[i][1]*verti2D[i][1]+verti2D[i][1]*verti2D[i+1][1]+verti2D[i+1][1]*verti2D[i+1][1])*var;
		Iyy += (verti2D[i][0]*verti2D[i][0]+verti2D[i][0]*verti2D[i+1][0]+verti2D[i+1][0]*verti2D[i+1][0])*var;
		Ixy += (verti2D[i][0]*verti2D[i+1][1]+2.0*verti2D[i][0]*verti2D[i][1]+
				2.0*verti2D[i+1][0]*verti2D[i+1][1]+verti2D[i+1][0]*verti2D[i][1])*var;
	}
	var  = (verti2D[vNum-1][0]*verti2D[0][1]-verti2D[0][0]*verti2D[vNum-1][1]);
	Ixx += (verti2D[vNum-1][1]*verti2D[vNum-1][1]+verti2D[vNum-1][1]*verti2D[0][1]+verti2D[0][1]*verti2D[0][1])*var;
	Iyy += (verti2D[vNum-1][0]*verti2D[vNum-1][0]+verti2D[vNum-1][0]*verti2D[0][0]+verti2D[0][0]*verti2D[0][0])*var;
	Ixy += (verti2D[vNum-1][0]*verti2D[0][1]+2.0*verti2D[vNum-1][0]*verti2D[vNum-1][1]+
			2.0*verti2D[0][0]*verti2D[0][1]+verti2D[0][0]*verti2D[vNum-1][1])*var;
	double alpha = std::atan(2.0*Ixy/(Iyy-Ixx+Tiny))/2.0;
	for (unsigned j=0; j<Dim; j++) {
		priDir[0][j] = e1[j]*std::cos(alpha)+e2[j]*std::sin(alpha);
		priDir[1][j] = -e1[j]*std::sin(alpha)+e2[j]*std::cos(alpha);
	}
	priDir[2] = e3;
	return priDir;
}

// Rotate a vector lying on xy-plane (x1,y1,0) normal to z-axis (0,0,1) to anyplane plane with normal vector n.
// The refernece point is origin (0,0,0)
Vec     Geometry::rotateVec                         (   const double        x1,
                                                        const double        y1,
                                                        const UniVec        n)
{
    Vec     outVec;
    if ((n[0]==0.0)&&(n[1]==0.0)) {
        outVec[0] = x1;
        outVec[1] = y1;
        outVec[2] = 0.0;
        return outVec;
    }
    double  c = (n[2]-1)/(n[0]*n[0]+n[1]*n[1]);
    outVec[0] = (1+n[0]*n[0]*c)*x1+ n[0]*n[1]*c*y1;
    outVec[1] = n[0]*n[1]*c*x1 + (1+n[1]*n[1]*c)*y1;
    outVec[2] = -n[0]*x1 - n[1]*y1;
    return outVec;
}

std::array<double,3>     Geometry::getSecMomentOfArea2D  (   std::vector<std::array<double,2> >      &v2D)
{
    std::array<double,3>     I;             //I[0] = Ixx, I[1] = Iyy, I[2] = Ixy
    I[0] = 0.0;     I[1] = 0.0;     I[2] = 0.0;
    double C;
    for (unsigned i=0; i<v2D.size()-1; i++) {
        C = v2D[i][0]*v2D[i+1][1]-v2D[i+1][0]*v2D[i][1];
        I[0] += (v2D[i][1]*v2D[i][1]+v2D[i][1]*v2D[i+1][1]+v2D[i+1][1]*v2D[i+1][1])*C;
        I[1] += (v2D[i][0]*v2D[i][0]+v2D[i][0]*v2D[i+1][0]+v2D[i+1][0]*v2D[i+1][0])*C;
        I[2] += (v2D[i][0]*v2D[i+1][1]+2.0*v2D[i][0]*v2D[i][1]+
                2.0*v2D[i+1][0]*v2D[i+1][1]+v2D[i+1][0]*v2D[i][1])*C;
    }
    I[0] /= 12.0;   I[1] /= 12.0;   I[2] /= 24.0;
    I[0] = std::fabs(I[0]);
    I[1] = std::fabs(I[1]);
    I[2] = std::fabs(I[2]);
    return I;
}

bool	Geometry::getCentroid2D						(	const Surface		&v,
														Point				&centroid)			//vertices coordinates
{
	centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;
	unsigned vNum = v.size();
	double xy;
	double Axy = 0.0;
	double diff;
	std::array<double,Dim>	mean;
	std::array<bool,Dim>	isCoplanarXYZ;
	isCoplanarXYZ[0] = false;	isCoplanarXYZ[1] = false;	isCoplanarXYZ[2] = false;
	for (unsigned d=0; d<Dim; d++) {
		mean[d] = 0.0;
		for (unsigned i=0; i<vNum; i++) {
			mean[d] += v[i][d];
		}
		mean[d] /= vNum;
		diff = 0.0;
		for (unsigned i=0; i<vNum; i++) {
			diff += abs(mean[d] - v[i][d]);
		}
		diff /= vNum;
		if (std::fabs(diff) < 1.0e-5) {
			isCoplanarXYZ[d] = true;
		}
	}
	unsigned coplanarNum = 0;
	for (unsigned d=0; d<Dim; d++) {
		coplanarNum += isCoplanarXYZ[d];
	}
	if (coplanarNum==0) {
		return false;
	}
	if (coplanarNum>=2) {
		for (unsigned d=0; d<Dim; d++) {
			centroid[d] = mean[d];
		}
	} else if (coplanarNum==1) {
		tripOut << "coplanarNum==1" << std::endl;
		std::array<unsigned,2> xyz;
		for (unsigned d=0; d<Dim; d++) {
			if (isCoplanarXYZ[d]) {
				centroid[d] = mean[d];
				xyz[0] = (d+1)%3;
				xyz[1] = (d+2)%3;
				// For all vertices except last
				unsigned i;
				for (i=0; i<vNum-1; i++) {
					xy = v[i][xyz[0]]*v[i+1][xyz[1]]-v[i+1][xyz[0]]*v[i][xyz[1]];
					Axy +=xy;
					centroid[xyz[0]] += (v[i][xyz[0]]+v[i+1][xyz[0]])*xy;
					centroid[xyz[1]] += (v[i][xyz[1]]+v[i+1][xyz[1]])*xy;
				}
				// Do last vertex
				xy = v[i][xyz[0]]*v[0][xyz[1]]-v[0][xyz[0]]*v[i][xyz[1]];
				Axy +=xy;
				centroid[xyz[0]] += (v[i][xyz[0]]+v[0][xyz[0]])*xy;
				centroid[xyz[1]] += (v[i][xyz[1]]+v[0][xyz[1]])*xy;
				centroid[xyz[0]] /= 3.0*Axy+Tiny; centroid[xyz[1]] /= 3.0*Axy+Tiny;
			}
		}
		tripOut << centroid[0] << " " << centroid[1] << " " << centroid[2] << std::endl;
	}
	return true;
}


Clock::Clock(string clockName)
{
    name = clockName ;
    gettimeofday(&tv,NULL);
    wallStart=tv.tv_sec+tv.tv_usec/1e6;
    cpuStart = clock()/(double) CLOCKS_PER_SEC;
    if (name.size()!=0)
    {
        dualOut << name << " clock is set ! " <<endl;
    }
    else
    {
    	dualOut << "A clock is set ! " <<endl;
    }
}

void Clock::start(string clockName)
{
    name = clockName ;
    gettimeofday(&tv,NULL);
    wallStart=tv.tv_sec+tv.tv_usec/1e6;
    cpuStart = clock()/(double) CLOCKS_PER_SEC;
    if (name.size()!=0)
    {
        dualOut << name << " clock is set ! " <<endl;
    }
}

void Clock::get()
{
    unsigned s,min,hr;

    gettimeofday(&tv,NULL);
    wall = tv.tv_sec+tv.tv_usec/1e6 -wallStart;
    cpu = clock()/(double) CLOCKS_PER_SEC - cpuStart;

    hr = (unsigned) wall/3600;
    min = (unsigned) wall%3600/60;
    s = (unsigned) wall%60;
    if (name.size()!=0)
    {
        dualOut << "Clock " <<  name << ": ";
    }
    dualOut << "Time = " << wall << " s, or " << hr << " hr " << min << " min " << s << " s , CPU Time = " << cpu <<"s" <<endl;
}

double Clock::get(string taskName)
{
    unsigned s,min,hr;

    gettimeofday(&tv,NULL);
    wall = tv.tv_sec+tv.tv_usec/1e6 -wallStart;
    cpu = clock()/(double) CLOCKS_PER_SEC - cpuStart;

    hr = (unsigned) wall/3600;
    min = (unsigned) wall%3600/60;
    s = (unsigned) wall%60;

    if (name.size()!=0)
    {
        dualOut << "Clock " <<  name << " - " << taskName << ": ";
    }
    dualOut << "Time = " << wall << " s, or " << hr << " hr " << min << " min " << s << " s , CPU Time = " << cpu <<"s" <<endl;

    return wall;
}

double Clock::seed ()
{
    gettimeofday(&tv,NULL);
    return tv.tv_sec+tv.tv_usec/1e6;
}

string Clock::getDateAndTime()
{
    // current date/time based on current system
    time_t now = time(0);

    // convert now to string form
    return ctime(&now);
}


void Mat3Dvec::resize (unsigned nx, unsigned ny, unsigned nz)
{
    size_x = nx;
    size_y = ny;
    size_z = nz;
    std::cout << "(nx,ny,nz) = " << size_x << ',' << size_y << ',' << size_z << std::endl;

    vecNodeID.clear();
    vecNodeID.resize(size_x*size_y*size_z);
}

void Mat2Dvec::resize (unsigned nx, unsigned ny)
{
    size_x = nx;
    size_y = ny;

    vecNodeID.clear();
    vecNodeID.resize(size_x*size_y);
}

std::array<bool,6> getFaceFixity ()
{
    std::array<bool,6> face;
    for (unsigned d=0; d<6; d++) {
        face[d] = false;
    }
    if (InitStress[0]) {
        face[Face3] = true;
        face[Face5] = true;
    }
    if (InitStress[1]) {
        face[Face2] = true;
        face[Face4] = true;
    }
    if (InitStress[2]) {
        face[Face1] = true;
        face[Face6] = true;
    }
    for (unsigned d=0; d<6; d++) {
        bool isRestrain = false;
        for (unsigned i=0; i<6; i++) {
            isRestrain = isRestrain || BoundaryMat[d][i];
        }
        if (isRestrain) {
            face[d] = true;
        }
    }
    if (Nx<=4) {
        face[Face3] = false;
        face[Face5] = false;
    }
    if (Ny<=4) {
        face[Face2] = false;
        face[Face4] = false;
    }
    if (Nz<=4) {
        face[Face1] = false;
        face[Face6] = false;
    }
    if (ProblemID==0) {
        face[LoadFace] = true;
    }
    return face;
}
