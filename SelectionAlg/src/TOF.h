#ifndef TOF_HH
#define TOF_HH

#include <iostream>
#include <math.h>

#include <TVector3.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

struct TOFParams
{
    TVector3 vertex;
    TVector3 pmt;
    double z_interface;
    double n_LS;
    double n_Water;
};


class TOFCalculator
{
private:
    double c;
    double PMT_R;
	double LS_R;
    
    TOFParams params;

    double TOFFunction(const double* xy);
    double TOFMinimizer();
public:
    // TOFCalculator(TVector3 Vertex, TVector3 PMT, double interface); //mm
    TOFCalculator();
    ~TOFCalculator();

    void setVertex(TVector3 vtx){params.vertex = vtx;};
    void setPMTVertex(TVector3 PMTvtx){params.pmt = PMTvtx;};
    void setInterface(double interface){params.z_interface = interface;};
    double CalTOF();
};


#endif
