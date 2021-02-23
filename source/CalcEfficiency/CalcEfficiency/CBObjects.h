#ifndef COMBINED_OBJ
#define COMBINED_OBJ

#include <vector>
class CBObject
{
  public :
    CBObject() :
    isPassed(-1),
    pt(0), 
    eta(0),
    phi(0),
    roiNum(-1),
    roiSector(-1),
    idpt(0), 
    ideta(0),
    idphi(0)
    {};
    ~CBObject() {};

  public :
    int isPassed;
    double pt; 
    double eta;
    double phi;
    int roiNum;
    int roiSector;
    double idpt;
    double ideta;
    double idphi;
};

//-----------------------------------------------------//
//-----------------------------------------------------//
typedef std::vector<CBObject> CBObjects;
//-----------------------------------------------------//
//-----------------------------------------------------//

#endif //COMBINED_OBJ
