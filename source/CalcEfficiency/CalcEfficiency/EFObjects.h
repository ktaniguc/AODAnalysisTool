#ifndef EVENTFILTER_OBJ
#define EVENTFILTER_OBJ

#include <vector>
class EFObject
{
  public :
    EFObject() :
    dRef(0),
    isPassed(-1),
    pt(0), 
    eta(0),
    phi(0)
    {};
    ~EFObject() {};

  public :
    double dRef;
    int isPassed;
    double pt; 
    double eta;
    double phi;
};

//-----------------------------------------------------//
//-----------------------------------------------------//
typedef std::vector<EFObject> EFObjects;
//-----------------------------------------------------//
//-----------------------------------------------------//

#endif //EVENTFILTER_OBJ
