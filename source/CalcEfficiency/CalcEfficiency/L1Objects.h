#ifndef LVL1ROIS_OBJ
#define LVL1ROIS_OBJ

#include <vector>
#include <stdint.h>
class L1Object
{
  public :
    L1Object() : 
    dRl1(0),
    isPassed(-1),
    isPassedByEvent(false),
    isPassedByEvent_TDT(false),
    eta(0),
    phi(0),
    thrValue(0),
    roiNum(-1),
    thrNumber(0),
    roiSector(-1),
    roiWord(0),
    isMoreCandInRoI(false)
    {};
    ~L1Object() {};

  public : 
    double dRl1;
    int isPassed;
    bool isPassedByEvent;
    bool isPassedByEvent_TDT;
    double eta;
    double phi;
    double thrValue;
    int roiNum;
    int thrNumber;
    int roiSector;
    uint32_t roiWord;
    bool isMoreCandInRoI;
};

//-----------------------------------------------------//
//-----------------------------------------------------//
typedef std::vector<L1Object> L1Objects;
//-----------------------------------------------------//
//-----------------------------------------------------//

#endif //LVL1ROIS_OBJ
