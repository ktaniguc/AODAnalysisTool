#ifndef OFFLINE_OBJ
#define OFFLINE_OBJ

#include <vector>

class OfflineObject
{
  public :
    OfflineObject() :
    pt(-99999), 
    eta(-99999),
    phi(-99999),
    extEta(-99999),
    extPhi(-99999),
    tpextdR(-99999),
    d0(-99999),
    z0(-99999),
    segmentN(0)
    {
      for(int i_seg = 0; i_seg < SegmentMaxNumber; i_seg++){
        segmentX[i_seg]              = -99999.;
        segmentY[i_seg]              = -99999.;
        segmentZ[i_seg]              = -99999.;
        segmentPx[i_seg]             = -99999.;
        segmentPy[i_seg]             = -99999.;
        segmentPz[i_seg]             = -99999.;
        segmentChiSquared[i_seg]     = -99999.;
        segmentNumberDoF[i_seg]      = -99999.;
        segmentSector[i_seg]         = -99999.;
        segmentChamberIndex[i_seg]   = -99999.;
        segmentEtaIndex[i_seg]       = -99999.;
        segmentNPrecisionHits[i_seg] = -99999.;
        segmentNPhiLayers[i_seg]     = -99999.;
        segmentNTrigEtaLayers[i_seg] = -99999.;
      }
    };
    ~OfflineObject() {};

  public :
    static const int SegmentMaxNumber = 10;

    double pt; 
    double eta;
    double phi;
    double extEta;
    double extPhi;
    double tpextdR;
    double d0;
    double z0;
    int segmentN;
    double segmentX[SegmentMaxNumber];
    double segmentY[SegmentMaxNumber];
    double segmentZ[SegmentMaxNumber];
    double segmentPx[SegmentMaxNumber];
    double segmentPy[SegmentMaxNumber];
    double segmentPz[SegmentMaxNumber];
    double segmentChiSquared[SegmentMaxNumber];
    double segmentNumberDoF[SegmentMaxNumber];
    double segmentSector[SegmentMaxNumber];
    double segmentChamberIndex[SegmentMaxNumber];
    double segmentEtaIndex[SegmentMaxNumber];
    double segmentNPrecisionHits[SegmentMaxNumber];
    double segmentNPhiLayers[SegmentMaxNumber];
    double segmentNTrigEtaLayers[SegmentMaxNumber];

};

//-----------------------------------------------------//
//-----------------------------------------------------//
typedef std::vector<OfflineObject> OfflineObjects;
//-----------------------------------------------------//
//-----------------------------------------------------//

#endif //OFFLINE_OBJ
