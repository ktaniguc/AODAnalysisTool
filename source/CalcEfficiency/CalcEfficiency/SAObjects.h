#ifndef STANDALONE_OBJ
#define STANDALONE_OBJ

#include <vector>
#include <string>
#include <stdint.h>

class SAObject
{
  public :
    SAObject() { Clear(); };
    ~SAObject() {};

  public :
  void Clear()
  {
    isPassed = -1;
    pt = -99999; 
    eta = -99999;
    phi = -99999;
    etaBE = -99999;
    phiBE = -99999;
    etaMS = -99999;
    phiMS = -99999;
    roiNum = -1;
    roiSector = -1;
    sAddress = -2;
    tgcPt = -99999;
    ptBarrelRadius = -99999;
    ptBarrelSagitta = -99999;
    ptEndcapAlpha = -99999;
    ptEndcapBeta = -99999;
    ptEndcapRadius = -99999;
    ptCSC = -99999;
    roiEta = -99999;
    roiPhi = -99999;
    isRpcFailure = -99999;
    isTgcFailure = -99999;
    superPointR_BI = -99999;
    superPointR_BM = -99999;
    superPointR_BO = -99999;
    superPointR_EI = -99999;
    superPointR_EM = -99999;
    superPointR_EO = -99999;
    superPointR_EE = -99999;
    superPointR_CSC = -99999;
    superPointR_BEE = -99999;
    superPointR_BME = -99999;
    superPointZ_BI = -99999;
    superPointZ_BM = -99999;
    superPointZ_BO = -99999;
    superPointZ_EI = -99999;
    superPointZ_EM = -99999;
    superPointZ_EO = -99999;
    superPointZ_EE = -99999;
    superPointZ_CSC = -99999;
    superPointZ_BEE = -99999;
    superPointZ_BME = -99999;
    superPointSlope_BI = -99999;
    superPointSlope_BM = -99999;
    superPointSlope_BO = -99999;
    superPointSlope_EI = -99999;
    superPointSlope_EM = -99999;
    superPointSlope_EO = -99999;
    superPointSlope_EE = -99999;
    superPointSlope_CSC = -99999;
    superPointSlope_BEE = -99999;
    superPointSlope_BME = -99999;
    superPointIntercept_BI = -99999;
    superPointIntercept_BM = -99999;
    superPointIntercept_BO = -99999;
    superPointIntercept_EI = -99999;
    superPointIntercept_EM = -99999;
    superPointIntercept_EO = -99999;
    superPointIntercept_EE = -99999;
    superPointIntercept_CSC = -99999;
    superPointIntercept_BEE = -99999;
    superPointIntercept_BME = -99999;
    superPointChi2_BI = -99999;
    superPointChi2_BM = -99999;
    superPointChi2_BO = -99999;
    superPointChi2_EI = -99999;
    superPointChi2_EM = -99999;
    superPointChi2_EO = -99999;
    superPointChi2_EE = -99999;
    superPointChi2_CSC = -99999;
    superPointChi2_BEE = -99999;
    superPointChi2_BME = -99999;
    rpcHitX.clear();
    rpcHitY.clear();
    rpcHitZ.clear();
    rpcHitR.clear();
    rpcHitEta.clear();
    rpcHitPhi.clear();
    rpcHitLayer.clear();
    rpcHitStationName.clear();
    rpcHitMeasPhi.clear();
    tgcHitZ.clear();
    tgcHitR.clear();
    tgcHitEta.clear();
    tgcHitPhi.clear();
    tgcHitWidth.clear();
    tgcHitStationNum.clear();
    tgcHitIsStrip.clear();
    tgcHitBCTag.clear();
    tgcHitInRoad.clear();

    mdtHitIsOutlier.clear();
    mdtHitChamber.clear();
    mdtHitR.clear();
    mdtHitZ.clear();
    mdtHitPhi.clear();
    mdtHitResidual.clear();

    roadAw.clear();
    roadBw.clear();
    zMin.clear();
    zMax.clear();
    rMin.clear();
    rMax.clear();
    etaMin.clear();
    etaMax.clear();
    stgcClusterR.clear();
    stgcClusterZ.clear();
    stgcClusterEta.clear();
    stgcClusterPhi.clear();
    stgcClusterResidualR.clear();
    stgcClusterResidualPhi.clear();
    stgcClusterStationEta.clear();
    stgcClusterStationPhi.clear();
    stgcClusterStationName.clear();
    stgcClusterType.clear();
    stgcClusterIsOutlier.clear();
    stgcClusterLayer.clear();
    mmClusterR.clear();
    mmClusterZ.clear();
    mmClusterEta.clear();
    mmClusterPhi.clear();
    mmClusterResidualR.clear();
    mmClusterResidualPhi.clear();
    mmClusterStationEta.clear();
    mmClusterStationPhi.clear();
    mmClusterStationName.clear();
    mmClusterIsOutlier.clear();
    mmClusterLayer.clear();
  };

  public :
    int isPassed;
    double pt; 
    double eta;
    double phi;
    double etaBE;
    double phiBE;
    double etaMS;
    double phiMS;
    int roiNum;
    int roiSector;
    int sAddress;  // -1 : endcap, else : barrel
    double tgcPt;
    double ptBarrelRadius;
    double ptBarrelSagitta;
    double ptEndcapAlpha;
    double ptEndcapBeta;
    double ptEndcapRadius;
    double ptCSC;
    float roiEta;
    float roiPhi;
    int isRpcFailure;
    int isTgcFailure;
    //the measured radious of the muon in one particular super point
    double superPointR_BI;
    double superPointR_BM;
    double superPointR_BO;
    double superPointR_EI;
    double superPointR_EM;
    double superPointR_EO;
    double superPointR_EE;
    double superPointR_CSC;
    double superPointR_BEE;
    double superPointR_BME;
    //the measured Z position of the muon in one particular super point
    double superPointZ_BI;
    double superPointZ_BM;
    double superPointZ_BO;
    double superPointZ_EI;
    double superPointZ_EM;
    double superPointZ_EO;
    double superPointZ_EE;
    double superPointZ_CSC;
    double superPointZ_BEE;
    double superPointZ_BME;
    //the measured slope of the muon in one particular super point
    double superPointSlope_BI;
    double superPointSlope_BM;
    double superPointSlope_BO;
    double superPointSlope_EI;
    double superPointSlope_EM;
    double superPointSlope_EO;
    double superPointSlope_EE;
    double superPointSlope_CSC;
    double superPointSlope_BEE;
    double superPointSlope_BME;
    //the measured intercept of the muon in one particular super point
    double superPointIntercept_BI;
    double superPointIntercept_BM;
    double superPointIntercept_BO;
    double superPointIntercept_EI;
    double superPointIntercept_EM;
    double superPointIntercept_EO;
    double superPointIntercept_EE;
    double superPointIntercept_CSC;
    double superPointIntercept_BEE;
    double superPointIntercept_BME;
    //the chi2 of the fit in one particular super point
    double superPointChi2_BI;
    double superPointChi2_BM;
    double superPointChi2_BO;
    double superPointChi2_EI;
    double superPointChi2_EM;
    double superPointChi2_EO;
    double superPointChi2_EE;
    double superPointChi2_CSC;
    double superPointChi2_BEE;
    double superPointChi2_BME;
    //RPC hits
    std::vector < float > rpcHitX;
    std::vector < float > rpcHitY;
    std::vector < float > rpcHitZ;
    std::vector < float > rpcHitR;
    std::vector < float > rpcHitEta;
    std::vector < float > rpcHitPhi;
    std::vector < uint32_t > rpcHitLayer;
    std::vector < std::string > rpcHitStationName;
    std::vector < uint32_t > rpcHitMeasPhi;
    //TGC hits
    std::vector < float > tgcHitZ;
    std::vector < float > tgcHitR;
    std::vector < float > tgcHitEta;
    std::vector < float > tgcHitPhi;
    std::vector < float > tgcHitWidth;
    std::vector < int > tgcHitStationNum;
    std::vector < bool > tgcHitIsStrip;
    std::vector < int > tgcHitBCTag;
    std::vector < bool > tgcHitInRoad;
    //MDT hits
    std::vector < int > mdtHitIsOutlier;
    std::vector < int > mdtHitChamber;
    std::vector < float > mdtHitR;
    std::vector < float > mdtHitZ;
    std::vector < float > mdtHitPhi;
    std::vector < float > mdtHitResidual;

    std::vector < float > roadAw;
    std::vector < float > roadBw;
    std::vector < float > zMin;
    std::vector < float > zMax;
    std::vector < float > rMin;
    std::vector < float > rMax;
    std::vector < float > etaMin;
    std::vector < float > etaMax;
    //

    //NSW
    std::vector < float > stgcClusterR;
    std::vector < float > stgcClusterZ;
    std::vector < float > stgcClusterEta;
    std::vector < float > stgcClusterPhi;
    std::vector < float > stgcClusterResidualR;
    std::vector < float > stgcClusterResidualPhi;
    std::vector < int > stgcClusterStationEta;
    std::vector < int > stgcClusterStationPhi;
    std::vector < int > stgcClusterStationName;
    std::vector < int > stgcClusterType;
    std::vector < int > stgcClusterIsOutlier;
    std::vector < unsigned int > stgcClusterLayer;
    std::vector < float > mmClusterR;
    std::vector < float > mmClusterZ;
    std::vector < float > mmClusterEta;
    std::vector < float > mmClusterPhi;
    std::vector < float > mmClusterResidualR;
    std::vector < float > mmClusterResidualPhi;
    std::vector < int > mmClusterStationEta;
    std::vector < int > mmClusterStationPhi;
    std::vector < int > mmClusterStationName;
    std::vector < int > mmClusterIsOutlier;
    std::vector < unsigned int > mmClusterLayer;
    //

};

//-----------------------------------------------------//
//-----------------------------------------------------//
typedef std::vector<SAObject> SAObjects;
//-----------------------------------------------------//
//-----------------------------------------------------//

#endif //STANDALONE_OBJ
