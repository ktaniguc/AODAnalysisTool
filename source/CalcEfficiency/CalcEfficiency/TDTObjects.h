#ifndef TDT___OBJ
#define TDT___OBJ

#include <vector>
#include <string>
#include <stdint.h>
#include "SAObjects.h"
#include "CBObjects.h"

class TDTObject
{
  public :
    TDTObject() { Clear(); };
    ~TDTObject() {};

  public :
  void Clear()
  {
    isPassedChain = false;
    isPassedL1_evt = false;
    isPassedSA_evt = false;
    isPassedCB_evt = false;
    isPassedSAIO_evt = false;
    isPassedCBIO_evt = false;
    isPassedEF_evt = false;
    isPassedL1.clear();
    isPassedSA.clear();
    isPassedCB.clear();
    isPassedSAIO.clear();
    isPassedCBIO.clear();
    isPassedEF.clear();
    L1RoINumber.clear();
    L1RoISector.clear();
    SARoINumber.clear();
    SARoISector.clear();
    CBRoINumber.clear();
    CBRoISector.clear();
    SAIORoINumber.clear();
    SAIORoISector.clear();
    CBIORoINumber.clear();
    CBIORoISector.clear();
    m_l1.clear();
    m_sa.clear();
    m_saio.clear();
    m_cb.clear();
    m_cbio.clear();
    m_ef.clear();
  };

  public :
    std::string trigChainName;
    bool isPassedChain;
    bool isPassedL1_evt;
    bool isPassedSA_evt;
    bool isPassedCB_evt;
    bool isPassedSAIO_evt;
    bool isPassedCBIO_evt;
    bool isPassedEF_evt;
    std::vector<bool> isPassedL1;
    std::vector<bool> isPassedSA;
    std::vector<bool> isPassedCB;
    std::vector<bool> isPassedSAIO;
    std::vector<bool> isPassedCBIO;
    std::vector<bool> isPassedEF;
    std::vector<int> L1RoINumber;
    std::vector<int> L1RoISector;
    std::vector<int> SARoINumber;
    std::vector<int> SARoISector;
    std::vector<int> CBRoINumber;
    std::vector<int> CBRoISector;
    std::vector<int> SAIORoINumber;
    std::vector<int> SAIORoISector;
    std::vector<int> CBIORoINumber;
    std::vector<int> CBIORoISector;
    L1Objects m_l1;
    SAObjects m_sa;
    CBObjects m_cb;
    SAObjects m_saio;
    CBObjects m_cbio;
    EFObjects m_ef;
};

//-----------------------------------------------------//
//-----------------------------------------------------//
typedef std::vector<TDTObject> TDTObjects;
//-----------------------------------------------------//
//-----------------------------------------------------//

#endif //TDT___OBJ
