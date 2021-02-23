#ifndef REQUIREMENTS_OBJ
#define REQUIREMENTS_OBJ

class Requirement
{
  public :
    Requirement() : 
    reqdRl1(0),
    CBmatchingdR(0),
    EFmatchingdR(0),
    invMass(0),
    tpdPhi(0)
    {};
    ~Requirement() {};

  public : 
    double reqdRl1;
    double CBmatchingdR;
    double EFmatchingdR;
    double invMass;
    double tpdPhi;
};

//-----------------------------------------------------//
//-----------------------------------------------------//
typedef std::vector<Requirement> Requirements;
//-----------------------------------------------------//
//-----------------------------------------------------//

#endif //REQUIREMENTS_OBJ
