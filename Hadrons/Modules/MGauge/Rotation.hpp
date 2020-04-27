#ifndef Hadrons_MGauge_Rotation_hpp_
#define Hadrons_MGauge_Rotation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Rotation                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class RotationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RotationPar,
                                    double,       charge,
                                    unsigned int, tA,
                                    unsigned int, tB,
                                    unsigned int, muMin,
                                    unsigned int, muMax);
};

template <typename GImpl>
class TRotation: public Module<RotationPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    typedef PhotonR::GaugeField     rotField; // borrowing for the time being
public:
    // constructor
    TRotation(const std::string name);
    // destructor
    virtual ~TRotation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasT_{false};
    std::string tName_;
};

MODULE_REGISTER_TMP(Rotation, TRotation<GIMPL>, MGauge);

/******************************************************************************
 *                 TRotation implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TRotation<GImpl>::TRotation(const std::string name)
: Module<RotationPar>(name)
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TRotation<GImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TRotation<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TRotation<GImpl>::setup(void)
{
    envCreateLat(rotField, getName());
    envTmpLat(rotField, "tmp");
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envTmpLat(LatticeComplex, "eta");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TRotation<GImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Rotating gauge links with charge " << par().charge
        << " for t = " << par().tA << std::endl;
    }
    else
    {
        LOG(Message) << "Rotating gauge links with charge " << par().charge
        << " for tA = " << par().tA << " tB = " << par().tB << std::endl;
    }
    
    auto    &theta = envGet(rotField, getName());
    envGetTmp(rotField, tmp);
    auto    &t = envGet(Lattice<iScalar<vInteger>>, tName_);
    envGetTmp(LatticeComplex, eta);

    if (!hasT_)
    {
        LatticeCoordinate(t, Tp);
        hasT_ = true;
    }
    eta = 1.*par().charge;
    eta = where((t >= par().tA) and (t <= par().tB), eta, 0.*eta);
    theta = 1.;
    tmp   = 0.;
    theta = theta*eta;

    for(unsigned int mu=par().muMin;mu<=par().muMax;mu++)
    {
        auto U = PeekIndex<LorentzIndex>(theta, mu);
        PokeIndex<LorentzIndex>(tmp, U, mu);
    }
    theta = tmp;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_Rotation_hpp_
