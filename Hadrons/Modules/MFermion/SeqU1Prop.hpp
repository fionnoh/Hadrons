#ifndef Hadrons_MFermion_SeqU1Prop_hpp_
#define Hadrons_MFermion_SeqU1Prop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SeqU1Prop                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class SeqU1PropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqU1PropPar,
                                    std::string,  qRotP,
                                    std::string,  qRotM,
                                    double,       charge);
};

template <typename FImpl>
class TSeqU1Prop: public Module<SeqU1PropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSeqU1Prop(const std::string name);
    // destructor
    virtual ~TSeqU1Prop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SeqU1Prop, TSeqU1Prop<FIMPL>, MFermion);
MODULE_REGISTER_TMP(ZSeqU1Prop, TSeqU1Prop<ZFIMPL>, MFermion);

/******************************************************************************
 *                 TSeqU1Prop implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqU1Prop<FImpl>::TSeqU1Prop(const std::string name)
: Module<SeqU1PropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqU1Prop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qRotP, par().qRotM};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqU1Prop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqU1Prop<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqU1Prop<FImpl>::execute(void)
{
    auto   &qRotP  = envGet(PropagatorField, par().qRotP);
    auto   &qRotM  = envGet(PropagatorField, par().qRotM);
    double charge2 = 2.0*par().charge;
    auto   &prop   = envGet(PropagatorField, getName());

    Complex i(0.0,1.0);
    prop = (i/charge2)*(qRotP - qRotM);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_SeqU1Prop_hpp_
