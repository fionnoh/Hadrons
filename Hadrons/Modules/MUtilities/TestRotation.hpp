#ifndef Hadrons_MUtilities_TestRotation_hpp_
#define Hadrons_MUtilities_TestRotation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         TestRotation                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class TestRotationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TestRotationPar,
                                    std::string,  qRotP,
                                    std::string,  qRotM,
                                    std::string,  qSeq,
                                    double,       charge,
                                    std::string,  origin,
                                    unsigned int, mu,
                                    unsigned int, tJ);
};

template <typename FImpl>
class TTestRotation: public Module<TestRotationPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TTestRotation(const std::string name);
    // destructor
    virtual ~TTestRotation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TestRotation, TTestRotation<FIMPL>, MUtilities);

/******************************************************************************
 *                 TTestRotation implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTestRotation<FImpl>::TTestRotation(const std::string name)
: Module<TestRotationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTestRotation<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qRotP, par().qRotM, par().qSeq};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TTestRotation<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTestRotation<FImpl>::setup(void)
{
    envTmpLat(PropagatorField, "deriv");
    envTmpLat(LatticeComplex, "c");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTestRotation<FImpl>::execute(void)
{
    auto   &qRotP  = envGet(PropagatorField, par().qRotP);
    auto   &qRotM  = envGet(PropagatorField, par().qRotM);
    auto   &qSeq   = envGet(PropagatorField, par().qSeq);
    double charge2 = 2.0*par().charge;
    unsigned int mu = par().mu;
    Complex i(0.0,1.0);
    int nt = env().getDim(Tp);

    SitePropagator        qSite;
    Complex               seq_S, seq_V, rot_S, rot_V, diff_VS;
    std::vector<int>      siteCoord = strToVec<int>(par().origin);
    std::vector<TComplex> seq_buf, rot_buf;

    LOG(Message) << "Comparing sequential insertion to U(1) rotation for "
                 << "charge = " << par().charge
                 << ", mu = " << mu
                 << ", tJ = " << par().tJ << std::endl; 

    LOG(Message) << "====================================" << std::endl; 
    LOG(Message) << "Peeking site " << par().origin << std::endl; 
    peekSite(qSite, qSeq, siteCoord);
    seq_S = trace(qSite);
    seq_V = trace(qSite*Gamma::gmu[mu]);
    LOG(Message) << std::setprecision(14) << "Seq S  = " << abs(seq_S) << std::endl;
    LOG(Message) << "Seq V  = " << abs(seq_V) << std::endl;

    envGetTmp(PropagatorField, deriv);
    deriv = (i/charge2)*(qRotP - qRotM);
    peekSite(qSite, deriv, siteCoord);
    rot_S = trace(qSite);
    rot_V = trace(qSite*Gamma::gmu[mu]);
    LOG(Message) << "Rot S  = " << abs(rot_S) << std::endl;
    LOG(Message) << "Rot V  = " << abs(rot_V) << std::endl;

    rot_S -= seq_S;
    rot_V -= seq_V;
    LOG(Message) << "Diff S  = " << abs(rot_S) << std::endl;
    LOG(Message) << "Diff V  = " << abs(rot_V) << std::endl;
    LOG(Message) << "====================================" << std::endl; 


    envGetTmp(LatticeComplex, c);
    c = trace(qSeq);
    sliceSum(c, seq_buf, Tp);

    c = trace(deriv);
    sliceSum(c, rot_buf, Tp);

    for (int t = 0; t < nt; ++t)
    {
        seq_S = TensorRemove(seq_buf[t]);
        rot_S = TensorRemove(rot_buf[t]);
        diff_VS = seq_S - rot_S;
        LOG(Message) << "t = " << t
                     << ", Seq S  = " << abs(seq_S)
                     << ", Rot S = " << abs(rot_S)
                     << ", diff = " << abs(diff_VS) << std::endl;
    }

    LOG(Message) << "====================================" << std::endl; 
    c = trace(qSeq*Gamma::gmu[mu]);
    sliceSum(c, seq_buf, Tp);
    c = trace(deriv*Gamma::gmu[mu]);
    sliceSum(c, rot_buf, Tp);


    for (int t = 0; t < nt; ++t)
    {
        seq_V = TensorRemove(seq_buf[t]);
        rot_V = TensorRemove(rot_buf[t]);
        diff_VS = seq_V - rot_V;
        LOG(Message) << "t = " << t
                     << ", Seq V  = " << abs(seq_V)
                     << ", Rot V = " << abs(rot_V)
                     << ", diff = " << abs(diff_VS) << std::endl;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_TestRotation_hpp_
