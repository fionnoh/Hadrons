/*
 * Test_QED.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;


    // run setup ///////////////////////////////////////////////////////////////
    Application         application;
    std::vector<double> residual = {1e-16, 1e-15, 1e-14, 1e-13,
                                    1e-12, 1e-11, 1e-10, 1e-09}; 
    std::vector<int> resSize     = {16, 15, 14, 13,
                                    12, 11, 10, 9}; 
    std::vector<double> charge   = {1e-10, 1e-09, 1e-08, 1e-07,
                                    1e-06, 1e-05, 1e-04, 1e-03};
    std::vector<int> charSize    = {10, 9, 8, 7,
                                    6, 5, 4, 3};
    double              mass     = 0.3;
    unsigned int        mu       = 3;
    unsigned int        tJ       = 0;

    unsigned int  nt    = GridDefaultLatt()[Tp];

    std::vector<ComplexD> omega;
    omega.push_back(std::complex<double>(6.863249884465925e-02,-5.506585308274019e-02));
    omega.push_back(std::complex<double>(6.863249884465925e-02,5.506585308274019e-02));
    omega.push_back(std::complex<double>(9.901366519626265e-02,-0.000000000000000e+00));
    omega.push_back(std::complex<double>(1.260742995029118e-01,-0.000000000000000e+00));
    omega.push_back(std::complex<double>(2.113790261902896e-01,-0.000000000000000e+00));
    omega.push_back(std::complex<double>(3.419850204537295e-01,-0.000000000000000e+00));
    omega.push_back(std::complex<double>(5.423524091567911e-01,-0.000000000000000e+00));
    omega.push_back(std::complex<double>(8.309511666859551e-01,-0.000000000000000e+00));
    omega.push_back(std::complex<double>(1.182313183893475e+00,-0.000000000000000e+00));
    omega.push_back(std::complex<double>(1.458064389850479e+00,-0.000000000000000e+00));

    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
    globalPar.trajCounter.step  = 20;
    globalPar.runId             = "test";
    application.setPar(globalPar);

    // gauge field
    application.createModule<MGauge::Random>("gauge");
    // pt source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 1";
    std::string twist    = "0. 0. 0. 0.";

    // action
    MAction::ZMobiusDWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls    = 10;
    actionPar.M5    = 1.0;
    actionPar.b    = 1.5;
    actionPar.c    = 0.5;
    actionPar.mass  = mass;
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    actionPar.omega = omega;
    std::string actionZero = "ZMob";
    application.createModule<MAction::ZMobiusDWF>(actionZero, actionPar);


    for (unsigned int r = 0; r < residual.size(); ++r)
    {            
        // solvers
        MSolver::ZRBPrecCG::Par solverPar;
        solverPar.action       = actionZero;
        solverPar.residual     = residual[r];
        solverPar.maxIteration = 100000;
        std::string cgZero = "cgZero_"
                             + std::to_string(resSize[r]) + actionZero;
        application.createModule<MSolver::ZRBPrecCG>(cgZero, solverPar);

        MFermion::ZGaugeProp::Par quarkPar;
        quarkPar.solver = cgZero;
        quarkPar.source = "pt";
        std::string propZero = "propZero_" + cgZero;
        application.createModule<MFermion::ZGaugeProp>(propZero, quarkPar);

        // seq sources with conserved vector insertion
        MSource::ZSeqConserved::Par seqPar;
        seqPar.q         = propZero + "_5d";
        seqPar.source    = "pt";
        seqPar.action    = actionZero;
        seqPar.tA        = tJ;
        seqPar.tB        = tJ;
        seqPar.curr_type = Current::Vector;
        seqPar.mu_min    = mu;
        seqPar.mu_max    = mu;
        seqPar.mom       = "0. 0. 0. 0.";
        seqPar.photon    = "";
        std::string seqSrc = "seqSrc_" + std::to_string(resSize[r]);
        application.createModule<MSource::ZSeqConserved>(seqSrc, seqPar);
        // seq propagator with conserved vector insertion
        MFermion::ZGaugeProp::Par seqPropPar;
        seqPropPar.solver = cgZero;
        seqPropPar.source = seqSrc;
        std::string seqProp = "seqProp_" + std::to_string(resSize[r]);
        application.createModule<MFermion::ZGaugeProp>(seqProp, seqPropPar);
    }

    for (unsigned int c = 0; c < charge.size(); ++c)
    {
        // gauges
        MGauge::Rotation::Par rotPar;
        rotPar.charge = charge[c];
        rotPar.tA     = tJ;
        rotPar.tB     = tJ;
        rotPar.muMin  = mu;
        rotPar.muMax  = mu;
        std::string rotPlusName = "rotPlus" + std::to_string(charSize[c]);
        application.createModule<MGauge::Rotation>(rotPlusName, rotPar);

        // MGauge::Rotation::Par rotMinusPar;
        rotPar.charge = -charge[c];
        std::string rotMinusName = "rotMinus" + std::to_string(charSize[c]);
        application.createModule<MGauge::Rotation>(rotMinusName, rotPar);


        MGauge::Electrify::Par electrifyPar;
        electrifyPar.gauge   = "gauge";
        electrifyPar.emField = rotPlusName;
        electrifyPar.e       = 1;
        electrifyPar.charge  = 1;
        std::string gaugePlusName = "gaugeRotPlus" + std::to_string(charSize[c]);
        application.createModule<MGauge::Electrify>(gaugePlusName, electrifyPar);

        electrifyPar.emField = rotMinusName;
        std::string gaugeMinusName = "gaugeRotMinus" + std::to_string(charSize[c]);
        application.createModule<MGauge::Electrify>(gaugeMinusName, electrifyPar);

        // actions
        actionPar.gauge = gaugePlusName;
        std::string actionPlus = "ZMobPlus" + std::to_string(charSize[c]);
        application.createModule<MAction::ZMobiusDWF>(actionPlus, actionPar);

        actionPar.gauge = gaugeMinusName;
        std::string actionMinus = "ZMobMinus" + std::to_string(charSize[c]);
        application.createModule<MAction::ZMobiusDWF>(actionMinus, actionPar);

        for (unsigned int r = 0; r < residual.size(); ++r)
        {            
            // solvers
            std::string cgZero = "cgZero_"
                                 + std::to_string(resSize[r]) + actionZero;

            MSolver::ZRBPrecCG::Par solverPar;
            solverPar.residual     = residual[r];
            solverPar.maxIteration = 100000;
            solverPar.action       = actionPlus;
            std::string cgPlus = "cgPlus_"
                                 + std::to_string(resSize[r]) + actionPlus;
            application.createModule<MSolver::ZRBPrecCG>(cgPlus, solverPar);

            solverPar.action = actionMinus;
            std::string cgMinus = "cgMinus_"
                                  + std::to_string(resSize[r]) + actionMinus;
            application.createModule<MSolver::ZRBPrecCG>(cgMinus, solverPar);
            
            // propagators
            std::string propZero = "propZero_" + cgZero;

            MFermion::ZGaugeProp::Par quarkPar;
            quarkPar.source = "pt";
            quarkPar.solver = cgPlus;
            std::string propPlus = "propPlus_" + cgPlus;
            application.createModule<MFermion::ZGaugeProp>(propPlus, quarkPar);

            quarkPar.solver = cgMinus;
            std::string propMinus = "propMinus_" + cgMinus;
            application.createModule<MFermion::ZGaugeProp>(propMinus, quarkPar);

            std::string seqProp = "seqProp_" + std::to_string(resSize[r]);

            // test
            MUtilities::TestRotation::Par testRotPar;
            testRotPar.qRotP  = propPlus;
            testRotPar.qRotM  = propMinus;
            testRotPar.qSeq   = seqProp;
            testRotPar.charge = charge[c];
            testRotPar.origin = "0 0 0 0";
            testRotPar.mu     = mu;
            testRotPar.tJ     = tJ;
            std::string testRot = "testRot_Residual_" + std::to_string(resSize[r])
                                  + "_Charge_" + std::to_string(charSize[c]);
            application.createModule<MUtilities::TestRotation>(testRot, testRotPar);
        }
    }



    
    // execution
    application.saveParameterFile("U1.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
