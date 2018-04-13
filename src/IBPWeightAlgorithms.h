/*////////////////////////////////////////////////////////////////////////////
 Copyright 2018 Istituto Nazionale di Fisica Nucleare

 Licensed under the EUPL, Version 1.2 or - as soon they will be approved by
 the European Commission - subsequent versions of the EUPL (the "Licence").
 You may not use this work except in compliance with the Licence.

 You may obtain a copy of the Licence at:

 https://joinup.ec.europa.eu/software/page/eupl

 Unless required by applicable law or agreed to in writing, software
 distributed under the Licence is distributed on an "AS IS" basis, WITHOUT
 WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 Licence for the specific language governing permissions and limitations under
 the Licence.
////////////////////////////////////////////////////////////////////////////*/


#ifndef IBPWEIGHTALGORITHMS_H
#define IBPWEIGHTALGORITHMS_H


#include <Core/StaticInterface.h>
#include <Math/Dense.h>


using namespace uLib;


template <class EventT>
class IBPWeightAlgorithm_PW {

public:
    static void evaluate(EventT *evc, Scalarf nominalp) {
        float pw_A = 1.42857;
        float pw_epsilon = 4.0;
        float Xres = 0;
        for (int j = evc->elements.size(); j --> 0;) //BACKWARD
        {
            float L = evc->elements[j].Wij(0,0);

            evc->elements[j].pw = pw_A * (nominalp) *
                    sqrt(pw_epsilon/(Xres + pw_epsilon));
            //        Xres += (evc->elements[j].voxel->Value * 1.E6 < 2.5) ?
            //                    L * evc->elements[j].voxel->Value * 40000 :
            //                    L * 2.5 * 0.04;
            Xres += L * evc->elements[j].voxel->Value * 40000;
        }
    }

};


template <class EventT>
class IBPWeightAlgorithm_SW {

public:
    static void evaluate(EventT *evc, Scalarf nominalp) {
        float scale = 0.8;   //this is experimental and empirical!
        float alpha = 0.436; // this is the distribution peak angle, 25 degrees
        // note: albeit NOT efficient, this procedure is fundamental in testing the algorithm capability
        TODO("Optimize PW functions");
        float b1 = 13.52; // Iron Values
        float b2 = 319.9;
        float c1 = 3.73E-4;
        float c2 = 2.55E-2;
        float d = 2.33;
        float e = 1.56;

        float sw_epsilon_1 = scale * sqrt(b1 + c1 * alpha * alpha);
        float sw_epsilon_2 = scale * sqrt(b2 + c2 * alpha * alpha);
        float sw_A = d + e * cos(alpha) * cos(alpha);

        float Xres = 0;
        Matrix4f Wij;
        for (int j = evc->elements.size(); j --> 0;) //BACKWARD
        {
            Wij = evc->elements[j].Wij;

            evc->elements[j].pw = (nominalp) * sqrt( sw_A /
                    (Xres/sw_epsilon_1 + pow(Xres/sw_epsilon_2,2) + 1) );
        /*
            Xres += (evc->elements[j].voxel->density * 1.E6 < 2.5) ?
                    Wij[0] * evc->elements[j].voxel->density * 40000 :
                    Wij[0] * 2.5 * 0.04;
            */
            Xres += Wij(0,0) * evc->elements[j].voxel->Value * 40000; //<- previous version of
                                          // p_weight: bugged but enhancing!
        }
    }


};


template <class EventT>
class IBAnalyzerAlgorithm_CW {
public:
    static void evaluate(EventT *evc, Scalarf nominalp){
        float cw_A = 1.42857;
        float cw_epsilon = 50;
        float X0_tot = 0;
        Matrix4f Wij;
        for (int j = 0; j < evc->elements.size(); ++j)
        {
            Wij = evc->elements[j].Wij;
            /*X0_tot += (evc->elements[j].voxel->density * 1.E6 < 2.5 ) ?
                     Wij[0] * evc->elements[j].voxel->density * 40000 :
                     Wij[0] * 2.5 * 0.04;
             */
            X0_tot += Wij(0,0) * evc->elements[j].voxel->Value * 40000;
        }

        float _pw = cw_A * (nominalp) *
                sqrt( cw_epsilon / ( X0_tot + cw_epsilon ));
        for (int i = 0; i < evc->elements.size(); ++i )
        {
            evc->elements[i].pw = _pw;
        }
    }
};




#endif // IBPWEIGHTALGORITHMS_H
