/*//////////////////////////////////////////////////////////////////////////////
// CMT Cosmic Muon Tomography project //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  Copyright (c) 2014, Universita' degli Studi di Padova, INFN sez. di Padova

  Coordinators: Prof. Gianni Zumerle < gianni.zumerle@pd.infn.it >
                Paolo Checchia       < paolo.checchia@pd.infn.it >

  Authors: Andrea Rigoni Garola < andrea.rigoni@pd.infn.it >
           Matteo Furlan        < nuright@gmail.com >
           Sara Vanini          < sara.vanini@pd.infn.it >

  All rights reserved
  ------------------------------------------------------------------

  This file can not be copied and/or distributed without the express
  permission of  Prof. Gianni Zumerle  < gianni.zumerle@pd.infn.it >

//////////////////////////////////////////////////////////////////////////////*/




#include "TFile.h"
#include "TTree.h"

#include "IBMinimizationVariablesEvaluator.h"
#include "IBNormalPlaneMinimizationVariablesEvaluator.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "testing-prototype.h"

#include "IBVoxRaytracer.h"
#include "IBVoxCollectionCap.h"

#include "Core/Options.h"

#include "Vtk/uLibVtkViewer.h"
#include "Vtk/vtkContainerBox.h"
#include "Vtk/vtkMuonScatter.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace uLib;





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// main


namespace {
struct Config : Options {

    // members //
    bool vtk;
    struct Img { Vector3f size,pos,spacing; } img;
    struct Var { bool free; float alpha; } eval;
    struct Mu { HLine3f in,out; } mu;

    Config() : Options("Config:") {
        add_options()
                ("help","geet the help")
                ("vtk",&vtk,true,"show vtk viewer")
                ("img.size",&img.size,Vector3f(100,100,100),"img size boudaries")
                ("img.position",&img.pos,Vector3f(0,0,0),"img position")
                ("img.spacing",&img.spacing,Vector3f(5,5,5),"img spacing")
                ("eval.free",&eval.free,false,"enable free rotation" )
                ("eval.alpha",&eval.alpha,(float)0,"rotation")
                ("mu.in",&mu.in, {HPoint3f(50,50,110),HVector3f(0,0,-1)} , "line in" )
                ("mu.out",&mu.out, {HPoint3f(110,110,-10),HVector3f(1,1,-1)} , "line in" )
                ;
    }
};
} // namespace




int main(int argc, char **argv) {
    BEGIN_TESTING(IBMinVar2);

    Config p;
    p.parse_command_line(argc,argv);

    Vector3i dim = (p.img.size.array() / p.img.spacing.array()).cast<int>();
    IBVoxCollection voxels(dim);

    voxels.SetSpacing(p.img.spacing);
    voxels.SetPosition(p.img.pos);

    IBVoxel zero = {0,0,0};
    voxels.InitVoxels(zero);

    MuonScatter mu;
    mu.LineIn() = p.mu.in;
    mu.LineOut() = p.mu.out;

    VoxRaytracer tracer(voxels);
    //    IBMinimizationVariablesEvaluator *eval = IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    IBNormalPlaneMinimizationVariablesEvaluator *eval = new IBNormalPlaneMinimizationVariablesEvaluator();
    eval->$$.use_free_rotation = p.eval.free;
    eval->$$.alphaXZ = p.eval.alpha;
    eval->setRaytracer(&tracer);

    bool res = eval->evaluate(mu);
    TEST1( res );
    Vector4f data = eval->getDataVector();

    if(res) {
        std::cout << "data: " << data.transpose() << "\n";
    }
    else {
        std::cerr << "eval failed\n";
    }

    {
        float dist = Vector2f(data(1),data(3)).norm();
        // compute test dist (point to line distance) //
        Vector3f a = mu.LineIn().origin.head(3);
        Vector3f n = mu.LineIn().direction.head(3);
        n = n.array() / n.norm();

        HPoint3f p;
        TEST1( tracer.GetExitPoint(mu.LineOut(),p) );

        Vector3f ap = a-p.head(3);
        float prj = ap.transpose() * n;
        float test_dist = (ap - prj*n).norm();

        Eigen::ParametrizedLine<float,3> eline(a,n);
        float test_dist2 = eline.distance(p.head(3));

        std::cout << "dists: " << dist << " vs " << test_dist << " vs " << test_dist2 << "\n";
        TEST1( fabs(dist-test_dist) < 0.001 );
        TEST1( fabs(dist-test_dist2) < 0.001 );
    }



    if(p.vtk) {
        Vtk::Viewer viewer;
        Vtk::vtkMuonScatter v_mu(mu);
        viewer.AddPuppet(v_mu);
        ContainerBox box;
        box.SetSize(p.img.size);
        box.SetPosition(voxels.GetPosition());
        Vtk::vtkContainerBox v_box(box);
        viewer.AddPuppet(v_box);
        viewer.Start();
    }

    delete eval;
    END_TESTING
}
