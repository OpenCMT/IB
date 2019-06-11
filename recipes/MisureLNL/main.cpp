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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
// Sara Vanini
// Analysis code for LNL data image reconstruction and detector rotation
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////// */


#include <iostream>
#include <fstream>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>
#include <pcl/common/common.h>
#include <pcl/common/angles.h>
#include <pcl/common/transforms.h>
#include <pcl/point_cloud.h>
#include <pcl/registration/transformation_estimation_svd.h>

#include <TFile.h>
#include <TTree.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>

#include "Core/Options.h"
#include "IBPocaEvaluator.h"
#include "IBAnalyzerEM.h"
#include "IBAnalyzerTrackCount.h"
#include "IBAnalyzerEMAlgorithm.h"
#include "IBVoxRaytracer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBNormalPlaneMinimizationVariablesEvaluator.h"
#include "IBVoxCollection.h"
#include "IBMuonError.h"
#include "Detectors/MuonScatter.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeLNLdataReader.h"
#include "IBVoxFilters.h"
#include "IBAnalyzerWPoca.h"
#include "IBAnalyzerPoca.h"

#include "IBAnalyzerWTrackLengths.h"
#include "IBMAPUpdateDensityAlgorithms.h"

#include "IBAnalyzerEMAlgorithmSGA.h"
#include "IBAnalyzerEMAlgorithmMGA.h"

#include <Math/Transform.h>

#include "IBMuonCollection.h"
#include "IB.h"

using namespace uLib;

////////////////////////////////////////////////////////////////////////////////
/////////////////////// RECIPE PARAMETERS ////////////////
////////////////////////////////////////////////////////////////////////////////
namespace {
static struct Parameters : Options
{
    std::string execute;

    char  *file_in;
    char  *file_out;
    float minutes;
    float start;
    float voxsize;
    bool dump;

    float momentum;
    float sijcut;
    Vector3i iterations;
    std::string analyzers;

    struct detTransform{
        std::string datafile;
        float min;
        float start;
        bool transform;
        bool findMatrix;
        bool voxelReferenceSystem;
        bool dump;
        std::string cloud1;
        std::string cloud2;
        Vector3f rotation;
        Vector3f translation;
    } detector;

    struct initImage {
        float lambda_air;
        bool maskVoxels;
        float maskThreshold;
        bool axis;
        int sliceAxis;
        bool reweight;
        HPoint3f B_block;
        HPoint3f E_block;
    } image;

    struct analysisOptions {
        int sets;
        int nslices;
        float sliceThickness;
        HVector3f sliceTolerance;
        float sliceBegin;
    } analysis;


    Parameters(const char *hello = "Program options usage: .exe file_in.root file_out 0 0 voxsize [options]") : Options(hello) {
        add_options()
                ("help", "printout help")

                // GENERAL //
                ("execute",       &execute,std::string("reconstruction"),"reconstruction clip test_init_image *_analysis")
                ("momentum",      &momentum,  (float)0.95,  "momentum [GeV]")
                ("sijcut",        &sijcut,    (float)30.,   "SijCut value")
                ("iterations",    &iterations, Vector3i(1, 1000, 5000), "iterations [start drop end]")
                ("analyzers",    &analyzers, std::string("no analyzer selected"),"selected analyzers class names")
                ("dump",    &dump,  (bool)0,"Dump muon root file")

                // DETECTOR //
                ("detector.datafile",&detector.datafile,std::string("datafile.txt"),"Muon data file from detector")
                ("detector.min",&detector.min,(float)0,"Daq time requested [min]")
                ("detector.start",&detector.start,(float)0,"Daq time starting minutes")
                ("detector.transform",    &detector.transform,  (bool)0,"Roto-translate muon in coordinate system given by rotation+translation")
                ("detector.findMatrix",    &detector.findMatrix,  (bool)0,"Find transformation matrix from input point clouds")
                ("detector.voxelReferenceSystem",    &detector.voxelReferenceSystem,  (bool)0,"Use voxel reference system, otherwise Paolo's coord")
                ("detector.dump",    &detector.dump,  (bool)0, "Matrix and transormation goodness dump on file")
                ("detector.rotation",    &detector.rotation, Vector3f(0.,0.,0.), "Simulated chambers rotation")
                ("detector.translation",    &detector.translation, Vector3f(0.,0.,0.), "Simulated chambers traslation")
                ("detector.cloud1",&detector.cloud1,std::string("file1.txt"),"Start reference system point set txt file")
                ("detector.cloud2",&detector.cloud2,std::string("file2.txt"),"End reference system point set txt file")

                // INITIMAGE //
                ("image.lambda_air",  &image.lambda_air,  (float)0.07, "initial density for air")
                ("image.maskVoxels",    &image.maskVoxels,  (bool)0,"Mask voxels under maskThreshold")
                ("image.maskThreshold", &image.maskThreshold, (float)0.07,"Mask threshold")
                ("image.axis",   &image.axis,  (bool)1,"0= no slice axis, tolerance applied in XYZ")
                ("image.sliceAxis",&image.sliceAxis,  (int)2,"slice axis : 2=along Z, 1= along Y, 0= along X")
                ("image.reweight", &image.reweight,  (bool)1,"slice density: reweight boundary voxels")
                ("image.B_block",    &image.B_block,    HPoint3f(0,0,0),"begin point Block")
                ("image.E_block",    &image.E_block,    HPoint3f(0,0,0),"end point Block")

                // ANALYSIS OPTIONS //
                ("analysis.sets",     &analysis.sets,  (int)1,"number of vtk sets to analyse")
                ("analysis.nslices",  &analysis.nslices,  (int)4,"number of slices to analyse")
                ("analysis.sliceThickness",  &analysis.sliceThickness,  (float)10.,"thickness of slices [cm]")
                ("analysis.sliceTolerance",  &analysis.sliceTolerance,  HVector3f(0,0,0),"tolerance slices [cm]")
                ("analysis.sliceBegin",  &analysis.sliceBegin,  (float)-5.,"first slice begin position [cm]")
                ;
    }
} p;   // <-- INSTANCE //
} // namespace

////////////////////////////////////////////////////////////////////////////////
////////////// FUNCTIONS ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void setInitImage(IBVoxCollection &img, const char * vtkout){

    std::cout << "*** Init lambda values (*1.E-6)" << "\nair: " << p.image.lambda_air << std::endl;

    /// fill all voxels with air
    IBVoxel air{p.image.lambda_air*1.E-6,0,0};
    img.InitLambda(air);

    /// test
    img.ExportToVtk(vtkout);

    return;
}

////////////////////////////////////////////////////////////////////////////////
void createVoxelMask(IBVoxCollection &img, const char * vtkout, float density){

    std::cout << "Create voxel mask with density threshold " << density << std::endl;

    /// mask the image
    img = img.maskImage(density);

    /// testing.... export image
    img.ExportToVtk(vtkout);

    return;
}

////////////////////////////////////////////////////////////////////////////////
std::string getFileName(std::string name)
{
    if(name.find_last_of(".") != std::string::npos)
        name = name.substr(0,name.find_last_of("."));
    if(name.find_last_of("/") != std::string::npos)
        name = name.substr(name.find_last_of("/")+1,std::string::npos);
    if(name.find_last_of("*") != std::string::npos)
        name = name.substr(0,name.find_last_of("*"));

    //std::cout << "stripped file name = " << name << std::endl;
    return name;
}

////////////////////////////////////////////////////////////////////////////////
float getDensitySubImg(const IBVoxCollection &img, HPoint3f B, HPoint3f E) {
    float density = 0;

    for(int i=0; i< img.GetDims().prod(); ++i) {
        // voxel position in world reference system
        Vector4f Plocal = Vector4f(0,0,img.UnMap(i).cast<float>()(2),1);
        Vector4f P = img.GetWorldPoint(Plocal);
        float PB = (P(2)+img.GetSpacing()(2))-B(2);
        float PE = E(2)-(P(2)+img.GetSpacing()(2));

        float weight = 1.;
        if(PB>0 && PB<img.GetSpacing()(2)){
            weight = PB/img.GetSpacing()(2);
            //std::cout << "PB " << PB << ", weight B" << weight << std::endl;
        }
        if(PE>0 && PE<img.GetSpacing()(2)){
            weight = PE/img.GetSpacing()(2);
            //std::cout  << "PE " << PE << ", weight E" << weight << std::endl;
        }
        if(img.At(i).Value > 0)
            density += img.At(i).Value * img.GetSpacing().prod() * weight;
    }
    return density;
}

////////////////////////////////////////////////////////////////////////////////
void clipObject(const char *file, HPoint3f B, HPoint3f E, float thr, HPoint3f& B_obj, HPoint3f& E_obj, const char * vtkname){

    IBVoxCollection img;
    if( !img.ImportFromVtk(file) ){
        std::cout << "ATTENTION : error opening image from file..." << std::endl;
        return;
    }

    /// clip a first large image to define boundaries
    img = img.clipImage(B,E);

    /// clip with density threshold to refine
    img = img.clipImage(thr*1.E-6);

    img.ExportToVtk(vtkname);

    B_obj = img.GetPosition().homogeneous();
    Vector3f v = Vector3f(img.GetDims().cast<float>().cwiseProduct(img.GetSpacing()));
    HVector3f vp = HVector3f(v);
    E_obj = B_obj + vp;

    std::cout << "Object clip position: " << B_obj.transpose() << "\n";
    std::cout << "Object clip dimension : " << v.transpose() << "\n";

    return;
}

////////////////////////////////////////////////////////////////////////////////
float getTotDensityImg(const IBVoxCollection &img) {
    float density = 0;
    for(int i=0; i< img.GetDims().prod(); ++i) {
        if(img.At(i).Value > 0){
            std::cout << "Adding Vox " << i << ", lambda " << img.At(i).Value << std::endl;
            density += img.At(i).Value * img.GetSpacing().prod();
        }
    }
    return density;
}

////////////////////////////////////////////////////////////////////////////////
float objectDensity(const char *file_in, HPoint3f B, HPoint3f E) {

  std::cout << "calculating object density " << std::endl;

    /// clip object
    IBVoxCollection img;
    if( !img.ImportFromVtk(file_in) )
        return 0;

    // returns in 1/m units
    img = img.LambdaToInvLrad(3.);

    HVector3f slice_tolerance = HVector3f(0.,0.,0.);

    IBVoxCollection slice;
    HPoint3f p1 = HPoint3f(B - slice_tolerance);
    HPoint3f p2 = HPoint3f(E + slice_tolerance);
    //std::cout << "\nobjectDensity : Object B " << p1 << ", E " << p2 << std::endl;

    std::cout << "In bounds?" << std::endl;
    if(!(img.IsInsideBounds(p1)) ||!(img.IsInsideBounds(p2)) ){
        std::cout << "ATTENTION requested point outside voxel collection !!! Aborting..." << std::endl;
        return 0;
    }
    std::cout << "Is in bounds" << std::endl;
    std::cout << "Clipping image..." << std::endl;
    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    slice = img.clipImage(p1,p2);
    std::cout << "Clipped" << std::endl;
    // object volume
    float vol_slice = 20.*20.*10.;

    std::cout << "Getting total density" << std::endl;
    float density = getTotDensityImg(slice) / vol_slice * 1E6;

    ///testing....
    Vector3i Nvox = slice.GetDims();
    std::cout << "Object B " << B << ", E " << E << std::endl;
    std::cout << "\n Object tot density " << getTotDensityImg(slice)
              << " n vox " <<  Nvox.prod()
              << "\n Volume " << vol_slice
              << "\n density " << density << std::endl;
    std::cout << "Num voxels X " << Nvox(0) << ", Y " << Nvox(1) << ", Z "
              << Nvox(2) << ", total " << Nvox.prod() << std::endl;
//    char file[200];
//    std::string run = getFileName(file_in);
//    sprintf(file,"test_object_%s.vtk",run.c_str());
//    slice.ExportToVtk(file);

//    std::cout << "\nTotal object density:" << density << std::endl;
//    std::cout << "Object B " << B << ", E " << E << std::endl;
//    std::cout << "Object volume " << vol_slice << std::endl;

    return density;
}

////////////////////////////////////////////////////////////////////////////////
Vector<Vector2f> objectVarSliceDensities(const char *file_in, HPoint3f B, HPoint3f E, int nslices) {

    /// clip object
    IBVoxCollection img;
    if( !img.ImportFromVtk(file_in) )
        return 0;

    // returns in 1/m units
    img = img.LambdaToInvLrad(3.);

    /// init density vector
    Vector<Vector2f> densities;

    /// compare with total object density and save it in vector
    float totDensity = objectDensity(file_in, B, E);
    //std::cout << "\n\nTotal object density:" << totDensity << std::endl;
    //std::cout << "\n Object B " << B << ", E " << E << std::endl;

    /// define slices
    // view
    int ixz = p.image.sliceAxis;
    // slice
    float slice_thickness = p.analysis.sliceThickness;
    HVector3f slice_tolerance = p.analysis.sliceTolerance;
    if(p.image.axis)
      slice_tolerance(ixz) = 0;

    // begin end point with tolerance
    B(ixz) +=  p.analysis.sliceBegin;
    float slice0Start = B(ixz);
    E(ixz) = B(ixz) + slice_thickness;

    float position = B(ixz) + (0.5 * slice_thickness) - slice0Start;
    float vol_slice;

    IBVoxCollection slice;
    for(int i=0; i<nslices; ++i) {
        //std::cout << "\n slice " << i << " B " << B << ", E " << E << std::endl;
        HPoint3f p1 = HPoint3f(B - slice_tolerance);
        HPoint3f p2 = HPoint3f(E + slice_tolerance);
        //std::cout << "with tolerance: B " << p1 << ", E " << p2 << std::endl;

        if(!(img.IsInsideBounds(p1)) ||!(img.IsInsideBounds(p2)) ){
                std::cout << "ATTENTION requested point outside voxel collection !!! Aborting..." << std::endl;
                return densities;
        }
        slice = img.clipImage(p1,p2);
//        char file[200];
//        std::string run = getFileName(file_in);
//        sprintf(file,"test_slice_%i_%s.vtk",i,run.c_str());
//        slice.ExportToVtk(file);
       if(p.image.sliceAxis==1)
                vol_slice = slice_thickness * 20. * 20.;
            else
                vol_slice = slice_thickness * 20. * 10.;

        float density = getTotDensityImg(slice) / vol_slice * 1E6;
        densities.push_back(Vector2f(position,density));

        B(ixz) += slice_thickness;
        E(ixz) += slice_thickness;

        ///testing....
//        std::cout << "\n Slice " << i << " density " << getTotDensityImg(slice)
//                  << " volume " << vol_slice << " density " << density << std::endl;
//        Vector3i Nvox = slice.GetDims();
//        std::cout << "Num voxels X " << Nvox(0) << ", Y " << Nvox(1) << ", Z "
//                  << Nvox(2) << ", total " << Nvox.prod() << std::endl;

//        Vector3f Ncm = Vector3f(slice.GetDims().cast<float>().cwiseProduct(slice.GetSpacing()));
//        std::cout << "Num cm X " << Ncm(0) << ", Y " << Ncm(1) << ", Z "
//                  << Ncm(2) << ", total " << Ncm.prod() << std::endl;
    }

    densities.push_back(Vector2f(0,totDensity));

    return densities;
}

////////////////////////////////////////////////////////////////////////////////
Vector<Vector2f> objectSliceDensities(const char *file_in, HPoint3f B, HPoint3f E, int nslices) {

    /// this function computes density running over voxels, voxel dim=slice thickness
  std::cout << "Importing from vtk" << std::endl;

    /// clip object
    IBVoxCollection img;
    if( !img.ImportFromVtk(file_in) )
        return 0;

    std::cout << "Done" << std::endl;

    // returns in 1/m units
    img = img.LambdaToInvLrad(3.);

    /// init density vector
    Vector<Vector2f> densities;

    std::cout << "Generating region densities" << std::endl;
    /// compare with total object density and save it in vector
    float totDensity = objectDensity(file_in, B, E);
    float totPos = E(2)/2. + B(2)/2. + p.analysis.sliceBegin + 10.;

    /// 20150331 carota mockup, III region densities
    /// ATTENTION: don't define region on voxel boundaries, it will add voxel of boundary
    float tol = 2.5 + 0.1;
    float start = B(2) + p.analysis.sliceBegin + 12.5;
    // region 1
    HPoint3f Br1 = HPoint3f(B(0),B(1),start + tol);
    HPoint3f Er1 = HPoint3f(E(0),E(1),start + 20 - tol);
    float IregionDensity = objectDensity(file_in, Br1, Er1);
    // region 2
    HPoint3f Br2 = HPoint3f(B(0),B(1),start + 20 + tol);
    HPoint3f Er2 = HPoint3f(E(0),E(1),start + 30 - tol);
    float IIregionDensity = objectDensity(file_in, Br2, Er2);
    // region 3
    HPoint3f Br3 = HPoint3f(B(0),B(1),start + 30 + tol);
    HPoint3f Er3 = HPoint3f(E(0),E(1),start + 47.5 - tol);
    float IIIregionDensity = objectDensity(file_in, Br3, Er3);

    std::cout << "Done" << std::endl;

    /// define slices
    // view
    int ixz = p.image.sliceAxis;
    // slice
    float slice_thickness = p.analysis.sliceThickness;
    HVector3f slice_tolerance = p.analysis.sliceTolerance;
    if(p.image.axis)
      slice_tolerance(ixz) = 0;

    // compute carota center
    HPoint3f C = HPoint3f((B+E)*0.5);
    std::cout << "\n carota " << " B " << B << ", E " << E << ", C " << C << std::endl;

    // begin end point with tolerance
    B(ixz) +=  p.analysis.sliceBegin;
    float slice0Start = B(ixz);
    E(ixz) = B(ixz) + slice_thickness;

    IBVoxCollection slice;
    for(int i=0; i<nslices; ++i) {
        std::cout << "\n slice " << i << " B " << B << ", E " << E << std::endl;
        HPoint3f p1 = HPoint3f(B - slice_tolerance);
        HPoint3f p2 = HPoint3f(E + slice_tolerance);
        std::cout << "with tolerance: B " << p1 << ", E " << p2 << std::endl;

        if(!(img.IsInsideBounds(p1)) ||!(img.IsInsideBounds(p2)) ){
                std::cout << "ATTENTION requested point outside voxel collection !!! Aborting..." << std::endl;
                return densities;
        }
        // find voxels
	std::cout << "Clipping..." << std::endl;
        Vector3i vB = img.Find(p1);
        Vector3i vE = img.Find(p2);
        vE(ixz) = vB(ixz);
        slice = img.clipImage(vB,vE);
//        char file[200];
//        std::string run = getFileName(file_in);
//        sprintf(file,"test_slice_%i_%s.vtk",i,run.c_str());
//        slice.ExportToVtk(file);

	std::cout << "clipped" << std::endl;

        float position = B(ixz) + (0.5 * slice_thickness) - slice0Start;
        float vol_slice;

       vol_slice = slice_thickness * 20. * 10.;

        float density = getTotDensityImg(slice) / vol_slice * 1E6;
        densities.push_back(Vector2f(position,density));
        std::cout << "density = " << density << std::endl;

        B(ixz) += slice_thickness;
        E(ixz) += slice_thickness;

//        ///testing....
//        Vector3i Nvox = slice.GetDims();
//        std::cout << "\n Slice " << i << " n vox " <<  Nvox.prod()
//                  << " density " << getTotDensityImg(slice)
//                  << " volume " << vol_slice << " density " << density << std::endl;
//        std::cout << "Num voxels X " << Nvox(0) << ", Y " << Nvox(1) << ", Z "
//                  << Nvox(2) << ", total " << Nvox.prod() << std::endl;
//        Vector3f Ncm = Vector3f(slice.GetDims().cast<float>().cwiseProduct(slice.GetSpacing()));
//        std::cout << "Num cm X " << Ncm(0) << ", Y " << Ncm(1) << ", Z "
//                  << Ncm(2) << ", total " << Ncm.prod() << std::endl;
    }

    densities.push_back(Vector2f(Er1(2)/2. + Br1(2)/2. - slice0Start,IregionDensity));
    densities.push_back(Vector2f(Er2(2)/2. + Br2(2)/2. - slice0Start,IIregionDensity));
    densities.push_back(Vector2f(Er3(2)/2. + Br3(2)/2. - slice0Start,IIIregionDensity));
    densities.push_back(Vector2f(totPos - slice0Start,totDensity));
//    std::cout << "Region I "   << Er1(2)/2. + Br1(2)/2. - slice0Start << " density " << IregionDensity << std::endl;
//    std::cout << "Region II "  << Er2(2)/2. + Br2(2)/2. - slice0Start << " density " << IIregionDensity << std::endl;
//    std::cout << "Region III " << Er3(2)/2. + Br3(2)/2. - slice0Start << " density " << IIIregionDensity << std::endl;
//    std::cout << "Total sample " << totPos - slice0Start << " density " << totDensity << std::endl;

    return densities;
}

////////////////////////////////////////////////////////////////////////////////
void goClips(const char *file){

    /// clip block
    std::cout << "Block clip:" << std::endl;
    const char * nameout = "Block_clip.vtk";
    clipObject(file,HPoint3f(-120,-177,-90),HPoint3f(-70,-135,30),p.image.lambda_air,p.image.B_block,p.image.E_block,nameout);

    std::string filename = "init_image_" + getFileName(file) + ".config";
    std::ofstream of(filename);

    of << "[image]\n";
    of << "B_block = " << p.image.B_block.transpose() << "\n";
    of << "E_block = " << p.image.E_block.transpose() << "\n";

    return;
}

////////////////////////////////////////////////////////////////////////////////
static std::string systemPipe(std::string cmd) {
    FILE *pipe = popen(cmd.c_str(),"r");
    if(pipe) {
        std::stringstream ss;
        int c = fgetc (pipe);
        while (c != EOF) {
            if(c == '\n') ss << ' ';
            else ss << (char)c;
            c = fgetc (pipe);
        }
        return ss.str();
    }
    else return std::string("");
}

////////////////////////////////////////////////////////////////////////////////
TGraphErrors * createTGraph(const Vector<Vector2f> &v, const Vector<float> &e, int mcolor, int mstyle=7) {
    if(mcolor==10)
        mcolor+=1;

    /// check that the error are flat over density values
    //return verifyErrors(v,e,3);

    // 20150129 assume the same error in every point, computed as the mean of the relative errors
    float err_rel_mean = 0;
    for(int i=0 ; i<v.size(); ++i){
        float err_rel = e[i]/v[i](1);
        err_rel_mean += err_rel;
    }
    err_rel_mean /= v.size();


    TGraphErrors *gr = new TGraphErrors();
    for(int i=0 ; i<v.size(); ++i) {
        gr->SetPoint(i,v[i](0),v[i](1));
        //gr->SetPointError(i,0,e[i]);
        float e = err_rel_mean*v[i](1);
        gr->SetPointError(i,0,e);
    }
    gr->SetMarkerColor(mcolor);
    gr->SetLineColor(mcolor);
    gr->SetLineWidth(2.4);
    gr->SetMarkerStyle(7);
    gr->SetMarkerSize(1.7);
    //mstyle++;
    mcolor++;
    return gr;
}

////////////////////////////////////////////////////////////////////////////////
void objectSliceAnalisis(int argc, char **argv, TGraph *&gr_mean){
    /// open output file
    std::ofstream of;
    of.open("MeanRms.txt", std::ofstream::app);

    /// loop on files sets and initialize variables
    std::vector<std::string> files;
    for(int is=0; is<p.analysis.sets; is++){
      std::string cmdout = systemPipe(std::string("ls ") + std::string(argv[1+is]));
      //std::cout << "cmd: " << cmdout << "\n";
      std::istringstream ss(cmdout);

      for (std::string each; std::getline(ss,each,' ');
           files.push_back(each));

      of << std::string(argv[1+is]);
      for(int isl=0; isl<p.analysis.nslices; isl++)
          of << ";";
      of << "\n";
      for(int isl=0; isl<p.analysis.nslices; isl++)
          of << "; slice "<< isl+1;
      of << "\n";
    }
    std::cout << "There are " << files.size() << " data samples" << std::endl;


      /// density MEAN
      Vector<Vector2f> densities_mean(p.analysis.nslices+4);
      Vector<float> densities_rms(p.analysis.nslices+4);
      for(int iv=0; iv<p.analysis.nslices+4; iv++){
          densities_mean[iv]=Vector2f(0,0);
          densities_rms[iv]=0;
      }

      /// define begin end point of object
      HPoint3f B = p.image.B_block;
      HPoint3f E = p.image.E_block;

      std::cout << "Begin = " << B << std::endl;
      std::cout << "End = " << E << std::endl;

      int c=0;
      foreach (std::string &str, files) {
          if(c<100){
              //if(c==0)
              std::cout << "processing file: " << str << "\n";
              Vector<Vector2f> densities = objectSliceDensities(str.c_str(),B,E,p.analysis.nslices);
              for(int iv=0; iv<p.analysis.nslices+4; iv++){
                  densities_mean[iv](0) =  densities[iv](0);
                  densities_mean[iv](1) +=  densities[iv](1);
              }
              //mg_set->Add(CreateTGraph(densities,densities_rms,c+2));
              c++;
          }
      }
      for(int iv=0; iv<p.analysis.nslices+4; iv++)
        densities_mean[iv](1) = densities_mean[iv](1) / c;

      // density RMS
      c=0;
      foreach (std::string &str, files) {
          if(c<100){
              Vector<Vector2f> densities = objectSliceDensities(str.c_str(),B,E,p.analysis.nslices);
              for(int iv=0; iv<densities.size(); iv++)
                  densities_rms[iv] += pow(densities[iv](1)-densities_mean[iv](1),2.);
              c++;
          }
      }
      for(int iv=0; iv<p.analysis.nslices+4; iv++)
        densities_rms[iv] = sqrt(densities_rms[iv]/(c-1)) / sqrt(c); //standard deviation from mean

      gr_mean = createTGraph(densities_mean,densities_rms,2);

      std::cout << "Mean " << densities_mean;
      std::cout << "Rms " << densities_rms << "\n";
      of << "mean";
      for(int iv=0; iv<p.analysis.nslices+4; iv++)
        of << ";" << densities_mean[iv](1);
      of << "\n";
      of << "rms";
      for(int iv=0; iv<p.analysis.nslices+4; iv++)
        of << ";" << densities_rms[iv];
      of << "\n\n";

//        if(mg_set) {
//          //TFile *file = new TFile(argv[2],"RECREATE");
//          TCanvas *c = new TCanvas("c","c",1200,600);
//          c->Divide(2,1);
//          c->cd(1);
//          mg_set->Draw("alp");
//          mg_set->GetXaxis()->SetTitle("slices");
//          mg_set->GetYaxis()->SetTitle("#lambda density in slice");
//          mg_set->GetYaxis()->SetTitleOffset(1.5);
//          mg_set->SetTitle("Slices Densities");

//          c->cd(2);
//          gr_mean->Draw("alp");
//          gr_mean->GetXaxis()->SetTitle("slices");
//          gr_mean->GetYaxis()->SetTitle("#lambda mean in slice");
//          gr_mean->GetYaxis()->SetTitleOffset(1.5);
//          gr_mean->SetTitle("Slices Mean Densities");

//          TPaveText *pt = new TPaveText(.7,.7,.95,.9,"brNDC");
//          pt->SetFillColor(kWhite);
//          pt->SetTextAlign(12);
//          pt->SetBorderSize(1);
//          pt->AddText("#color[2]{Run 2161}");
//          pt->AddText("Vox size 5 cm");
//          pt->AddText("Slices 10 cm");
//          pt->AddText("DAQ Time 60 min");
//          pt->Draw();

//          c->Draw();
//          c->Write();
//          c->SaveAs("mg_2161_img_vox5_60min.png");
//          //file->Close();
//      }
//  }
    of.close();
}

//////////////////////////////////////////////////////////////////////////////////
bool find_analyzer(IBAnalyzer *an) {
    return p.analyzers.find(an->type_name(),0,strlen(an->type_name())) != std::string::npos;
}

//////////////////////////////////////////////////////////////////////////////////
HPoint3f voxelToCoord(Vector3i vox, Vector3f spacing, HVector3f B_collection){

    Vector3f v = Vector3f(vox.cast<float>().cwiseProduct(spacing));
    HVector3f vp = HVector3f(v);
    HPoint3f B_vox = HPoint3f(B_collection + vp);

    std::cout << "Voxel " << vox << "-------> " << B_vox << std::endl;

    return B_vox;
}

////////////////////////////////////////////////////////////////////////////////
pcl::PointCloud<pcl::PointXYZ> csvToPCD(std::string filename){

    std::cout << "\nReading file  " << filename << std::endl;

    // reference system compatibility from Paolo's to IB
    // scale factor form m->cm, origin on first voxel B_vertex
    float scale = 100.;
    pcl::PointXYZ origin;
    origin.x = -113.;
    origin.y = -183.;
    origin.z = -101.;

    std::ifstream file(filename);
    std::string line;
    std::string col;

    // header
    std::getline(file, line, ',');
    int N = atoi(line.c_str());
    //std::cout << "N=" << N << std::endl;

    // spacing
    float voxSize = 0.;
    HVector3f B_collection;
    if(p.detector.voxelReferenceSystem){
        std::getline(file, line, ',');
        voxSize = atof(line.c_str());

        std::getline(file, line, ',');
        float Bx = atof(line.c_str());
        std::getline(file, line, ',');
        float By = atof(line.c_str());
        std::getline(file, line, ',');
        float Bz = atof(line.c_str());
        B_collection = HVector3f(Bx,By,Bz);
    }

    //Fill Point Cloud Data
    pcl::PointCloud<pcl::PointXYZ> cloud;
    cloud.width = N;
    cloud.height = 1;
    cloud.is_dense = false;
    cloud.points.resize(cloud.width * cloud.height);

    // data set
    std::getline(file, line);
    int i=0;

    while ( std::getline(file, line) && i<N) {

        std::istringstream csvStream(line);

        std::getline(csvStream, col, ',');
        float X = atof(col.c_str());
        std::getline(csvStream, col, ',');
        float Y = atof(col.c_str());
        std::getline(csvStream, col, ',');
        float Z = atof(col.c_str());

        if(p.detector.voxelReferenceSystem){
            Vector3f spacing(voxSize, voxSize, voxSize);

            HPoint3f p = voxelToCoord(Vector3i(X,Y,Z),spacing,B_collection);
            cloud.points[i].x = p[0];
            cloud.points[i].y = p[1];
            cloud.points[i].z = p[2];
        } else {
            cloud.points[i].x = X * scale + origin.x;
            cloud.points[i].y = Y * scale + origin.y;
            cloud.points[i].z = Z * scale + origin.z;
        }

        //std::cout << "Point " << i << "\t x=" << cloud.points[i].x << ", y=" << cloud.points[i].y << ", z=" << cloud.points[i].z << std::endl;
        i++;
    }

//    filename += ".pcd";
//    pcl::io::savePCDFileASCII(filename,cloud);

    return cloud;
}

////////////////////////////////////////////////////////////////////////////////
int doIterations(const char *file_in,
                  const char *file_out,
                  float min,
                  float start,
                  float vox_size) {

    IB::Version::PrintSelf(std::cout);

    ////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// READER /////
    TFile* f = new TFile (file_in);
    IBMuonEventTTreeReader *reader = IBMuonEventTTreeReader::New(f);
    reader->setTFile(f);
    reader->setMomentum(p.momentum);
    reader->setAcquisitionTime(min);
    reader->setStartTime(start);
    std::cout << "\n--------- READER : open file with " << reader->getNumberOfEvents() << " events!\n";

    /// error2
    IBMuonError sigma(6.02, 7.07);
    sigma.crossChamberErrorCorrection(true);
    reader->setError(sigma);

    /// chambers alignment
    {
        IBMuonEventTTreeLNLdataReader *rd = ((IBMuonEventTTreeLNLdataReader *)reader);
        Eigen::Affine3f tr(Matrix4f::Identity());
//        rd->setAlignment(tr.matrix());
        rd->setAlignmentFromData(60);
        std::cout << "READER ALIGNMENT: \n" << rd->getAlignment() << "\n\n";
    }

    /////////////////////////////////// ACQUISITION /////
    std::cout << "\n--------- READER : loading muons" << std::endl;
    IBMuonCollection muons;
    int tot=0;
    int ev = reader->getNumberOfEvents();

    for( int i=0; i<80; ++i) std::cout << "-"; std::cout << "\r";
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            muons.AddMuon(mu);
            tot++;
        }
        if(tot++%(ev/80) == 0) std::cout << "o" << std::flush;
    }
    std::cout << "\n";


    /////////////////////////////////// MUONS SELF ALIGNMENT /////
    muons.PerformMuonSelfAlignment();

    std::cout << "\n Initial muon collection : " << std::endl;
    muons.PrintSelf(std::cout);

    /////////////////////////////////// DETECTOR  ROTO-TRANSLATION /////
    if(p.detector.transform){
        std::cout << "\n--------- DETECTOR  ROTO-TRANSLATION: Muons coordinate system transformation ---------\n ";

        /// load detector file
        TFile* f = new TFile (p.detector.datafile.c_str());
        reader->setTFile(f);
        reader->setStartTime(p.detector.start);
        reader->setAcquisitionTime(p.detector.min);

        /// acquisition
        IBMuonCollection rototrans_muons;
        int tot=0;
        int ev = reader->getNumberOfEvents();

        for( int i=0; i<80; ++i) std::cout << "-"; std::cout << "\r";
        for (int i=0; i<ev; i++) {
            MuonScatter mu;
            if(reader->readNext(&mu)) {
                rototrans_muons.AddMuon(mu);
                tot++;
            }
            if(tot++%(ev/80) == 0) std::cout << "o" << std::flush;
        }
        std::cout << "\n";

        /// self alignment
        rototrans_muons.PerformMuonSelfAlignment();

        std::cout << "\n To be rotated muon collection : " << std::endl;
        rototrans_muons.PrintSelf(std::cout);


        if(p.detector.findMatrix){
             // load point sets
            pcl::PointCloud<pcl::PointXYZ> cloud1 = csvToPCD(p.detector.cloud1);
            pcl::PointCloud<pcl::PointXYZ> cloud2 = csvToPCD(p.detector.cloud2);

            // trasformation estimation: find the rigid transformation matrix that can be used to transform one to other
            pcl::registration::TransformationEstimationSVD<pcl::PointXYZ,pcl::PointXYZ> TESVD;
            pcl::registration::TransformationEstimationSVD<pcl::PointXYZ,pcl::PointXYZ>::Matrix4 transformation;
            TESVD.estimateRigidTransformation (cloud1,cloud2,transformation);

            std::cout << "\n\nThe Estimated Rotation-Translation matrix is : \n" << std::endl;
            printf ("\n");
            printf ("    | %6.3f %6.3f %6.3f | \n", transformation (0,0), transformation (0,1), transformation (0,2));
            printf ("R = | %6.3f %6.3f %6.3f | \n", transformation (1,0), transformation (1,1), transformation (1,2));
            printf ("    | %6.3f %6.3f %6.3f | \n", transformation (2,0), transformation (2,1), transformation (2,2));
            printf ("\n");
            printf ("t = < %0.3f, %0.3f, %0.3f >\n", transformation (0,3), transformation (1,3), transformation (2,3));

            // rotate muon collection reference system
            Eigen::Matrix4f matrix = Eigen::Matrix4f::Identity();
            for(int ir=0; ir<3; ir++){
                matrix(ir,3) = transformation (ir,3);
                for(int ic=0; ic<3; ic++)
                    matrix(ir,ic) =  transformation(ir,ic);
            }

            /// Evaluate transformation "goodness" via Sum(T*point - point) / Npoints
            // Executing the transformation from cloud 1
            pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud (new pcl::PointCloud<pcl::PointXYZ> ());
            pcl::transformPointCloud (cloud1, *transformed_cloud, matrix);

            float sumDiscrep = 0.;
            for(int ip=0; ip< cloud2.size(); ip++){
                float Dx =  transformed_cloud->points[ip].x-cloud2.points[ip].x;
                float Dy = transformed_cloud->points[ip].y-cloud2.points[ip].y;
                float Dz =  transformed_cloud->points[ip].z-cloud2.points[ip].z;

                Vector3f discrep(Dx,Dy,Dz);
                float norm = discrep.norm();
                sumDiscrep += norm;
                //std::cout << "Point " << ip << "\t Dx=" << Dx << ", Dy=" << Dy << ", Dz=" << Dz << "----> " << norm << " - Sum norm " << sumDiscrep << std::endl;
            }
            sumDiscrep /= cloud2.size();
            std::cout << "\n -----> Mean |T*point - point| = " << sumDiscrep << " [cm]" << std::endl;

            if(p.detector.dump){
                char file[200];
                sprintf(file,"%s_transformation_analysis.txt",p.file_out);
                std::ofstream ofs;
                ofs.open (file, std::ofstream::out | std::ofstream::app);

                ofs << "\n--------- DETECTOR  ROTO-TRANSLATION: Muons coordinate system transformation --------- ";
                ofs << "\nComputing transformation " << "\n  * from points in file " << p.detector.cloud1 << "\n  * to points in file " << p.detector.cloud2 << std::endl;
                ofs <<  "\nThe Estimated Rotation-Translation matrix is : \n" << matrix << std::endl;
                ofs <<  "\n -----> Mean |T*point - point| = " << sumDiscrep << " [cm]" << std::endl;

                ofs.close();

                return 0;
            }

            /// transform muons
            rototrans_muons.dataRotoTranslation(matrix);
        } else {
            std::cout << "Rotation : " << p.detector.rotation.transpose() << "\n";
            std::cout << "Translation : " << p.detector.translation.transpose() << "\n\n";

            rototrans_muons.dataRotoTranslation(p.detector.rotation,p.detector.translation);
        }

        std::cout << "\n Rotated muon collection : " << std::endl;
        rototrans_muons.PrintSelf(std::cout);

        /// add muons to previous collection
        muons.AddCollection(rototrans_muons);
    }

    std::cout << "\n Total muon collection : " << std::endl;
    muons.PrintSelf(std::cout);

    /////////////////////////////////// MUON COLLECTION DUMP /////
    if(p.dump){
        char rootfile[100];
        sprintf(rootfile, "%s.root",p.file_out);
        //muons.DumpTTree(rootfile);
        muons.DumpSimpleTree(rootfile);
    }


    ////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// VOX COLLECTION /////
//    // all volume between chambers
//    Vector3f vox_bounding( 300, 183, 240  ); // centered bounding size //
//    Vector3f vox_pos(     -150,-183,-120  );
//    // only selected volume for carota analisys
//    Vector3f vox_bounding( 300, 63,  240   );
//    Vector3f vox_pos(     -150, -183,-120  );
    // volume for ispra block
    Vector3f vox_bounding( 226, 183, 202  ); // centered bounding size //
    Vector3f vox_pos(     -113,-183,-101  );

    IBVoxCollection voxels(Vector3i(vox_bounding(0)/vox_size,
                                    vox_bounding(1)/vox_size,
                                    vox_bounding(2)/vox_size));

    IBVoxel zero{0,0,0};
    IBVoxel air{p.image.lambda_air*1.E-6,0,0};
    voxels.InitLambda(air);
    voxels.SetSpacing (Vector3f(vox_size,
                                vox_size,
                                vox_size));
    voxels.SetPosition(vox_pos);

    HPoint3f B_vox = HPoint3f(voxels.GetPosition().homogeneous());
    Vector3f v = Vector3f(voxels.GetDims().cast<float>().cwiseProduct(voxels.GetSpacing()));
    HVector3f vp = HVector3f(v);
    HPoint3f E_vox = HPoint3f(B_vox + vp);

    std::cout << "Voxel collection boundaries " << B_vox << ", " << E_vox << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// ALGORITHMS /////
    // poca //
    IBPocaEvaluator* processor =
            IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);
    IBNormalPlaneMinimizationVariablesEvaluator* minimizator =
            new IBNormalPlaneMinimizationVariablesEvaluator();
    minimizator->SetAlphaXZ(0); // X
    minimizator->setRaytracer(tracer);

    // ML Algorithm //
    IBAnalyzerEMAlgorithmSGA *ml_algorithm = new IBAnalyzerEMAlgorithmSGA_PX;
    //IBAnalyzerEMAlgorithmSGA *ml_algorithm = new IBAnalyzerEMAlgorithmSGA_P;

    // analyzer //
    IBAnalyzerEM* aem = new IBAnalyzerEM(voxels);
    aem->SetMLAlgorithm(ml_algorithm);
    aem->SetPocaAlgorithm(processor);
    aem->SetRayAlgorithm(tracer);
    aem->SetVarAlgorithm(minimizator);

    IBAnalyzerPoca *anpc = new IBAnalyzerPoca;
    anpc->SetVoxCollection(&voxels);
    anpc->SetPocaAlgorithm(processor);

    IBAnalyzerTrackCount *antrk = new IBAnalyzerTrackCount();
    antrk->SetVoxCollection(&voxels);
    antrk->SetPocaAlgorithm(processor);
    antrk->SetRayAlgorithm(tracer);


    /////////////////////////////////// TRACK ANALYZER /////
    if(find_analyzer(antrk))
    {
        std::cout << "\n--------- TRACK ANALYZER....\n";
        antrk->SetMuonCollection(&muons);
        voxels.InitLambda(zero);
        antrk->Run(1,1);
        char file[200];
        sprintf(file,"%s_%s.vtk",p.file_out,antrk->type_name());
        voxels.ExportToVtk(file,1);

        delete antrk;
    }

    /////////////////////////////////// POCA ANALYZER /////
    if(find_analyzer(anpc))
    {
        std::cout << "\n--------- POCA ANALYZER....\n";
        anpc->SetMuonCollection(&muons);
        voxels.InitLambda(zero);
        anpc->Run(1,1);
        char file[200];
        sprintf(file,"%s_%s.vtk",p.file_out,anpc->type_name());
        voxels.ExportToVtk(file,1);

        delete anpc;
    }

    /////////////////////////////////// EM ANALYZER /////
    if(find_analyzer(aem)){
        aem->SetMuonCollection(&muons);

        /////////////////////////////////// DUMP /////
        if(p.dump){
            char rootfile[100];
            sprintf(rootfile, "%s_aem.root",p.file_out);
            aem->dumpEventsTTree(rootfile);
        }

        /////////////////////////////////// SIJ CUTS /////
        std::cout << "\n--------- CUTS...." << std::endl;
        if(p.sijcut > 0.){
            std::cout << "Sij cut =  " << p.sijcut << std::endl;
            int count_precut = aem->Size();
            std::cout << "Sij cut -> muons precut: " << count_precut << std::endl;
            voxels.InitLambda(air);
            aem->SijCut(p.sijcut);
            int count_postcut = aem->Size();
            std::cout << " ---> cut: " << count_precut - count_postcut << "  (" << float(count_precut - count_postcut)/float(count_precut)*100.  << "%) \n";
        }

        /////////////////////////////////// VOXEL MASK /////
        if(p.image.maskVoxels){
            std::cout << "\n--------- VOXEL MASK...." << std::endl;
            voxels.InitLambda(air);

            /// flag voxels with density under threshold with -1*value
            char vtkoutfile[100];
            sprintf(vtkoutfile, "test_init_image_mask%f_%s.vtk",p.image.maskThreshold,getFileName(p.file_out).c_str());
            std::cout << "Create voxel mask with density threshold " << p.image.maskThreshold << std::endl;

            /// mask the voxels between the bottom layer and the upper chmaber
            HPoint3f B_mask_up = HPoint3f(-113,-60,-101);
            HPoint3f E_mask_up = HPoint3f(113,-3,101);
            voxels = voxels.maskImage(B_mask_up,E_mask_up,-0.07);

            /// mask the voxels between the lower chamber and the reference blocks
            HPoint3f B_mask_dw = HPoint3f(-113,-180.5,-101);
            HPoint3f E_mask_dw = HPoint3f(113,-163,101);
            voxels = voxels.maskImage(B_mask_dw,E_mask_dw,-0.07);

            /// testing.... export image
            voxels.ExportToVtk(vtkoutfile);

            aem->filterEventsVoxelMask();
        }

        /////////////////////////////////// EM ITERATIONS /////
        std::cout << "\n--------- ITERATIONS...." << std::endl;
        voxels.InitLambda(air);
        char file[100];
        for (int i=p.iterations(0); i<=p.iterations(2)/p.iterations(1); ++i) {
          std::cout << "Running iteration " << i << std::endl;
            aem->Run(p.iterations(1),1);
                int n = i*p.iterations(1);
                sprintf(file, "%s_%05i.vtk",file_out,n);
                voxels.ExportToVtk(file,0);
        }
        delete aem;
    }

    delete minimizator;
    delete processor;
    delete tracer;
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////// MAIN /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

    //usage: .exe file_in file_out 0 0 voxsize [options]";
    std::string config_file("cmt.config");
    std::string init_image_file("NULL");

    p.add_options()
            ("config",&config_file,"set config file")
            ("init",&init_image_file,"image init file")
            ;

    p.parse_command_line(argc,argv);
    p.parse_config_file(config_file);
    if(init_image_file!=std::string("NULL"))
      p.parse_config_file(init_image_file);

    p.parse_command_line(argc,argv);

    /// open file
    const char * file = (char *)argv[1];
    p.file_in  = (char *)argv[1];
    p.file_out = (char *)argv[2];

    ////////////////////////////////////////////////////////////////////////////////
    if(p.execute == "reconstruction") {
        std::cout << "Executing reconstruction on file " << file << "..." << std::endl;

        ///image reconstruction algorithm
        p.minutes  = atof(argv[3]);
        p.start    = atof(argv[4]);
        p.voxsize  = atof(argv[5]);

        doIterations(p.file_in,
                      p.file_out,
                      p.minutes,
                      p.start,
                      p.voxsize
                      );
    }
    ////////////////////////////////////////////////////////////////////////////////
    else if(p.execute == "clip") {

        std::cout << "Executing clips.... on file " << file << "..." << std::endl;

        ///clip image
        goClips(file);

        return 1;

    }
    ////////////////////////////////////////////////////////////////////////////////
    else if(p.execute == "mask") {

        std::cout << "Executing mask voxels on file " << p.file_in << "..." << std::endl;

        ///set init image testing
        IBVoxCollection img;
        if( !img.ImportFromVtk(p.file_in) )
            return 0;

        ///mask image
        createVoxelMask(img,p.file_out,p.image.maskThreshold);

        return 1;

    }
    ////////////////////////////////////////////////////////////////////////////////
    else if(p.execute == "test_init_image") {

        std::cout << "Testing init image on file " << file << "..." << std::endl;

        ///set init image testing
        IBVoxCollection img;
        if( !img.ImportFromVtk(file) )
            return 0;

        char outfile[100];
        sprintf(outfile, "test_init_image_%s.vtk",getFileName(p.file_out).c_str());
        setInitImage(img,outfile);

        return 1;
    }
    ////////////////////////////////////////////////////////////////////////////////
    else if(p.execute == "block_analysis") {
        std::cout << "Executing analysis on file " << file << "..." << std::endl;

        TApplication a("a", 0, 0);

        ///Slice analisis
        char gr_name[100];
        std::string run = getFileName(file);

        run = run.substr(0,4);
        sprintf(gr_name,"mg_%s",run.c_str());
        //sprintf(gr_name,"mg",run.c_str());
        TGraph * gr_mean = new TGraph(gr_name,"Density mean for different samples");

        std::cout << " Doing slice analysis " << std::endl;
        objectSliceAnalisis(argc,argv,gr_mean);
        std::cout << " Done slice analysis " << std::endl;

        gr_mean->SetMinimum(0);
        gr_mean->SetMaximum(14.);
        gr_mean->GetXaxis()->SetTitle("Block slices");
        gr_mean->GetYaxis()->SetTitle("#lambda mean density in block");
        gr_mean->GetYaxis()->SetTitleOffset(1.5);
        gr_mean->SetTitle("Density mean for different samples");
        gr_mean->SetName(gr_name);

        /// draw canvas
//        TCanvas *cm = new TCanvas("cm","cm",900,600);
//        cm->cd(1);
//        //gr_mean->Draw("AB1");
//        gr_mean->Draw("ALP");

//        TPaveText *pt = new TPaveText(.7,.7,.95,.9,"brNDC");
//        pt->SetFillColor(kWhite);
//        pt->SetTextAlign(12);
//        pt->SetBorderSize(1);
//        pt->AddText("#bf{Run 2161: slice analysis}");
//        pt->AddText("#bf{Image reconstruction}");
//        //pt->AddText("#bf{Voxel reweight}");
//        pt->AddText("#bf{NO Sij guess}");
//        pt->AddText("#bf{Sij cut N=30}");
//        pt->AddText("#bf{4000 iterations}");
//        pt->AddText("#bf{Init image: Fe=10, carota=2, air=0.07}");
//        //            pt->AddText("#color[1]{vox 2.5 cm, 60 min}");
//        //            pt->AddText("#color[2]{vox 2.5 cm, 120 min}");
//        //            pt->AddText("#color[3]{vox 2.5 cm, 240 min}");
//        pt->AddText("#color[2]{vox 2.5 cm, 480 min}");
//        //            pt->AddText("#color[4]{vox 5 cm, 60 min}");
//        //            pt->AddText("#color[5]{vox 5 cm, 120 min}");
//        //            pt->AddText("#color[6]{vox 5 cm, 240 min}");
//        pt->AddText("#color[3]{vox 5 cm, 480 min}");
//        pt->Draw();
//        a.Run(kTRUE);

        /// dump on file
       std::cout << "Writing output..." << std::endl;

        TFile *fout = new TFile("BlockOutput.root","UPDATE");
        fout->cd();
        gr_mean->Write();
        fout->Close();

        return 1;
    }

    return 0;
}






