#include <TMultiGraph.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>

#include "IBROC.h"



using namespace uLib;

namespace ROCPlot {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ROC ROOT GRAPH
enum EnumROCStyle {
    StyleROC = 0,
    StyleFAMP
};

class ROCGraph : public IBROC {
    typedef IBROC BaseClass;
public:

    ROCGraph(EnumROCStyle style = StyleFAMP) :
        m_style(style)
    {
        m_awo.SetLineWidth(2);
        m_awo.SetMarkerSize(0);
        m_owa.SetLineWidth(3);
        m_owa.SetMarkerSize(0);
    }


    virtual void Update() {
        this->Draw();
    }

    void SetStyle(EnumROCStyle style) {
        this->m_style = style;
        this->Update();
    }

    uLibRefMacro(awo,TGraph)
    uLibRefMacro(owa,TGraph)

protected:
    int Draw() {
        m_awo.Set(this->size());
        m_owa.Set(this->size());

        m_awo.Clear();
        m_owa.Clear();

        switch(m_style) {
        case ROCPlot::StyleROC:
            m_awo.SetLineColor(kRed);
            m_owa.SetLineColor(kRed);
            for(int i=0; i<this->size(); i++) {
                m_awo.SetPoint(i,this->at(i).Awo(),100 - this->at(i).Owa());
                m_owa.SetPoint(i,100 - this->at(i).Owa(),this->at(i).Awo());
            }
            break;
        case ROCPlot::StyleFAMP:
            m_awo.SetLineColor(kRed);
            m_owa.SetLineColor(kBlue);
            for(int i=0; i<this->size(); i++) {
                m_awo.SetPoint(i,this->at(i).X(),this->at(i).Awo());
                m_owa.SetPoint(i,this->at(i).X(),this->at(i).Owa());
            }
            break;
        }
        return this->size();
    }



private:
    // members //
    TGraph m_awo, m_owa;
    EnumROCStyle m_style;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ROC ROOT MULTI GRAPH

class ROCMultiGraph : public TMultiGraph {

    typedef TMultiGraph BaseClass;
public:
    ROCMultiGraph(EnumROCStyle style = StyleFAMP) :
        m_min(0), m_max(100), m_style(style)
    {
        m_legend=new TLegend(0.685,0.12,0.88,0.22);
        m_legend->SetTextFont(72);
        m_legend->SetTextSize(0.03);
        m_legend->SetFillColor(0);
        m_legend->SetCornerRadius(2);
        m_legend->SetMargin(0.3);
        m_legend->SetBorderSize(0);
    }

    ~ROCMultiGraph() {
        delete m_legend;
    }

    void AddROC( ROCGraph rocg, int line_style = 1 ) {
        rocg.awo().SetLineStyle(line_style);
        rocg.owa().SetLineStyle(line_style);
        this->m_rocs.push_back(rocg);
    }

    inline ROCGraph & operator[] (int i) {
        return m_rocs[i];
    }

    void SetRange(float min,float max) {
        m_min = min;
        m_max = max;
        this->Draw();
    }

    void SetStyle(EnumROCStyle style) {
        m_style = style;
        if(m_rocs.size() > 0)
            this->Draw();
    }

    void Draw()  {
        m_legend->Clear();
        if(m_rocs.size() == 1 && m_style == StyleFAMP)
        {
            ROCGraph & roc = m_rocs.at(0);
            m_legend->AddEntry(&roc.awo(),roc.m_labels[1].c_str(),"l");
            m_legend->AddEntry(&roc.owa(),roc.m_labels[2].c_str(),"l");
        }
        else
        {
            for (int i=0; i< m_rocs.size(); ++i) {
                m_legend->AddEntry(&m_rocs.at(i).awo(), m_rocs.at(i).m_labels[0].c_str(), "l");
            }
        }

        // add graphs //
        this->Clear();

        if(m_style == StyleFAMP) {
            for(int i=0; i< m_rocs.size(); ++i)
            {
                m_rocs.at(i).SetStyle(m_style);
                this->Add(&m_rocs.at(i).awo());
                this->Add(&m_rocs.at(i).owa());
            }
            BaseClass::Draw("apl");
            this->GetXaxis()->SetLimits(m_min,m_max);
            //            this->GetYaxis()->SetRangeUser(0.,100.);
            this->GetXaxis()->SetTitle("Normalized density");
            this->GetYaxis()->SetTitle("Rate (%)");
            m_legend->Draw();
        }
        else {
            for(int i=0; i< m_rocs.size(); ++i)
            {
                m_rocs.at(i).SetStyle(m_style);
                this->Add(&m_rocs.at(i).awo());
            }
            BaseClass::Draw("apl");
            this->GetXaxis()->SetRangeUser(0,100);
            this->GetYaxis()->SetRangeUser(0,100);
            this->GetXaxis()->SetTitle("False alarms (%)");
            this->GetYaxis()->SetTitle("Detection efficiency (%)");
            if(m_rocs.size() > 1)
                m_legend->Draw();
        }
        this->GetXaxis()->SetTitleOffset(1.3);
        this->GetYaxis()->SetTitleOffset(1.6);
        this->GetXaxis()->SetTitleFont(42);
        this->GetYaxis()->SetTitleFont(42);
        this->GetXaxis()->SetTitleSize(0.03);
        this->GetYaxis()->SetTitleSize(0.03);
        this->GetXaxis()->SetLabelFont(42);
        this->GetYaxis()->SetLabelFont(42);
        this->GetXaxis()->SetLabelSize(0.03);
        this->GetYaxis()->SetLabelSize(0.03);

    }

private:
    TLegend * m_legend;
    Vector<ROCGraph> m_rocs;
    float m_min,m_max;
    EnumROCStyle m_style;
};



} // ROCPlot










////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// MAIN

main(int argc, char * argv[]){

    // .................................................................... //
    struct Parameters {
        char * file0;
        char * file1;
        char * file_out;
        int type;
        char * file0_caption;
        char * file1_caption;
        float x_min;
        float x_max;
    } p ={
        NULL,
        NULL,
        (char *)"rocplot.eps",
        1,
        NULL,
        NULL,
        NAN,
        NAN
    };

    if(argc>=2) {
        p.file0=argv[1];
        std::cout << "file1 = " << p.file0 << "\n";
    }
    if(argc>=3 && std::string(argv[2]) != "NULL" ) {
        p.file1=argv[2];
        std::cout << "file2 = " << p.file1 << "\n";
    }
    if(argc>=4) {
        p.file_out=argv[3];
    }
    if(argc>=5) {
        p.type = atoi(argv[4]);
    }
    if(argc>=7) {
        p.file0_caption = argv[5];
        p.file1_caption = argv[6];
    }
    if(argc>=9){
        p.x_min = atof(argv[7]);
        p.x_max = atof(argv[8]);
    }
    // .................................................................... //

    gStyle->SetPaperSize(30,30);  //default

    ROCPlot::ROCMultiGraph mg;
    mg.SetStyle( (ROCPlot::EnumROCStyle)p.type  );

    ROCPlot::ROCGraph roc1;
    roc1.read_roc_with_header(p.file0);
    if(p.file0_caption)
        roc1.m_labels[0] = std::string(p.file0_caption);
    if(p.file1)
        mg.AddROC(roc1,7);
    else
       mg.AddROC(roc1);

    float range_max = 0;
    float range_min = 0;
    if(isnan(p.x_max))
        range_max = roc1.GetRange()(1);
    if(isnan(p.x_min))
        range_min = roc1.GetRange()(0);



    if(p.file1)
    {
        ROCPlot::ROCGraph roc2;
        roc2.read_roc_with_header(p.file1);
        if(p.file1_caption)
            roc2.m_labels[0] = std::string(p.file1_caption);
        mg.AddROC(roc2);
        if(isnan(p.x_max))
            range_max = range_max > roc2.GetRange()(1) ? range_max : roc2.GetRange()(1);
        if(isnan(p.x_min))
            range_min = range_min < roc2.GetRange()(1) ? range_min : roc2.GetRange()(1);
    }

    if(isnan(p.x_max))
        p.x_max = range_max + range_max/5;
    if(isnan(p.x_min))
        p.x_min = range_min - range_min/5;


    TCanvas * c0 = new TCanvas("c0","multigraph L3",200,10,700,500);
    c0->SetFillColor(0);
    c0->SetGrid();
    c0->SetBorderMode(0);

    mg.SetRange(p.x_min,p.x_max);
    mg.Draw();

    c0->Print(p.file_out);
}
