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


#include <string>
#include <iostream>
#include <iomanip>

#include <TMultiGraph.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TAxis.h>
#include <TBox.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <Math/QuantFunc.h>
#include <Math/QuantFuncMathCore.h>

#include <TMath.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <TArrow.h>
#include <TMarker.h>
#include <TText.h>
#include <TLatex.h>
//#include <TArrowEditor.h>
#include <TPaveText.h>


#include "IBROC.h"

#include "Core/Options.h"

using namespace uLib;


namespace ROCPlot {



class ROCHQuote {
public:

    ROCHQuote() :
        m_p1(0,0), m_p2(0,0),
        m_distance(0),
        m_marker(23)
    {}


    ROCHQuote(const ROCHQuote &copy)
    {
        m_p1 = copy.m_p1;
        m_p2 = copy.m_p2;
        m_distance = copy.m_distance;
        m_marker = copy.m_marker;
        m_Label = copy.m_Label;
        m_LineAttr = copy.m_LineAttr;
    }

    ~ROCHQuote() {
        this->Clear();
    }

    // compare for sorting ... second point comparison
    static bool Compare(const ROCHQuote &q1, const ROCHQuote &q2) {
        return q1.m_p2(0) < q2.m_p2(0);
    }

    void Draw() {

        this->Clear();
        //        m_m1 = TMarker(m_p1(0),m_p1(1),m_marker);
        //        m_m2 = TMarker(m_p2(0),m_p2(1),m_marker);

        float x12 = fabs(m_p1(0) - m_p2(0));

        m_Arrow1 = TArrow(m_p1(0),m_p1(1),m_p1(0),m_p1(1)+m_distance,0.01,"<|");
        m_Arrow2 = TArrow(m_p2(0),m_p2(1),m_p2(0),m_p2(1)+m_distance,0.01,"<|");
        m_Line   = TLine(m_p1(0),m_p1(1)+m_distance,m_p2(0),m_p2(1)+m_distance);
        m_Hline  = TLine(std::min(m_p1(0),m_p2(0))-x12*0.7,m_p1(1),std::max(m_p1(0),m_p2(0))+x12*0.7,m_p2(1));

        m_LineAttr.Copy(m_Arrow1);
        m_LineAttr.Copy(m_Arrow2);
        m_LineAttr.Copy(m_Line);
        m_LineAttr.Copy(m_Hline);


        Vector2f tp;
        tp = (m_p1 + m_p2)/2;
        tp(0) = fmin(m_p1(0),m_p2(0)) - (tp(0)-m_p1(0))/10;
        tp(1) += m_distance + m_LabelOffset;


        // LEAK !!!!
        TLatex *Text = TLatex().DrawLatex(tp(0),tp(1),m_Label.c_str());
        Text->SetTextSize(0.05);


        m_Pave.AddText(m_Label.c_str());

        m_Pave.SetX1(tp(0));

        m_Pave.SetY1(tp(1));

        //        m_box = TBox(tp(0),tp(1),tp(0)+gPad->PadtoX(ts(0)),tp(1)+gPad->PadtoY(ts(1)));
        //        m_box.SetFillColor(kWhite);
        //        m_box.SetFillStyle(1001);
        //        m_box.Draw();
        //        Text->Draw();

        m_Arrow1.Draw();
        m_Arrow2.Draw();
        m_Line.Draw();
        m_Hline.Draw();
    }

    void Clear() {
        m_Arrow1.Clear();
        m_Arrow2.Clear();
        m_Line.Clear();
    }

    uLibRefMacro(p1,Vector2f)
    uLibRefMacro(p2,Vector2f)
    uLibSetMacro(distance,Double_t)
    uLibRefMacro(LineAttr,TAttLine)
    uLibRefMacro(Label,std::string)
    uLibSetMacro(LabelOffset,Double_t)

private:
    Double_t m_distance;
    Vector2f m_p1,m_p2;
    std::string m_Label;

    TPaveText m_Pave;
    TLine     m_Line,m_Hline;
    TArrow    m_Arrow1,m_Arrow2;
    TAttLine  m_LineAttr;
    Int_t     m_marker;

    Double_t  m_LabelOffset;

    TBox m_box;
    Vector<TObject *> props;
};











////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ROC ROOT GRAPH
enum EnumROCStyle {
    StyleROC = 0,
    StyleNFPN = 1,
    StyleFPN = 2
};

class ROCGraph : public IBROC {
    typedef IBROC BaseClass;
public:

    ROCGraph(EnumROCStyle style = StyleNFPN, bool optimism = 0) :
        m_style(style),
        m_Optimism(optimism)
    {
        m_FPF.SetLineWidth(2);
        m_FPF.SetMarkerSize(0);
        m_FNF.SetLineWidth(3);
        m_FNF.SetMarkerSize(0);
    }


    virtual void Update() {
        this->Draw();
    }

    void SetStyle(EnumROCStyle style) {
        this->m_style = style;
        this->Update();
    }

    uLibRefMacro(FPF,TGraph);
    uLibRefMacro(FNF,TGraph);
    uLibRefMacro(FPF_1sigma,TGraph);
    uLibRefMacro(FPF_2sigma,TGraph);
    uLibRefMacro(FNF_1sigma,TGraph);
    uLibRefMacro(FNF_2sigma,TGraph);

    float Awo(int i) {
        return fabs(this->at(i).Awo()/100 - m_Optimism);
    }

    float Owa(int i) {
        return fabs(this->at(i).Owa()/100 - m_Optimism);
    }

//protected:
public:
    /**
     * Get Poisson confidence interval fo clevel percentile
     * see: PDG (2010) Statistics 33.3.2.5 Poisson data
     */
    Vector2f get_Poisson_interval(Scalarf clevel, int n, int N) {
        Scalarf  alpha = 1 - clevel/100;
        Vector2f nu(0,0); // low - high //
        nu <<
              ROOT::Math::chisquared_quantile(alpha/2,   2*n    ),
                ROOT::Math::chisquared_quantile(1.-alpha/2, 2*(n+1));
        return nu/2/N;
    }

    /**
     * Get Bernulli confidence interval fo clevel percentile
     * see: PDG (2010) Statistics 33.3.2.5 Poisson data
     */
    Vector2f get_Bernulli_interval(Scalarf clevel, int n, int N) {
        Scalarf  alpha = 1 - clevel/100;
        Vector2f nu(0,0); // low - high //
        float Ff1, Ff2;
        Ff1 =     n * ROOT::Math::fdistribution_quantile( alpha/2,   2*n,     2*(N-n+1) );
        Ff2 = (n+1) * ROOT::Math::fdistribution_quantile( 1.-alpha/2, 2*(n+1), 2*(N-n)+0.00001   );
        nu <<   Ff1 / (N - n + 1 + Ff1),
                Ff2 / (N - n + Ff2);
        return nu;
    }


    int Draw() {
        m_FPF.Set(this->size());
        m_FNF.Set(this->size());
        m_FPF_1sigma.Set(2*this->size());
        m_FPF_2sigma.Set(2*this->size());
        m_FNF_1sigma.Set(2*this->size());
        m_FNF_2sigma.Set(2*this->size());

        m_FPF.Clear();
        m_FNF.Clear();
        m_FPF_1sigma.Clear();
        m_FPF_2sigma.Clear();
        m_FNF_1sigma.Clear();
        m_FNF_2sigma.Clear();

        float lf = 1;
        switch(m_style) {
        case ROCPlot::StyleROC:
            m_FPF.SetLineColor(kRed);
            m_FNF.SetLineColor(kRed);
            for(int i=0; i<this->size(); i++) {
                m_FPF.SetPoint(i,this->at(i).Awo()/100 ,1 - this->at(i).Owa()/100);
                m_FNF.SetPoint(i,1 - this->at(i).Owa()/100,this->at(i).Awo()/100); // not used
            }
            break;
        case ROCPlot::StyleFPN:
        case ROCPlot::StyleNFPN:
            if(m_style == ROCPlot::StyleFPN)
                lf = s_Shultz2Lambda;
            m_FPF.SetLineColor(kRed);
            m_FNF.SetLineColor(kBlue);
            m_FPF.SetLineWidth(2);
            m_FNF.SetLineWidth(2);

            m_FPF_1sigma.SetFillStyle(1001);
            m_FPF_2sigma.SetFillStyle(1001);
            m_FNF_1sigma.SetFillStyle(1001);
            m_FNF_2sigma.SetFillStyle(1001);
            m_FPF_1sigma.SetFillColor(kGreen-6);
            m_FPF_2sigma.SetFillColor(kYellow);
            m_FNF_1sigma.SetFillColor(kGreen-6);
            m_FNF_2sigma.SetFillColor(kYellow);
            for(int i=0; i<this->size(); i++) {
                m_FPF.SetPoint(i,this->at(i).X()*lf, Awo(i) );
                m_FNF.SetPoint(i,this->at(i).X()*lf, Owa(i) );
                if(this->Samples()(0) > 0 && this->Samples()(1) > 0)
                {
                    Vector2f FPF_interval, FNF_interval;
                    // 1 sigma //
                    FPF_interval = get_Bernulli_interval(68.3, Awo(i)*this->Samples()(0), this->Samples()(0) );
                    m_FPF_1sigma.SetPoint(i,this->at(i).X()*lf,FPF_interval[0]);
                    m_FPF_1sigma.SetPoint(2*this->size()-i-1,this->at(i).X()*lf,FPF_interval[1]);

                    FNF_interval = get_Bernulli_interval(68.3, Owa(i)*this->Samples()(1), this->Samples()(1) );
                    m_FNF_1sigma.SetPoint(i,this->at(i).X()*lf,FNF_interval[0]);
                    m_FNF_1sigma.SetPoint(2*this->size()-i-1,this->at(i).X()*lf,FNF_interval[1]);

                    // 2 sigma //
                    FPF_interval = get_Bernulli_interval(95, Awo(i)*this->Samples()(0), this->Samples()(0) );
                    m_FPF_2sigma.SetPoint(i,this->at(i).X()*lf,FPF_interval[0]);
                    m_FPF_2sigma.SetPoint(2*this->size()-i-1,this->at(i).X()*lf,FPF_interval[1]);

                    FNF_interval = get_Bernulli_interval(95, Owa(i)*this->Samples()(1), this->Samples()(1) );
                    m_FNF_2sigma.SetPoint(i,this->at(i).X()*lf,FNF_interval[0]);
                    m_FNF_2sigma.SetPoint(2*this->size()-i-1,this->at(i).X()*lf,FNF_interval[1]);

                }
            }
            break;
        }

        return this->size();
    }

    Vector2f GetPoint(const TGraph &g, int i) {
        double p1,p2;
        g.GetPoint(i,p1,p2);
        return Vector2f(p1,p2);
    }



    TGraph GetSubgraph(const TGraph &graph, bool uplow=0) {
        TGraph subgraph;
        int size = graph.GetN()/2;
        subgraph.Set(size);
        for(int i=0; i<size; ++i) {
            subgraph.SetPoint(i,graph.GetX()[size*uplow+i],graph.GetY()[size*uplow+i]);
        }
        return subgraph;
    }

    Vector3f GetRatio(Scalarf y) {
        this->Update();

        int sgn = (1-2*m_Optimism);

        float x_min = GetIntersection(m_FPF,y,-sgn);
        float x_max = GetIntersection(m_FNF,y,+sgn);

        // Error bands starts from lower boundary //
        // upper boundary has inverted direction  //
        TGraph sgr;
        sgr = GetSubgraph(m_FPF_1sigma, m_Optimism);
        float x_min_l = GetIntersection(sgr, y, -1);
        sgr = GetSubgraph(m_FPF_1sigma, !m_Optimism);
        float x_min_r = GetIntersection(sgr, y, +1);
        sgr = GetSubgraph(m_FNF_1sigma, !m_Optimism);
        float x_max_l = GetIntersection(sgr, y, -1);
        sgr = GetSubgraph(m_FNF_1sigma, m_Optimism);
        float x_max_r = GetIntersection(sgr, y, +1);

        //        std::cout << " intercepts  -> " << Vector2f(x_min,x_max).transpose()     << "\n";
        //        std::cout << " bandA min L/R -> " << Vector2f(x_min_l,x_min_r).transpose() << "\n";
        //        std::cout << " bandA max L/R -> " << Vector2f(x_max_l,x_max_r).transpose() << "\n";

        x_max_l -= x_max;
        x_max_r -= x_max;
        x_min_l -= x_min;
        x_min_r -= x_min;

        //        std::cout << " intercepts  -> " << Vector2f(x_min,x_max).transpose()     << "\n";
        //        std::cout << " err min L/R -> " << Vector2f(x_min_l,x_min_r).transpose() << "\n";
        //        std::cout << " err max L/R -> " << Vector2f(x_max_l,x_max_r).transpose() << "\n";

        float ratio = x_max/x_min;                     
        float rup = ratio * sqrt(pow(x_max_r/x_max,2) + pow(x_min_l/x_min,2) );
        float rdown = ratio * sqrt(pow(x_max_l/x_max,2) + pow(x_min_r/x_min,2) );

        std::cout << "RATIO: " << ratio << " +" << rup << " -" << rdown << "\n";

        return Vector3f(ratio,rup,rdown);
    }

    float GetIntersection(const TGraph &graph, Double_t y, int a){
        Double_t * vx = graph.GetX();
        Double_t * vy = graph.GetY();
        float dy = vy[0] - y;
        float x=0;
        int i=1;   while( vy[i]*a <= y*a && i<graph.GetN() ) { dy=vy[i]-y; i++; }
        int j = i-1; while( vy[j] == y && j > 0 ) j--;
        if(j<i-1) {
            // plateau: get midpt //
            for(int k=j; k<i; ++k) x += vx[k];
            x /= i-j;
        }
        else {
            // get x between 2 points interp //
            Vector2f p1,p2;
            p1 = GetPoint(graph, i            );
            p2 = GetPoint(graph, i+a*sign(dy) );
            x = (p2(0)-p1(0))/(p2(1)-p1(1))*(y-p1(1))+p1(0);
        }
        return x;
    }

    void SetOptimism(bool val) {
        this->m_Optimism = val;
        this->Update();
    }

private:
    // members //
    TGraph               m_FPF, m_FNF;
    TGraph               m_FPF_1sigma, m_FPF_2sigma, m_FNF_1sigma, m_FNF_2sigma;
    EnumROCStyle         m_style;
    bool                 m_Optimism;
    static constexpr float s_Shultz2Lambda = 0.04;
};












////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ROC ROOT MULTI GRAPH

class ROCMultiGraph : public TMultiGraph {

    typedef TMultiGraph BaseClass;
public:
    ROCMultiGraph(EnumROCStyle style = StyleFPN) :
        m_xmin(0),
        m_xmax(100),
        m_ymin(0),
        m_ymax(1),
        m_legendpos(0.685,0.12), //0.685,0.12,0.88,0.22
        m_ShowLegend(1),
        m_style(style),
        m_RLevel(0),
        m_DisplayBands(0),
        m_Optimism(0)
    {
        m_legend=new TLegend(m_legendpos(0),     m_legendpos(1),
                             m_legendpos(0)+0.20,m_legendpos(1)+0.10);
        m_legend->SetTextFont(72);
        m_legend->SetTextSize(0.03);
        m_legend->SetFillColor(0);
        m_legend->SetCornerRadius(2);
        m_legend->SetMargin(0.2);
        m_legend->SetBorderSize(0);

        m_FreeLabel.SetTextSize(0.06);

        //m_ErrorBands << 1.5,4;
    }

    ~ROCMultiGraph() {
        delete m_legend;
    }

    void AddROC(ROCGraph rocg, int line_style = 1 ) {
        rocg.SetStyle(m_style);
        rocg.SetOptimism(m_Optimism);
        rocg.FPF().SetLineStyle(line_style);
        rocg.FNF().SetLineStyle(line_style);
        this->m_rocs.push_back(rocg);        
        this->m_quotes.push_back(ROCHQuote());
    }

    inline ROCGraph & operator[] (int i) {
        return m_rocs[i];
    }

    void SetXRange(float min,float max) {
        if(!isnan(min))
            m_xmin = min;
        if(!isnan(max))
            m_xmax = max;
    }

    void SetYRange(float min, float max) {
        if(!isnan(min))
            m_ymin = min;
        if(!isnan(max))
            m_ymax = max;
    }

    void SetStyle(EnumROCStyle style) {
        m_style = style;
    }

    void Update() {
        if(m_rocs.size() > 0) {
            gPad->Update();
            this->Draw();
        }
    }

    void SetRatioLevel(float clevel) {
        m_RLevel = clevel;
    }

    void Draw()  {
        // Legend //
        m_legend->Clear();
        if(m_TimeLabel != std::string("")) {
            m_legend->SetHeader(m_TimeLabel.c_str());
        }
        if(m_rocs.size() == 1 && (m_style == StyleNFPN || m_style == StyleFPN) )
        {
            ROCGraph & roc = m_rocs.at(0);
            m_legend->AddEntry(&roc.FPF(),"FPF","l");
            m_legend->AddEntry(&roc.FNF(),"FNF","l");
        }
        else
        {
            for (int i=0; i< m_rocs.size(); ++i) {
                char lratio[200] = "";
                std::string rlabel = m_rocs.at(i).m_labels[0] + std::string(lratio);
                m_legend->AddEntry(&m_rocs.at(i).FPF(), rlabel.c_str(), "l");
            }
        }
        m_legend->SetTextSize(0.04);

        // clear all plots //
        this->Clear();

        // draw plots //
        switch(m_style) {
        case ROCPlot::StyleROC:
            DrawROC();
            // set graphics options //
            this->GetXaxis()->SetLimits(0.,1.);
            this->SetMinimum(0);
            this->SetMaximum(1);
            this->GetXaxis()->SetNdivisions(5);
            this->GetYaxis()->SetNdivisions(5);
            this->GetXaxis()->SetTitle("1 - Specificity     (FPF)");
            this->GetYaxis()->SetTitle("Sensitivity      (TPF)");
            BaseClass::Draw("a");
            if(m_rocs.size() > 1 && m_ShowLegend)
                m_legend->Draw();
            m_FreeLabel.Draw();
            break;
        case ROCPlot::StyleNFPN:
            DrawFPN();
            // set graphics options //
            this->SetMinimum((Double_t)m_ymin);
            this->SetMaximum((Double_t)m_ymax);
            this->GetYaxis()->SetLimits(m_ymin,m_ymax);
            this->GetYaxis()->SetNdivisions(5);
            this->GetXaxis()->SetLimits(m_xmin,m_xmax);
            this->GetXaxis()->SetNdivisions(10);
            this->GetXaxis()->SetTitle("Normalized ssd threshold");
            if(m_Optimism)
                this->GetYaxis()->SetTitle("Sensitivity (TPF),   \nSpecificity (TNF)");
            else
                this->GetYaxis()->SetTitle("FPF, FNF");
            BaseClass::Draw("a");
            for(int i=0;i<m_quotes.size();++i) m_quotes[i].Draw();
            if(ShowLegend()) m_legend->Draw();
            m_FreeLabel.Draw();
            break;
        case ROCPlot::StyleFPN:
            DrawFPN();
            // set graphics options //
            this->SetMinimum((Double_t)m_ymin);
            this->SetMaximum((Double_t)m_ymax);
            this->GetYaxis()->SetLimits(m_ymin,m_ymax);
            this->GetYaxis()->SetNdivisions(5);
            this->GetXaxis()->SetLimits(m_xmin,m_xmax);
            this->GetXaxis()->SetNdivisions(10);
            this->GetXaxis()->SetTitle("ssd threshold [rad^{2}/cm]");
            if(m_Optimism)
                this->GetYaxis()->SetTitle("Sensitivity (TPF),   \nSpecificity (TNF)");
            else
                this->GetYaxis()->SetTitle("FPF, FNF");
            BaseClass::Draw("a");
            for(int i=0;i<m_quotes.size();++i) m_quotes[i].Draw();
            if(ShowLegend()) m_legend->Draw();
            m_FreeLabel.Draw();
            break;
        }

        // set generic graphics options //
        this->GetXaxis()->SetTitleOffset(1);
        this->GetYaxis()->SetTitleOffset(1);
        this->GetXaxis()->SetTitleFont(42);
        this->GetYaxis()->SetTitleFont(42);
        this->GetXaxis()->SetTitleSize(0.05);
        this->GetYaxis()->SetTitleSize(0.05);
        this->GetXaxis()->SetLabelFont(42);
        this->GetYaxis()->SetLabelFont(42);
        this->GetXaxis()->SetLabelSize(0.05);
        this->GetYaxis()->SetLabelSize(0.05);

    }



    uLibRefMacro (ErrorBands,Vector<float>)
    uLibRefMacro (TimeLabel,std::string)
    uLibRefMacro (ShowLegend,bool)

    void SetTitle(const char *title) {
        BaseClass::SetTitle(title);
        this->Update();
    }

    uLibGetMacro (Optimism,bool)
    void SetOptimism(bool opt) {
        if( this->m_Optimism != opt )
        {
            this->m_Optimism = opt;
            for(int i=0; i<m_rocs.size(); i++)
            {
                ROCGraph &roc = m_rocs.at(i);
                roc.SetOptimism(opt);
            }
            this->Update();
        }
    }

    void SetDisplayBands(int b) { this->m_DisplayBands = b; }

    void SetFreeLabel(const char *text, float pad_x = 0.5, float pad_y = 0.85) {
        if(m_Optimism) pad_y = 1-pad_y;
        m_FreeLabel.SetText(0,0,text);
        UInt_t w,h;
        m_FreeLabel.GetBoundingBox(w,h);
        m_FreeLabel.SetNDC();
        Vector2f size;
        size << 1.*w/gPad->GetWw(), 1.*h/gPad->GetWh();
        m_FreeLabel.SetX(pad_x-size(0)/2);
        m_FreeLabel.SetY(pad_y-size(1)/2);
        this->Update();
    }



private:

    void DrawFPN() {

        // set Style //
        for(int i=0; i< m_rocs.size(); ++i)
        {
            ROCGraph &roc = m_rocs.at(i);
            roc.SetStyle(m_style);
        }

        // Confidence intervals //
        if(m_DisplayBands > 1)
        for(int i=0; i< m_rocs.size(); ++i)
        {
            ROCGraph &roc = m_rocs.at(i);
            this->Add(&roc.FPF_2sigma(),"f");
            this->Add(&roc.FNF_2sigma(),"f");
        }
        if(m_DisplayBands > 0)
        for(int i=0; i< m_rocs.size(); ++i)
        {
            ROCGraph &roc = m_rocs.at(i);
            this->Add(&roc.FPF_1sigma(),"f");
            this->Add(&roc.FNF_1sigma(),"f");
        }

        // Graphs //
        for(int i=0; i< m_rocs.size(); ++i)
        {
            ROCGraph &roc = m_rocs.at(i);
            this->Add(&roc.FPF(),"L");
            this->Add(&roc.FNF(),"L");
        }

        // RATIO //
        if(m_RLevel > 0) {
            for(int i=0; i< m_rocs.size(); ++i)
            {
                Vector2f p1,p2;
                p1 << m_rocs[i].GetIntersection(m_rocs[i].FPF(),m_RLevel,-(1-2*m_Optimism)),m_RLevel;
                p2 << m_rocs[i].GetIntersection(m_rocs[i].FNF(),m_RLevel,+(1-2*m_Optimism)),m_RLevel;
                m_quotes[i].p1() = p1;
                m_quotes[i].p2() = p2;
                ROCGraph &g = m_rocs[i];
                TAttLine att(kBlack, g.FPF().GetLineStyle(), 1);
                m_quotes[i].LineAttr() = att;


                std::cout << "INIZIO loop: \n";
                float val = 0.002;
                Vector3f ratio(0,0,0);
                for(int j=-5; j<=5; j++) {
                    std::cout << "val: " << m_RLevel + j*val << "-> ";
                    ratio += g.GetRatio(m_RLevel + j*val);
                }

                g.SetOptimism( !m_Optimism );
                std::cout << "INIZIO loop ribaltato: \n";
                for(int j=-5; j<=5; j++) {
                    std::cout << "val: " << (1-m_RLevel) + j*val << "-> ";
                    ratio += g.GetRatio((1-m_RLevel) + j*val);
                }

                ratio /= 22;
                std::cout << "RATIO: " << ratio(0) << " +" << ratio(1) << " -" << ratio(2) << "\n";

                g.SetOptimism( m_Optimism );



                char lratio[200] = "";
                sprintf(lratio,"R = %.2f^{+%.2f}_{-%.2f}",
                        ratio(0),ratio(1),ratio(2));
                m_quotes[i].Label() = std::string(lratio);
            }

            // Set Quote y position //
            std::sort(m_quotes.begin(),m_quotes.end(),ROCHQuote::Compare);
            for(int i=0; i< m_rocs.size(); ++i)
            {
                m_quotes[i].Setdistance( (m_ymax*(!m_Optimism) - m_RLevel)/(m_quotes.size()+1) * (i+1) );
                m_quotes[i].SetLabelOffset( 0.08 * (1-2*m_Optimism) );
            }
        }        

        BaseClass::Draw("a");
    }






    void DrawROC() {
        for(int i=0; i< m_rocs.size(); ++i)
        {
            m_rocs.at(i).SetStyle(m_style);
            this->Add(&m_rocs.at(i).FPF(),"L");
        }
        BaseClass::Draw("al");
    }



private:
    std::string m_TimeLabel;
    TLatex m_FreeLabel;
    TLegend * m_legend;
    bool m_ShowLegend;
    Vector<ROCGraph> m_rocs;
    Vector<ROCHQuote> m_quotes;
    int m_DisplayBands;
    Vector<float> m_ErrorBands; // error band in yband percent //
    float m_xmin,m_xmax;
    float m_ymin,m_ymax;
    Vector2f m_legendpos;
    float m_RLevel;
    EnumROCStyle m_style;
    bool m_Optimism;
};



} // ROCPlot











/**
 * Print Canvas is a trick to get better looking graphs:
 * Actually it does not work properly:
 * see.. http://www-pnp.physics.ox.ac.uk/~rodriges/psfrag/#faq
 */
void print_canvas(TCanvas* c, const char* filename)
{
  TString tmpname(Form("/tmp/psfrag_tmp_%d.eps", gSystem->GetPid()));
  c->Print(tmpname);
  gSystem->Exec(Form("make_psfrags.pl -d -s 1.5 -f %d %s %s",
                     gStyle->GetTitleFont(),
             tmpname.Data(),
             filename));
  //gSystem->Unlink(tmpname);
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// MAIN
struct Parameters {
    char * file0;
    char * file1;
    char * file_out;
    int type;
    std::string file0_caption;
    std::string file1_caption;
    float x_min;
    float x_max;
    float y_min;
    float y_max;
    int error_bands;
    float clevel;
    unsigned int samples;
    std::string header;
    bool showleg;
    std::string hlegend;
    std::string title;
    bool optimism;
};

main(int argc, char * argv[]){

    // .................................................................... //
    Parameters p = {
        NULL,
        NULL,
        (char *)"rocplot.eps",
        1,
        "no label",
        "no label",
        NAN,
        NAN,
        NAN,
        NAN,
        0,
        0,
        400,
        "",
        1,
        "",
        "",
        true
    };

    Options opt("usage: IB_ROCPlot file1 file2 outfile [options]");
    opt.add_options()
            ("help", "get this help and exit")
            ("type",   &p.type,          "Plot type: 0=ROC, 1=FPN, 2=FPN-shultzconverted")
            ("label1", &p.file0_caption, "label of first plot")
            ("label2", &p.file1_caption, "label of second plot")
            ("xmin",   &p.x_min,         "min of x plot thresholds")
            ("xmax",   &p.x_max,         "max of x plot thresholds")
            ("ymin",   &p.y_min,         "min of y plot thresholds")
            ("ymax",   &p.y_max,         "max of y plot thresholds")
            ("bands",  &p.error_bands,   "adds error bands 1 or 2 sigma")
            ("clevel", &p.clevel,        "set level (%) to get ratio plot")
            ("samples",&p.samples,       "set fixed number of samples for rocs")
            ("header" ,&p.header,        "header label")
            ("showleg",&p.showleg,       "show legeng")
            ("hlegend",&p.hlegend,    "legend header label")
            ("title"  ,&p.title,      "set graph title")
            ("optimism", &p.optimism, "see FPN(=false pos and neg)-->false or TNP(=true neg and pos)-->true");


    opt.parse_command_line(argc,argv);


    if(argc>=4) {
        p.file0=argv[1];
        std::cout << "file1 = " << p.file0 << "\n";
        if(std::string(argv[2]) != "NULL" ) {
            p.file1=argv[2];
            std::cout << "file2 = " << p.file1 << "\n";
        }
        p.file_out=argv[3];
    }
    else {
        std::cerr << "error parsing command line: use --help \n";
        exit(1);
    }

    // .................................................................... //



    //    gStyle->SetPaperSize(30,30);     // default
    gStyle->SetTitleFont(60, "XY"); // don't ask (psfrag trick)

    ROCPlot::ROCMultiGraph mg;
    mg.SetStyle( (ROCPlot::EnumROCStyle)p.type  );
    mg.SetOptimism( p.optimism );

    ROCPlot::ROCGraph roc1,roc2;

    // ROC GRAPH 1 //
    roc1.read_csv(p.file0);
    if(roc1.Samples() == Vector2i(0,0))
        roc1.Samples() = Vector2i(p.samples,p.samples); // force samples //
    if(p.file0_caption.c_str())
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


    // ROC GRAPH 2 //
    if(p.file1)
    {
        roc2.read_csv(p.file1);
        if(roc2.Samples() == Vector2i(0,0))
            roc2.Samples() = Vector2i(p.samples,p.samples); // force samples //
        if(p.file1_caption.c_str())
            roc2.m_labels[0] = std::string(p.file1_caption);
        mg.AddROC(roc2);
        if(isnan(p.x_max))
            range_max = range_max > roc2.GetRange()(1) ? range_max : roc2.GetRange()(1);
        if(isnan(p.x_min))
            range_min = range_min < roc2.GetRange()(1) ? range_min : roc2.GetRange()(1);
    }

    // CANVAS //
    TCanvas * c0 = new TCanvas("c0","multigraph L3",200,10,700,500);
    c0->SetFillColor(0);
    c0->SetGrid();
    c0->SetBorderMode(0);

    // RANGE //
    if(isnan(p.x_max))
        p.x_max = range_max + range_max/5;
    if(isnan(p.x_min))
        p.x_min = range_min - range_min/5;
    mg.SetXRange(p.x_min,p.x_max);
    mg.SetYRange(p.y_min,p.y_max);
    mg.Update();

    // TITLE //
    if(p.title != "") {
        std::cout << "setting title: " << p.title << "\n";
        mg.SetTitle(p.title.c_str());
        gPad->Modified();
    }

    // HEADER LABEL //
    if(p.header != std::string("") ) {
        mg.SetFreeLabel(p.header.c_str());
    }


    // LEGEND //
    mg.ShowLegend() = p.showleg;
    if(p.hlegend != std::string("") ) {
        mg.TimeLabel() = p.hlegend;
    }


    // BANDS //
    if(p.error_bands) {
        mg.ErrorBands() << 6./4 + 1./100 , 17./4 + 1./100;
        mg.SetDisplayBands(p.error_bands);
    }


    if(p.clevel > 0) {
        mg.SetRatioLevel(p.clevel/100);

        std::cout << "********** ROC 1 *************\n";
        roc1.SetStyle(ROCPlot::StyleNFPN);
        roc1.GetRatio(p.clevel/100);

        if(roc2.size() > 0) {
            std::cout << "********** ROC 2 *************\n";
            roc2.SetStyle(ROCPlot::StyleNFPN);
            roc2.GetRatio(p.clevel/100);
        }
    }

    // DRAW //
    gPad->Update();
    mg.Draw();

    // SAVE //
    c0->Print(p.file_out,"eps");

    // BEAUTIFY //
    //    print_canvas(c0, p.file_out); // better looking graph (not working)
}
