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


#include <TTree.h>
#include <TFile.h>
#include "IBMuonEventTTreeLNLdataReader.h"
#include <Core/Macros.h>

using namespace uLib;


IBMuonEventTTreeLNLdataReader::IBMuonEventTTreeLNLdataReader() {
    m_tree         = NULL;
    m_total_events = 0;
    m_max_event    = 0;
    m_pos          = 0;
    m_hitX         = 6;
    m_hitZ         = 4;
    m_momentum     = 0.f;
    m_integrity    = true;
    m_align = Matrix4f::Identity();
#ifndef NDEBUG
    m_out = new TFile("evDistro.root","RECREATE");
    //m_out = out;
    m_dumpster = new TTree("ev","ev");
    //m_dumpster = dumpster;

    m_dumpster->Branch("x_up", &m_xu, "xu");
    m_dumpster->Branch("z_up", &m_zu, "zu");
    m_dumpster->Branch("p_up", &m_pu, "pu");
    m_dumpster->Branch("t_up", &m_tu, "tu");
    m_dumpster->Branch("x_down", &m_xd, "xd");
    m_dumpster->Branch("z_down", &m_zd, "zd");
    m_dumpster->Branch("p_down", &m_pd, "pd");
    m_dumpster->Branch("t_down", &m_td, "td");
#endif // NDEBUG
    m_align.translate(Vector3f(-0.96,-182.9,-0.3)); // DEFAULT RAW ALIGNMENT
}

IBMuonEventTTreeLNLdataReader::~IBMuonEventTTreeLNLdataReader()
{
#ifndef NDEBUG
    m_out->cd();
    m_dumpster->Write();
    m_out->Close();
#endif
}

void IBMuonEventTTreeLNLdataReader::setTTree(TTree *tree)
{
    init(tree);
}

void IBMuonEventTTreeLNLdataReader::setTFile(TFile* file)
{
    init(file);
}

void IBMuonEventTTreeLNLdataReader::setHitCuts(int nx_cut, int nz_cut)
{
    m_hitX = nx_cut;
    m_hitZ = nz_cut;
}

void IBMuonEventTTreeLNLdataReader::setMomentum(Scalarf p)
{
    m_momentum = p;
}

void IBMuonEventTTreeLNLdataReader::selectionCode(short code)
{
}

void IBMuonEventTTreeLNLdataReader::setError(IBMuonError &e)
{
    m_error = &e;
}

void IBMuonEventTTreeLNLdataReader::setAcquisitionTime(float min)
{
    std::cout << "\n---------------" << std::endl;
    m_total_events = min * events_per_minute();
    if(m_total_events + m_pos > m_max_event) {
        printf("ATTENTION : Requested events (%i) from start event (%i) are not available.... processing only %i events\n", m_total_events, m_pos, m_max_event - m_pos);
         m_total_events = m_max_event - m_pos;
    }
    std::cout << "Processing " << m_total_events << " events from event n. " <<  m_pos << " to event n. " << m_total_events + m_pos << "\n" << std::endl;
}

unsigned long IBMuonEventTTreeLNLdataReader::getNumberOfEvents()
{
    return m_total_events;
}

unsigned long IBMuonEventTTreeLNLdataReader::getCurrentPosition()
{
    return m_pos;
}

bool IBMuonEventTTreeLNLdataReader::readNext(MuonScatter *event)
{
    // Acquire Event
    m_integrity = true;
    if (likely(m_pos <= m_max_event)) {
        m_tree->GetEntry(m_pos);
    } else m_integrity = false;
    m_pos++;

    GetMuonScatter(event);
    return m_integrity;
}


void IBMuonEventTTreeLNLdataReader::setStartTime(float min)
{
    std::cout << "\n---------------" << std::endl;
    m_pos = (unsigned long)(min * events_per_minute());
    if (unlikely(m_pos>m_max_event)) {
        printf("\nATTENTION : Requested start time event number %i greater than total events %i.\nAborting...\n", m_pos, m_max_event);
        exit(0);
    }
    std::cout << "Set start time " << min << " minutes i.e. event n. " << m_pos << std::endl;
}

void IBMuonEventTTreeLNLdataReader::setAlignmentFromData(float min)
{
    Matrix4f M = Matrix4f::Zero();
    Vector4f K = Vector4f::Zero();
    Vector2f w(1,1);
    MuonScatter mu;

    unsigned long ev;
    if(min == 0)
        ev = m_max_event;
    else {
        ev = (unsigned long)(min * events_per_minute());
    }

    int tot=0;
    std::cout << "Aligning Reader with data.. \n";
    for( int i=0; i<80; ++i) std::cout << "-"; std::cout << "\r";
    for (unsigned long i=0; i<ev; i++) {
        m_integrity = true;
        m_tree->GetEntry(i);
        GetMuonScatter(&mu);

        if(m_integrity) {
            float &x1 = mu.LineIn().origin(0);
            float &x2 = mu.LineOut().origin(0);
            float &z1 = mu.LineIn().origin(2);
            float &z2 = mu.LineOut().origin(2);
            float &tphi1 = mu.LineIn().direction(0);
            float &tphi2 = mu.LineOut().direction(0);
            float &tthe1 = mu.LineIn().direction(2);
            float &tthe2 = mu.LineOut().direction(2);

            Matrix4f Mi;
            Mi <<
                  w(0) * tphi1*tphi1 + w(1) * tthe1*tthe1, w(1) * x2 * tthe1 - w(0) * z2 * tphi1, w(0) * tphi1, w(1) * tthe1,
                  w(1) * x2 * tthe1 - w(0) * z2 * tphi1,   w(0) * z2 * z2 + w(1) * x2 *x2,        -w(0) * z2,   w(1) * x2,
                  w(0) * tphi1,                            -w(0) * z2,                             w(0),         0,
                  w(1) * tthe1,                             w(1) * x2,                               0 ,         w(1);

            M += Mi * 2;

            Vector4f Ki = Vector4f::Zero();
            Ki <<
                  w(0) * (x1 - x2) * tphi1 + w(1) * (z1 - z2) * tthe1,
                  w(0) * (x2 - x1) * z2  + w(1) * (z1 - z2) * x2,
                  w(0) * (x1 - x2),
                  w(1) * (z1 - z2);
            K += Ki * 2;

            tot++;
        }
        if(tot++%(ev/80) == 0) std::cout << "o" << std::flush;
    }
    std::cout << "\n";
    Vector4f X = M.inverse() * K;
    std::cout << "ALIGNMENT FOUND [H, Phi, x0, z0]: " << X.transpose() << "\n";

    { // APPLY TO PIMPL TRANSFORM //
        Eigen::Affine3f &tr = m_align;
        tr.rotate(Eigen::AngleAxisf(-X(1),Vector3f(0,1,0)));
        tr.translation() += Vector3f(X(2),X(0)-tr.translation()(1),X(3));
    }

}

void IBMuonEventTTreeLNLdataReader::setAlignment(Matrix4f align)
{
    m_align = align;
}

Matrix4f IBMuonEventTTreeLNLdataReader::getAlignment()
{
    return m_align.matrix();
}



void IBMuonEventTTreeLNLdataReader::init(TFile * file) {
    // init variablers
    m_pos = 0;
    m_max_event = 0;
    m_total_events = 0;

    if (file->IsZombie()) {
        printf("Requested file not found!\nAborting...\n");
        // TODO remove exit(0) from library
        exit(0);
    }

    // init trees
    TTree* t = (TTree*)file->Get("RADMU");
    if (t) {
        init(t);
    } else {
        printf("Requested TTree not found in file! Maybe wrong TTree name?\nAborting...\n");
        // TODO remove exit(0) from library
        exit(0);
    }
}

void IBMuonEventTTreeLNLdataReader::init(TTree * tree) {
    m_tree = tree;
    m_max_event = m_tree->GetEntries();
    //std::cout << "Requested event " << m_total_events+m_pos << " from max of " << m_max_event << std::endl;
    if (m_total_events == 0.f) m_total_events = m_max_event;
    if(m_total_events+m_pos>m_max_event) {
        printf("Requested time interval rejected at TTree initialization.\nAborting...\n");
        // TODO remove exit(0) from library
        exit(0);
    }
    tree->SetBranchAddress("EVENT",      &m_buffer.event);
    tree->SetBranchAddress("SEG_ns_glo", &m_buffer.iseg);
    tree->SetBranchAddress("SEG_sx_glo",  m_buffer.sx);
    tree->SetBranchAddress("SEG_ss_glo",  m_buffer.ss);
    tree->SetBranchAddress("SEG_sn_glo",  m_buffer.sn);
    tree->SetBranchAddress("SEG_ersx_glo",m_buffer.ersx);
    tree->SetBranchAddress("SEG_erss_glo",m_buffer.erss);
}

void IBMuonEventTTreeLNLdataReader::GetMuonScatter(MuonScatter * muon_event)
{
    if(unlikely(!m_integrity)) return;
    float *track_ptr;
    MuonTracks *traks = &m_track;
    MuonHits hit_count;
    hit_count.data = 0;

    //fill tracks and hits
    for (int is = 0; is < 4; is++) {
        if ((m_buffer.sn[is] > 1300) && (m_buffer.sn[is] < 1310)) {     // Phi up
            track_ptr = (float *)&traks->phi_up;
            hit_count.phi_up = m_buffer.sn[is] - 1300;
        }
        else if ((m_buffer.sn[is] < -1300) && (m_buffer.sn[is] > -1310)) { // Theta up
            track_ptr = (float *)&traks->theta_up;
            hit_count.theta_up = -m_buffer.sn[is] - 1300;
        }
        else if ((m_buffer.sn[is] > 2300) && (m_buffer.sn[is] < 2310)) { // Phi down
            track_ptr = (float *)&traks->phi_down;
            hit_count.phi_down = m_buffer.sn[is] - 2300;
        }
        else if ((m_buffer.sn[is] < -2300) && (m_buffer.sn[is] > -2310)) { // Theta down
            track_ptr = (float *)&traks->theta_down;
            hit_count.theta_down = -m_buffer.sn[is] - 2300;
        } else {
            //debug_debug("error reading muon track code");
            break;
        }
        *track_ptr++ = m_buffer.sx[is];
        *track_ptr++ = m_buffer.ss[is];
        *track_ptr++ = m_buffer.ersx[is];
        *track_ptr =   m_buffer.erss[is];
    }

    if (    hit_count.phi_up < m_hitX     ||
            hit_count.phi_up > 8          ||
            hit_count.phi_down < m_hitX   ||
            hit_count.phi_down > 8        ||
            hit_count.theta_up < m_hitZ   ||
            hit_count.theta_up > 4        ||
            hit_count.theta_down < m_hitZ ||
            hit_count.theta_down > 4       )
    {
        m_integrity = false;
        return;
    }

    muon_event->LineIn().origin     <<  m_track.phi_up.position,   0,  m_track.theta_up.position,   1;  // HARDCODED
    muon_event->LineIn().direction  << -m_track.phi_up.slope,     -1, -m_track.theta_up.slope,      0;
    muon_event->LineOut().origin    <<  m_track.phi_down.position, 0,  m_track.theta_down.position, 1; // AUTO
    muon_event->LineOut().direction << -m_track.phi_down.slope,   -1, -m_track.theta_down.slope,    0;


    { // ALIGNMENT //
        const Eigen::Affine3f &tr = m_align;
        muon_event->LineOut().origin = tr * muon_event->LineOut().origin;
        muon_event->LineOut().direction = tr * muon_event->LineOut().direction;
        muon_event->LineOut().direction /= fabs(muon_event->LineOut().direction(1)); // back to slopes //
    }

    muon_event->SetMomentum(m_momentum);
    muon_event->SetMomentumPrime(m_momentum);

    // HardCoded cuts on potition and delta slope
    if (fabs(muon_event->LineIn().origin(0))>154.                                      ||
        fabs(muon_event->LineIn().origin(2))>125.                                      ||
        fabs(muon_event->LineOut().origin(0))>154.                                     ||
        fabs(muon_event->LineOut().origin(2))>125.                                     ||
        fabs(muon_event->LineIn().direction(0)) > 1.4                                  ||
        fabs(muon_event->LineIn().direction(2)) > 1.26                                 ||
        fabs(muon_event->LineOut().direction(0)-muon_event->LineIn().direction(0))>0.5 ||
        fabs(muon_event->LineOut().direction(2)-muon_event->LineIn().direction(2))>0.5  )
    {
        m_integrity = false;
    }

    muon_event->ErrorIn().direction_error  = Vector4f::Zero();
    muon_event->ErrorOut().direction_error = Vector4f::Zero();

    m_error->evaluate(*muon_event, 2, 2);

#ifndef NDEBUG
    m_xu = muon_event->LineIn().origin(0);
    m_zu = muon_event->LineIn().origin(2);
    m_pu = muon_event->LineIn().direction(0);
    m_tu = muon_event->LineIn().direction(2);
    m_xd = muon_event->LineOut().origin(0);
    m_zd = muon_event->LineOut().origin(2);
    m_pd = muon_event->LineOut().direction(0);
    m_td = muon_event->LineOut().direction(2);
    m_dumpster->Fill();
#endif // NDEBUG
}




