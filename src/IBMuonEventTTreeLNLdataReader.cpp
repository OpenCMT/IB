#include <TTree.h>
#include <TFile.h>
#include "IBMuonEventTTreeLNLdataReader.h"

using namespace uLib;

class IBMuonEventTTreeLNLdataReaderPimpl
{
    struct MuonTrack_projection {
        float position;
        float slope;
        float position_error;
        float slope_error;
    };

    typedef struct
    {
        struct MuonTrack_projection phi_up;
        struct MuonTrack_projection phi_down;
        struct MuonTrack_projection theta_up;
        struct MuonTrack_projection theta_down;
    } MuonTracks;

    typedef union {
        struct {
            unsigned short int phi_up : 4;
            unsigned short int phi_down : 4;
            unsigned short int theta_up : 4;
            unsigned short int theta_down : 4;
        };
        unsigned short int data;
    } MuonHits;

#define SAFE_TTREE_BUFFER 50
    struct Demo_LNL_buffer {
        int event;
        int iseg;
        int sn[SAFE_TTREE_BUFFER];
        float sx[SAFE_TTREE_BUFFER];
        float ss[SAFE_TTREE_BUFFER];
        float ersx[SAFE_TTREE_BUFFER];
        float erss[SAFE_TTREE_BUFFER];
    };
#undef SAFE_TTREE_BUFFER

public:
    typedef struct Demo_LNL_buffer Buffer;

    IBMuonEventTTreeLNLdataReaderPimpl()
    {
        m_tree         = NULL;
        m_total_events = 0;
        m_max_event    = 0;
        m_pos          = 0;
        m_hitX         = 6;
        m_hitZ         = 4;
        m_momentum     = 0.f;
        m_integrity    = true;
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
    }

    void init(TFile * file)
    {
        if (file->IsZombie()) {
            printf("Requested file not found!\nAborting...\n");
            exit(0);
        }
        TTree* t = (TTree*)file->Get("RADMU");
        if (t) {
            init(t);
        } else {
            printf("Requested TTree not found in file! Maybe wrong TTree name?\nAborting...\n");
            exit(0);
        }
    }

    void init(TTree * tree)
    {
        m_tree = tree;
        m_max_event = m_tree->GetEntries();
        if (m_total_events == 0.f) m_total_events = m_max_event;
        if(m_total_events+m_pos>m_max_event) {
            printf("Requested time interval rejected at TTree initialization.\nAborting...\n");
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

    void AcquireEvent()
    {
        m_integrity = true;
        if (likely(m_pos <= m_max_event)) {
            m_tree->GetEntry(m_pos);
        } else m_integrity = false;
        m_pos++;
        return;
    }

    void GetMuonEvent(MuonScatter * muon_event)
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
        if (hit_count.phi_up < m_hitX || hit_count.phi_down < m_hitX ||
            hit_count.theta_up < m_hitZ || hit_count.theta_down < m_hitZ ) {
            m_integrity = false;
            return;
        }

        muon_event->LineIn().origin     << m_track.phi_up.position, 0, m_track.theta_up.position,1;  // HARDCODED
        muon_event->LineIn().direction  << -m_track.phi_up.slope, -1, -m_track.theta_up.slope,0;
        muon_event->LineOut().origin    << m_track.phi_down.position -1.004 , -183.43, m_track.theta_down.position -0.112 ,1; // HARDCODED
        muon_event->LineOut().direction << -m_track.phi_down.slope -0.0004, -1, -m_track.theta_down.slope -0.0027, 0;
        muon_event->SetMomentum(m_momentum);

        // HardCoded cuts on potition and delta slope
        if (fabs(muon_event->LineIn().origin(0))>308.                                      ||
            fabs(muon_event->LineIn().origin(2))>250.                                      ||
            fabs(muon_event->LineOut().origin(0))>308.                                     ||
            fabs(muon_event->LineOut().origin(2))>250.                                     ||
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


public:
#ifndef NDEBUG
    TFile*        m_out;
    TTree*        m_dumpster;
    float         m_xu, m_zu, m_xd, m_zd, m_pu, m_tu, m_pd, m_td;
#endif
    TTree*        m_tree;
    MuonTracks    m_track;
    Buffer        m_buffer;
    Scalarf       m_momentum;
    IBMuonError  *m_error;
    bool          m_integrity;
//    short         m_code;
    int           m_hitX;
    int           m_hitZ;
    unsigned long m_total_events;
    unsigned long m_max_event;
    unsigned long m_pos;

};

IBMuonEventTTreeLNLdataReader::IBMuonEventTTreeLNLdataReader() :
    d(new IBMuonEventTTreeLNLdataReaderPimpl) {}

IBMuonEventTTreeLNLdataReader::~IBMuonEventTTreeLNLdataReader()
{
#ifndef NDEBUG
    d->m_out->cd();
    d->m_dumpster->Write();
    d->m_out->Close();
#endif
    delete d;
}

void IBMuonEventTTreeLNLdataReader::setTTree(TTree *tree)
{
    d->init(tree);
}

void IBMuonEventTTreeLNLdataReader::setTFile(TFile* file)
{
    d->init(file);
}

void IBMuonEventTTreeLNLdataReader::setHitCuts(int nx_cut, int nz_cut)
{
    d->m_hitX = nx_cut;
    d->m_hitZ = nz_cut;
}

void IBMuonEventTTreeLNLdataReader::setMomentum(Scalarf p)
{
    d->m_momentum = p;
}

void IBMuonEventTTreeLNLdataReader::selectionCode(short code)
{
}

void IBMuonEventTTreeLNLdataReader::setError(IBMuonError &e)
{
    d->m_error = &e;
}

void IBMuonEventTTreeLNLdataReader::setAcquisitionTime(float min)
{
    d->m_total_events = min * 6.7E4;
    if(d->m_total_events+d->m_pos>d->m_max_event) {
        printf("Requested time interval rejected at acquisition time setting.\nAborting...\n");
        exit(0);
    }
    printf("Processing %.2f minutes (%i events)\n", min, d->m_total_events);
}

unsigned long IBMuonEventTTreeLNLdataReader::getNumberOfEvents()
{
    return d->m_total_events;
}

unsigned long IBMuonEventTTreeLNLdataReader::getCurrentPosition()
{
    return d->m_pos;
}

bool IBMuonEventTTreeLNLdataReader::readNext(MuonScatter *event)
{
    d->AcquireEvent();
    d->GetMuonEvent(event);
    return d->m_integrity;
}


void IBMuonEventTTreeLNLdataReader::setStartTime(float min)
{
    d->m_pos = min * 6.7E4;
    if (unlikely(d->m_total_events+d->m_pos>d->m_max_event)) {
        printf("Requested time interval rejected at start time setting\nAborting...\n");
        exit(0);
    }
    printf("Starting processing from %.2f minutes (event n. %i)\n", min, d->m_pos);
}
