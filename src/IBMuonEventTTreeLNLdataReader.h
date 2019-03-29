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



#ifndef IBMUONEVENTTTREELNLDATAREADER_H
#define IBMUONEVENTTTREELNLDATAREADER_H

#include "IBMuonEventTTreeReader.h"

class IBMuonEventTTreeLNLdataReader : public IBMuonEventTTreeReader
{
public:
    IBMuonEventTTreeLNLdataReader();
    ~IBMuonEventTTreeLNLdataReader();
    void setTTree(class TTree* tree);
    void setTFile(class TFile* file);
    void setHitCuts(int nx_cut, int nz_cut);
    void setMomentum(Scalarf p);
    void selectionCode(short code);
    void setError(IBMuonError &e);
    void setAcquisitionTime(float min);
    void setStartTime(float min);

    void setAlignmentFromData(float min = 0.0);
    void setAlignment(Matrix4f align);
    Matrix4f getAlignment();


    unsigned long getNumberOfEvents();
    unsigned long getCurrentPosition();

    bool readNext(uLib::MuonScatter *event);

private:

    typedef struct Demo_LNL_buffer Buffer;

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

    void init(TFile * file);
    void init(TTree * tree);
    void GetMuonScatter(MuonScatter * muon_event);
    static const float events_per_minute() { return 242.3912*60; } //converted to minutes

#ifndef NDEBUG
    TFile*        m_out;
    TTree*        m_dumpster;
    float         m_xu, m_zu, m_xd, m_zd, m_pu, m_tu, m_pd, m_td;
#endif
    TTree*        m_tree;
    MuonTracks    m_track;
    Scalarf       m_momentum;
    IBMuonError  *m_error;
    bool          m_integrity;
    int           m_hitX;
    int           m_hitZ;
    unsigned long m_total_events;
    unsigned long m_max_event;
    unsigned long m_pos;

    Eigen::Affine3f        m_align;
    struct Demo_LNL_buffer m_buffer;

};

#endif // IBMUONEVENTTTREELNLDATAREADER_H
