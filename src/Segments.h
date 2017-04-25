#ifndef MMAAPP_SEGMENTS_H
#define MMAAPP_SEGMENTS_H
#include "Info.h"
class Segments
{
public:
    Segments()
        :m_segids()
    {
    }
    virtual ~Segments() {}
    SupervoxelSet & segids() {return m_segids;}
private:
    SupervoxelSet m_segids;
};

class Axons : public Segments
{
public:
    Axons()
        :Segments()
        ,m_segment()
        ,m_freeends()
    {
    }
    virtual ~Axons() {}
    id_type segid(id_type freeend) const { return m_segment[freeend]; }
    QSet<id_type > & freeends(id_type segid) {return m_freeends[segid];}
private:
    SegmentDict m_segment;
    SupervoxelDict m_freeends;
};

class Dendrites : public Segments
{
public:
    Dendrites()
        :Segments()
        ,m_segment()
        ,m_freeends()
        ,m_shafts()
        ,m_anchors()
    {
    }
    virtual ~Dendrites() {}
    id_type segid(id_type freeend) { return m_segment[freeend]; }
    QSet<id_type > & freeends(id_type segid) {return m_freeends[segid];}
    id_type shaft(id_type segid) const { return m_shafts[segid]; }
    id_type anchor(id_type freeend) const { return m_anchors[freeend]; }
private:
    SegmentDict m_segment;
    SupervoxelDict m_freeends;
    SegmentDict m_shafts;
    SegmentDict m_anchors;
};

class Spines : public Segments
{
public:
    Spines()
        :Segments()
        ,m_segment()
        ,m_freeends()
        ,m_psd()
        ,m_anchors()
    {
    }
    virtual ~Spines() {}
    id_type segid(id_type freeend) const { return m_segment[freeend]; }
    QSet<id_type > & freeends(id_type segid) {return m_freeends[segid];}
    id_type psd(id_type segid) const { return m_psd[segid]; }
    id_type anchor(id_type freeend) const { return m_anchors[freeend]; }
private:
    SegmentDict m_segment;
    SupervoxelDict m_freeends;
    SegmentDict m_psd;
    SegmentDict m_anchors;
};

class Segmentation
{
public:
    Segmentation()
        :m_axons()
        ,m_dendrites()
        ,m_spines()
        ,m_smallSegments()
        ,m_processedSegments()
    {
    }
    ~Segmentation() {}
private:
    Axons m_axons;
    Dendrites m_dendrites;
    Spines m_spines;
    Segments m_smallSegments;
    Segments m_processedSegments;
};
#endif
