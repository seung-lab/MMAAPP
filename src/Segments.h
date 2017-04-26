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
    virtual void insertSegment(id_type segid) { m_segids.insert(segid); };
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
    virtual void insertFreeEnds(id_type segid, const SupervoxelSet & freeEnds);
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
    const SupervoxelSet & freeEnds(id_type segid) { return m_freeends[segid]; }
    const SupervoxelSet & shaft(id_type segid) { return m_shafts[segid]; }
    id_type anchor(id_type freeend) const { return m_anchors[freeend]; }
    virtual void insertFreeEnds(id_type segid, id_type anchor, const SupervoxelSet & shaft);
private:
    SegmentDict m_segment;
    SupervoxelDict m_freeends;
    SupervoxelDict m_shafts;
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
    enum SegmentType {
        Unknown = 0,
        Axon,
        Dendrite,
        Glial
    };
    Segmentation(SupervoxelInfo * svInfo, SegmentInfo * segInfo)
        :m_sizeThreshold(1000000)
        ,m_axons()
        ,m_dendrites()
        ,m_spines()
        ,m_smallSegments()
        ,m_processedSegments()
    {
        m_svInfo = svInfo;
        m_segInfo = segInfo;
    }
    ~Segmentation() {}
    void init();
    int segLength(id_type segid);
    size_type segSize(id_type segid);
    size_type sumSize(const SupervoxelSet & svList);
    QVector<value_type > segSem(id_type segid);
    QVector<value_type > sumSem(const SupervoxelSet & svList);
    id_type largestSupervoxel(id_type segid, size_type * maxSize);
    id_type maximumPSD(id_type segid, value_type * max);
    id_type findPSD(id_type segid);
    SegmentType classifySegment(id_type segid, SupervoxelSet & freeEnds);
    SegmentType checkSemantic(QVector<value_type > sem, size_type vol);
    bool isGlial(QVector<value_type > sem, size_type vol);
    bool isAxon(QVector<value_type > sem, size_type vol);
    bool isDendrite(QVector<value_type > sem, size_type vol);
    int checkConnectivity(const SupervoxelSet & svList, const SupervoxelSet & exclude);
    SupervoxelSet findEnds(const SupervoxelSet & svList, id_type seed, const SupervoxelSet & exclude, bool free = true);
    SupervoxelSet verifyFreeEnds(const SupervoxelSet & ends, const SupervoxelSet & svList);
    void processDendrite(id_type segid);
private:
    size_type m_sizeThreshold;
    Axons m_axons;
    Dendrites m_dendrites;
    Spines m_spines;
    Segments m_smallSegments;
    Segments m_processedSegments;
    SupervoxelInfo * m_svInfo;
    SegmentInfo * m_segInfo;
};
#endif
