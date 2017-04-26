#include "Segments.h"
#include "Info.h"
#include <QQueue>
#include <QtDebug>

void Axons::insertFreeEnds(id_type segid, const SupervoxelSet & freeEnds)
{
    m_freeends[segid] |= freeEnds;
    foreach (auto s, freeEnds) {
        m_segment[s] = segid;
    }
}

void Dendrites::insertFreeEnds(id_type segid, id_type anchor, const SupervoxelSet & freeEnds)
{
    m_freeends[segid] |= freeEnds;
    foreach (auto s, freeEnds) {
        m_segment[s] = segid;
        m_anchors[s] = anchor;
    }
}

size_type Segmentation::segSize(id_type segid)
{
    return sumSize(m_segInfo->supervoxelList(segid));
}

size_type Segmentation::sumSize(const SupervoxelSet & svList)
{
    size_type sizeSum = 0;
    foreach (auto s, svList) {
        sizeSum += m_svInfo->supervoxelSize(s);
    }
    return sizeSum;
}

QVector<value_type > Segmentation::segSem(id_type segid)
{
    return sumSem(m_segInfo->supervoxelList(segid));
}

QVector<value_type > Segmentation::sumSem(const SupervoxelSet & svList)
{
    QVector<value_type > semSum({0,0,0,0});
    foreach (auto s, svList) {
        const QVector<value_type > & sSem = m_svInfo->semanticInfo(s);
        for (int i = 0; i < 4; i++) {
            semSum[i] += sSem[i];
        }
    }
    return semSum;
}

int Segmentation::segLength(id_type segid)
{
    return m_segInfo->supervoxelList(segid).size();
}

id_type Segmentation::largestSupervoxel(id_type segid, size_type * max)
{
    id_type c = 0;
    size_type c_size = 0;
    SupervoxelSet & svList = m_segInfo->supervoxelList(segid);
    foreach (auto s, svList) {
        if (c_size < m_svInfo->supervoxelSize(s)) {
            c = s;
            c_size = m_svInfo->supervoxelSize(s);
        }
    }
    if (max != NULL) {
        (*max) = c_size;
    }
    return c;
}

id_type Segmentation::maximumPSD(id_type segid, value_type * max)
{
    id_type c = 0;
    value_type c_value = 0;
    SupervoxelSet & svList = m_segInfo->supervoxelList(segid);
    foreach (auto s, svList) {
        const QVector<value_type > & sSem = m_svInfo->semanticInfo(s);
        if (c_value < sSem[3]) {
            c = s;
            c_value = sSem[3];
        }
    }
    if (max != NULL) {
        (*max) = c_value;
    }
    return c;
}

id_type Segmentation::findPSD(id_type segid)
{
    auto sem = segSem(segid);
    if (sem[3] < 200) {
        return 0;
    }
    value_type psd_max = 0;
    id_type m = maximumPSD(segid, &psd_max);
    if (psd_max > 0.5 * sem[3]) {
        return m;
    }
    return 0;
}

bool Segmentation::isGlial(QVector<value_type > sem, size_type vol)
{
    if (sem[2] > 0.5*vol) {
        return true;
    }
    return false;
}

bool Segmentation::isDendrite(QVector<value_type > sem, size_type vol)
{
    if (sem[1] > 0.5*vol && sem[1] > sem[0]) {
        return true;
    }
    return false;
}

bool Segmentation::isAxon(QVector<value_type > sem, size_type vol)
{
    if (sem[0] > 0.5*vol && sem[0] > sem[1]) {
        return true;
    }
    return false;
}

Segmentation::SegmentType Segmentation::checkSemantic(QVector<value_type > sem, size_type vol)
{
    SegmentType seg_type = Segmentation::Unknown;
    if (isGlial(sem, vol)) {
        seg_type = Segmentation::Glial;
    } else if (isDendrite(sem, vol)) {
        seg_type = Segmentation::Dendrite;
    } else if (isAxon(sem, vol)) {
        seg_type = Segmentation::Axon;
    }
    return seg_type;
}

int Segmentation::checkConnectivity(const SupervoxelSet & svList, const SupervoxelSet & exclude)
{
    int cc = 0;
    auto visited = exclude;
    foreach(auto root, svList) {
        if (visited.contains(root)) {
            continue;
        }
        QQueue<id_type > queue;
        queue.enqueue(root);
        while (!queue.isEmpty()) {
            root = queue.dequeue();
            if (visited.contains(root)) {
                continue;
            }
            visited.insert(root);
            foreach(auto neighbour,  m_svInfo->neighbours(root)) {
                if (visited.contains(neighbour) || !svList.contains(neighbour)) {
                    continue;
                }
                queue.enqueue(neighbour);
            }
        }
        cc += 1;
    }
    return cc;
}

SupervoxelSet Segmentation::findEnds(const SupervoxelSet & svList, id_type seed, const SupervoxelSet & exclude, bool free)
{
    auto visited = exclude;
    QQueue<id_type > queue;
    queue.enqueue(seed);
    SupervoxelSet ends;
    SupervoxelSet children;
    int depth = 0;
    while (!queue.isEmpty()) {
        auto root = queue.dequeue();
        if (!visited.contains(root)) {
            visited.insert(root);
            bool atEnd = true;
            foreach(auto neighbour,  m_svInfo->neighbours(root)) {
                if (visited.contains(neighbour) || queue.contains(neighbour) || !svList.contains(neighbour)) {
                    continue;
                }
                children.insert(neighbour);
                atEnd = false;
            }
            if (root != seed && m_svInfo->semanticInfo(root)[3] > 500) {
                atEnd = false;
            }
            if (atEnd && !m_svInfo->atBoundary(root)) {
                ends.insert(root);
            }
        }
        if (queue.isEmpty() && !children.isEmpty()) {
            if (!free) {
                ends.clear();
            }
            foreach (auto c, children) {
                queue.enqueue(c);
                if (!free && !m_svInfo->atBoundary(c)) {
                    ends.insert(c);
                }
            }
            children.clear();
            depth += 1;
        }
    }
    if (ends.size() == 0 && depth <= 1 && !m_svInfo->atBoundary(seed)) {
        ends.insert(seed);
    }
    return ends;
}

SupervoxelSet Segmentation::verifyFreeEnds(const SupervoxelSet & ends, const SupervoxelSet & svList)
{
    SupervoxelSet free_ends;
    foreach(auto s, ends) {
        SupervoxelSet end_segments;
        bool near_boundary = false;
        if (m_svInfo->atBoundary(s)) {
            continue;
        }
        SupervoxelSet neighbours_of_neighbours;
        foreach(auto neighbour, m_svInfo->neighbours(s)) {
            if (svList.contains(neighbour)) {
                if (m_svInfo->atBoundary(neighbour)) {
                    near_boundary = true;
                    break;
                }
                neighbours_of_neighbours |= (SupervoxelSet::fromList(m_svInfo->neighbours(neighbour)) & svList);
                end_segments.insert(neighbour);
            }
        }
        if (near_boundary || neighbours_of_neighbours.size() > 7) {
            continue;
        }
        end_segments.insert(s);
        int cc = checkConnectivity(svList, end_segments);
        if (cc == 1) {
            free_ends.insert(s);
        }
    }
    return free_ends;
}

Segmentation::SegmentType Segmentation::classifySegment(id_type segid, SupervoxelSet & freeEnds)
{
    auto total_size = segSize(segid);
    auto total_sem = segSem(segid);
    size_type max_size = 0;
    auto largest_sv = largestSupervoxel(segid, &max_size);
    auto seg_type = checkSemantic(total_sem, total_size);

    int cc = checkConnectivity(m_segInfo->supervoxelList(segid), QSet<id_type >({largest_sv}));
    if (total_size > m_sizeThreshold && (cc > 2 && max_size > (0.5 * total_size))) {
        if (seg_type == Segmentation::Unknown) {
            seg_type = Segmentation::Dendrite;
        }
        return seg_type;
    }
    auto ends = findEnds(m_segInfo->supervoxelList(segid), largest_sv, SupervoxelSet());
    if (cc == 1) {
        ends.insert(largest_sv);
    }

    freeEnds |= verifyFreeEnds(ends, m_segInfo->supervoxelList(segid));

    return seg_type;
}

void Segmentation::processDendrite(id_type segid)
{
    SupervoxelSet shaft;
    id_type m = largestSupervoxel(segid, NULL);
    shaft.insert(m);
    SupervoxelSet branches = m_segInfo->supervoxelList(segid) - shaft;
    SupervoxelSet anchors = branches & SupervoxelSet::fromList(m_svInfo->neighbours(m));
    if (segid == 82124) {
        qDebug() << anchors;
    }
    foreach (auto a, anchors) {
        auto free_ends = findEnds(m_segInfo->supervoxelList(segid), a, shaft);
        m_dendrites.insertFreeEnds(segid, a, free_ends);
    }
    qDebug() << "Dendrite: segid:" << segid << "parts:" << segLength(segid) << "size:" << segSize(segid) << "free_ends:" << m_dendrites.freeEnds(segid).size();
    qDebug() << m_dendrites.freeEnds(segid);
}

void Segmentation::init()
{
    auto segids = m_segInfo->compositeSegments();
    qDebug() << segids.size() << "segments";
    foreach (auto a, segids) {
        auto size_a = segSize(a);
        auto sem_a = segSem(a);
        auto length_a = segLength(a);
        if (segLength(a) >= 5) {
            SupervoxelSet free_ends;
            SegmentType seg_type = classifySegment(a, free_ends);
            qDebug() << "segid:" << a << "parts:" << segLength(a) << "size:" << size_a << "free_ends:" << free_ends.size() << free_ends << seg_type;

            if ((seg_type == Axon || seg_type == Unknown) && !free_ends.isEmpty()) {
                m_axons.insertSegment(a);
                m_axons.insertFreeEnds(a, free_ends);
                qDebug() << "Axon: segid:" << a << "parts:" << segLength(a) << "size:" << size_a << "free_ends:" << free_ends.size() << free_ends ;
            }

            if (seg_type == Dendrite && size_a > m_sizeThreshold) {
                processDendrite(a);
            }

            if (length_a <= 30 && size_a < m_sizeThreshold) {
                m_smallSegments.insertSegment(a);
                continue;
            }

            if (seg_type == Glial && length_a > 30) {
                m_processedSegments.insertSegment(a);
                continue;
            }

        }
        //id_type shaft_a = 0;
        //id_type syn_a = 0;
        //auto maxSize_a = maxSize(a, &shaft_a);
        //auto maxSyn_a = maxSyn(a, &syn_a);
        //qDebug() << "segid:" << a << "length:" << segLength(a) << "segsize:" << size_a << "sem:" << sem_a;
        //qDebug() << "max supervoxel:" << shaft_a << maxSize_a << "synapse:" << syn_a << maxSyn_a;
    }
}

