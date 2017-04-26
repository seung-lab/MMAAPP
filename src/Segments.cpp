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

void Spines::insertFreeEnds(id_type segid, id_type anchor, const SupervoxelSet & freeEnds)
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

bool Segmentation::processDendrite(id_type segid)
{
    SupervoxelSet shaft;
    id_type m = largestSupervoxel(segid, NULL);
    shaft.insert(m);
    m_dendrites.insertSegment(segid, shaft);
    SupervoxelSet branches = m_segInfo->supervoxelList(segid) - shaft;
    SupervoxelSet anchors = branches & SupervoxelSet::fromList(m_svInfo->neighbours(m));
    foreach (auto a, anchors) {
        auto free_ends = findEnds(m_segInfo->supervoxelList(segid), a, shaft);
        m_dendrites.insertFreeEnds(segid, a, free_ends);
    }
    qDebug() << "Dendrite: segid:" << segid << "parts:" << segLength(segid) << "size:" << segSize(segid) << "free_ends:" << m_dendrites.freeEnds(segid).size();
    return true;
}

bool Segmentation::processSpine(id_type segid, const SupervoxelSet & freeEnds)
{
    auto p = findPSD(segid);
    if (p == 0) {
        return false;
    }
    m_spines.insertSegment(segid, p);
    auto free_ends = findEnds(m_segInfo->supervoxelList(segid), p, SupervoxelSet(), false) & freeEnds;
    m_spines.insertFreeEnds(segid, p, free_ends);
    qDebug() << "Spine: segid:" << segid << "parts:" << segLength(segid) << "size:" << segSize(segid) << "free_ends:" << m_spines.freeEnds(segid).size();
    return true;
}

void Segmentation::init()
{
    auto segids = m_segInfo->compositeSegments();
    qDebug() << segids.size() << "segments";
    foreach (auto a, segids) {
        auto size_a = segSize(a);
        auto sem_a = segSem(a);
        auto length_a = segLength(a);
        if (length_a >= 5) {
            SupervoxelSet free_ends;
            SegmentType seg_type = classifySegment(a, free_ends);
            qDebug() << "segid:" << a << "parts:" << length_a << "size:" << size_a << "free_ends:" << free_ends.size() << free_ends << seg_type;

            if (seg_type == Glial && length_a > 30) {
                m_processedSegments.insertSegment(a);
                continue;
            }

            if (length_a <= 30 && size_a < m_sizeThreshold) {
                m_smallSegments.insertSegment(a);
            }

            if ((seg_type == Axon || seg_type == Unknown) && !free_ends.isEmpty()) {
                m_axons.insertSegment(a);
                m_axons.insertFreeEnds(a, free_ends);
                qDebug() << "Axon: segid:" << a << "parts:" << segLength(a) << "size:" << size_a << "free_ends:" << free_ends.size() << free_ends ;
                continue;
            }

            if (seg_type == Dendrite && size_a > m_sizeThreshold) {
                if (processDendrite(a)) {
                    continue;
                }
            }

            if (length_a < 100 && size_a < m_sizeThreshold) {
                if (processSpine(a, free_ends)) {
                    continue;
                }
            }
        } else if (length_a < 5 && size_a < m_sizeThreshold) {
                m_smallSegments.insertSegment(a);
                auto p = findPSD(a);
                if (p != 0) {
                    m_spines.insertSegment(a, p);
                    m_spines.insertFreeEnds(a, p, m_segInfo->supervoxelList(a));
                }
        }

        foreach (auto a, m_segInfo->allSegments() - segids) {
                m_smallSegments.insertSegment(a);
                if (m_svInfo->semanticInfo(a)[3] > 200) {
                    m_spines.insertSegment(a, a);
                    m_spines.insertFreeEnds(a, a, SupervoxelSet({a}));
                }

        }
    }
}

bool Segmentation::checkSegEdge(const MeanPlusEdge * edge)
{
    const SupervoxelSet & seg1 = m_segInfo->supervoxelList(edge->p1);
    const SupervoxelSet & seg2 = m_segInfo->supervoxelList(edge->p2);
    SupervoxelSet neighbour1;
    SupervoxelSet neighbour2;
    if (seg1.contains(edge->v1)) {
        neighbour1 = SupervoxelSet::fromList(m_svInfo->neighbours(edge->v1));
        neighbour2 = SupervoxelSet::fromList(m_svInfo->neighbours(edge->v2));
    } else {
        neighbour1 = SupervoxelSet::fromList(m_svInfo->neighbours(edge->v2));
        neighbour2 = SupervoxelSet::fromList(m_svInfo->neighbours(edge->v1));
    }
    if  ((edge->aff/edge->area) < m_reliableMeanAffinity) {
        return false;
    }

    int len1 = (seg2 & neighbour1).size();
    int len2 = (seg1 & neighbour2).size();

    if (len1 == 1 && len2 == 1) {
        if (edge->num > 1000) {
            return false;
        }
    } else if (len1 > 1 && len2 > 1) {
        return false;
    }
    return true;
}

void Segmentation::matchAxons(SupervoxelDict & mergeGraph)
{
    QSet<const MeanPlusEdge *> visitedEdges;
    SupervoxelSet processed;
    const SupervoxelSet & axons = m_axons.segids();
    qDebug() << "match axons:" << axons.size() << "candidates";
    auto all_free_ends = SupervoxelSet::fromList(m_axons.allFreeEnds());
    foreach (auto a, axons) {
        auto matches = SupervoxelSet::fromList(m_segInfo->neighbours(a)) & axons;
        foreach (auto b, matches) {
            const MeanPlusEdge * seg_edge = m_segInfo->edge(a,b);
            if (visitedEdges.contains(seg_edge)) {
                continue;
            }
            visitedEdges.insert(seg_edge);
            if (!all_free_ends.contains(seg_edge->v1) || m_processedSegments.segids().contains(seg_edge->v1)) {
                continue;
            }
            if (!all_free_ends.contains(seg_edge->v2) || m_processedSegments.segids().contains(seg_edge->v2)) {
                continue;
            }
            if (checkSegEdge(seg_edge)) {
                qDebug() << qMin(a,b) << qMax(a,b);
                //FIXME, this is wrong
                processed << seg_edge->v1 << seg_edge->v2;
                mergeGraph[a].insert(b);
                mergeGraph[b].insert(a);
            }
        }
    }
    foreach (auto a, processed) {
        m_processedSegments.insertSegment(a);
    }
}

value_type Segmentation::checkSupervoxelEdges(const SupervoxelSet & set_a, const SupervoxelSet & set_b)
{
    value_type max_aff = 0;
    foreach(auto a, set_a) {
        foreach(auto b, (set_b & SupervoxelSet::fromList(m_svInfo->neighbours(a)))) {
            const MeanPlusEdge * edge = m_svInfo->edge(a,b);
            if (edge->num > 1000) {
                continue;
            }
            if (set_b.size() > 5 && m_svInfo->semanticInfo(b)[3] > 500) {
                continue;
            }
            if (edge->aff/edge->area > 0.6 || (edge->aff/edge->area > m_reliableMeanAffinity && edge->num > 20 && edge->num < 500)) {
                if (max_aff < (edge->aff/edge->area)) {
                    max_aff = edge->aff/edge->area;
                }
            }
        }
    }
    return max_aff;
}

void Segmentation::attachSmallSegments(SupervoxelDict & mergeGraph)
{
    QSet<const MeanPlusEdge *> visitedEdges;
    SupervoxelSet dend_candidates = m_dendrites.segids();
    const SupervoxelSet & small_segs = m_smallSegments.segids();
    SupervoxelSet attached;
    bool keep_going = true;
    while (keep_going) {
        keep_going = false;
        SupervoxelSet new_candidates;
        foreach (auto b, (small_segs - attached)) {
            auto set_b = m_segInfo->supervoxelList(b);
            if (set_b.isEmpty()) {
                set_b.insert(b);
            }
            if (m_processedSegments.segids().contains(b)) {
                continue;
            }
            auto sem_b = sumSem(set_b);
            auto size_b = sumSize(set_b);
            if (set_b.size() > 10 && !isDendrite(sem_b, size_b)) {
                continue;
            }
            id_type max_a = 0;
            value_type max_aff = 0;
            foreach (auto a, (SupervoxelSet::fromList(m_segInfo->neighbours(b)) & dend_candidates)) {
                SupervoxelSet & set_a = m_segInfo->supervoxelList(a);
                if (set_a.isEmpty()) {
                    set_a.insert(a);
                }
                auto tmp = checkSupervoxelEdges(set_a, set_b);
                if (tmp > m_reliableMeanAffinity && tmp > max_aff) {
                    max_aff = tmp;
                    max_a = a;
                }
            }
            if (max_aff > m_reliableMeanAffinity) {
                new_candidates.insert(b);
                attached.insert(b);
                mergeGraph[max_a].insert(b);
                mergeGraph[b].insert(max_a);
                qDebug() << "Merge: " << b << max_a;
                SupervoxelSet & set_max_a = m_segInfo->supervoxelList(max_a);
                set_max_a.insert(b);
                keep_going = true;
            }
        }
        dend_candidates = new_candidates;
    }
    //return attached;
}

void Segmentation::postProcess()
{
    SupervoxelDict mergeGraph;
    matchAxons(mergeGraph);
    attachSmallSegments(mergeGraph);
}
