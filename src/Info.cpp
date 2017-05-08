#include "Info.h"
#include <QString>
#include <QFile>
#include <QTextStream>
#include <QtDebug>
#include <QtConcurrent>

std::vector<MeanPlusEdge *> generate_edges(const QVector<QByteArray *> lines)
{
    std::vector<MeanPlusEdge * > edges;
    edges.reserve(lines.size());
    foreach (auto line, lines) {
        edges.push_back(new MeanPlusEdge(line));
    }
    return edges;
}

MeanPlusEdge::MeanPlusEdge(const QByteArray * line)
    :p1(0),p2(0),sum(0),num(0),v1(0),v2(0),aff(0),area(0)
{
    QList<QByteArray> data = line->trimmed().split(' ');
    p1 = data[0].toLongLong();
    p2 = data[1].toLongLong();
    sum = data[2].toDouble();
    num = data[3].toDouble();
    v1 = data[4].toLongLong();
    v2 = data[5].toLongLong();
    aff = data[6].toDouble();
    area = data[7].toDouble();
}

SupervoxelInfo::SupervoxelInfo()
    :m_maxSegId(0)
    ,m_regionGraph()
    ,m_boundingBoxes()
    ,m_supervoxelSizes()
    ,m_semanticInfo()
    ,m_segmentDict()
    ,m_boundarySupervoxels()
{
    loadSupervoxelInfo();
}

void SupervoxelInfo::loadSupervoxelInfo()
{
    readRegionGraph("rg_volume.in");
    readSupervoxelSizes("sv_volume.in");
    readBoundingBoxes("bbox_volume.in");
    readSemanticInfo("sem_volume.in");
    readMST("test_mst.in");
    findBoundaries(QVector<coord_type > ({1,1,1,2048,2048,256}));
}

void SupervoxelInfo::readBoundingBoxes(const QString & filename)
{
    qDebug() << "Reading the bounding boxes";
    QFile inputFile(filename);
    if (!inputFile.open(QIODevice::ReadOnly))
    {
        return;
    }
    m_boundingBoxes.reserve(m_maxSegId+1);
    for (unsigned int i = 0; i <= m_maxSegId; i++) {
        m_boundingBoxes.append(QVector<coord_type>());
    }
    while (!inputFile.atEnd()) {
        auto data = inputFile.readLine().trimmed().split(' ');
        auto segid = data[0].toLongLong();
        for (int i = 0; i < 6; i++){
            m_boundingBoxes[segid].append(data[i+1].toInt());
        }
    }
    inputFile.close();
}

void SupervoxelInfo::readSemanticInfo(const QString & filename)
{
    qDebug() << "Reading the semantic information";
    QFile inputFile(filename);
    if (!inputFile.open(QIODevice::ReadOnly))
    {
        return;
    }
    m_semanticInfo.reserve(m_maxSegId+1);
    for (unsigned int i = 0; i <= m_maxSegId; i++) {
        m_semanticInfo.append(QVector<value_type>(4,0));
    }
    while (!inputFile.atEnd()) {
        auto data = inputFile.readLine().trimmed().split(' ');
        auto segid = data[0].toLongLong();
        for (int i = 0; i < 4; i++){
            m_semanticInfo[segid][i] = data[i+1].toDouble();
        }
    }
    inputFile.close();
}

void SupervoxelInfo::readSupervoxelSizes(const QString & filename)
{
    qDebug() << "Reading the supervoxel sizes";
    QFile inputFile(filename);
    if (!inputFile.open(QIODevice::ReadOnly))
    {
        return;
    }
    m_supervoxelSizes.reserve(m_maxSegId+1);
    for (unsigned int i = 0; i <= m_maxSegId; i++) {
        m_supervoxelSizes.append(0);
    }
    while (!inputFile.atEnd()) {
        auto data = inputFile.readLine().trimmed().split(' ');
        m_supervoxelSizes[data[0].toLongLong()] = data[1].toLongLong();
    }
    inputFile.close();
}

void SupervoxelInfo::findBoundaries(const QVector<coord_type > & boundingBox)
{
    for (unsigned int i = 1; i <= m_maxSegId; i++) {
        const QVector<coord_type > & bbox = m_boundingBoxes[i];
        if (bbox.size() != 6) {
            continue;
        }
        for (int j = 0; j < 3; j++) {
            if (bbox[j] <= boundingBox[j]) {
                m_boundarySupervoxels.insert(i);
                break;
            }
            if (bbox[j+3] >= boundingBox[j+3]) {
                m_boundarySupervoxels.insert(i);
                break;
            }
        }
    }
}

void SupervoxelInfo::readMST(const QString & filename)
{
    qDebug() << "Reading the MST";
    QFile inputFile(filename);
    if (!inputFile.open(QIODevice::ReadOnly))
    {
        return;
    }
    m_segmentDict.reserve(m_maxSegId+1);
    for (unsigned int i = 0; i <= m_maxSegId; i++) {
        m_segmentDict.append(i);
    }
    while (!inputFile.atEnd()) {
        auto data = inputFile.readLine().trimmed().split(' ');
        auto child1 = data[0].toLongLong();
        auto child2 = data[1].toLongLong();
        auto parent = data[2].toLongLong();
        if (parent == child1)
        {
            m_segmentDict[child2] = parent;
        } else {
            m_segmentDict[child1] = parent;
        }
    }
    inputFile.close();
}

void SupervoxelInfo::readRegionGraph(const QString & filename)
{
    qDebug() << "Reading the region graph";
    QFile inputFile(filename);
    if (!inputFile.open(QIODevice::ReadOnly))
    {
        return;
    }
    QList<QByteArray> metaData = inputFile.readLine().trimmed().split(' ');
    m_maxSegId = metaData[0].toLongLong();
    uint64_t num_edges = metaData[2].toLongLong();
    int batch_size = 1000000;
    QVector<QByteArray *> rg_entries;
    QVector<QFuture<std::vector<MeanPlusEdge *> > > futures;
    m_regionGraph.reserve(m_maxSegId+1);
    for (unsigned int i = 0; i <= m_maxSegId; i++) {
        m_regionGraph.append(QHash<id_type, MeanPlusEdge * >());
    }
    for (unsigned int i = 0; i < num_edges; i++) {
        rg_entries << new QByteArray(inputFile.readLine());
        if ((i > 0) && (i % batch_size == 0)) {
            futures << QtConcurrent::run(generate_edges, QVector<QByteArray *>(rg_entries));
            rg_entries.clear();
        }
    }
    if (rg_entries.size() > 0) {
        futures << QtConcurrent::run(generate_edges, QVector<QByteArray *>(rg_entries));
    }
    foreach (auto f, futures) {
        auto edges = f.result();
        foreach (auto edge, edges) {
            m_regionGraph[edge->p1][edge->p2] = edge;
            m_regionGraph[edge->p2][edge->p1] = edge;
        }
    }
    inputFile.close();
}

void SupervoxelInfo::agglomerate(SupervoxelDict & supervoxelDict)
{
    qDebug() << "agglomerating";
    for (unsigned int c = 1; c <= m_maxSegId; c++) {
        auto p = m_segmentDict[c];
        if (p == c) {
            continue;
        }
        QSet<id_type > clst;
        clst.insert(c);
        while (m_segmentDict[p] != p) {
            clst.insert(p);
            p = m_segmentDict[p];
        }
        foreach(auto x, clst) {
            m_segmentDict[x] = p;
        }
        if (supervoxelDict.contains(p)) {
            supervoxelDict[p] |= clst;
        } else {
            clst.insert(p);
            supervoxelDict[p] = clst;
        }
    }
}

SegmentInfo::SegmentInfo()
    :m_regionGraph()
    ,m_supervoxelDict()
{
    readRegionGraph("new_rg.in");
}

void SegmentInfo::readRegionGraph(const QString & filename)
{
    qDebug() << "Reading the region graph after agglomeration";
    QFile inputFile(filename);
    if (!inputFile.open(QIODevice::ReadOnly))
    {
        return;
    }
    QList<QByteArray> metaData = inputFile.readLine().trimmed().split(' ');
    uint64_t num_edges = metaData[2].toLongLong();

    for (unsigned int i = 0; i < num_edges; i++) {
        QByteArray * line = new QByteArray(inputFile.readLine());
        MeanPlusEdge * edge = new MeanPlusEdge(line);
        if (!m_regionGraph.contains(edge->p1)) {
            m_regionGraph[edge->p1] = QHash<id_type, MeanPlusEdge * >();
        }
        if (!m_regionGraph.contains(edge->p2)) {
            m_regionGraph[edge->p2] = QHash<id_type, MeanPlusEdge * >();
        }
        //qDebug() << "p1: " << edge->p1 << "p2: " << edge->p2;
        m_regionGraph[edge->p1][edge->p2] = edge;
        m_regionGraph[edge->p2][edge->p1] = edge;
    }
    inputFile.close();
}
