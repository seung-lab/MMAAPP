#include "Segmentation.h"
#include <QString>
#include <QFile>
#include <QTextStream>
#include <QtDebug>

MeanPlusEdge::MeanPlusEdge(const QByteArray & line)
    :p1(0),p2(0),sum(0),num(0),v1(0),v2(0),aff(0),area(0)
{
    QList<QByteArray> data = line.split(' ');
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
}

void SupervoxelInfo::readRegionGraph(const QString & filename)
{
    QFile inputFile(filename);
    if (!inputFile.open(QIODevice::ReadOnly))
    {
        return;
    }
    QList<QByteArray> metaData = inputFile.readLine().trimmed().split(' ');
    qDebug() << metaData;
    m_maxSegId = metaData[0].toLongLong();
    uint64_t num_edges = metaData[2].toLongLong();

    m_regionGraph.reserve(m_maxSegId+1);
    for (unsigned int i = 0; i <= m_maxSegId; i++) {
        m_regionGraph.append(QHash<id_type, MeanPlusEdge * >());
    }
    for (unsigned int i = 0; i < num_edges; i++) {
        QByteArray line = inputFile.readLine().trimmed();
        MeanPlusEdge * edge = new MeanPlusEdge(line);
        //qDebug() << "p1: " << edge->p1 << "p2: " << edge->p2;
        m_regionGraph[edge->p1][edge->p2] = edge;
        m_regionGraph[edge->p2][edge->p1] = edge;
    }
    inputFile.close();
}
