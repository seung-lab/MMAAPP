#include <QByteArray>
#include <QHash>
#include <QSet>
#include <QVector>
#include <cstdint>

typedef uint64_t id_type;
typedef uint64_t size_type;
typedef int32_t coord_type;
typedef double value_type;
typedef double count_type;

class MeanPlusEdge
{
public:
    MeanPlusEdge(const QByteArray & line);
    id_type p1;
    id_type p2;
    value_type sum;
    count_type num;
    id_type v1;
    id_type v2;
    value_type aff;
    count_type area;
};

typedef QList<QHash<id_type, MeanPlusEdge * > > RegionGraphArray;
typedef QList<QVector<coord_type > > BoundingBoxes;
typedef QVector<size_type > SupervoxelSizes;
typedef QVector<QVector<value_type > > SemanticInfo;
typedef QVector<id_type > SegmentArray;
typedef QHash<id_type, QHash<id_type, MeanPlusEdge * > > RegionGraph;
typedef QHash<id_type, QSet<id_type > > SupervoxelDict;
typedef QSet<id_type > SupervoxelSet;

class SupervoxelInfo
{
public:
    SupervoxelInfo();
    ~SupervoxelInfo();
    id_type maxSegId() { return m_maxSegId; };
    void loadSupervoxelInfo();
private:
    void readRegionGraph(const QString & filename);
    id_type m_maxSegId;
    RegionGraphArray m_regionGraph;
    BoundingBoxes m_boundingBoxes;
    SupervoxelSizes m_supervoxelSizes;
    SemanticInfo m_semanticInfo;
    SegmentArray m_segmentDict;
    SupervoxelSet m_boundarySupervoxels;
};
