#ifndef MMAAPP_OUTPUT_H
#define MMAAPP_OUTPUT_H
#include "Info.h"
#include <QHash>

class OutputRegionGraph
{
public:
    OutputRegionGraph(RegionGraph & regionGraph)
        :m_edges()
        ,m_regionGraph(regionGraph)
    {
    }
    ~OutputRegionGraph();
    void appendEdges(const SupervoxelDict & mergeGraph, value_type threshold);
    void updateRegionGraph(id_type maxSegId);
private:
    QString printEdge(const MeanPlusEdge * edge, value_type affinity);
    QHash<const MeanPlusEdge *, value_type>  m_edges;
    RegionGraph & m_regionGraph;
};
#endif
