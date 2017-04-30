#include "Output.h"
#include <QQueue>
#include <QFile>
#include <QStringList>
#include <QtDebug>

void OutputRegionGraph::appendEdges(const SupervoxelDict & mergeGraph, value_type threshold)
{
    SupervoxelSet visited;
    qDebug() << "Merge" << mergeGraph.keys().size() << "segments at threshold:" << threshold ;
    foreach (auto p, mergeGraph.keys()) {
        QVector<id_type> segGroup;
        if (visited.contains(p)) {
            continue;
        }
        QQueue<id_type> queue;
        queue.enqueue(p);
        segGroup.append(p);
        while (!queue.isEmpty()) {
            auto root = queue.dequeue();
            if (visited.contains(root)) {
                continue;
            }
            visited.insert(root);
            foreach(auto neighbour, mergeGraph[root]) {
                if (!visited.contains(neighbour)) {
                    queue.enqueue(neighbour);
                    segGroup.append(neighbour);
                }
            }
        }
        for (int i = 0; i < segGroup.size(); i++) {
            for (int j = i+1; j < segGroup.size(); j++) {
                if (m_regionGraph[segGroup[i]].contains(segGroup[j])) {
                    m_edges[m_regionGraph[segGroup[i]][segGroup[j]]] = threshold;
                }
            }
        }
    }
}

QString OutputRegionGraph::printEdge(const MeanPlusEdge * edge, value_type meanAffinity)
{
    auto area = edge->num;
    auto sum_aff = edge->sum;
    if (meanAffinity > 0) {
        sum_aff = area*meanAffinity;
    } else {
        sum_aff *= 0.9995;
    }
    return QString("%1 %2 %3 %4 %5 %6 %7 %8\n").arg(edge->p1).arg(edge->p2).arg(sum_aff).arg(area).arg(edge->v1).arg(edge->v2).arg(edge->aff).arg(edge->area);
}

void OutputRegionGraph::updateRegionGraph(id_type maxSegId)
{
    QSet<const MeanPlusEdge *> visited;
    QStringList region_graph_output;
    foreach (auto edge, m_edges.keys()) {
        region_graph_output << printEdge(edge, m_edges[edge]);
        visited.insert(edge);
    }
    foreach (auto a, m_regionGraph.keys()) {
        foreach (auto b, m_regionGraph[a].keys()) {
            auto edge = m_regionGraph[a][b];
            if (visited.contains(edge)) {
                continue;
            }
            region_graph_output << printEdge(edge, -1);
            visited.insert(edge);
        }
    }
    QFile outputFile("axon.in");
    outputFile.open(QIODevice::WriteOnly | QIODevice::Text);
    if(!outputFile.isOpen()){
        qDebug() << "- Error, unable to open" << "axon.in" << "for output";
        return;
    }
    QTextStream outputStream(&outputFile);
    outputStream << maxSegId << " " << maxSegId+1 << " "<< region_graph_output.size() << "\n";
    foreach (auto s, region_graph_output) {
        outputStream << s;
    }
}
