#include "Info.h"
#include "Segments.h"
#include <QtDebug>

int main(int argc, char* argv[])
{
    qDebug() << "Start loading";
    SupervoxelInfo * svInfo = new SupervoxelInfo();
    SegmentInfo * segInfo = new SegmentInfo();
    svInfo->agglomerate(segInfo->supervoxelDict());
    Segmentation * segmentation = new Segmentation(svInfo, segInfo);
    segmentation->init();
    qDebug() << "Finish loading";
    segmentation->postProcess();
}
