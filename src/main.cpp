#include "Info.h"
#include "Segments.h"
#include <iostream>

int main(int argc, char* argv[])
{
    std::cout << "Start loading" << std::endl;
    SupervoxelInfo * svInfo = new SupervoxelInfo();
    SegmentInfo * segInfo = new SegmentInfo();
    svInfo->agglomerate(segInfo->supervoxelDict());
    Segmentation * segmentation = new Segmentation(svInfo, segInfo);
    segmentation->init();
    std::cout << "Finish loading" << std::endl;
    segmentation->postProcess();
}
