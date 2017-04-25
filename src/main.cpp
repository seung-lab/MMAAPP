#include "Info.h"
#include "Segments.h"
#include <iostream>

int main(int argc, char* argv[])
{
    std::cout << "Start loading" << std::endl;
    SupervoxelInfo * svInfo = new SupervoxelInfo();
    SegmentInfo * segInfo = new SegmentInfo();
    auto svDict = segInfo->supervoxelDict();
    svInfo->agglomerate(svDict);
    Segmentation * segmentation = new Segmentation();
    std::cout << "Finish loading" << std::endl;
}
