#include "Segmentation.h"
#include <iostream>

int main(int argc, char* argv[])
{
    std::cout << "Start loading" << std::endl;
    SupervoxelInfo * svInfo = new SupervoxelInfo();
    SegmentInfo * segInfo = new SegmentInfo();
    auto svDict = segInfo->supervoxelDict();
    svInfo->agglomerate(svDict);
    std::cout << "Finish loading" << std::endl;
}
