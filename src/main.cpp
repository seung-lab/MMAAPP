#include "Info.h"
#include "Segments.h"
#include <QCommandLineParser>
#include <QtDebug>

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName("MMAAPP");
    QCoreApplication::setApplicationVersion("0.1");
    QCommandLineParser parser;
    parser.setApplicationDescription("Mostly Mean Affinity Agglomeration and Post-Processing");
    parser.addHelpOption();
    parser.addVersionOption();

    QCommandLineOption agglomerationThresholdOption(QStringList() << "a" << "agglomerationThreshold", "Agglomeration Threshold", "0.27");
    parser.addOption(agglomerationThresholdOption);

    QCommandLineOption postprocessThresholdOption(QStringList() << "p" << "postprocessThreshold", "Post Process Threshold", "0.17");
    parser.addOption(postprocessThresholdOption);

    parser.process(app);

    double agglomerationThreshold = 0.27;
    double postprocessThreshold = 0.17;

    if (parser.isSet(agglomerationThresholdOption)) {
        qDebug() << "agglomeration threshold" << parser.value(agglomerationThresholdOption);
        agglomerationThreshold = parser.value(agglomerationThresholdOption).toDouble();
    }

    if (parser.isSet(postprocessThresholdOption)) {
        qDebug() << "post process threshold" << parser.value(postprocessThresholdOption);
        postprocessThreshold = parser.value(postprocessThresholdOption).toDouble();
    }

    qDebug() << "Start loading";
    SupervoxelInfo * svInfo = new SupervoxelInfo();
    SegmentInfo * segInfo = new SegmentInfo();
    svInfo->agglomerate(segInfo->supervoxelDict());
    Segmentation * segmentation = new Segmentation(svInfo, segInfo, agglomerationThreshold, postprocessThreshold);
    segmentation->init();
    qDebug() << "Finish loading";
    segmentation->postProcess();
}
