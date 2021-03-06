#include <QCoreApplication>
#include <QtCore>
#include <QDebug>
#include "storagemodel.h"

int bearRiverAnalysis();
int test();
int run();
int validate();
int logan();
int pataha();
int testStats();
int runXYZ();
int soil();
int soilFixer();
int soilFixerLoop();
int bridgeCreekGW();
int addStorage();

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QDateTime startTime = QDateTime::currentDateTime();

    //test();
    //validate();
    //logan();
    //testStats();
    //runXYZ();
    //soil();
    //soilFixer();
    //soilFixerLoop();
    //bridgeCreekGW();
    //bearRiverAnalysis();
    addStorage();

    QDateTime endTime = QDateTime::currentDateTime();

    qDebug()<<"done"<<startTime.secsTo(endTime)<<startTime.secsTo(endTime)/60.0;

    return a.exec();
}

int bearRiverAnalysis()
{
    // Run suface storage for 0.05, 0.1, 0.25, 0.5, 0.75, and 1.0 of capacity
    QString dirPath = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010204/HUC12";
    QDir dir(dirPath);
    QFileInfoList fiList = dir.entryInfoList();
    const char *statPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/csv/volume_estimates2.tif";
    qDebug()<<"starting loop";

    for (int i=2; i<fiList.length(); i++)
    {
        qDebug()<<fiList[i].baseName();
        std::string basePath = dirPath.toStdString() + "/" + fiList[i].baseName().toStdString();
        QString runDir = "/03_out_05";
        QString qBasePath = dirPath + fiList[i].baseName();
        if (!QDir(qBasePath + runDir).exists())
        {
            QDir().mkdir(qBasePath + runDir);
        }
        std::string bratPath = basePath + "/01_shpIn/brat_cap_20170224.shp";
        std::string outDir = basePath + runDir.toStdString();
        std::string demPath = basePath + "/02_rasIn/dem_vbfac.tif";
        std::string filPath = basePath + "/02_rasIn/fil_vbfac.tif";
        std::string facPath = basePath + "/02_rasIn/fac_01km_vbfac.tif";
        std::string fdirPath = basePath + "/02_rasIn/fdird_vbfac.tif";
        //qDebug()<<"names set"<<bratPath.toStdString().c_str()<<outDir.toStdString().c_str()<<demPath.toStdString().c_str()<<facPath.toStdString().c_str()<<fdirPath.toStdString().c_str();

        Raster raster;
        raster.setNoData(demPath.c_str(), -9999.0, 50, 5000);
        //raster.setNoData(filPath.c_str(), -9999.0, 50, 5000);
        raster.setNoData(fdirPath.c_str(), 0, 1, 200);
        raster.setNoData(facPath.c_str(), -1, 1, 2);
        //qDebug()<<"running"<<i+1<<"of"<<fiList.length()<<basePath;

        //value of 3 = use recursive HAND algorithm (greatly increases speed)
        StorageModel model(bratPath.c_str(), outDir.c_str(), demPath.c_str(), fdirPath.c_str(), facPath.c_str(), 0.05, 3, statPath);

        //value of 3 = the number of dams placed on a reach is randomly selected from dam complex distribution, after all reaches processed additional complexes added if below desired capacity
        model.run(3);
        //qDebug()<<basePath<<"done";
    }

    return 0;
}

int test()
{

    const char *bratPath = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010101/HUC12/160101010801/01_shpIn/brat_cap_20170224.shp";
    const char *demPath = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010101/HUC12/160101010801/02_rasIn/dem_vbfac.tif";
    const char *filPath = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010101/HUC12/160101010801/02_rasIn/fil_vbfac.tif";
    const char *fdirPath = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010101/HUC12/160101010801/02_rasIn/fdird_vbfac.tif";
    const char *facPath = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010101/HUC12/160101010801/02_rasIn/fac_01km_vbfac.tif";
    const char *outDir = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010101/HUC12/160101010801/03_out";
    const char *damsPath = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010101/HUC12/160101010801/01_shpIn/dams_brat_join.shp";
    const char *csvPath = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010101/HUC12/160101010801/03_out/comparison.csv";
    const char *statPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/csv/volume_estimates2.tif";

//    StorageModel model(bratPath, outDir, filPath, fdirPath, facPath, 0.5, 3);
//    model.run(2);

    //value of 3 = use recursive HAND algorithm (greatly increases speed)
    StorageModel model(bratPath, outDir, demPath, fdirPath, facPath, 1.0, 3, statPath);

    //value of 3 = the number of dams placed on a reach is randomly selected from dam complex distribution, after all reaches processed additional complexes added if below desired capacity
    model.run(3);

//    const char *soil = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/Soil/UT/gssurgo_ut10m_utm12.tif";
//    const char *points = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/Soil/UT/grid_points/ut_points.txt";
//    Raster raster;
//    raster.toXYZ(soil, points);



    return 0;
}

int run()
{
    const char *shpIn = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn/BRAT_TempleFk_WS.shp";
    const char *shpInDir = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn";
    const char *shpInLogan = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn/lc_brat.shp";
    const char *shpOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut";
    const char *vbtm = "TempleFork_VB";
    const char *vbtm200buf = "TempleFork_200Buf";
    const char *demIn10m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/templefk_10m_ws.tif";
    const char *demIn10m_clip = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/templefk_10m_vb.tif";
    const char *fdir10m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fdir10m.tif";
    const char *fac10m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fac_2500.tif";
    const char *fil10m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fil10m.tif";
    const char *fdir10m_clip = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fdir10m_clip.tif";
    const char *fac10m_clip = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fac_2500_clip.tif";
    const char *fil10mclip = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fil10m_clip.tif";
    const char *demIn1m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fme450000.tif";
    const char *fil1m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fil_1m.tif";
    const char *fac1m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fac_200000.tif";
    const char *fdir1m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fdir1m.tif";
    const char *demIn1m_clip = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fme450000_clip.tif";
    const char *fil1m_clip = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fil1m_clip.tif";
    const char *fac1m_clip = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fac_200000_clip.tif";
    const char *fdir1m_clip = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fdir1m_clip.tif";
    const char *demInLogan = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/lc_fil_vb_utm12.tif";
    const char *fdirLogan = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/lc_fdir_vb_utm12.tif";
    const char *facLogan = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/lc_fac2000_vb_utm12.tif";
    const char *exDams = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/ValidationLayers/Dams_BRAT_join5_UTM12N_pondArea.shp";
    const char *csvOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut/comparison10m_dist.csv";
    const char *heightOut10m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut/hand10m.tif";
    const char *test10m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut/test.tif";

    StorageModel model(shpIn, shpOut, demIn10m, fdir10m, fac10m, 0.5, 1);
    model.run(1);

    //StorageModel model(shpIn, shpOut, demIn10m_clip, fdir10m_clip, fac10m_clip, 0.5);
    //model.run();
    //model.runFromPoints(exDams, csvOut);



    return 0;
}

int validate()
{
    //Ogden NF 1m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_1m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_1m/02_rasIn/fil1m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_1m/02_rasIn/fdir1m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_1m/02_rasIn/fac1m_500000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_1m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_1m/01_shpIn/dams_brat_join.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_1m/03_out/comparison.csv";

    //Ogden NF 10m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_10m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_10m/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_10m/02_rasIn/fdir10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_10m/02_rasIn/fac10m_5000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_10m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_10m/01_shpIn/dams_brat_join.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/OgdenNF_10m/03_out/comparison.csv";

    //Ogden NF 10m inputs - preserve dam location
//    const char *bratPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/OgdenNF_10m_DamLoc/01_shpIn/brat.shp";
//    const char *demPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/OgdenNF_10m_DamLoc/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/OgdenNF_10m_DamLoc/02_rasIn/fdir10m_vb.tif";
//    const char *facPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/OgdenNF_10m_DamLoc/02_rasIn/fac10m_5000_vb.tif";
//    const char *outDir = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/OgdenNF_10m_DamLoc/03_out";
//    const char *damsPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/OgdenNF_10m_DamLoc/01_shpIn/dams_brat_join.shp";
//    const char *csvPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/OgdenNF_10m_DamLoc/03_out/comparison.csv";

    //Santa Clara 1m inputs
//    const char *bratPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/SantaClara_1m/01_shpIn/brat.shp";
//    const char *demPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/SantaClara_1m/02_rasIn/fil1m_vb.tif";
//    const char *fdirPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/SantaClara_1m/02_rasIn/fdir1m_vb.tif";
//    const char *facPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/SantaClara_1m/02_rasIn/fac1m_4000000_vb.tif";
//    const char *outDir = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/SantaClara_1m/03_out";
//    const char *damsPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/SantaClara_1m/01_shpIn/dams_20160629.shp";
//    const char *csvPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/SantaClara_1m/03_out/comparison.csv";

    //Farmington 3m inputs
//    const char *bratPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Farmington_3m/01_shpIn/brat.shp";
//    const char *demPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Farmington_3m/02_rasIn/fil3m_vb.tif";
//    const char *fdirPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Farmington_3m/02_rasIn/fdir3m_vb.tif";
//    const char *facPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Farmington_3m/02_rasIn/fac3m_10000_vb.tif";
//    const char *outDir = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Farmington_3m/03_out";
//    const char *damsPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Farmington_3m/01_shpIn/dams_brat_join.shp";
//    const char *csvPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Farmington_3m/03_out/comparison.csv";

    //Logan HUC 8 inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/02_rasIn/fdir10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/02_rasIn/fac10m_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/01_shpIn/Dams_LoganHUC8_JoinBRAT.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/03_out/comparison.csv";

    //Logan HUC 8 5m inputs
//    const char *bratPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8_5m/01_shpIn/brat.shp";
//    const char *demPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8_5m/02_rasIn/fil_vb.tif";
//    const char *fdirPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8_5m/02_rasIn/fdir_vb.tif";
//    const char *facPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8_5m/02_rasIn/fac_100000_vbq.tif";
//    const char *outDir = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8_5m/03_out";
//    const char *damsPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8_5m/01_shpIn/Dams_LoganHUC8_JoinBRAT.shp";
//    const char *csvPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8_5m/03_out/comparison.csv";

    //Duchesne HUC 8 inputs
//    const char *bratPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Duchesne_HUC8/01_shpIn/brat.shp";
//    const char *demPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Duchesne_HUC8/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Duchesne_HUC8/02_rasIn/fdir10m_vb.tif";
//    const char *facPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Duchesne_HUC8/02_rasIn/fac10m_1700_vb.tif";
//    const char *outDir = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Duchesne_HUC8/03_out";
//    const char *damsPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Duchesne_HUC8/01_shpIn/Dams_BRAT_join_20160624.shp";
//    const char *csvPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Duchesne_HUC8/03_out/comparison.csv";

      //Temple Fork 10m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/01_shpIn/BRAT_TempleFk_WS.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/02_rasIn/fdir10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/02_rasIn/fac_2500_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/01_shpIn/DamArea_BRAT_joined.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/03_out/comparison.csv";

//    //Temple Fork 10m inputs volume comparison
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/01_shpIn/BRAT_TempleFk_WS.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/02_rasIn/fdird10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/02_rasIn/fac_2500_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/01_shpIn/input.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/03_out/comparison.csv";

    //Temple Fork 10m inputs, MODFLOW test
//    const char *bratPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/TempleFork_10m_MODFLOW/01_shpIn/BRAT_TempleFk_WS.shp";
//    const char *demPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/TempleFork_10m_MODFLOW/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/TempleFork_10m_MODFLOW/02_rasIn/fdir10m_vb.tif";
//    const char *facPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/TempleFork_10m_MODFLOW/02_rasIn/fac_1300_vb.tif";
//    const char *outDir = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/TempleFork_10m_MODFLOW/03_out";
//    const char *damsPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/TempleFork_10m_MODFLOW/01_shpIn/DamArea_BRAT_joined.shp";
//    const char *csvPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/TempleFork_10m_MODFLOW/03_out/comparison.csv";

    //Temple Fork 1m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m/01_shpIn/BRAT_TempleFk_WS.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m/02_rasIn/fil1m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m/02_rasIn/fdir1m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m/02_rasIn/fac_200000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m/01_shpIn/DamArea_BRAT_joined.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m/03_out/comparison.csv";

    //Temple Fork 1m inputs volume comparison
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m_vol/01_shpIn/BRAT_TempleFk_WS.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m_vol/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m_vol/02_rasIn/fdird1m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m_vol/02_rasIn/fac_200000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m_vol/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m_vol/01_shpIn/input.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m_vol/03_out/comparison.csv";

    //Price (Gordon) 2m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_2m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_2m/02_rasIn/fil_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_2m/02_rasIn/fdir_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_2m/02_rasIn/fac_100000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_2m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_2m/01_shpIn/dams_join";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_2m/03_out/comparison.csv";

    //Price (Gordon) 2m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_2m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_2m/02_rasIn/fil_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_2m/02_rasIn/fdir_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_2m/02_rasIn/fac_1000000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_2m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_2m/01_shpIn/dam_loc";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_2m/03_out/comparison.csv";

    //Price (Gordon) 5m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_5m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_5m/02_rasIn/fil_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_5m/02_rasIn/fdir_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_5m/02_rasIn/fac_100000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_5m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_5m/01_shpIn/dam_loc";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Price_Gordon_5m/03_out/comparison.csv";

    //Heber 2m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Heber_2m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Heber_2m/02_rasIn/fil_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Heber_2m/02_rasIn/fdir_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Heber_2m/02_rasIn/fac_60000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Heber_2m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Heber_2m/01_shpIn/dams_join";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Heber_2m/03_out/comparison.csv";

    //Test recursive HAND algorithm
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/01_shpIn/BRAT_TempleFk_WS.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/02_rasIn/fdird10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/02_rasIn/fac_2500_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/01_shpIn/DamArea_BRAT_joined.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/03_out/comparison.csv";

    //Recursive with Temple Fork 1m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/01_shpIn/BRAT_TempleFk_WS.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/02_rasIn/fdird1m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/02_rasIn/fac_200000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/01_shpIn/DamArea_BRAT_joined.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/03_out/comparison.csv";

    //Recursive with Ogden NF 1m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/02_rasIn/fdird1m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/02_rasIn/fac1m_500000.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/01_shpIn/dams_brat_join.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/03_out/comparison.csv";

    //Recursive with Ogden NF 10m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/02_rasIn/fdird10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/02_rasIn/fac10m_5000.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/01_shpIn/dams_brat_join.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/03_out/comparison.csv";

    //Recursive with Logan HUC 8 inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan100/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan100/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan100/02_rasIn/fdir10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan100/02_rasIn/fac10m_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan100/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan100/01_shpIn/dams_20160629.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan100/03_out/comparison.csv";

    //Beaver Creek 10m Inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/05_HUC12/BeaverCreek/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/05_HUC12/BeaverCreek/02_rasIn/fil.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/05_HUC12/BeaverCreek/02_rasIn/fdir.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/05_HUC12/BeaverCreek/02_rasIn/fac.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/05_HUC12/BeaverCreek/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/05_HUC12/BeaverCreek/01_shpIn/Dams_LoganHUC8_JoinBRAT.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/05_HUC12/BeaverCreek/03_out/comparison.csv";

    //Bridge Creek 10m inputs volume comparison
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/01_shpIn/BRAT_bcHUC10.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/02_rasIn/fdird_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/02_rasIn/fac_50000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/01_shpIn/points_snap.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/03_out/comparison.csv";

    //Bridge Creek 1m inputs volume comparison
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/01_shpIn/BRAT_bcHUC10.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/02_rasIn/fdird_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/02_rasIn/fac_1mil_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/01_shpIn/points_snap.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/03_out/comparison.csv";

    //Curtis Creek 10m inputs (whole Logan HUC 8) volume comparison
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/02_rasIn/fdird10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/02_rasIn/fac10m_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/01_shpIn/points_snap.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/03_out/comparison.csv";

//    //Curtis Creek 10m for groundwater validation
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030203/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030203/02_rasIn/fil_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030203/02_rasIn/fdir.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030203/02_rasIn/fac_5.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030203/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030203/01_shpIn/points_snap.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030203/03_out/comparison.csv";

    //Crash when running this setup in batch mode, figuring out problems
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/1601020302/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030308/02_rasIn/fil.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030308/02_rasIn/fdir.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030308/02_rasIn/fac_5.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030308/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030308/01_shpIn/DamArea_BRAT_joined.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/160102030308/03_out/comparison.csv";

    //Test algorithm limiting pond size
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/02_rasIn/fdird_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/02_rasIn/fac_1km_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/01_shpIn/dams_20160629.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203/03_out/comparison.csv";

    //Test recursive HAND algorithm with Laptop
    const char *bratPath = "C:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/01_shpIn/BRAT_TempleFk_WS.shp";
    const char *demPath = "C:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/02_rasIn/dem_vb.tif";
    const char *fdirPath = "C:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/02_rasIn/fdird10m_vb.tif";
    const char *facPath = "C:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/02_rasIn/fac_2500_vb.tif";
    const char *outDir = "C:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/03_out";
    const char *damsPath = "C:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/01_shpIn/DamArea_BRAT_joined.shp";
    const char *csvPath = "C:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/03_out/comparison.csv";


//    Raster raster;
//    raster.setNoData(demPath, -9999.0, 50, 5000);
//    raster.setNoData(fdirPath, 0, 1, 200);
//    raster.setNoData(facPath, -1, 1, 2);

    const char *txtPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/csv/volume_estimates2.txt";
    const char *statPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/csv/volume_estimates2.tif";

    Raster raster;
    qDebug()<<"raster initialized";
    raster.fromXYZ(statPath, txtPath, 79, 79, -9999.0);
    qDebug()<<"raster done";

/*
 * ******************************************************************************************
 * ************************* RUN TYPES WHEN USING BRAT (Initialization types)****************
 * ******************************************************************************************
 * type = 1: vector based (pond determined with maximum extent polygon)
 * type = 2: raster based (pond determined with flow direction algebra)
 * type = 3: raster based, backwards and recursive flow direction algorithm (faster)
 * ******************************************************************************************
 */

    //Initialize surface storage model with statistical correction
    //StorageModel model(bratPath, outDir, demPath, fdirPath, facPath, 1.0, 3, statPath);

    //Initialize surface storage model without statistical correction
//    StorageModel model(bratPath, outDir, demPath, fdirPath, facPath, 0.5, 3);

 /*
 * *****************************************************************************************************
 * ************************* RUN TYPES WHEN USING EXISTING POINTS***************************************
 * *****************************************************************************************************
 * type = 1 (default): use existing dam locations (copy created)
 * type = 2: existing dam locations, maintain locations (do not move to flow accumulation)
 * type = 3: use existing dam points with heights
 * type = 4: use existing dam points with heights, maintain locations (do not move to flow accumulation)
 * *****************************************************************************************************
 */

    /*
    * *****************************************************************************************************
    * ************************* RUN TYPES FOR POINT PLACEMENT**********************************************
    * *****************************************************************************************************
    * type = 1 (default): distribute dams evenly on each reach according to BRAT capacity
    * type = 2: add a dam complex to best habitats first, then continue to lesser habitats
    * type = 3: same as type 2 but updated to get actual 100% capacity
    * *****************************************************************************************************
    */

    // Run surface storage model
    //model.runFromPoints(damsPath, csvPath);
    //model.runFromPoints(damsPath, csvPath, 1);
    //model.runFromPointsWithHeights(damsPath,csvPath, 3);
    //Value of 3 is for full 100% capacity scenario while still placing dams in complexes
    //model.run(3);

    return 0;
}

int logan()
{
    QString dirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010203";
    QDir dir(dirPath);
    QFileInfoList fiList = dir.entryInfoList();
    qDebug()<<"starting loop";

    for (int i=2; i<fiList.length(); i++)
    {
        qDebug()<<fiList[i].baseName();
        std::string basePath = dirPath.toStdString() + "/" + fiList[i].baseName().toStdString();
        QString qBasePath = dirPath + fiList[i].baseName();
        if (!QDir(qBasePath + "/03_out").exists())
        {
            QDir().mkdir(qBasePath + "/03_out");
        }
        std::string bratPath = basePath + "/01_shpIn/brat.shp";
        std::string outDir = basePath + "/03_out";
        std::string demPath = basePath + "/02_rasIn/fil.tif";
        std::string facPath = basePath + "/02_rasIn/fac_5.tif";
        std::string fdirPath = basePath + "/02_rasIn/fdir.tif";
        //qDebug()<<"names set"<<bratPath.toStdString().c_str()<<outDir.toStdString().c_str()<<demPath.toStdString().c_str()<<facPath.toStdString().c_str()<<fdirPath.toStdString().c_str();

        Raster raster;
        raster.setNoData(demPath.c_str(), -9999.0, 50, 5000);
        raster.setNoData(fdirPath.c_str(), 0, 1, 200);
        raster.setNoData(facPath.c_str(), -1, 1, 2);
        //qDebug()<<"running"<<i+1<<"of"<<fiList.length()<<basePath;
        StorageModel model(bratPath.c_str(), outDir.c_str(), demPath.c_str(), fdirPath.c_str(), facPath.c_str(), 0.5, 3);
        model.run(2);
        //qDebug()<<basePath<<"done";
    }
}

int pataha()
{
    const char *shpIn = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelIn/Pataha_BRAT.shp";
    const char *shpOut = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelOut";
    const char *demIn = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelIn/pataha_10m.tif";
    const char *fdir = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelIn/fdir_10m.tif";
    const char *fac = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelIn/fac_10m.tif";

    StorageModel model(shpIn, shpOut, demIn, fdir, fac, 0.5, 1);
    model.run(1);

    return 0;
}

int testStats()
{
    QString fn = "C:/Users/khafe/Desktop/teststat.txt";
    Statistics normDist(Random::randomSeries(100000, RDT_norm, 0.93, 0.17), RDT_norm);
    QVector<double> data = normDist.getData();

    QFile file(fn);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outStream(&file);

    for (int i=0; i<data.length(); i++)
    {
        outStream<<data[i]<<"\n";
    }

    file.close();

    return 0;
}

int runXYZ()
{
    const char *elev = "E:/etal/Projects/NonLoc/BeaverModeling/08_MODFLOW/tf_10m/01_Inputs/fil10m.tif";
    const char *elev_xyz = "E:/etal/Projects/NonLoc/BeaverModeling/08_MODFLOW/tf_10m/01_Inputs/dem.xyz";
    const char *head = "E:/etal/Projects/NonLoc/BeaverModeling/08_MODFLOW/tf_10m/01_Inputs/head.tif";
    const char *head_xyz = "E:/etal/Projects/NonLoc/BeaverModeling/08_MODFLOW/tf_10m/01_Inputs/head.xyz";

    Raster raster;
    raster.toXYZ(elev, elev_xyz);
    raster.toXYZ(head, head_xyz);

    return 0;
}

int soil()
{
    qDebug()<<"setting up";
    const char *mukeyUT = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/Soil/gSSURGO/Utah/TIF/ut_mapunit_10m.tif";
    const char *ksatUT = "C:/etal/Projects/USA/gSSURGO/Utah/ut_ksat_10m.tif";
    const char *kvUT = "C:/etal/Projects/USA/gSSURGO/Utah/ut_kv_10m.tif";
    const char *porUT = "C:/etal/Projects/USA/gSSURGO/Utah/ut_por_10m.tif";
    const char *ksatUTout = "F:/01_etal/GIS_Data/USA/Soil/gSSURGO/Utah/ut_ksat_10m_fix.tif";
    const char *kvUTout = "F:/01_etal/GIS_Data/USA/Soil/gSSURGO/Utah/ut_kv_10m_fix.tif";
    const char *porUTout = "F:/01_etal/GIS_Data/USA/Soil/gSSURGO/Utah/ut_por_10m_fix.tif";

    const char *mukeyID = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/Soil/gSSURGO/Idaho/TIF/id_mapunit_10m.tif";
    const char *ksatID = "C:/etal/Projects/USA/gSSURGO/Idaho/id_ksat_10m.tif";
    const char *kvID = "C:/etal/Projects/USA/gSSURGO/Idaho/id_kv_10m.tif";
    const char *porID = "C:/etal/Projects/USA/gSSURGO/Idaho/id_por_10m.tif";
    const char *ksatIDout = "F:/01_etal/GIS_Data/USA/Soil/gSSURGO/Idaho/id_ksat_10m_fix.tif";
    const char *kvIDout = "F:/01_etal/GIS_Data/USA/Soil/gSSURGO/Idaho/id_kv_10m_fix.tif";
    const char *porIDout = "F:/01_etal/GIS_Data/USA/Soil/gSSURGO/Idaho/id_por_10m_fix.tif";

    const char *mukeyOR = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/Soil/gSSURGO/Oregon/TIF/or_mapunit_10m.tif";
    const char *ksatOR = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/Soil/gSSURGO/Oregon/TIF/or_ksat_10m";
    const char *kvOR = "C:/etal/Projects/USA/gSSURGO/Oregon/or_kv_10m.tif";
    const char *porOR = "C:/etal/Projects/USA/gSSURGO/Oregon/or_por_10m.tif";
    const char *ksatORout = "C:/etal/Projects/USA/gSSURGO/Oregon/or_ksat_10m_fix.tif";
    const char *kvORout = "C:/etal/Projects/USA/gSSURGO/Oregon/or_kv_10m_fix.tif";
    const char *porORout = "C:/etal/Projects/USA/gSSURGO/Oregon/or_por_10m_fix.tif";

    qDebug()<<"start";
    Raster raster;
    raster.adjustSoil(mukeyID, porID, porIDout);
    qDebug()<<"done";
}

int soilFixer()
{
    //Bear River
    const char *dem = "G:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010204/02_rasIn/dem_vbfac.tif";
    const char *huc8 = "G:/01_etal/GIS_Data/USA/NHD/Utah/Watersheds/BearRiverNHD/HUC8.tif";
    const char *huc12 = "G:/01_etal/GIS_Data/USA/NHD/Utah/Watersheds/BearRiverNHD/HUC12.tif";
    const char *soil = "G:/01_etal/GIS_Data/USA/Soil/SSURGO/Utah/BearRiver/EntireDrianage/kv_clip.tif";
    const char *out = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010204/02_rasIn/kv_vbfac.tif";

    double kvmax = 229.0;
    double ksatmax = 345.0;
    double pormax = 100.0;

    //Bridge Creek
//    const char *dem = "F:/01_etal/GIS_Data/USA/DEM/NED_10m/Oregon/JohnDay/LJD_10m_vb.tif";
//    dem = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/02_rasIn/dem_vb.tif";
//    const char *huc8 = "F:/01_etal/GIS_Data/USA/NHD/USA/Oregon/JohnDay/LJD_HUC8_utm10.tif";
//    const char *huc12 = "F:/01_etal/GIS_Data/USA/NHD/USA/Oregon/JohnDay/LJD_HUC12_utm10.tif";
//    const char *soil = "F:/01_etal/GIS_Data/USA/Soil/gSSURGO/Oregon/JohnDay/LJD_ksat_vb.tif";
//    const char *out = "F:/01_etal/GIS_Data/USA/Soil/gSSURGO/Oregon/JohnDay/BridgeCreek/ksat_10m_vb.tif";

//    double kvmax = 422.0;
//    double ksatmax = 423.0;
//    double pormax = 100.0;

    Raster_BeaverPond fixer;
    qDebug()<<"starting soil raster function";
    fixer.soilRasterCreation(dem, huc8, huc12, soil, out, kvmax);
    //fixer.soilRasterCreation_table(dem, huc8, huc12, out, 100.0);
    qDebug()<<"done";
}

int soilFixerLoop()
{
    QStringList names;
    names<<"por"<<"fc"<<"ksat"<<"kv";
    QStringList hucs;
    hucs<<"16010101"<<"16010102"<<"16010201"<<"16010202"<<"16010203"<<"16010204";
    double max[4] = {100.0, 100.0, 345.0, 229.0};

    //const char *dem = "G:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/16010204/02_rasIn/dem_vbfac.tif";
    const char *huc8 = "G:/01_etal/GIS_Data/USA/NHD/Utah/Watersheds/BearRiverNHD/HUC8.tif";
    const char *huc12 = "G:/01_etal/GIS_Data/USA/NHD/Utah/Watersheds/BearRiverNHD/HUC12.tif";

    QString dembase = "G:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/";
    QString soilbase = "G:/01_etal/GIS_Data/USA/Soil/SSURGO/Utah/BearRiver/EntireDrianage/";
    QString outbase = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/";

    Raster_BeaverPond fixer;

    for (int i=0; i<hucs.length()-1; i++)
    {
        QString dem = dembase + hucs[i] + "/02_rasIn/dem_vbfac.tif";

        for (int j=0; j<names.length(); j++)
        {
            QString soil = soilbase + names[j] + "_clip.tif";
            QString out = outbase + hucs[i] + "/02_rasIn/" + names[j] + "_vbfac.tif";
            QFileInfo finfo(out);
            if (!finfo.exists())
            {
                fixer.soilRasterCreation(dem.toStdString().c_str(), huc8, huc12, soil.toStdString().c_str(), out.toStdString().c_str(), max[j]);
            }
            qDebug()<<hucs[i]<<names[j]<<"done"<<i*names.length()+j+1<<"of"<<names.length()*(hucs.length()-1);
        }
    }
}

int bridgeCreekGW()
{
    //Bridge Creek 10m for groundwater validation
    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BridgeCreek/10m/2015/HUC12/170702040303/01_shpIn/BRAT_bcHUC10.shp";
    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BridgeCreek/10m/2015/HUC12/170702040303/02_rasIn/fil_vb.tif";
    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BridgeCreek/10m/2015/HUC12/170702040303/02_rasIn/fdir_vb.tif";
    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BridgeCreek/10m/2015/HUC12/170702040303/02_rasIn/fac_50000_vb.tif";
    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BridgeCreek/10m/2015/HUC12/170702040303/03_out";
    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BridgeCreek/10m/2015/HUC12/170702040303/01_shpIn/points_snap_intact.shp";
    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BridgeCreek/10m/2015/HUC12/170702040303/03_out/comparison.csv";

    Raster raster;
    raster.setNoData(demPath, -9999.0, 50, 5000);
    raster.setNoData(fdirPath, 0, 1, 200);
    raster.setNoData(facPath, -1, 1, 2);

    const char *txtPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/csv/volume_estimates.txt";
    const char *statPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/csv/volume_estimates.tif";

//    Raster raster;
//    qDebug()<<"raster initialized";
//    raster.fromXYZ(statPath, txtPath, 79, 79, -9999.0);
//    qDebug()<<"raster done";

/*
 * ******************************************************************************************
 * ************************* RUN TYPES WHEN USING BRAT (Initialization types)****************
 * ******************************************************************************************
 * type = 1: vector based (pond determined with maximum extent polygon)
 * type = 2: raster based (pond determined with flow direction algebra)
 * type = 3: raster based, backwards and recursive flow direction algorithm (faster)
 * ******************************************************************************************
 */

    //Initialize surface storage model with statistical correction
    //StorageModel model(bratPath, outDir, demPath, fdirPath, facPath, 1.0, 3, statPath);
    //Initialize surface storage model without statistical correction
    StorageModel model(bratPath, outDir, demPath, fdirPath, facPath, 1.0, 3);

 /*
 * *****************************************************************************************************
 * ************************* RUN TYPES WHEN USING EXISTING POINTS***************************************
 * *****************************************************************************************************
 * type = 1 (default): use existing dam locations (copy created)
 * type = 2: existing dam locations, maintain locations (do not move to flow accumulation)
 * type = 3: use existing dam points with heights
 * type = 4: use existing dam points with heights, maintain locations (do not move to flow accumulation)
 * *****************************************************************************************************
 */

    /*
    * *****************************************************************************************************
    * ************************* RUN TYPES FOR POINT PLACEMENT**********************************************
    * *****************************************************************************************************
    * type = 1 (default): distribute dams evenly on each reach according to BRAT capacity
    * type = 2: add a dam complex to best habitats first, then continue to lesser habitats
    * *****************************************************************************************************
    */

    // Run surface storage model
    //model.runFromPoints(damsPath, csvPath);
    //model.runFromPoints(damsPath, csvPath, 3);
    model.runFromPointsWithHeights(damsPath,csvPath, 3);
    //model.run(2);

    return 0;
}

int addStorage()
{
    Raster raster;

    const char *gw = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/EntireBasin/out_05/hdch_mid_m3.tif";
    const char *sw = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/EntireBasin/out_05/depMid_m.tif";
    const char *out = "E:/konrad/Projects/Modeling/BeaverWaterStorage/wrk_Data/AnalysisRuns/BearRiverHUC8/EntireBasin/out_05/totalMid.tif";

    raster.add(gw, sw, out);

    return 0;
}
