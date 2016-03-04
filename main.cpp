#include <QCoreApplication>
#include <QtCore>
#include <QDebug>
#include "storagemodel.h"

int test();
int run();

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QDateTime startTime = QDateTime::currentDateTime();

    test();

    QDateTime endTime = QDateTime::currentDateTime();

    qDebug()<<"done"<<startTime.secsTo(endTime)<<startTime.secsTo(endTime)/60.0;

    return a.exec();
}

int test()
{
    const char *shpIn = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn/BRAT_TempleFk_WS.shp";
    const char *shpInLogan = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn/lc_brat.shp";
    const char *shpOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut";
    const char *demIn10m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/templefk_10m_ws.tif";
    const char *demIn1m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fme450000.tif";
    const char *demInLogan = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Raster/lc_dem_10m.tif";
    const char *exDams = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/ValidationLayers/Dams_BRAT_join5_UTM12N_pondArea.shp";
    const char *csvOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut/comparison10m_dist.csv";

    StorageModel model(shpIn, shpOut, demIn10m, 0.5);
    //model.run();
    model.runFromPoints(exDams, csvOut);

    return 0;
}

int run()
{
    GDALAllRegister();
    OGRRegisterAll();

//    const char *shpIn = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn";
//    const char *shpOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut";

//    const char *shpIn = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn";
//    const char *shpOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut";

//    const char *demIn = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/templefk_10m_ws.tif";
//    const char *depOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/ponddepth_10m.tif";
//    const char *freqOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_10m.tif";
//    const char *csv = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_10m.csv";

//    const char *demIn = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/templefk_10m_ws.tif";
//    const char *depOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/ponddepth_10m.tif";
//    const char *freqOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_10m.tif";
//    const char *csv = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_10m.csv";

//    const char *demIn = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fme450000.tif";
//    const char *depOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/ponddepth_1m.tif";
//    const char *freqOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_1m.tif";
//    const char *csv = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_1m.csv";

    const char *shpIn = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Shapefiles";
    const char *shpOut = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Results/Shapefile";

    const char *demIn = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Raster/lc_dem_10m.tif";
    const char *depOut = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Results/Raster/dep.tif";
    const char *freqOut = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Results/Raster/freqwet.tif";
    const char *csv = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Results/csv/freqwet.csv";


    const char *binRas = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Results/Raster/bin.tif";
    const char *regRas = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Results/Raster/reg.tif";

//    const char *binRas = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/bin.tif";
//    const char *regRas = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/reg.tif";

    const char *damLayerName = "ESRIDams_BRAT_join";
    const char *bratLayerName = "lc_brat";
    //const char *damLayerName = "Dams_ESRI";



    return 0;
}


