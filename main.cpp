#include <QCoreApplication>
#include <QtCore>
#include <QDebug>
#include "storagemodel.h"

int test();
int run();
int validate();
int pataha();
int testStats();

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QDateTime startTime = QDateTime::currentDateTime();

    //test();
    validate();
    //testStats();

    QDateTime endTime = QDateTime::currentDateTime();

    qDebug()<<"done"<<startTime.secsTo(endTime)<<startTime.secsTo(endTime)/60.0;

    return a.exec();
}

int test()
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
    const char *demInLogan = "C:/Users/khafe/Desktop/Classes/CEE_6400_Hydrology/FinalProject/Raster/lc_dem_10m.tif";
    const char *exDams = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/ValidationLayers/Dams_BRAT_join5_UTM12N_pondArea.shp";
    const char *csvOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut/comparison10m_dist.csv";
    const char *heightOut10m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut/hand10m.tif";
    const char *test10m = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut/test.tif";

//    StorageModel model(shpIn, shpOut, demIn10m_clip, fdir10m_clip, fac10m_clip, 0.5);
//    model.run();
//    model.runFromPoints(exDams, csvOut);

    const char *handfil = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/WSE_start.tif";
    const char *watsurf = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/dem_dep.tif";
    const char *handfdir = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/HAND_fdir.tif";
    const char *depmid = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/depMid.tif";
    const char *pondmid = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/depMidPond.tif";
    const char *handin = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/HAND_inMid.tif";
    const char *handnew = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/HAND_outMid.tif";
    const char *handchange = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/HAND_change.tif";
    const char *pondid = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/HAND_pondId.tif";
    const char *gwextend = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/gw_extend.tif";
    const char *gwdepth = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/20160606_10m_GW/gw_depth.tif";

    Raster_BeaverPond rasterBP;
    Raster raster;
    //raster.setNoData(fil1m, -9999, 0, 5000);
    qDebug()<<"creating hand in";
    rasterBP.createHANDInput(pondmid, fac10m_clip, handin);
    qDebug()<<"adding pond depth to dem";
    rasterBP.add(fil10m, depmid, watsurf);
    qDebug()<<"height above network";
    raster.heightAboveNetwork(fil10m, handfdir, fac10m_clip, handfil);
    rasterBP.heightAboveNetwork(watsurf, handfdir, handin, handnew, pondid);
    qDebug()<<"subtract";
    rasterBP.subtractHAND(handfil, handnew, handchange);
    qDebug()<<"done";

    //testing gw extension via D8 flow direction
    //rasterBP.groundwaterDepth(handfil, handnew, gwdepth);
    rasterBP.flowDownstream(handchange, fdir10m, fac10m, demIn10m, pondmid, gwdepth, gwextend);


    //Raster raster;
    //raster.extractByMask_CellCenters(demIn1m, demIn1m_clip, shpInDir, vbtm200buf);
    //raster.extractByMask_CellCenters(fdir1m, fdir1m_clip, shpInDir, vbtm);
    //raster.extractByMask_CellCenters(fac1m, fac1m_clip, shpInDir, vbtm);
    //raster.setNoData(fac1m_clip, -9999, 0, 2);
    //raster.heightAboveNetwork(fil10m, fdir10m, fac10m, heightOut10m);
    //raster.subtract(demIn10m, heightOut10m, test10m);

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

    StorageModel model(shpIn, shpOut, demIn10m, fdir10m, fac10m, 0.5);
    model.run();

    //StorageModel model(shpIn, shpOut, demIn10m_clip, fdir10m_clip, fac10m_clip, 0.5);
    //model.run();
    //model.runFromPoints(exDams, csvOut);



    return 0;
}

int validate()
{
    //Logan HUC 8 inputs
    const char *bratPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8/01_shpIn/brat.shp";
    const char *demPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8/02_rasIn/fil10m2.tif";
    const char *fdirPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8/02_rasIn/fdir10m2.tif";
    const char *facPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8/02_rasIn/fac10m_1700.tif";
    const char *outDir = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8/03_out";
    const char *damsPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8/01_shpIn/Dams_LoganHUC8_JoinBRAT.shp";
    const char *csvPath = "E:/etal/Projects/NonLoc/BeaverModeling/06_ValidationSurfaceStorage/Logan_HUC8/03_out/comparison.csv";

    StorageModel model(bratPath, outDir, demPath, fdirPath, facPath, 1.0);
    model.runFromPoints(damsPath, csvPath);

    return 0;
}

int pataha()
{
    const char *shpIn = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelIn/Pataha_BRAT.shp";
    const char *shpOut = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelOut";
    const char *demIn = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelIn/pataha_10m.tif";
    const char *fdir = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelIn/fdir_10m.tif";
    const char *fac = "E:/etal/Projects/NonLoc/BeaverModeling/03_Results/Pataha/Inputs/ModelIn/fac_10m.tif";

    StorageModel model(shpIn, shpOut, demIn, fdir, fac, 0.5);
    model.run();

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
