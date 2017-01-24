#include <QCoreApplication>
#include <QtCore>
#include <QDebug>
#include "storagemodel.h"

int test();
int run();
int validate();
int logan();
int pataha();
int testStats();
int runXYZ();
int soil();

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QDateTime startTime = QDateTime::currentDateTime();

    //test();
    validate();
    //logan();
    //testStats();
    //runXYZ();
    //soil();

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

//    Raster_BeaverPond rasterBP;
//    Raster raster;
//    //raster.setNoData(fil1m, -9999, 0, 5000);
//    qDebug()<<"creating hand in";
//    rasterBP.createHANDInput(pondmid, fac10m_clip, handin);
//    qDebug()<<"adding pond depth to dem";
//    rasterBP.add(fil10m, depmid, watsurf);
//    qDebug()<<"height above network";
//    raster.heightAboveNetwork(fil10m, handfdir, fac10m_clip, handfil);
//    rasterBP.heightAboveNetwork(watsurf, handfdir, handin, handnew, pondid);
//    qDebug()<<"subtract";
//    rasterBP.subtractHAND(handfil, handnew, handchange);
//    qDebug()<<"done";

//    //testing gw extension via D8 flow direction
//    //rasterBP.groundwaterDepth(handfil, handnew, gwdepth);
//    rasterBP.flowDownstream(handchange, fdir10m, fac10m, demIn10m, pondmid, gwdepth, gwextend);


    //Raster raster;
    //raster.extractByMask_CellCenters(demIn1m, demIn1m_clip, shpInDir, vbtm200buf);
    //raster.extractByMask_CellCenters(fdir1m, fdir1m_clip, shpInDir, vbtm);
    //raster.extractByMask_CellCenters(fac1m, fac1m_clip, shpInDir, vbtm);
    //raster.setNoData(fac1m_clip, -9999, 0, 2);
    //raster.heightAboveNetwork(fil10m, fdir10m, fac10m, heightOut10m);
    //raster.subtract(demIn10m, heightOut10m, test10m);

    const char *soil = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/Soil/UT/gssurgo_ut10m_utm12.tif";
    const char *points = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/Soil/UT/grid_points/ut_points.txt";
    Raster raster;
    raster.toXYZ(soil, points);



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
    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/01_shpIn/BRAT_TempleFk_WS.shp";
    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/02_rasIn/fil10m_vb.tif";
    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/02_rasIn/fdir10m_vb.tif";
    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/02_rasIn/fac_2500_vb.tif";
    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/03_out";
    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/01_shpIn/DamArea_BRAT_joined.shp";
    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m/03_out/comparison.csv";

//    //Temple Fork 10m inputs volume comparison
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/01_shpIn/BRAT_TempleFk_WS.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_10m_vol/02_rasIn/fdir10m_vb.tif";
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
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m_vol/02_rasIn/fil1m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/TempleFork_1m_vol/02_rasIn/fdir1m_vb.tif";
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
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/02_rasIn/fdir10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/02_rasIn/fac_2500_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/01_shpIn/DamArea_BRAT_joined.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTest/03_out/comparison.csv";

    //Recursive with Temple Fork 1m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/01_shpIn/BRAT_TempleFk_WS.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/02_rasIn/fdir1m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/02_rasIn/fac_200000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/01_shpIn/DamArea_BRAT_joined.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestTF1m/03_out/comparison.csv";

    //Recursive with Ogden NF 1m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/02_rasIn/fdir1m.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/02_rasIn/fac1m_500000.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/01_shpIn/dams_brat_join.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO1m/03_out/comparison.csv";

    //Recursive with Ogden NF 10m inputs
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/02_rasIn/dem_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/z_RecursiveTestNFO10m/02_rasIn/fdir10m.tif";
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
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/02_rasIn/fil_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/02_rasIn/fdir_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/02_rasIn/fac_50000_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/01_shpIn/points_snap.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_10m_vol/03_out/comparison.csv";

    //Bridge Creek 1m inputs volume comparison
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/01_shpIn/BRAT_bcHUC10.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/02_rasIn/fil_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/02_rasIn/fdir_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/02_rasIn/fac_1mil_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/01_shpIn/points_snap.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/BridgeCreek_1m_vol/03_out/comparison.csv";

    //Curtis Creek 10m inputs (whole Logan HUC 8) volume comparison
//    const char *bratPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/01_shpIn/brat.shp";
//    const char *demPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/02_rasIn/fil10m_vb.tif";
//    const char *fdirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/02_rasIn/fdir10m_vb.tif";
//    const char *facPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/02_rasIn/fac10m_vb.tif";
//    const char *outDir = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/03_out";
//    const char *damsPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/01_shpIn/points_snap.shp";
//    const char *csvPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/CurtisCreek_10m_vol/03_out/comparison.csv";

//    Raster raster;
//    raster.setNoData(demPath, -9999.0, 100, 4000);
//    raster.setNoData(fdirPath, 0, 1, 200);
//    raster.setNoData(facPath, -1, 1, 2);

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
    StorageModel model(bratPath, outDir, demPath, fdirPath, facPath, 0.5, 3);

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
    //model.runFromPoints(damsPath, csvPath, 1);
    //model.runFromPointsWithHeights(damsPath,csvPath, 3);
    model.run(2);

    return 0;
}

int logan()
{
    QString dirPath = "F:/01_etal/Projects/Modeling/BeaverWaterStorage/wrk_Data/ValidationRuns/Logan_HUC8/05_HUC12";
    QDir dir(dirPath);
    QFileInfoList fiList = dir.entryInfoList();
    qDebug()<<"starting loop";

    for (int i=2; i<fiList.length(); i++)
    {
        qDebug()<<fiList[i].baseName();
        std::string basePath = dirPath.toStdString() + "/" + fiList[i].baseName().toStdString();
        std::string bratPath = basePath + "/01_shpIn/brat.shp";
        std::string outDir = basePath + "/03_out";
        std::string demPath = basePath + "/02_rasIn/fil.tif";
        std::string facPath = basePath + "/02_rasIn/fac.tif";
        std::string fdirPath = basePath + "/02_rasIn/fdir.tif";
        //qDebug()<<"names set"<<bratPath.toStdString().c_str()<<outDir.toStdString().c_str()<<demPath.toStdString().c_str()<<facPath.toStdString().c_str()<<fdirPath.toStdString().c_str();

        //qDebug()<<"running"<<i+1<<"of"<<fiList.length()<<basePath;
        StorageModel model(bratPath.c_str(), outDir.c_str(), demPath.c_str(), fdirPath.c_str(), facPath.c_str(), 0.82, 2);
        model.run(1);
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
