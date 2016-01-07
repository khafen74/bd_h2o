#include "dampolygons.h"

DamPolygons::DamPolygons(DamPoints pts, double maxDamHt)
{
    setMaxDamHeight(maxDamHt);
    init(pts);
}

void DamPolygons::init(DamPoints pondPts)
{
    loadDriver();
    setOutDirPath(pondPts.getOutDirPath());
    setDemPath(pondPts.getDemPath());
    setMaxDistance(150.0);
    m_layerName = "PondPolygons";

    OGRDataSource *pOutDs;
    OGRLayer *pPoly_lyr, *pPts_lyr;

    pOutDs = m_pDriverShp->CreateDataSource(m_outDir, NULL);
    pPts_lyr = pOutDs->GetLayerByName(pondPts.getLayerName());
    pPoly_lyr = pOutDs->CreateLayer(m_layerName, pPts_lyr->GetSpatialRef(), wkbPolygon, NULL);

    createPondPolygons(pPts_lyr, pPoly_lyr);

    OGRDataSource::DestroyDataSource(pOutDs);
}

void DamPolygons::loadDriver()
{
    OGRRegisterAll();
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    m_pDriverShp = registrar->GetDriverByName("ESRI Shapefile");
}

void DamPolygons::createFields(OGRLayer *pLayer)
{
    OGRFieldDefn field("brat_ID", OFTInteger);
    pLayer->CreateField(&field);
    field.SetName("pond_ID");
    field.SetType(OFTInteger);
    pLayer->CreateField(&field);
    field.SetName("WSE_max");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
}

void DamPolygons::createPondPolygons(OGRLayer *pPts, OGRLayer *pPolys)
{
    double xCoords[5], yCoords[5];
    double azimuthCurrent, azimuthStart, distance, slope;
    int nDams = pPts->GetFeatureCount();

    for (int i=0; i<nDams; i++)
    {
        OGRFeature *pDamFeature;
        OGRPolygon polygon;
        OGRLinearRing ring;
        OGRPoint point, *damPoint;
        OGRFeature *pPolyFeature = OGRFeature::CreateFeature(pPolys->GetLayerDefn());

        pDamFeature = pPts->GetFeature(i);
        azimuthStart = pDamFeature->GetFieldAsDouble("az_us");
        slope = pDamFeature->GetFieldAsDouble("slope");
        distance = getUpstreamDistance(m_maxDamHt/slope);
        damPoint = (OGRPoint*) pDamFeature->GetGeometryRef();

        for (int j=0; j<5; j++)
        {
            azimuthCurrent = Geometry::addDegrees(azimuthStart, ANGLE_OFFSET[j]);
            Geometry::calcCoords(damPoint->getX(), damPoint->getY(), azimuthCurrent, distance, xCoords[j], yCoords[j]);

            point.setX(xCoords[j]);
            point.setY(yCoords[j]);

            ring.addPoint(&point);
        }

        ring.addPoint(xCoords[0], yCoords[0]);
        polygon.addRing(&ring);
        pPolyFeature->SetGeometry(&polygon);
        pPolys->CreateFeature(pPolyFeature);
        OGRFeature::DestroyFeature(pPolyFeature);
        OGRFeature::DestroyFeature(pDamFeature);
    }
}

double DamPolygons::getUpstreamDistance(double dist)
{
    if (dist > m_maxDist)
    {
        return m_maxDist;
    }
    else
    {
        return dist;
    }
}

void DamPolygons::setDemPath(const char *demPath)
{
    m_demPath = demPath;
}

void DamPolygons::setFieldValues(OGRFeature *pFeat, int bratID, int pondID, double maxWSE)
{
    pFeat->SetField("brat_ID", bratID);
    pFeat->SetField("pond_ID", pondID);
    pFeat->SetField("WSE_max", maxWSE);
}

void DamPolygons::setMaxDamHeight(double maxDamHt)
{
    m_maxDamHt = maxDamHt;
}

void DamPolygons::setMaxDistance(double distance)
{
    m_maxDist = distance;
}

void DamPolygons::setOutDirPath(const char *outDir)
{
    m_outDir = outDir;
}

