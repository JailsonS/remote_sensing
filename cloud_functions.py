
import ee

# cloud params
CLD_PRB_THRESH = 50
CLOUD_FILTER = 60
NIR_DRK_THRESH = 0.15
CLD_PRJ_DIST = 1
BUFFER = 50
SR_BAND_SCALE = 1e4


def addCloudBands(img):
    
    cldPrb = ee.Image(img.get('s2cloudless')).select('probability')
    isCloud = cldPrb.gt(CLD_PRB_THRESH).rename('clouds')

    return img.addBands(ee.Image([cldPrb, isCloud]))



def addShadowBands(img):
  
    notWater = img.select('SCL').neq(6)

    darkPixels = img.select('B8').lt(NIR_DRK_THRESH * SR_BAND_SCALE).multiply(notWater).rename('dark_pixels')

    # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    shadowAzimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')))

    cldProj = img.select('clouds').directionalDistanceTransform(shadowAzimuth, CLD_PRJ_DIST * 10)\
        .reproject(**{'crs': img.select(0).projection(), 'scale': 100})\
        .select('distance')\
        .mask()\
        .rename('cloud_transform')

    shadows = cldProj.multiply(darkPixels).rename('shadows')

    return img.addBands(ee.Image([darkPixels, cldProj, shadows]))


def addCloudShadowMask(img):

    imgCloud = addCloudBands(img)
    imgCloudShadow = addShadowBands(imgCloud)

    isCldShdw = imgCloudShadow.select('clouds').add(imgCloudShadow.select('shadows')).gt(0)

    isCldShdw_ = isCldShdw.focalMin(2).focalMax(BUFFER*2/20)\
        .reproject(**{'crs': img.select([0]).projection(), 'scale': 20})\
        .rename('cloudmask')
        
    return imgCloudShadow.addBands(isCldShdw_)


