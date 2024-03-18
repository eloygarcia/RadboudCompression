import sys
sys.path.append('/aux/')

global elasticParameters
global NHParameters
global TIParameters
global anisotropicDirection
global anisotropicRatio

anisotropicRatio = 0.5
def bulkModulus( Eyoung, nuPoisson):
    return (Eyoung/(3*(1-2*nuPoisson)))

def shearModulus( Eyoung, nuPoisson):
    return (Eyoung/(2*(1+nuPoisson)))

def firstLameParameter( Eyoung, nuPoisson):
    return ((Eyoung*nuPoisson)/((1+nuPoisson)*(1-2*nuPoisson)))

def secondLameParameter( Eyoung, nuPoisson):
    return (Eyoung/(2*(1+nuPoisson)))

def anisotropyTerm( Eyoung, anisotropyRatio):
    return (Eyoung*anisotropyRatio)

anisotropicDirection = [ 1, 0, 0]
# material_parameters.push_back(1491.64);
# material_parameters.push_back(147672.24); // Fatty
# material_parameters.push_back(1141760);
# material_parameters.push_back(1);
# material_parameters.push_back(0);
# material_parameters.push_back(0);

###
#elasticParameters = [[4.46*1000, 0.499], # adipose
#                      [15.1*1000, 0.499], #glandular
#                      [60.0*1000, 0.499]] # skin

elasticParameters = [[4.46, 0.499], # adipose
                      [15.1, 0.499], #glandular
                      [20.0, 0.499]] # skin
def computeNHparameters( elasticParameters: list):
    NHParameters = []
    for i in range(len(elasticParameters)):
        NHParameters.append([shearModulus(elasticParameters[i][0],elasticParameters[i][1]),
                             bulkModulus(elasticParameters[i][0],elasticParameters[i][1])])
    return NHParameters

NHParameters = computeNHparameters(elasticParameters)

def computeTIparameters(elasticParameters: list, anisotropicRatio: float, anisotropicDirection: list):
    TIParameters = []
    for i in range(len(elasticParameters)):
        TIParameters.append([ shearModulus(elasticParameters[i][0],elasticParameters[i][1]),
                              bulkModulus(elasticParameters[i][0],elasticParameters[i][1]),
                              anisotropyTerm(elasticParameters[i][0], anisotropicRatio)] + anisotropicDirection)
    return TIParameters

TIParameters = computeTIparameters(elasticParameters, anisotropicRatio, anisotropicDirection)