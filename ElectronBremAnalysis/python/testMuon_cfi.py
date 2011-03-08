import FWCore.ParameterSet.Config as cms

XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/materials.xml', 
        'SimG4CMS/MuonBremAnalysis/data/testMuon.xml'),
    rootNodeName = cms.string('testMuon:TestMuon')
)


