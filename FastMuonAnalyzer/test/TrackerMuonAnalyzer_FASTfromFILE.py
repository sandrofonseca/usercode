import FWCore.ParameterSet.Config as cms

process = cms.Process("test")

# The number of events to be processed.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) )

# Geometry:
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.DTGeometryBuilder.dtGeometry_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "MC_3XY_V21::All"
process.GlobalTag.globaltag = "START42_V9::All"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_1_2_FCH.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_13_1_IrU.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_103_1_NjV.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_73_1_iaA.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_49_1_bfO.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_3_1_RPM.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_82_1_s4Y.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_25_1_QBQ.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_24_1_REL.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_102_1_yvw.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_75_1_lsY.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_9_1_KYe.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_46_1_6pH.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_99_1_Yq7.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_81_1_ywj.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_76_1_Lqe.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_48_1_Eht.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_55_1_GFw.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_52_1_hU4.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_72_1_hfb.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_16_1_I6T.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_104_1_dnL.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_105_1_2w7.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_61_1_oXt.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_85_1_mo3.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_40_1_euA.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_94_1_SLm.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_56_1_gQd.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_5_1_ixi.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_39_1_KJT.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_2_1_HpP.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_4_1_aF6.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_45_1_fGH.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_106_1_9Br.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_26_1_269.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_108_1_jXC.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_44_1_YGk.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_38_1_1gp.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_93_1_SRu.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_50_1_b12.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_29_1_lu3.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_66_1_MKy.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_59_1_26J.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_107_1_VAO.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_86_1_Obe.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_79_1_lTM.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_31_1_Mbl.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_6_1_rAi.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_68_1_6WA.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_96_1_6Ah.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_65_1_sGC.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_67_1_X0F.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_36_1_DTt.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_95_1_suQ.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_80_1_57d.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_53_1_sGa.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_19_1_1rc.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_58_1_1L4.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_35_1_f4C.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_63_1_hCz.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_14_1_kDi.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_70_1_5IY.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_74_1_j4Z.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_18_1_Hne.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_42_1_HSX.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_43_1_rwV.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_10_1_pD9.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_37_1_WDG.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_91_1_7Py.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_110_1_lh8.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_109_1_nzv.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_89_1_YDe.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_101_1_Eaq.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_98_1_Aif.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_60_1_2jT.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_17_1_jzj.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_100_1_Lv4.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_87_1_nO3.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_78_1_5At.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_12_1_ehO.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_20_1_3t9.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_51_1_rC3.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_32_1_bDE.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_57_1_evY.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_22_1_pna.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_90_1_IYz.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_69_1_jxy.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_34_1_Ic1.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_15_1_GcK.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_30_1_MYI.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_92_1_ksW.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_28_1_KaP.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_47_1_XZV.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_8_1_jw1.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_62_1_G8C.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_23_1_dmN.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_54_1_AG8.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_77_1_Tpb.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_7_1_NOO.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_11_1_mqM.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_33_1_jeO.root',
## 'file:crab_FastSimStandard1000GeV_MuonBrem/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_97_1_nzt.root'


##FastMuon_Standard

## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_57_1_smF.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_67_1_Lnr.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_56_1_xwS.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_4_1_raU.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_105_1_HPv.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_76_1_zt5.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_86_1_RJO.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_70_1_jbE.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_43_1_E6z.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_93_1_g16.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_95_1_b9u.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_31_1_4rQ.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_10_2_QtX.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_92_1_liw.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_5_1_HXj.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_18_1_iR8.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_7_2_WVl.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_81_1_2vi.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_102_1_iPe.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_42_1_8DB.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_103_1_5LX.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_25_1_GGZ.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_88_1_V3h.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_68_1_8Zg.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_22_1_va9.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_50_2_5tu.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_107_1_tEI.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_109_1_sj0.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_23_1_zq6.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_34_1_GOa.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_106_1_HZj.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_85_2_Jis.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_24_1_vQK.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_15_1_DWh.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_38_1_h7k.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_17_1_Ixq.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_53_1_EpY.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_104_1_EiA.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_9_1_RK5.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_55_1_J4t.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_28_1_Zzg.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_80_2_avW.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_73_1_xWJ.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_87_1_5AA.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_36_1_jFP.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_74_1_iZE.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_47_1_Bj7.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_90_1_Vpe.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_45_2_QQA.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_65_1_76a.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_1_1_ZrL.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_21_1_Scv.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_37_4_1uH.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_12_1_Rat.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_98_1_xLX.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_8_3_AHy.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_27_2_Mme.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_91_1_Bqy.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_48_1_ZkB.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_64_2_rP4.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_79_2_Loh.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_66_1_IlA.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_77_2_Jne.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_16_1_rfN.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_49_1_pTO.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_89_1_WfK.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_30_1_P8V.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_72_1_f7a.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_71_1_h2x.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_32_2_TYT.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_82_1_Ao3.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_29_2_Zq7.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_52_2_P2X.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_54_2_SA7.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_3_3_RfU.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_58_2_oTE.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_6_1_3hJ.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_61_1_DM8.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_78_1_Wpk.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_35_4_0af.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_84_1_PcV.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_2_1_P9H.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_11_1_djC.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_13_1_9GD.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_69_1_Jiz.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_39_1_OoN.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_60_1_UJ7.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_46_1_SFg.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_101_2_pIy.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_75_2_DXt.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_59_2_kLL.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_63_1_xfd.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_26_2_OVF.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_51_1_KTi.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_110_1_FQM.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_19_2_wWd.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_20_1_okX.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_94_2_dkg.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_41_2_NR4.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_44_2_uiY.root',
## 'file:crab_FastSimStandard1000GeV_standard/res/SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_83_1_e8T.root'


###Full Sim



## Single mu, pt = 1000 GeV
#
#     'rfio:/castor/cern.ch/user/a/aperrott/FastSim/CMSSW_3_5_8_SingleMuPt1000_cfi_py_GEN_FASTSIM_VALIDATION.root'
'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_default_perrotta_test.root'
#'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_default_TrueMuonBrem.root'
#'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_MuonBremTrue.root'
#'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION.root'    
# 'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case1.root'
#  'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case2.root'   
##  'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case3.root'
##  'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case4.root'
##  'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case5.root'
##  'file:SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case6.root'

    )
)

process.TrackerMuonAnalyzer = cms.EDAnalyzer("TrackerMuonAnalyzer",
     rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_default_perrotta.root'),
 #   rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_default_TrueMuonBrem.root'),
   # rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE.root'),
##     rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_Case1.root'),
##     rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_Case2.root'),
##      rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_Case3.root'),
##      rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_Case4.root'),
##      rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_Case5.root'),
##      rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_Case6.root'),
                                             
   # rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_MuonBrem.root'),
   
    FastSim = cms.untracked.bool(True),
    IgnoreMissingCollections = cms.untracked.bool(True),
    RootTree = cms.untracked.bool(False),
    Debug = cms.untracked.bool(False),
    useGen = cms.untracked.bool(False),
    useReco = cms.untracked.bool(False),
    useEff = cms.untracked.bool(False),
    test_Track= cms.untracked.bool(True),

    labelSimVertex = cms.InputTag("famosSimHits"),
    labelSimTrack = cms.InputTag("famosSimHits","MuonSimTracks"),
    dtSimHits = cms.InputTag("MuonSimHits","MuonDTHits"),
    cscSimHits = cms.InputTag("MuonSimHits","MuonCSCHits"),
    rpcSimHits = cms.InputTag("MuonSimHits","MuonRPCHits"),

  ##   #Tracker Labels
##     labelSimHitsTOB = cms.InputTag("famosSimHits","TrackerHits"),
##     labelSimHitsTIB = cms.InputTag("famosSimHits","TrackerHits"),   
                                             

    # VERIFICA SEMPRE TUTTE QUESTE LABELS!
    labelRecoTrackFast = cms.InputTag("generalTracks"),
    labelL1MuonFast = cms.InputTag("l1extraParticles"),
#    labelL1MuonFast = cms.InputTag("MISSING"),
    labelL2MuonFast = cms.InputTag("hltL2Muons"),
    labelL3MuonFast = cms.InputTag("hltL3Muons"),
#    labelL3MuonFast = cms.InputTag("MISSING"),
    labelGlobalMuonFast = cms.InputTag("muons"),
    labelGlobalMuonTrackFast = cms.InputTag("globalMuons"),
    labelStandAloneMuonTrackFast = cms.InputTag("standAloneMuons"),
#    labelStandAloneMuonTrackFast = cms.InputTag("standAloneMuons","UpdatedAtVtx"),

)

process.evtContent = cms.EDAnalyzer("EventContentAnalyzer")

process.p1 = cms.Path(process.TrackerMuonAnalyzer)

## To write out events (not need: FastSimulation _is_ fast!)
#process.o1 = cms.OutputModule(
#    "PoolOutputModule",
#    fileName = cms.untracked.string("NewMuonAnalyzerOutputFile.root"),
##    outputCommands = cms.untracked.vstring("keep *",
##                                           "drop *_mix_*_*",
##                                          "drop *_*_*_RECO",
##                                           "drop *_*_*_HLT")
#    )
#process.outpath = cms.EndPath(process.o1)

# Keep the logging output to a nice level #
# process.Timing =  cms.Service("Timing")
process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.destinations = cms.untracked.vstring("pyDetailedInfo.txt")

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )
