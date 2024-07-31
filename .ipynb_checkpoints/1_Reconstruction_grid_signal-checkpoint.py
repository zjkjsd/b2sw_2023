# +
import basf2 as b2
import modularAnalysis as ma
from variables import variables as vm
import variables.collections as vc
import variables.utils as vu
import vertex as vx
import stdPi0s


b2.B2INFO(f"Prepending analysis GT: {ma.getAnalysisGlobaltag()}")
b2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())
#b2.conditions.prepend_globaltag('chargedpidmva_rel6_v5')

# Define the path
main_path = b2.Path()

input_file = ['']
output_file = 'grid.root'
# import sys
# input_file = sys.argv[1]
# output_file = sys.argv[2]

ma.inputMdstList(environmentType='default', filelist=input_file, path=main_path)

########## Tracks and Clusters ###########

# Hadron tracks: https://confluence.desy.de/display/BI/Hadron+ID+Performance
goodTrack = '[abs(dz)<4] and [dr<2] and thetaInCDCAcceptance and [nCDCHits>20]'

# Neutral performance: https://confluence.desy.de/display/BI/Neutrals+Performance
vm.addAlias('goodphoton','passesCut([clusterE>0.05] and [abs(clusterTiming)<formula(2*clusterErrorTiming)] and \
[abs(clusterTiming)<200] and inECLAcceptance and [beamBackgroundSuppression>0.5] and [fakePhotonSuppression>0.8])')


########## Final State Particles ###########
ma.fillParticleList('K+:myk', cut=goodTrack + ' and [kaonID_noSVD>0.9]', path=main_path)

stdPi0s.stdPi0s(listtype='eff60_May2020',beamBackgroundMVAWeight="MC15ri",
                fakePhotonMVAWeight="MC15ri",path=main_path)
ma.cutAndCopyLists(outputListName='pi0:mypi0', inputListNames='pi0:eff60_May2020',
                   cut='[0.115<M<0.15] and [daughter(0,goodphoton)] and [daughter(1,goodphoton)]', path=main_path)

vm.addAlias('daughterAngle','daughterAngle(0, 1)')
vm.addAlias('daughterAngle_CMS','useCMSFrame(daughterAngle(0, 1))')

daughterDiff_vars = [f'daughterAngle{frame}' for frame in ['','_CMS']]
for var in ['phi','theta','p']:
    vm.addAlias(f'daughterDiff_{var}',f'daughterDiffOf(0,1,{var})')
    vm.addAlias(f'daughterDiff_CMS_{var}',f'useCMSFrame(daughterDiffOf(0,1,{var}))')
    #vm.addAlias(f'daughterDiff_rest_{var}',f'useRestFrame(daughterDiffOf(0, 1, {var}))')
    daughterDiff_vars += [f'daughterDiff{frame}_{var}' for frame in ['','_CMS']]

# for grid jobs
# ma.applyCuts("K+:myk", "[binaryPID_noSVD(321, 211)>0.9] and [pionIDNN<0.1]", path=main_path) # take off for systematics calculation
# ma.applyCuts("pi0:mypi0", '[daughterAngle<1] and [abs(daughterDiff_phi)<1] and [abs(daughterDiff_theta)<0.5]', path=main_path)

########## Event Kinematics/Shape ##########
ma.buildEventKinematics(fillWithMostLikely=True,path=main_path)
ma.buildEventShape(cleoCones=False, collisionAxis=False, foxWolfram=True, 
                   harmonicMoments=False, jets=False, sphericity=False, path=main_path)
ESVariables = ['foxWolframR1', 'foxWolframR4', 'thrustAxisCosTheta']


########### Reconstruct B ##########
ma.reconstructDecay('B+:Kpi0 -> K+:myk pi0:mypi0', cut='[Mbc>5.26] and [abs(deltaE)<0.5]', path=main_path)
# required for the tagV
vx.treeFit('B+:Kpi0', conf_level=0.00, updateAllDaughters=False, massConstraint=[], ipConstraint=True, path=main_path)

DCSV_momentum = ['p', 'pErr', 'phi', 'phiErr', 'theta', 'thetaErr']
DCSV_cluster = ['clusterNHits', 'clusterTiming', 'clusterE9E21', 'clusterE', 'clusterReg']
DCSV_track = ['nCDCHits', 'pValue', 'charge','electronID_noSVD', 'kaonID_noSVD',
              'muonID_noSVD', 'pionID_noSVD', 'protonID_noSVD', 'deuteronID_noSVD']

# for grid jobs
# ma.applyCuts("B+:Kpi0", "[Mbc>5.25] and [daughterAngle_CMS>2.96] and [abs(daughterDiff_CMS_phi)>2.9] and \
# [abs(daughterDiff_CMS_p)<0.5]", path=main_path)

########## Best Candidate Selection ###########
ma.rankByHighest("B+:Kpi0", variable="Mbc", numBest=1, path=main_path)

########### MC Truth Matching ##########
ma.matchMCTruth('B+:Kpi0', path=main_path)
ma.applyCuts("B+:Kpi0", "isSignal==1", path=main_path)

########### Build the ROE ##########
ma.fillParticleList('gamma:all', '', path=main_path)
ma.getBeamBackgroundProbability('gamma:all','MC15ri', path=main_path)
ma.getFakePhotonProbability('gamma:all','MC15ri', path=main_path)
ma.buildRestOfEvent('B+:Kpi0', fillWithMostLikely=True, path=main_path)

good_track = '[abs(dz)<20] and [dr<10] and thetaInCDCAcceptance and [nCDCHits>0] and [pt>0.075]' 
good_gamma = f'[clusterE>0.05] and [clusterNHits>1.5] and [abs(clusterTiming)<formula(2*clusterErrorTiming)] and \
[abs(clusterTiming)<200] and [beamBackgroundSuppression>0.2] and [fakePhotonSuppression>0.2]'
roe_mask = ('my_mask', good_track, good_gamma)
ma.appendROEMasks('B+:Kpi0', [roe_mask], path=main_path)


# creates V0 particle lists and uses V0 candidates to update/optimize the Rest Of Event
ma.updateROEUsingV0Lists('B+:Kpi0', mask_names='my_mask', default_cleanup=True, selection_cuts=None,
                         apply_mass_fit=True, fitter='treefit', path=main_path)

# save roe tracks and clusters to the ntuple
roe_path = b2.Path()
deadEndPath = b2.Path()
ma.signalSideParticleFilter('B+:Kpi0', '', roe_path, deadEndPath)

roe_cut = '[isInRestOfEvent==1] and passesROEMask(my_mask)'
roe_tracks = ('pi+:roe',roe_cut)
roe_gammas = ('gamma:roe', roe_cut)
ma.fillParticleLists([roe_tracks, roe_gammas], path=roe_path)
ma.variablesToNtuple("pi+:roe", DCSV_momentum + DCSV_track + ['mcPDG','isSignal'], 
                     filename=output_file, treename='roe_tracks', path=roe_path)
ma.variablesToNtuple("gamma:roe", DCSV_momentum + DCSV_cluster + ['mcPDG','isSignal'], 
                     filename=output_file, treename='roe_clusters', path=roe_path)
main_path.for_each('RestOfEvent', 'RestOfEvents', roe_path)

vm.addAlias('roeP','roeP(my_mask)')
vm.addAlias('roePTheta','roePTheta(my_mask)')
vm.addAlias('roeCharge','roeCharge(my_mask)')
DCSV_roe = ['roeP', 'roePTheta', 'roeCharge']

# fit the tag-side B vertex
vx.TagV('B+:Kpi0',confidenceLevel=0.0,trackFindingType='standard_PXD',MCassociation='breco',constraintType='tube', 
        reqPXDHits=0, maskName='my_mask', fitAlgorithm='KFit', kFitReqReducedChi2=5.0, path=main_path)
vm.addAlias('TagVReChi2','formula(TagVChi2/TagVNDF)')
vm.addAlias('TagVReChi2IP','formula(TagVChi2IP/TagVNDF)')

TVVariables = ['DeltaZ',    'DeltaZErr',    'TagVReChi2',    'TagVReChi2IP',
               'TagVx',     'TagVxErr',     'TagVy',         'TagVyErr',
               'TagVz',     'TagVzErr',     'TagVNTracks']


########## Continuum Suppression #########
ma.buildContinuumSuppression(list_name="B+:Kpi0", roe_mask="my_mask", path=main_path)

vm.addAlias('KSFWV_et','KSFWVariables(et)') #correlates with p_D + p_l
vm.addAlias('KSFWV_mm2','KSFWVariables(mm2)') #correlates with mm2
vm.addAlias('KSFWV_hso00','KSFWVariables(hso00)')
vm.addAlias('KSFWV_hso01','KSFWVariables(hso01)')
vm.addAlias('KSFWV_hso02','KSFWVariables(hso02)')
vm.addAlias('KSFWV_hso03','KSFWVariables(hso03)')
vm.addAlias('KSFWV_hso04','KSFWVariables(hso04)')
vm.addAlias('KSFWV_hso10','KSFWVariables(hso10)')
vm.addAlias('KSFWV_hso12','KSFWVariables(hso12)')
vm.addAlias('KSFWV_hso14','KSFWVariables(hso14)')
vm.addAlias('KSFWV_hso20','KSFWVariables(hso20)') #correlates with mm2
vm.addAlias('KSFWV_hso22','KSFWVariables(hso22)')
vm.addAlias('KSFWV_hso24','KSFWVariables(hso24)')
vm.addAlias('KSFWV_hoo0','KSFWVariables(hoo0)')
vm.addAlias('KSFWV_hoo1','KSFWVariables(hoo1)')
vm.addAlias('KSFWV_hoo2','KSFWVariables(hoo2)')
vm.addAlias('KSFWV_hoo3','KSFWVariables(hoo3)')
vm.addAlias('KSFWV_hoo4','KSFWVariables(hoo4)')
vm.addAlias('CC1','CleoConeCS(1)')
vm.addAlias('CC2','CleoConeCS(2)')
vm.addAlias('CC3','CleoConeCS(3)')
vm.addAlias('CC4','CleoConeCS(4)')
vm.addAlias('CC5','CleoConeCS(5)')
vm.addAlias('CC6','CleoConeCS(6)')
vm.addAlias('CC7','CleoConeCS(7)')
vm.addAlias('CC8','CleoConeCS(8)')
vm.addAlias('CC9','CleoConeCS(9)')

CSVariables = ['isContinuumEvent',    "R2",            "thrustBm",
               "thrustOm",            "cosTBTO",       "cosTBz",
               'KSFWV_et',            'KSFWV_mm2',     'KSFWV_hso00',
               'KSFWV_hso01',         'KSFWV_hso02',   'KSFWV_hso03',
               'KSFWV_hso04',         'KSFWV_hso10',   'KSFWV_hso12',
               'KSFWV_hso14',         'KSFWV_hso20',   'KSFWV_hso22',
               'KSFWV_hso24',         'KSFWV_hoo0',    'KSFWV_hoo1',
               'KSFWV_hoo2',          'KSFWV_hoo3',    'KSFWV_hoo4',
               "CC1",                 "CC2",           "CC3",
               "CC4",                 "CC5",           "CC6",
               "CC7",                 "CC8",           "CC9", ]


########## Save Ntuples ##########
# ma.printMCParticles(onlyPrimaries=False, maxLevel=-1, path=main_path, suppressPrint=True,
#                     showProperties=False,showMomenta=False,showVertices=False,showStatus=False)

cms_kinematics = vu.create_aliases(vc.kinematics, "useCMSFrame({variable})", "CMS")
cms_mc_kinematics = vu.create_aliases(vc.mc_kinematics, "useCMSFrame({variable})", "CMS")
vm.addAlias('CMS_cosTheta', 'useCMSFrame(cosTheta)')

b_vars = vu.create_aliases_for_selected(
    list_of_variables= vc.deltae_mbc + CSVariables + ESVariables + TVVariables + daughterDiff_vars
    + DCSV_momentum + DCSV_roe + ['mcErrors','mcPDG','isSignal','CMS_cosTheta'],
    decay_string='^B+:Kpi0 -> K+:myk pi0:mypi0',
    prefix=['B'])

k_vars = vu.create_aliases_for_selected(
    list_of_variables= DCSV_momentum + DCSV_track + ['genMotherPDG','mcPDG','isSignal'],
    decay_string='B+:Kpi0 -> ^K+:myk pi0:mypi0',
    prefix=['K'])

pi0_vars = vu.create_aliases_for_selected(
    list_of_variables= cms_kinematics + ['genMotherPDG','mcErrors', 'mcPDG','isSignal','M'],
    decay_string='B+:Kpi0 -> K+:myk ^pi0:mypi0',
    prefix=['pi'])

gamma_vars = vu.create_aliases_for_selected(
    list_of_variables= DCSV_momentum + DCSV_cluster
    + ['genMotherPDG','mcPDG','isSignal','beamBackgroundSuppression', 'fakePhotonSuppression'],
    decay_string='B+:Kpi0 -> K+:myk [pi0:mypi0 -> ^gamma:eff60_May2020 ^gamma:eff60_May2020]',
    prefix=['g1','g2'])

candidate_vars = ['Ecms'] + b_vars + k_vars + pi0_vars + gamma_vars

ma.variablesToNtuple('B+:Kpi0', candidate_vars, basketsize=20000000,
                     filename=output_file, treename='Bsig', path=main_path)

b2.process(path=main_path)

