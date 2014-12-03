#! /usr/bin/env python
import os
import re
import sys
import subprocess

from optparse import OptionParser

parser = OptionParser()

parser.add_option('--sqlite_file', metavar='F', type='string', action='store',
                  dest='sqlite_file',
                  help='Input sqlite file')

parser.add_option('--version', metavar='F', type='string', action='store',
                  dest='version',
                  help='Version')

parser.add_option('--prep', metavar='F', action='store_true',
                  dest='prep',
                  help='Upload to prep area')

(options, args) = parser.parse_args()

#******************   template file  **********************************
if options.prep :
	templateFile = open('templateForDropbox_PREP.txt', 'r')
else :
	templateFile = open('templateForDropbox_PRODUCTION.txt', 'r')
fileContents = templateFile.read(-1)
print '--------------- TEMPLATE :  -----------------'
print fileContents
p1 = re.compile(r'TAGNAME')
p2 = re.compile(r'PRODNAME')

#******************   definitions  **********************************

##version = '2014Apr29'
version = options.version

payloads = [
#    "RecoTauTag_againstMuonMVAv1",
#     "RecoTauTag_againstMuonMVAv1_WPeff98_0",
#     "RecoTauTag_againstMuonMVAv1_WPeff99_0",
#     "RecoTauTag_againstMuonMVAv1_WPeff99_5",
#     "RecoTauTag_againstMuonMVAv1_mvaOutput_normalization",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff99",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff79",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff85",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff91",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff96",
    # "RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff99",
    "RecoTauTag_tauIdMVAnewDMwLTv1",
    "RecoTauTag_tauIdMVAnewDMwLTv1_WPEff40",
    "RecoTauTag_tauIdMVAnewDMwLTv1_WPEff50",
    "RecoTauTag_tauIdMVAnewDMwLTv1_WPEff60",
    "RecoTauTag_tauIdMVAnewDMwLTv1_WPEff70",
    "RecoTauTag_tauIdMVAnewDMwLTv1_WPEff80",
    "RecoTauTag_tauIdMVAnewDMwLTv1_WPEff90",
    "RecoTauTag_tauIdMVAnewDMwLTv1_mvaOutput_normalization",
    "RecoTauTag_tauIdMVAnewDMwoLTv1",
    "RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff40",
    "RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff50",
    "RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff60",
    "RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff70",
    "RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff80",
    "RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff90",
    "RecoTauTag_tauIdMVAnewDMwoLTv1_mvaOutput_normalization",
    "RecoTauTag_tauIdMVAoldDMwLTv1",
    "RecoTauTag_tauIdMVAoldDMwLTv1_WPEff40",
    "RecoTauTag_tauIdMVAoldDMwLTv1_WPEff50",
    "RecoTauTag_tauIdMVAoldDMwLTv1_WPEff60",
    "RecoTauTag_tauIdMVAoldDMwLTv1_WPEff70",
    "RecoTauTag_tauIdMVAoldDMwLTv1_WPEff80",
    "RecoTauTag_tauIdMVAoldDMwLTv1_WPEff90",
    "RecoTauTag_tauIdMVAoldDMwLTv1_mvaOutput_normalization",
    "RecoTauTag_tauIdMVAoldDMwoLTv1",
    "RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff40",
    "RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff50",
    "RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff60",
    "RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff70",
    "RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff80",
    "RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff90",
    "RecoTauTag_tauIdMVAoldDMwoLTv1_mvaOutput_normalization"
]

#*********************************************************************

files = []

# add all MVA configs plus WPs
for payload in payloads:

    s1 = payload
    s2 = payload
    k1 = p1.sub( s1, fileContents )
    k2 = p2.sub( s2, k1 )
    k2outfile = s2 + '.txt'
    print '--------------------------------------'
    print 'ORCOFF File for payload : ' + s2
    print 'Written to ' + k2outfile
    FILE = open(k2outfile,"w")
    FILE.write(k2)       
    FILE.close() #shake (or save) before use! 
    files.append( k2outfile )
    
for ifile in files :
    if options.prep :
        append = '_test'
    else :
         append = ''
    s = "./dropBoxOffline" + append + ".sh "+options.sqlite_file+" " + ifile
    print s
    subprocess.call([ "./dropBoxOffline" + append + ".sh", options.sqlite_file, ifile ])
  
