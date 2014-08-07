NOTE: These directions have been updated to work with the version of DropBox made 
	  available as of November 2012.

NOTE: You must have an account with the AlCa group in order to use these directions.



1) Create the SQLlite file. 

cd RecoTauTag/TauAnalysisTools/test/
cmsRun writeMVAs_cfg.py
 
This will print out a great deal of information for debugging purposes.


2) Check the payloads that are present in the created file, as of this writing, "Jec11_V12", but change
to the latest tag as you've done above. 

cd uploadToDropbox
cmscond_list_iov -c sqlite_file:RecoTauTag_MVAs_2014Apr29.db -a

NOTE: replace RecoTauTag_MVAs_2014Apr29.db by the name of the SQLlite file that you created in step 1) 


3) Then check each of the payloads individually. There is a script to help you out called "testAllIOVs.py".

python testAllIOVs.py

Example output:
...
	===============================================================
	Tag: RecoTauTag_againstMuonMVAv1_WPeff99_0
	===============================================================
	OID: 0001-00000062
	Scope: Unknown
	Description:  
	TimeType: runnumber
	Since                 Till                  Payload OID    Payload Class       
	--------------------  --------------------  -------------  --------------------
	                   1  18446744073709551615  0002-00000051  PhysicsTGraphPayload

	Total # of payload objects: 1
	===============================================================
	Tag: RecoTauTag_againstMuonMVAv1_WPeff99_5
	===============================================================
	OID: 0001-00000061
	Scope: Unknown
	Description:  
	TimeType: runnumber
	Since                 Till                  Payload OID    Payload Class       
	--------------------  --------------------  -------------  --------------------
	                   1  18446744073709551615  0002-00000050  PhysicsTGraphPayload

	Total # of payload objects: 1
	=================================================================
	Tag: RecoTauTag_againstMuonMVAv1_mvaOutput_normalization
	=================================================================
	OID: 0001-00000064
	Scope: Unknown
	Description:  
	TimeType: runnumber
	Since                 Till                  Payload OID    Payload Class         
	--------------------  --------------------  -------------  ----------------------
	                   1  18446744073709551615  0003-00000000  PhysicsTFormulaPayload

	Total # of payload objects: 1
	===========================================================
	Tag: RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL
	===========================================================
	OID: 0001-00000008
	Scope: Unknown
	Description:  
	TimeType: runnumber
	Since                 Till                  Payload OID    Payload Class
	--------------------  --------------------  -------------  -------------
	                   1  18446744073709551615  0000-00000008      GBRForest

	Total # of payload objects: 1
...

