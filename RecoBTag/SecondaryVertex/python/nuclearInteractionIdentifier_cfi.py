import FWCore.ParameterSet.Config as cms


nuclearInteractionCandIdentifier = cms.EDProducer("NuclearInteractionCandidateIdentifier",
      primaryVertices  = cms.InputTag("offlinePrimaryVertices"),
      secondaryVertices = cms.InputTag("inclusiveCandidateSecondaryVertices"),
      selection = cms.PSet(
		# maxZ = cms.double(29.)				# maximum Z of SV to be identified as NI
		# position = cms.vdouble(0., 15.)			# vector of doubles, even entries setting the lower limits for NI id, odd entries the upper limits
		# minMass = cms.double()			# minimum and maximum vertex mass for NI id
		# maxMass = cms.double()			
		# minNtracks = cms.int32()			# minimum and maximum vertex number of daughter for NI id
		# maxNtracks = cms.int32()
		# minNctau = cms.double()			# minimum Nctau = flightDistance2D/(pt/mass*0.05) for NI id (checking for compatibility with B hadron flight length)
		# distToNI = cms.double()			# if another SV (not flagged as NI) is closer than distToNI to a NI, it will also be flagged as one
		)
)
