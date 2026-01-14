# python run.py --input-correlations /sps/juno/mlecocq/RUN.10270.JUNODAQ.Physics.ds-2.global_trigger.20250922075923.2563_J25.5.0.rtraw \
# --input /sps/juno/mlecocq/RUN.10270.JUNODAQ.Physics.ds-2.global_trigger.20250922075923.2563_J25.5.0.esd \
# --recEDMPath /Event/CdVertexRecMixedPhase


# ----------- SNIPERTOPLAINTREE -----------
# python ${SNIPERTOPLAINTREEROOT}/share/run.py --input-correlations /sps/juno/mlecocq/Data/Physics/RUN.9789.JUNODAQ.Physics.ds-2.global_trigger.20250830102924.314_J25.5.0.rtraw \
# --input /sps/juno/mlecocq/Data/Physics/RUN.9789.JUNODAQ.Physics.ds-2.global_trigger.20250830102924.314_J25.5.0.esd \
# --recEDMPath /Event/CdVertexRecOMILREC

# ----------- SELECTION ALG -----------

# 1 IBD here
time python run.py --input-correlations /sps/juno/mlecocq/Data/Physics/RUN.9789.JUNODAQ.Physics.ds-2.global_trigger.20250830102924.314_J25.5.0.rtraw \
--input /sps/juno/mlecocq/Data/Physics/RUN.9789.JUNODAQ.Physics.ds-2.global_trigger.20250830102924.314_J25.5.0.esd \
--output RUN9789_IBD.root --recEDMPath /Event/CdVertexRecOMILREC
