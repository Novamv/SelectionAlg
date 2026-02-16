# ----------- SELECTION ALG -----------
time python ${SELECTIONALGROOT}/share/run.py --input /sps/juno/mlecocq/Data/Physics/RUN.9848.JUNODAQ.Physics.ds-2.global_trigger.20250902183819.121_J25.5.0.esd \
--output IBD.root --recEDMPath /Event/CdVertexRecOMILREC_JVtx

# ----------- SNIPERTOPLAINTREE -----------
# time python ${SNIPERTOPLAINTREEROOT}/share/run.py --input /sps/juno/mlecocq/Data/Physics/RUN.9848.JUNODAQ.Physics.ds-2.global_trigger.20250902181850.015_J25.5.0.esd --output RUN_9848_plain.root \
# --recEDMPath /Event/CdVertexRecOMILREC_JVtx --saveMuon