import MDSplus as mds
import os 
calib_shot = 161800
shot = 163114
comment = '???'
os.environ['cermain_path']='atlas.gat.com::'
T = mds.Tree('cermain',shot)
node = T.getNode('.calibration.comments')
node.putData('hello')

T.getNode('calib_shot').put(calib_shot)
T.getNode('.calibration:comments').put(comments)

for i in range(len(chords)):
    mds_path = '.calibration.tangential.channel{:02d}:fiducial'
    #wc = get the fiducial somehow?
    #Check there is no error, then write the data
    error = False
    if not error:
        T.getNode('calib_shot').put(wc)


  #     IF ~wc.ierr THEN BEGIN
  #         PRINT,'Writing '+mds_tag
  #         MDSPUT,mds_tag,'$',wc.fiducial

  # FOR i=0,N_ELEMENTS(chords)-1 DO BEGIN
  #     mds_channel = 'channel'+STRMID(chords[i],1,2)
  #     mds_tag = '\top.calibration.tangential.'+mds_channel+':fiducial'
  #     wc = DALPHA_GET_CERFIT_WAVECAL(shot,chords[i])

#MDSPUT, '\top:calib_shot', '$', LONG(shot)
#MDSPUT, '\top.calibration:comments', '$', com

#server=MDSplus.Connection('atlas.gat.com')
#server=mds.Connection('atlas.gat.com')
#server.openTree('cermain',163114)
#server.put('\top.calibration:comments','$', 'this is a test')
#Brian's way of doing it
  # MDSCONNECT,'atlas'
  # MDSOPEN,'cermain',shot

  # ;; For text and numeric nodes
  # ;; mdsput, NODE_NAME, '$', NODE_CONTENT
  # MDSPUT, '\top:calib_shot', '$', LONG(shot)
  # MDSPUT, '\top.calibration:comments', '$', com

  # FOR i=0,N_ELEMENTS(chords)-1 DO BEGIN
  #     mds_channel = 'channel'+STRMID(chords[i],1,2)
  #     mds_tag = '\top.calibration.tangential.'+mds_channel+':fiducial'
  #     wc = DALPHA_GET_CERFIT_WAVECAL(shot,chords[i])
  #     IF ~wc.ierr THEN BEGIN
  #         PRINT,'Writing '+mds_tag
  #         MDSPUT,mds_tag,'$',wc.fiducial
  #     ENDIF
  # ENDFOR

  # ;; Close tree
  # MDSCLOSE


  # ;; Now open and print out fiducials
  
  # MDSOPEN,'cermain',shot
  # FOR i=0,N_ELEMENTS(chords)-1 DO BEGIN
  #     mds_channel = 'channel'+STRMID(chords[i],1,2)
  #     mds_tag = '\top.calibration.tangential.'+mds_channel+':fiducial'
  #     fid = MDSVALUE(mds_tag,STATUS=status)
  #     IF status EQ 1 THEN $
  #       PRINT,'Fiducial: '+STRTRIM(fid,2) $
  #     ELSE $
  #       PRINT,'No fiducial for '+mds_channel
      
  # ENDFOR
  # MDSCLOSE
