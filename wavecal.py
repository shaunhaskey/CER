#!/task/imd/anaconda/bin/python
import pexpect, time, glob, os
import pidly
import sys

help_txt=''
start_shot = int(sys.argv[1])
end_shot = int(sys.argv[2])

#Set defaults
run_wavecal=False
write_mds=False

if len(sys.argv)>=3:
    run_wavecal=True if (sys.argv[3].lower()=='true') else False
else:
    run_wavecal = True

if len(sys.argv)>=4:
    write_mds=True if (sys.argv[4].lower()=='true') else False
else:
    write_mds = False

#Location of cerfit setup files
setup_dir='/u/grierson/cerfit/wavecal/setup/'

#write_mds = False
#run_wavecal = True

#chord = 'm01'
#shot = 163114

#start_shot = 163180
#end_shot = start_shot + 3
#start_shot = 162705
#end_shot = 162705+260

#start_shot = 162705+260
#end_shot = 163240
#start_time = time.time()

shots = range(start_shot,end_shot+1)
print 'shots to process:', shots
print 'run_wavecal:{}, write_mds:{}'.format(run_wavecal, write_mds)

#Generate all of the commands for cerfit
cmds = []
if run_wavecal:
    HOME = os.environ['HOME']
    wavecal_dir = '{}/cerfit/wavecal/'.format(HOME)
    os.system('mkdir -p {}'.format(wavecal_dir))
    os.chdir(wavecal_dir)
    for shot in shots:
        for num in range(1,33):
            chord  = 'm{:02d}'.format(num)
            setup_file = -1
            print chord
            a = glob.glob('{setup_dir}/in_{chord}*'.format(setup_dir=setup_dir, chord=chord))
            for i in a:
                i_split = i.split('.')[0].split('_')
                if (shot>=int(i_split[2])) and (shot<=int(i_split[3])):
                    setup_file = i
                    print shot, i
                    break
            if setup_file!=-1:
                print 'Using Wavecal Setup File: ', setup_file
                new_cmds = ['use {}\n'.format(setup_file),
                            'shot={}\n'.format(shot),
                            'tplot=none\n',
                            'go\n']
                            # 'type\n',
                            # 'verbose=y\n'
                            # 'go\n',]
                #print new_cmds
                for i_tmp in new_cmds:cmds.append(i_tmp)
            else:
                print 'setup file not found....'
if run_wavecal:
    cmds.append('exit\n')
    cmd_prompt = '\*'
    child = pexpect.spawn('cerfit')
    print 'Started CERFIT child process'
    for cmd in cmds:
        returned = child.expect(cmd_prompt)
        print 'sending ', cmd, returned
        child.sendline(cmd)
    count = 0
    #Finished with cerfit, try to close it
    returned = child.expect(cmd_prompt)
    child.sendline('exit\n')
    time.sleep(4)
    #If still alive kill the process
    while child.isalive():
        print 'Killing cerfit process, attempt {}'.format(count)
        time.sleep(1)
        child.kill(9)
        count+=1
        #Give up after some time....
        if count>=10:
            print 'Unable to kill cerfit, giving up'
            break
    print 'Give the file system some time....'
    time.sleep(5)

def move_files(shot, wavecal_dir):
    wavecal_dir = '{}/cerfit/wavecal/{}'.format(HOME, shot)
    cmd = 'mv w{}* {}'.format(shot, wavecal_dir)
    print cmd
    os.system(cmd)

if run_wavecal:
    #Create the directory structure 
    print 'Creating the directory structure in HOME/cerfit/wavecal/'
    for shot in shots:
        wavecal_dir = '{}/cerfit/wavecal/{}'.format(HOME, shot)
        os.system('mkdir -p {}'.format(wavecal_dir))
    print 'Give the file system some time....'
    time.sleep(5)
    #Move the files into it
    print 'Move files from HOME/cerfit/wavecal/ to HOME/cerfit/wavecal/shot/'
    for shot in shots:
        move_files(shot, wavecal_dir)
    #Should run a check to see if the files have been moved incase the filesystem is extra slow
    #a = glob.glob('{setup_dir}/in_{chord}*'.format(setup_dir=setup_dir, chord=chord))
    #if l
if write_mds:
    print 'Putting the wavecal data into the tree'
    #Start and setup IDL
    idl = pidly.IDL('/usr/local/bin/idl')
    idl('addanon')
    idl('@/u/kaplan/cerview/cerview.idl')
    idl("get_chord,149661,'t01'")
    idl('.compile /u/grierson/idlpros/add/addbrian')
    idl('.compile /u/grierson/idlpros/add/addnovi')
    idl('addbrian')
    idl('addnovi')
    idl('.compile bst_dispersion')
    #Make the list of channels up to m32
    channels = '[{}]'.format(','.join(["'m{:02d}'".format(i) for i in range(1,33)]))
    for shot in shots:
        print '='*10, shot, '='*10
        cmd = "DALPHA_MDSPLUS_PUT_FIDUCIALS,{},{}".format(shot, channels)
        print cmd
        idl(cmd)
    idl.close()
print 'Finished'
