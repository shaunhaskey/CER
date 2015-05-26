import pexpect, time
chord = 't01'
shot = '158672'
list_of_cmds = []
cmds = ["use in_{}.dat,shot={}\n".format(chord, shot),
        "write in_{}.dat\n".format(chord),]
#        "exit\n"]
cmd_prompt = '\*'
child = pexpect.spawn('cerfitw')
print 'expecting'

for cmd in cmds:
    i = child.expect(cmd_prompt)
    print 'sending ', cmd, i
    child.sendline(cmd)
with file('t01/30lt/tssub.dat') as filehandle:
    cmds = filehandle.readlines()

for cmd in cmds:
    i = child.expect(cmd_prompt)
    print 'sending ', cmd, i
    child.sendline(cmd)

print 'finished expecting'

