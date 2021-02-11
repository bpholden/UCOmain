import subprocess

MACHINELST = ["dresden","mainz","hamburg","warsaw"]


def generic(machine,services,command):
    if command in ['start','restart','stop','status'] and machine in MACHINELST:
        servicestr = "apf %s " % (command)
        servicestr += " ".join(services)
        cmdlst = ["ssh","-f",machine, servicestr]
        try:
            p = subprocess.check_output(cmdlst)
            return(True,p)
        except Exception, e:
            return(False,e)
    else:
        s = "%s is not in the list of correct machines" % (machine)
        return(False,s)

def restart(machine,services):
    rv, s = generic(machine,services,'restart')
    return rv, s 

def status(machine,services):
    rv, s = generic(machine,services,'status')
    return rv, s 

