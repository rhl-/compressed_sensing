
import resource, platform

def memory_usage():
    # # return the memory usage in gigabytes on OS-X
    osx = ('darwin' in  platform.platform().lower())
    linux = ('linux' in platform.platform().lower())
    if( osx):
        return (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/float(1e9))
    elif( linux):
        return (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/float(1e6))
    print "ERROR: PLAT %s NOT DETERMINED PRINTING RU_MAXRSS"%(platform.platform().lower())
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss 
