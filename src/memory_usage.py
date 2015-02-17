
import resource

def memory_usage():
    # # return the memory usage in GB
    # import psutil, os
    # process = psutil.Process(os.getpid())
    # return process.get_memory_info()[0] / float(10 ** 9)

    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / float(10 ** 9)
