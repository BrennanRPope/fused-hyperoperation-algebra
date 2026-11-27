import resource, runpy, sys
sys.argv = sys.argv[1:]
try:
    runpy.run_path(sys.argv[0], run_name='__main__')
finally:
    print('peak RSS MB: %.0f' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
