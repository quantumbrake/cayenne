from ..utils_cython import sumfunc


def direct_cython(arg1, arg2):
    print(sumfunc(arg1, arg2))
