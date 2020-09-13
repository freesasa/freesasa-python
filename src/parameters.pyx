from cfreesasa cimport *

## Used to specify the algorithm by Shrake & Rupley
ShrakeRupley = 'ShrakeRupley'

## Used to specify the algorithm by Lee & Richards
LeeRichards = 'LeeRichards'

cdef class Parameters:
    """
    Stores parameter values to be used by calculation.

    Default parameters are ::

        Parameters.defaultParameters = {
            'algorithm'    : LeeRichards,
            'probe-radius' : freesasa_default_parameters.probe_radius,
            'n-points'     : freesasa_default_parameters.shrake_rupley_n_points,
            'n-slices'     : freesasa_default_parameters.lee_richards_n_slices,
            'n-threads'    : freesasa_default_parameters.n_threads
        }

    Attributes:
        defaultParamers (dict): The default parameters
    """


    cdef freesasa_parameters _c_param

    defaultParameters = {
        'algorithm'    : LeeRichards,
        'probe-radius' : freesasa_default_parameters.probe_radius,
        'n-points'     : freesasa_default_parameters.shrake_rupley_n_points,
        'n-slices'     : freesasa_default_parameters.lee_richards_n_slices,
        'n-threads'    : freesasa_default_parameters.n_threads
    }

    def __init__(self,param=None):
        """
        Initializes Parameters object.

        Args:
              param (dict): optional argument to specify parameter-values,
                  see :py:attr:`.Parameters.defaultParameters`.

        Raises:
              AssertionError: Invalid parameter values supplied
        """
        self._c_param = freesasa_default_parameters
        if param != None:
            if 'algorithm' in param:    self.setAlgorithm(param['algorithm'])
            if 'probe-radius' in param: self.setProbeRadius(param['probe-radius'])
            if 'n-points' in param:     self.setNPoints(param['n-points'])
            if 'n-slices' in param:     self.setNSlices(param['n-slices'])
            if 'n-threads' in param:    self.setNThreads(param['n-threads'])
            unknownKeys = []
            for key in param:
                if not key in self.defaultParameters:
                    unknownKeys.append(key)
                if len(unknownKeys) > 0:
                    raise AssertionError('Key(s): ',unknownKeys,', unknown')

    def setAlgorithm(self,alg):
        """
        Set algorithm.

        Args:
             alg (str): algorithm name, only allowed values are
             :py:data:`freesasa.ShrakeRupley` and :py:data:`freesasa.LeeRichards`

        Raises:
             AssertionError: unknown algorithm specified
        """
        if alg == ShrakeRupley:
            self._c_param.alg = FREESASA_SHRAKE_RUPLEY
        elif alg == LeeRichards:
            self._c_param.alg = FREESASA_LEE_RICHARDS
        else:
            raise AssertionError("Algorithm '%s' is unknown" % alg)

    def algorithm(self):
        """
        Get algorithm.

        Returns:
            str: Name of algorithm
        """
        if self._c_param.alg == FREESASA_SHRAKE_RUPLEY:
            return ShrakeRupley
        if self._c_param.alg == FREESASA_LEE_RICHARDS:
            return LeeRichards
        raise Exception("No algorithm specified, shouldn't be possible")

    def setProbeRadius(self,r):
        """
        Set probe radius.

        Args:
            r (float): probe radius in Å (>= 0)

        Raises:
            AssertionError: r < 0
        """
        assert(r >= 0)
        self._c_param.probe_radius = r

    def probeRadius(self):
        """
        Get probe radius.

        Returns:
             float: Probe radius in Å
        """
        return self._c_param.probe_radius

    def setNPoints(self,n):
        """
        Set number of test points in Shrake & Rupley algorithm.

        Args:
            n (int): Number of points (> 0).

        Raises:
            AssertionError: n <= 0.
        """
        assert(n > 0)
        self._c_param.shrake_rupley_n_points = n

    def nPoints(self):
        """
        Get number of test points in Shrake & Rupley algorithm.

        Returns:
            int: Number of points.
        """
        return self._c_param.shrake_rupley_n_points

    def setNSlices(self,n):
        """
        Set the number of slices per atom in Lee & Richards algorithm.

        Args:
            n (int): Number of slices (> 0)

        Raises:
            AssertionError: n <= 0
        """
        assert(n> 0)
        self._c_param.lee_richards_n_slices = n

    def nSlices(self):
        """
        Get the number of slices per atom in Lee & Richards algorithm.

        Returns:
            int: Number of slices.
        """
        return self._c_param.lee_richards_n_slices

    def setNThreads(self,n):
        """
        Set the number of threads to use in calculations.

        Args:
            n (int): Number of points (> 0)

        Raises:
            AssertionError: n <= 0
        """
        assert(n>0)
        self._c_param.n_threads = n

    def nThreads(self):
        """
        Get the number of threads to use in calculations.

        Returns:
            int: Number of threads.
        """
        return self._c_param.n_threads

    # not pretty, but only way I've found to pass pointers around
    def _get_address(self, size_t ptr2ptr):
        cdef freesasa_parameters **p = <freesasa_parameters**> ptr2ptr
        p[0] = &self._c_param