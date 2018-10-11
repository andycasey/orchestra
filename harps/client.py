



def ingest_observation(observation):

    # request it from ESO.

    # download it. unpack it.

    # store it.

    # update database.

    # return it.

    raise None




class Observation(object):

    @property
    def ra(self):
        return None

    @property
    def dec(self):
        return None

    @property
    def instrument(self):
        return None


    @lazy_ingest
    def get_spectrum(self, wavelength_lower=None, wavelength_upper=None):

        # if this exists in the orchestra, great. otherwise we have to get it
        # from eso, ingest it, and deliver it to you.

        pass


    @lazy_ingest
    def get_metadata(self, keys):
        pass


    @lazy_ingest
    def get_calibration_frames(self, keys):
        pass


    @property
    def in_orchestra(self):
        """ is this observation already ingested into orchestra? """
        return False

    @property
    def calibrations_in_orchestra(self):
        return False






class Harps(object):

    def __init__(self, credentials):

        self.credentials = credentials


    def query(self, object_name):

        # return many observations.
        # yield these from many pages?

        pass



    def map_reduce(self, callable, filter=None):
        """
        perform work on spectra (without having to download spectra)
        """

        yield callable(observation)





def callable(observation):

    if observation.object_name not in ("HD122563", "HD 122563"):
        return None


    spectrum = observation.get_spectrum()
    rv = observation.get_metadata("RV")

    return (spectrum.sum(), rv)
