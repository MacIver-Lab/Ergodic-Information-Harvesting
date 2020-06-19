import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm, rice


class ObservationModel(object):
    @staticmethod
    def Gaussian(sigma, num_samples=101):
        sampGrid = np.linspace(0.0, 1.0, num_samples)
        mm = norm.pdf(sampGrid, 0.5, sigma)
        mm /= mm.max()
        return interp1d(sampGrid, mm)

    @staticmethod
    def Rice(sigma=0.2, num_samples=101):
        sampGrid = np.linspace(0.0, 1.0, num_samples)
        mm = rice.pdf(sampGrid, 0.0, 0.0, sigma)
        mm /= mm.max()
        return interp1d(sampGrid, mm)


class BaseTarget(object):
    def __init__(self, initial_position):
        self.initial_position = initial_position

    def measure(self, sensor_position):
        offset = sensor_position - self.initial_position + 0.5
        offset = np.clip(offset, 0.0, 1.0)
        return self.observationModel(offset)


class Distractor(BaseTarget):
    def __init__(self, initial_position, sigma=0.2):
        super(Distractor, self).__init__(initial_position)
        self.observationModel = ObservationModel.Rice(sigma)


class RealTarget(BaseTarget):
    def __init__(self, initial_position, sigma=0.06):
        super(RealTarget, self).__init__(initial_position)
        self.observationModel = ObservationModel.Gaussian(sigma)
