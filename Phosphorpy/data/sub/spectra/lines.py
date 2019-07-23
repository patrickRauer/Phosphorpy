from astropy import units as u


class SpectraLine:

    name = None
    lambda0 = None
    label = None

    def __init__(self, name, lambda0, label=None):
        self.name = name
        self.lambda0 = lambda0
        if label is None:
            label = name
        self.label = label

    def get_lambda0(self, unit=None):
        if isinstance(self.lambda0, u.Unit):
            return self.lambda0.to(unit)
        else:
            return (self.lambda0*u.AA).to(unit)
