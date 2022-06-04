from .slateratom import SlaterAtom, SlaterAtomSettings
from .sktwocnt import Sktwocnt, SktwocntSettings

__all__ = [ "ONECENTER_CALCULATORS", "ONECENTER_CALCULATOR_SETTINGS",
            "TWOCENTER_CALCULATORS", "TWOCENTER_CALCULATOR_SETTINGS" ]


class RegisteredCalculator:

    def __init__(self, settings, calculator):
        self.settings = settings
        self.calculator = calculator


ONECENTER_CALCULATOR_SETTINGS = {
    "slateratom": SlaterAtomSettings
}

ONECENTER_CALCULATORS = {
    RegisteredCalculator(SlaterAtomSettings, SlaterAtom)
}

TWOCENTER_CALCULATOR_SETTINGS = {
    "sktwocnt": SktwocntSettings
}

TWOCENTER_CALCULATORS = {
    RegisteredCalculator(SktwocntSettings, Sktwocnt)
}