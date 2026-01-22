import numpy as np
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms
from matplotlib.ticker import FuncFormatter

class AsymScale(mscale.ScaleBase):
    name = 'asym'

    def __init__(self, axis, **kwargs):
        mscale.ScaleBase.__init__(self)
        self.a = kwargs.get("a", 1)

    def get_transform(self):
        return self.AsymTrans(self.a)

    def set_default_locators_and_formatters(self, axis):
        # possibly, set a different locator and formatter here.
        fmt = lambda x,pos: "{}".format(np.abs(x))
        axis.set_major_formatter(FuncFormatter(fmt))

    class AsymTrans(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, a):
            mtransforms.Transform.__init__(self)
            self.a = a

        def transform_non_affine(self, x):
            return (x >= 0)*x + (x < 0)*x*self.a

        def inverted(self):
            return AsymScale.InvertedAsymTrans(self.a)

    class InvertedAsymTrans(AsymTrans):

        def transform_non_affine(self, x):
            return (x >= 0)*x + (x < 0)*x/self.a
        def inverted(self):
            return AsymScale.AsymTrans(self.a)