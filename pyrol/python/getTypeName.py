"""Wrappers around the Binder-generated ROL classes, and lookup by name.

Every class in pyrol.pyrol.ROL is replaced by a subclass created by
withTrackingConstructor. Instances of these subclasses keep references to
the PyROL objects passed to their constructor and to their add* methods:
C++ holds those objects only through Teuchos::RCP, so without these
references a Python-defined Objective or Vector could be garbage collected
while ROL still uses it.
"""
from pyrol.pyrol import ROL
import inspect
import sys

# All wrapped classes; arguments of these types are kept alive by track().
trackedTypes = []


def track(self, *args, **kwargs):
    """Keep references on self to the PyROL objects among the arguments."""
    tracked = tuple(trackedTypes)
    for arg in (*args, *kwargs.values()):
        if isinstance(arg, tracked):
            self._tracked_constructor_args.append(arg)


def tracking_method(cls_method):
    """Wrap a method so that self keeps references to its PyROL arguments."""
    def new_method(self, *args, **kwargs):
        track(self, *args, **kwargs)
        return cls_method(self, *args, **kwargs)
    return new_method


def _decorator23(value_method):
    """Let value(u, z, tol) of a SimOpt objective also be called as
    value(x, tol), where x is a Vector_SimOpt holding u and z."""
    def value(*args):
        if len(args) == 2:
            x, tol = args
            return value_method(x.get_1(), x.get_2(), tol)
        if len(args) == 3:
            u, z, tol = args
            return value_method(u, z, tol)
        raise TypeError("Unexpected number of arguments")
    return value


def _decorator34(value_method):
    """Let value(c, u, z, tol) of a SimOpt constraint also be called as
    value(c, x, tol), where x is a Vector_SimOpt holding u and z."""
    def value(*args):
        if len(args) == 3:
            c, x, tol = args
            return value_method(c, x.get_1(), x.get_2(), tol)
        if len(args) == 4:
            c, u, z, tol = args
            return value_method(c, u, z, tol)
        raise TypeError("Unexpected number of arguments")
    return value


def withTrackingConstructor(cls_name, cls):
    """Return a subclass of cls whose instances keep their PyROL arguments
    alive (see the module docstring)."""
    class newCls(cls):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

            self._tracked_constructor_args = []
            track(self, *args, **kwargs)

            # Also track the arguments of the add* methods, e.g. Problem.addConstraint.
            for name in dir(cls):
                if name.startswith('add') and callable(getattr(cls, name)):
                    setattr(newCls, name, tracking_method(getattr(cls, name)))

            if cls_name.startswith('Objective_SimOpt'):
                self.value = _decorator23(self.value)
            elif cls_name.startswith('Constraint_SimOpt'):
                self.value = _decorator34(self.value)

    return newCls


# Wrap the classes of the ROL module in place, and record all of its members.
ROL_members = {}
for cls_name, cls_obj in inspect.getmembers(sys.modules['pyrol.pyrol.ROL']):
    if inspect.isclass(cls_obj):
        cls_obj = withTrackingConstructor(cls_name, cls_obj)
        trackedTypes.append(cls_obj)
        setattr(sys.modules['pyrol.pyrol.ROL'], cls_name, cls_obj)
    ROL_members[cls_name] = (cls_obj, inspect.isclass(cls_obj))

# The wrapped classes, sorted by name.
ROL_classes = [(cls_name, cls_obj)
               for cls_name, cls_obj in inspect.getmembers(sys.modules['pyrol.pyrol.ROL'])
               if inspect.isclass(cls_obj)]


def getDefaultScalarType():
    return "double"


def _find_class(prefix, class_name):
    """Return the first class in ROL_classes whose lowercase name starts with
    prefix; otherwise warn about class_name and return None."""
    for name, cls in ROL_classes:
        if name.lower().startswith(prefix):
            return cls
    print('Warning: Unknown type "{}", the function returns None.'.format(class_name))
    return None


def getTypeName(class_name, scalar_type=getDefaultScalarType()):
    """Return the ROL class for a class name and scalar type, e.g.
    getTypeName('Objective') -> Objective_double_t (case-insensitive)."""
    return _find_class(class_name.lower() + "_" + scalar_type.lower(), class_name)


def getTypeName2(class_name):
    """Like getTypeName, but matching class_name alone."""
    return _find_class(class_name.lower(), class_name)
