All component types are are defined by 5 things. first off are 2 integers that determine
the number of terminals the component has, and how many parameters describe its behaviour:

const int <component_type>_terminals_count;
const int <component_type>_parameters_count;

As this is a non-linear circuit simulator, a component has an entirely custom curve of current
vs voltage. This function uses the voltage at the terminals to determine the currents entering
those terminals:

void <component_type>_currentCurve(const double *parameters, const double *v, double timestep, double *i);

Make note of the current direction being Into the component, so for example a resistor would have one
terminal current be positive where its voltage is higher, and the other terminal current would be negative.
All components in the model are assumed to obey kirchoff's current law.

In order to perform the multidimensional newtons method, the jacobian of the terminal currents with
respect to voltage must be calculated for each component. j is arranged so that the first n elements
are the derivative of current w.r.t the 0'th voltage, then the next n are w.r.t the 1'st voltage, then the
2'nd, etc.

void <component_type>_jacobian(const double *parameters, const double *v, double timestep, double *j);

In order to enable time-domain simulation, we store the currents and voltages of the previous timestep
in the parameter space of the component for the next timestep.

void <component_type>_updateState(double *parameters, const double *v, double timestep, const double *i);
