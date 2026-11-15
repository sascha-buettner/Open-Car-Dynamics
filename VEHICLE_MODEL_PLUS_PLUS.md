# Physical model configuration

The standalone double-track model and MF5.2 tire retain the generic public parameter values from `main`. They can be constructed and evaluated without an application parameter file. Newly added contact, belt, suspension, MF6.2, and recovery mechanisms are disabled unless explicitly selected. Their unset values are inert placeholders and do not represent a calibrated vehicle.

Application-specific vehicle, tire, and track parameters must be supplied by the application from its private configuration. No private parameter set is embedded in this repository.
