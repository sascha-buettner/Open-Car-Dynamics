// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include "ocd_tire_models_cpp/mf_52.hpp"
namespace tam::ocd::tire_models
{
// The original class name is retained for binary/source compatibility. FITTYP
// selects the implemented 52/62 branch; a name alone never establishes coverage.
using MagicFormula = MF52;
}  // namespace tam::ocd::tire_models
