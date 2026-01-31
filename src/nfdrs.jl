"""
    National Fire Danger Rating System (NFDRS) Fuel Moisture Models

Implementation of the NFDRS equations from:
Cohen, Jack D.; Deeming, John E. "The National Fire-Danger Rating System: basic equations."
Gen. Tech. Rep. PSW-82. Berkeley, CA: Pacific Southwest Forest and Range Experiment Station,
Forest Service, U.S. Department of Agriculture; 1985. 16 p.

This module implements:
- Equilibrium Moisture Content (EMC) calculations (Eq. 1a, 1b, 1c)
- Dead fuel moisture models (1-hr, 10-hr, 100-hr, 1000-hr timelag)
- Live fuel moisture models (herbaceous and woody)
- Fuel loading transfer (Eq. 5-8)
- Fire characteristics (Spread Component, Energy Release Component, Burning Index)
- Fire occurrence indexes (Ignition Component, Human-Caused and Lightning-Caused Fire Occurrence)
- Fuel model database (Appendix, 20 models A-U excluding M)

## Unit Convention

This implementation uses SI units throughout for inputs and outputs:
- Temperature: K (Kelvin)
- Fuel loading: kg/m²
- Fuel bed depth: m (meters)
- Surface-area-to-volume ratio: m⁻¹
- Heat of combustion: J/kg
- Wind speed: m/s
- Rate of spread: m/s
- Moisture content: dimensionless fraction (0-1)
- Density: kg/m³
- Precipitation duration: s (seconds)
- Latitude: rad (radians)

The original NFDRS equations use imperial units. The empirical coefficients have been
converted to SI units where possible, and unit reference constants are used for proper
dimensional analysis with ModelingToolkit.
"""

# t and D are imported in the main WildlandFire module

# =============================================================================
# Unit Conversion Constants (for fuel model data - these are numerical factors)
# =============================================================================

# Note: These are used for static fuel model data conversion only, not for symbolic computation.
# For symbolic computations, we use @constants with proper unit annotations.
const _FT_TO_M = 0.3048
const _M_TO_FT = 1.0 / _FT_TO_M
const _MI_TO_M = 1609.344
const _MPH_TO_MS = 0.44704
const _MS_TO_FPM = 196.8504
const _LB_PER_FT3_TO_KG_PER_M3 = 16.01846
const _LB_PER_FT2_TO_KG_PER_M2 = 4.88243
const _KG_PER_M2_TO_LB_PER_FT2 = 1.0 / _LB_PER_FT2_TO_KG_PER_M2
const _BTU_PER_LB_TO_J_PER_KG = 2326.0
const _CAL_PER_G_TO_J_PER_KG = 4184.0
const _TONS_PER_ACRE_TO_KG_PER_M2 = 0.2241702

# =============================================================================
# Equilibrium Moisture Content (EMC)
# =============================================================================

"""
    EquilibriumMoistureContent(; name=:EquilibriumMoistureContent)

Calculate equilibrium moisture content (EMC) from temperature and relative humidity.

Based on regression equations developed by Simard (1968) from Wood Handbook tables.
EMC is expressed as fraction moisture content (dimensionless).

Implements Equations 1a, 1b, 1c from Cohen & Deeming (1985), page 1.

# Parameters
- `TEMP`: Dry bulb temperature (K)
- `RH`: Relative humidity (dimensionless, 0-1)

# Variables
- `EMC`: Equilibrium moisture content (dimensionless)
"""
@component function EquilibriumMoistureContent(; name=:EquilibriumMoistureContent)
    @constants begin
        # Unit reference constants for proper dimensional analysis
        K_ref = 1.0, [description = "Reference temperature for non-dimensionalization", unit = u"K"]
        K_to_F_offset = 459.67, [description = "Kelvin to Fahrenheit offset (dimensionless)", unit = u"1"]
    end

    @parameters begin
        TEMP, [description = "Dry bulb temperature", unit = u"K"]
        RH, [description = "Relative humidity (dimensionless)", unit = u"1"]
    end

    @variables begin
        EMC(t), [description = "Equilibrium moisture content (dimensionless)", unit = u"1"]
    end

    # Convert temperature to Fahrenheit for the empirical equations
    # Original equations use °F and RH as percentage
    # Non-dimensionalize TEMP first, then convert: °F = K × 9/5 - 459.67
    TEMP_F = (TEMP / K_ref) * 9 / 5 - K_to_F_offset
    RH_pct = RH * 100

    eqs = [
        # EMC equations based on RH ranges (Eq. 1a, 1b, 1c)
        EMC ~ ifelse(RH < 0.10,
            # Eq. 1a: RH < 10%
            (0.03229 + 0.281073 * RH_pct - 0.000578 * TEMP_F * RH_pct) / 100,
            ifelse(RH < 0.50,
                # Eq. 1b: 10% <= RH < 50%
                (2.22749 + 0.160107 * RH_pct - 0.014784 * TEMP_F) / 100,
                # Eq. 1c: RH >= 50%
                (21.0606 + 0.005565 * RH_pct^2 - 0.00035 * RH_pct * TEMP_F - 0.483199 * RH_pct) / 100
            )
        )
    ]

    return System(eqs, t; name)
end

# =============================================================================
# 1-Hour Timelag Fuel Moisture
# =============================================================================

"""
    OneHourFuelMoisture(; name=:OneHourFuelMoisture)

Calculate 1-hour timelag fuel moisture content.

The response of 1-hour timelag fuels to environmental changes is so rapid that only
the potential moisture content (equivalent to EMC at fuel-atmosphere interface) is required.

From Cohen & Deeming (1985), page 3:
- When fuel moisture sticks are not used: MC1 = 1.03 × EMCPRM
- When fuel moisture sticks are used: MC1 = (4.0 × EMCPRM + MC10) / 5.0
- If raining at observation time: MC1 = 35%

# Parameters
- `EMCPRM`: EMC at fuel-atmosphere interface (dimensionless)
- `MC10`: 10-hour fuel moisture content (dimensionless), used when fuel sticks are used
- `use_fuel_sticks`: Flag (1 if fuel sticks used, 0 otherwise) (dimensionless)
- `is_raining`: Flag (1 if raining, 0 otherwise) (dimensionless)

# Variables
- `MC1`: 1-hour fuel moisture content (dimensionless)
"""
@component function OneHourFuelMoisture(; name=:OneHourFuelMoisture)
    @parameters begin
        EMCPRM, [description = "EMC at fuel-atmosphere interface (dimensionless)", unit = u"1"]
        MC10, [description = "10-hour fuel moisture content (dimensionless)", unit = u"1"]
        use_fuel_sticks, [description = "Flag: 1 if fuel sticks used, 0 otherwise (dimensionless)", unit = u"1"]
        is_raining, [description = "Flag: 1 if raining, 0 otherwise (dimensionless)", unit = u"1"]
    end

    @variables begin
        MC1(t), [description = "1-hour fuel moisture content (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Page 3-4: MC1 calculation based on conditions
        MC1 ~ ifelse(is_raining > 0.5,
            0.35,  # 35% if raining (page 4)
            ifelse(use_fuel_sticks > 0.5,
                (4.0 * EMCPRM + MC10) / 5.0,  # With fuel sticks (page 3)
                1.03 * EMCPRM  # Without fuel sticks (page 3)
            )
        )
    ]

    return System(eqs, t; name)
end

# =============================================================================
# 10-Hour Timelag Fuel Moisture
# =============================================================================

"""
    TenHourFuelMoisture(; name=:TenHourFuelMoisture)

Calculate 10-hour timelag fuel moisture content.

From Cohen & Deeming (1985), page 4:
- Without fuel sticks: MC10 = 1.28 × EMCPRM
- With fuel sticks: MC10 = AA×CC + BB×CC×(WT - 100.0)
  where age correction factors are based on Haines and Frost (1978):
  - AA = 0.5 × AGE / 30.0
  - BB = 1.0 + (0.02 × AGE / 30.0)
  - CC = CLIMAT / 4.0

# Parameters
- `EMCPRM`: EMC at fuel-atmosphere interface (dimensionless)
- `use_fuel_sticks`: Flag (1 if fuel sticks used, 0 otherwise) (dimensionless)
- `WT`: Weight of fuel sticks (kg)
- `AGE`: Days since sticks were set out (dimensionless)
- `CLIMAT`: NFDRS climate class 1-4 (dimensionless)

# Variables
- `MC10`: 10-hour fuel moisture content (dimensionless)
"""
@component function TenHourFuelMoisture(; name=:TenHourFuelMoisture)
    @constants begin
        kg_ref = 1.0, [description = "Reference mass for non-dimensionalization", unit = u"kg"]
        kg_to_g = 1000.0, [description = "Kilograms to grams conversion (dimensionless)", unit = u"1"]
        WT_ref_g = 100.0, [description = "Reference fuel stick weight (dimensionless)", unit = u"1"]
    end

    @parameters begin
        EMCPRM, [description = "EMC at fuel-atmosphere interface (dimensionless)", unit = u"1"]
        use_fuel_sticks, [description = "Flag: 1 if fuel sticks used, 0 otherwise (dimensionless)", unit = u"1"]
        WT = 0.1, [description = "Weight of fuel sticks", unit = u"kg"]
        AGE = 0.0, [description = "Days since sticks were set out (dimensionless)", unit = u"1"]
        CLIMAT = 2.0, [description = "NFDRS climate class 1-4 (dimensionless)", unit = u"1"]
    end

    @variables begin
        MC10(t), [description = "10-hour fuel moisture content (dimensionless)", unit = u"1"]
    end

    # Age correction factors (page 4)
    AA = 0.5 * AGE / 30.0
    BB = 1.0 + (0.02 * AGE / 30.0)
    CC = CLIMAT / 4.0

    # Weight difference: non-dimensionalize WT, convert to grams, subtract reference 100g
    WT_grams = (WT / kg_ref) * kg_to_g
    WT_diff = WT_grams - WT_ref_g

    eqs = [
        # Page 4: MC10 calculation
        MC10 ~ ifelse(use_fuel_sticks > 0.5,
            # With fuel sticks (age corrected)
            (AA * CC + BB * CC * WT_diff) / 100.0,  # Convert to fraction
            # Without fuel sticks
            1.28 * EMCPRM
        )
    ]

    return System(eqs, t; name)
end

# =============================================================================
# 100-Hour Timelag Fuel Moisture
# =============================================================================

"""
    HundredHourFuelMoisture(; name=:HundredHourFuelMoisture)

Calculate 100-hour timelag fuel moisture content using daily discrete update.

Uses weighted 24-hour average EMC (EMCBAR) based on hours of daylight.
The moisture content responds to the boundary condition with a response coefficient.

From Cohen & Deeming (1985), page 4-5:
- EMCBAR = (DAYLIT × EMCMIN + (24 - DAYLIT) × EMCMAX) / 24
- BNDRYH = ((24 - PPTDUR) × EMCBAR + PPTDUR × (0.5×PPTDUR + 41)) / 24
- MC100 = YMC100 + (BNDRYH - YMC100) × (1.0 - 0.87 × exp(-0.24))

# Parameters
- `EMCMIN`: EMC at max temp/min RH (dimensionless)
- `EMCMAX`: EMC at min temp/max RH (dimensionless)
- `PPTDUR`: Duration of precipitation (s)
- `LAT`: Station latitude (rad)
- `JDATE`: Julian day of year 1-366 (dimensionless)
- `YMC100`: Previous day's MC100 value (dimensionless)

# Variables
- `MC100`: 100-hour fuel moisture content (dimensionless)
- `EMCBAR`: Weighted 24-hour average EMC (dimensionless)
- `BNDRYH`: Weighted 24-hour average boundary condition (dimensionless)
"""
@component function HundredHourFuelMoisture(; name=:HundredHourFuelMoisture)
    @constants begin
        s_ref = 1.0, [description = "Reference time for non-dimensionalization", unit = u"s"]
        rad_ref = 1.0, [description = "Reference angle for non-dimensionalization", unit = u"rad"]
        hr_to_s = 3600.0, [description = "Hours to seconds conversion (dimensionless)", unit = u"1"]
        response_coef = 1.0 - 0.87 * exp(-0.24), [description = "Moisture response coefficient (dimensionless)", unit = u"1"]
    end

    @parameters begin
        EMCMIN, [description = "EMC at max temp/min RH (dimensionless)", unit = u"1"]
        EMCMAX, [description = "EMC at min temp/max RH (dimensionless)", unit = u"1"]
        PPTDUR = 0.0, [description = "Duration of precipitation", unit = u"s"]
        LAT = 0.6981, [description = "Station latitude", unit = u"rad"]
        JDATE = 180.0, [description = "Julian day of year 1-366 (dimensionless)", unit = u"1"]
        YMC100, [description = "Previous day's MC100 value (dimensionless)", unit = u"1"]
    end

    @variables begin
        MC100(t), [description = "100-hour fuel moisture content (dimensionless)", unit = u"1"]
        EMCBAR(t), [description = "Weighted 24-hour average EMC (dimensionless)", unit = u"1"]
        BNDRYH(t), [description = "Weighted 24-hour average boundary condition (dimensionless)", unit = u"1"]
    end

    # Convert precipitation duration from seconds to hours for equation (non-dimensionalize first)
    PPTDUR_hr = (PPTDUR / s_ref) / hr_to_s

    # Daylength calculation (page 4-5) - non-dimensionalize LAT for trig functions
    LAT_nd = LAT / rad_ref
    DECL = 0.41008 * sin((JDATE - 82) * 0.01745)
    tan_product = tan(LAT_nd) * tan(DECL)
    tan_product_clamped = max(-0.99, min(0.99, tan_product))
    DAYLIT = 24.0 * (1.0 - acos(tan_product_clamped) / 3.1416)

    eqs = [
        # Weighted 24-hour average EMC (page 5)
        EMCBAR ~ (DAYLIT * EMCMIN + (24.0 - DAYLIT) * EMCMAX) / 24.0,

        # Weighted 24-hour average boundary condition (page 5)
        BNDRYH ~ ((24.0 - PPTDUR_hr) * EMCBAR + PPTDUR_hr * (0.5 * PPTDUR_hr + 41.0) / 100.0) / 24.0,

        # 100-hour fuel moisture daily update equation (page 5)
        MC100 ~ YMC100 + (BNDRYH - YMC100) * response_coef
    ]

    return System(eqs, t; name)
end

# =============================================================================
# 1000-Hour Timelag Fuel Moisture
# =============================================================================

"""
    ThousandHourFuelMoisture(; name=:ThousandHourFuelMoisture)

Calculate 1000-hour timelag fuel moisture content using daily discrete update.

Uses 7-day running average of boundary conditions (BDYBAR).
From Cohen & Deeming (1985), page 5:
- BNDRYT = ((24 - PPTDUR) × EMCBAR + PPTDUR × (2.7×PPTDUR + 76)) / 24
- MC1000 = PM1000 + (BDYBAR - PM1000) × (1.0 - 0.82 × exp(-0.168))

Note: In the original NFDRS, BNDRYT is calculated daily but BDYBAR and MC1000
are calculated only every seventh day. This implementation calculates daily
with the running average provided as a parameter.

# Parameters
- `EMCBAR`: Weighted 24-hour average EMC (dimensionless)
- `PPTDUR`: Duration of precipitation (s)
- `BDYBAR`: 7-day running average boundary condition (dimensionless)
- `PM1000`: Previous MC1000 value (dimensionless) - MC1000 from 7 days ago in original NFDRS

# Variables
- `MC1000`: 1000-hour fuel moisture content (dimensionless)
- `BNDRYT`: Weighted 24-hour average boundary condition for 1000-hr (dimensionless)
"""
@component function ThousandHourFuelMoisture(; name=:ThousandHourFuelMoisture)
    @constants begin
        s_ref = 1.0, [description = "Reference time for non-dimensionalization", unit = u"s"]
        hr_to_s = 3600.0, [description = "Hours to seconds conversion (dimensionless)", unit = u"1"]
        response_coef = 1.0 - 0.82 * exp(-0.168), [description = "Moisture response coefficient (dimensionless)", unit = u"1"]
    end

    @parameters begin
        EMCBAR, [description = "Weighted 24-hour average EMC (dimensionless)", unit = u"1"]
        PPTDUR = 0.0, [description = "Duration of precipitation", unit = u"s"]
        BDYBAR, [description = "7-day running average boundary condition (dimensionless)", unit = u"1"]
        PM1000, [description = "Previous MC1000 value (dimensionless)", unit = u"1"]
    end

    @variables begin
        MC1000(t), [description = "1000-hour fuel moisture content (dimensionless)", unit = u"1"]
        BNDRYT(t), [description = "Weighted 24-hour average boundary condition (dimensionless)", unit = u"1"]
    end

    # Convert precipitation duration from seconds to hours for equation (non-dimensionalize first)
    PPTDUR_hr = (PPTDUR / s_ref) / hr_to_s

    eqs = [
        # Weighted 24-hour average boundary condition for 1000-hr (page 5)
        BNDRYT ~ ((24.0 - PPTDUR_hr) * EMCBAR + PPTDUR_hr * (2.7 * PPTDUR_hr + 76.0) / 100.0) / 24.0,

        # 1000-hour fuel moisture daily update equation (page 5)
        MC1000 ~ PM1000 + (BDYBAR - PM1000) * response_coef
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Herbaceous Fuel Moisture
# =============================================================================

"""
    HerbaceousFuelMoisture(; name=:HerbaceousFuelMoisture)

Calculate herbaceous fuel moisture content (MCHERB).

The model tracks five stages: pregreen, greenup, green, transition, and cured/frozen.
MCHERB is a function of X1000 (modified 1000-hr moisture), herbaceous type (annual/perennial),
and NFDRS climate class.

From Cohen & Deeming (1985), pages 6-8:
- Greenup: MCHERB = MCHRBI + (MCHRBP - MCHRBI) × GREN
- Green: MCHERB = HERBGA + HERBGB × X1000
- Transition (annuals): MCHERB = ANNTA + ANNTB × X1000
- Transition (perennials): MCHERB = PERTA + PERTB × X1000

Climate-dependent parameters (page 7):
| Climate | HERBGA | HERBGB | ANNTA  | ANNTB | PERTA | PERTB |
|---------|--------|--------|--------|-------|-------|-------|
| 1       | -70.0  | 12.8   | -150.5 | 18.4  | 11.2  | 7.4   |
| 2       | -100.0 | 14.0   | -187.7 | 19.6  | -10.3 | 8.3   |
| 3       | -137.5 | 15.5   | -245.2 | 22.0  | -42.7 | 9.8   |
| 4       | -185.0 | 17.4   | -305.2 | 24.3  | -93.5 | 12.2  |

# Parameters
- `X1000`: Modified 1000-hr moisture for herbaceous model (dimensionless)
- `MC1`: 1-hour fuel moisture content (dimensionless)
- `CLIMAT`: NFDRS climate class 1-4 (dimensionless)
- `is_annual`: Flag: 1 if annual, 0 if perennial (dimensionless)
- `GRNDAY`: Days since greenup started (dimensionless)
- `is_greenup`: Flag: 1 if in greenup, 0 otherwise (dimensionless)
- `is_cured`: Flag: 1 if cured/frozen, 0 otherwise (dimensionless)

# Variables
- `MCHERB`: Herbaceous fuel moisture content (dimensionless)
"""
@component function HerbaceousFuelMoisture(; name=:HerbaceousFuelMoisture)
    @parameters begin
        X1000, [description = "Modified 1000-hr moisture (dimensionless)", unit = u"1"]
        MC1, [description = "1-hour fuel moisture content (dimensionless)", unit = u"1"]
        CLIMAT = 2.0, [description = "NFDRS climate class 1-4 (dimensionless)", unit = u"1"]
        is_annual = 0.0, [description = "Flag: 1 if annual, 0 if perennial (dimensionless)", unit = u"1"]
        GRNDAY = 0.0, [description = "Days since greenup started (dimensionless)", unit = u"1"]
        is_greenup = 0.0, [description = "Flag: 1 if in greenup, 0 otherwise (dimensionless)", unit = u"1"]
        is_cured = 0.0, [description = "Flag: 1 if cured/frozen, 0 otherwise (dimensionless)", unit = u"1"]
    end

    @variables begin
        MCHERB(t), [description = "Herbaceous fuel moisture content (dimensionless)", unit = u"1"]
    end

    # Climate class dependent parameters (page 7) - X1000 as percent for these equations
    X1000_pct = X1000 * 100

    # Greenup parameters (page 7)
    HERBGA = ifelse(CLIMAT < 1.5, -70.0,
             ifelse(CLIMAT < 2.5, -100.0,
             ifelse(CLIMAT < 3.5, -137.5, -185.0)))
    HERBGB = ifelse(CLIMAT < 1.5, 12.8,
             ifelse(CLIMAT < 2.5, 14.0,
             ifelse(CLIMAT < 3.5, 15.5, 17.4)))

    # Transition parameters for annuals (page 8)
    ANNTA = ifelse(CLIMAT < 1.5, -150.5,
            ifelse(CLIMAT < 2.5, -187.7,
            ifelse(CLIMAT < 3.5, -245.2, -305.2)))
    ANNTB = ifelse(CLIMAT < 1.5, 18.4,
            ifelse(CLIMAT < 2.5, 19.6,
            ifelse(CLIMAT < 3.5, 22.0, 24.3)))

    # Transition parameters for perennials (page 8)
    PERTA = ifelse(CLIMAT < 1.5, 11.2,
            ifelse(CLIMAT < 2.5, -10.3,
            ifelse(CLIMAT < 3.5, -42.7, -93.5)))
    PERTB = ifelse(CLIMAT < 1.5, 7.4,
            ifelse(CLIMAT < 2.5, 8.3,
            ifelse(CLIMAT < 3.5, 9.8, 12.2)))

    # Greenup fraction (Eq. 4, page 7)
    GREN = min(1.0, GRNDAY / (7.0 * CLIMAT))

    # Potential herbaceous moisture (as fraction)
    MCHRBP = (HERBGA + HERBGB * X1000_pct) / 100.0

    eqs = [
        # MCHERB calculation based on stage (pages 6-8)
        MCHERB ~ ifelse(is_cured > 0.5,
            # Cured/frozen stage (page 8)
            ifelse(is_annual > 0.5, MC1, max(0.30, min(1.50, (PERTA + PERTB * X1000_pct) / 100.0))),
            ifelse(is_greenup > 0.5,
                # Greenup stage (page 7): phase from 0.30 to MCHRBP
                0.30 + (max(0.30, MCHRBP) - 0.30) * GREN,
                # Green or transition stage
                ifelse(MCHRBP > 1.20,
                    # Green stage (page 7): MCHERB > 120%
                    min(2.50, MCHRBP),
                    # Transition stage (page 8): 30% < MCHERB < 120%
                    ifelse(is_annual > 0.5,
                        max(0.30, min(1.50, (ANNTA + ANNTB * X1000_pct) / 100.0)),
                        max(0.30, min(1.50, (PERTA + PERTB * X1000_pct) / 100.0))
                    )
                )
            )
        )
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Woody (Shrub) Fuel Moisture
# =============================================================================

"""
    WoodyFuelMoisture(; name=:WoodyFuelMoisture)

Calculate woody shrub fuel moisture content (MCWOOD).

The model has four stages: pregreen, greenup, green, and frozen.
MCWOOD is a function of MC1000 and NFDRS climate class.

From Cohen & Deeming (1985), pages 8-9:
- Pregreen: MCWOOD = PREGRN
- Greenup: MCWOOD = MCWODI + (MCWODP - MCWODI) × GREN
- Green: MCWOOD = WOODGA + WOODGB × MC1000

Climate-dependent parameters (page 8):
| Climate | PREGRN | WOODGA | WOODGB |
|---------|--------|--------|--------|
| 1       | 50%    | 12.5   | 7.5    |
| 2       | 60%    | -5.0   | 8.2    |
| 3       | 70%    | -22.5  | 8.9    |
| 4       | 80%    | -45.0  | 9.8    |

# Parameters
- `MC1000`: 1000-hour fuel moisture content (dimensionless)
- `CLIMAT`: NFDRS climate class 1-4 (dimensionless)
- `GRNDAY`: Days since greenup started (dimensionless)
- `is_greenup`: Flag: 1 if in greenup, 0 otherwise (dimensionless)
- `is_frozen`: Flag: 1 if frozen/dormant, 0 otherwise (dimensionless)

# Variables
- `MCWOOD`: Woody fuel moisture content (dimensionless)
"""
@component function WoodyFuelMoisture(; name=:WoodyFuelMoisture)
    @parameters begin
        MC1000, [description = "1000-hour fuel moisture content (dimensionless)", unit = u"1"]
        CLIMAT = 2.0, [description = "NFDRS climate class 1-4 (dimensionless)", unit = u"1"]
        GRNDAY = 0.0, [description = "Days since greenup started (dimensionless)", unit = u"1"]
        is_greenup = 0.0, [description = "Flag: 1 if in greenup, 0 otherwise (dimensionless)", unit = u"1"]
        is_frozen = 0.0, [description = "Flag: 1 if frozen/dormant, 0 otherwise (dimensionless)", unit = u"1"]
    end

    @variables begin
        MCWOOD(t), [description = "Woody fuel moisture content (dimensionless)", unit = u"1"]
    end

    # MC1000 as percent for equations
    MC1000_pct = MC1000 * 100

    # PREGRN values by climate class (page 8)
    PREGRN = ifelse(CLIMAT < 1.5, 0.50,
             ifelse(CLIMAT < 2.5, 0.60,
             ifelse(CLIMAT < 3.5, 0.70, 0.80)))

    # WOODGA and WOODGB by climate class (page 8)
    WOODGA = ifelse(CLIMAT < 1.5, 12.5,
             ifelse(CLIMAT < 2.5, -5.0,
             ifelse(CLIMAT < 3.5, -22.5, -45.0)))
    WOODGB = ifelse(CLIMAT < 1.5, 7.5,
             ifelse(CLIMAT < 2.5, 8.2,
             ifelse(CLIMAT < 3.5, 8.9, 9.8)))

    # Greenup fraction (page 9)
    GREN = min(1.0, GRNDAY / (7.0 * CLIMAT))

    # Potential woody moisture (as fraction)
    MCWODP = (WOODGA + WOODGB * MC1000_pct) / 100.0

    eqs = [
        # MCWOOD calculation (pages 8-9)
        MCWOOD ~ ifelse(is_frozen > 0.5,
            # Frozen/dormant stage (page 9)
            PREGRN,
            ifelse(is_greenup > 0.5,
                # Greenup stage (page 9)
                PREGRN + (max(PREGRN, MCWODP) - PREGRN) * GREN,
                # Green stage (page 9)
                max(PREGRN, min(2.0, MCWODP))
            )
        )
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Fuel Loading Transfer (Herbaceous to 1-hour)
# =============================================================================

"""
    FuelLoadingTransfer(; name=:FuelLoadingTransfer)

Calculate the transfer of herbaceous fuel loading to 1-hour dead fuel class
based on herbaceous moisture content (curing).

Implements Equations 5, 6, 7, 8 from Cohen & Deeming (1985), page 7:
- Eq. 5: FCTCUR = 1.33 - 0.0111 × MCHERB (where MCHERB in percent)
- Eq. 6: WHERBC = FCTCUR × WHERB
- Eq. 7: W1P = W1 + WHERBC
- Eq. 8: WHERBP = WHERB - WHERBC

# Parameters
- `MCHERB`: Herbaceous fuel moisture content (dimensionless)
- `W1`: Base 1-hour fuel loading (kg/m²)
- `WHERB`: Herbaceous fuel loading (kg/m²)

# Variables
- `FCTCUR`: Fraction transferred to 1-hour class (dimensionless)
- `WHERBC`: Herbaceous loading transferred (kg/m²)
- `W1P`: Effective 1-hour fuel loading (kg/m²)
- `WHERBP`: Remaining herbaceous loading (kg/m²)
"""
@component function FuelLoadingTransfer(; name=:FuelLoadingTransfer)
    @parameters begin
        MCHERB, [description = "Herbaceous fuel moisture content (dimensionless)", unit = u"1"]
        W1, [description = "Base 1-hour fuel loading", unit = u"kg/m^2"]
        WHERB, [description = "Herbaceous fuel loading", unit = u"kg/m^2"]
    end

    @variables begin
        FCTCUR(t), [description = "Fraction transferred to 1-hour class (dimensionless)", unit = u"1"]
        WHERBC(t), [description = "Herbaceous loading transferred", unit = u"kg/m^2"]
        W1P(t), [description = "Effective 1-hour fuel loading", unit = u"kg/m^2"]
        WHERBP(t), [description = "Remaining herbaceous loading", unit = u"kg/m^2"]
    end

    # MCHERB as percent for equation
    MCHERB_pct = MCHERB * 100

    eqs = [
        # Eq. 5 (page 7): Fraction transferred
        FCTCUR ~ max(0.0, min(1.0, 1.33 - 0.0111 * MCHERB_pct)),

        # Eq. 6 (page 7): Amount transferred
        WHERBC ~ FCTCUR * WHERB,

        # Eq. 7 (page 7): Effective 1-hour loading
        W1P ~ W1 + WHERBC,

        # Eq. 8 (page 7): Remaining herbaceous loading
        WHERBP ~ WHERB - WHERBC
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Spread Component (SC)
# =============================================================================

"""
    SpreadComponent(; name=:SpreadComponent)

Calculate the NFDRS Spread Component (SC) based on Rothermel's fire spread model.

The SC is the forward rate of spread of the flaming front in m/s.
Implements equations from pages 9-11 of Cohen & Deeming (1985).

Based on Rothermel (1972) as modified by Albini (1976).

# Parameters
Fuel loadings (kg/m²):
- `W1P`: 1-hour fuel loading (including transferred herbaceous)
- `W10`: 10-hour fuel loading
- `W100`: 100-hour fuel loading
- `WHERBP`: Remaining herbaceous loading
- `WWOOD`: Woody fuel loading

Surface-area-to-volume ratios (m⁻¹):
- `SG1`, `SG10`, `SG100`: Dead fuel SAV ratios
- `SGHERB`, `SGWOOD`: Live fuel SAV ratios

Other parameters:
- `MC1`, `MC10`, `MC100`: Dead fuel moisture contents (dimensionless)
- `MCHERB`, `MCWOOD`: Live fuel moisture contents (dimensionless)
- `MXD`: Dead fuel moisture of extinction (percent)
- `HD`, `HL`: Heat of combustion (J/kg) for dead and live fuels
- `DEPTH`: Fuel bed depth (m)
- `WS`: 20-ft windspeed (m/s)
- `WNDFC`: Wind reduction factor (dimensionless)
- `slope_class`: NFDRS slope class 1-5 (dimensionless)
- `fuels_wet`: Flag: 1 if fuels wet or snow-covered (dimensionless)

# Variables
- `SC`: Spread component (m/s)
- `ROS`: Rate of spread (m/s)
- `IR`: Reaction intensity (W/m²)
"""
@component function SpreadComponent(; name=:SpreadComponent)
    @constants begin
        # Unit reference constants for non-dimensionalization
        W_ref = 1.0, [description = "Reference fuel loading for non-dimensionalization", unit = u"kg/m^2"]
        SG_ref = 1.0, [description = "Reference SAV ratio for non-dimensionalization", unit = u"m^-1"]
        H_ref = 1.0, [description = "Reference heat content for non-dimensionalization", unit = u"J/kg"]
        D_ref = 1.0, [description = "Reference depth for non-dimensionalization", unit = u"m"]
        V_ref = 1.0, [description = "Reference velocity for non-dimensionalization", unit = u"m/s"]
        V_zero = 0.0, [description = "Zero velocity for conditional expressions", unit = u"m/s"]
        IR_ref = 1.0, [description = "Reference reaction intensity for output", unit = u"W/m^2"]

        # Unit conversion factors (dimensionless numerical factors)
        ft_to_m = 0.3048, [description = "Feet to meters conversion (dimensionless)", unit = u"1"]
        m_to_ft = 3.28084, [description = "Meters to feet conversion (dimensionless)", unit = u"1"]
        ms_to_fpm = 196.8504, [description = "m/s to ft/min conversion (dimensionless)", unit = u"1"]
        fpm_to_ms = 0.00508, [description = "ft/min to m/s conversion (dimensionless)", unit = u"1"]
        kg_m2_to_lb_ft2 = 0.204816, [description = "kg/m² to lb/ft² conversion (dimensionless)", unit = u"1"]
        J_kg_to_Btu_lb = 0.000429923, [description = "J/kg to Btu/lb conversion (dimensionless)", unit = u"1"]

        # Particle densities (imperial for Rothermel equations: 32 lb/ft³)
        RHOD_imp = 32.0, [description = "Dead fuel particle density, lb/ft³ (dimensionless)", unit = u"1"]
        RHOL_imp = 32.0, [description = "Live fuel particle density, lb/ft³ (dimensionless)", unit = u"1"]

        # Mineral contents (page 9, dimensionless)
        STD = 0.0555, [description = "Dead fuel total mineral content (dimensionless)", unit = u"1"]
        STL = 0.0555, [description = "Live fuel total mineral content (dimensionless)", unit = u"1"]
        SD = 0.01, [description = "Dead fuel silica-free mineral content (dimensionless)", unit = u"1"]
        SL = 0.01, [description = "Live fuel silica-free mineral content (dimensionless)", unit = u"1"]

        # Reaction intensity conversion (Btu/ft²/min to W/m²)
        IR_SI_factor = 189.27, [description = "IR conversion factor (dimensionless)", unit = u"1"]
    end

    @parameters begin
        # Fuel loadings (kg/m²)
        W1P, [description = "1-hour fuel loading including transferred herb", unit = u"kg/m^2"]
        W10, [description = "10-hour fuel loading", unit = u"kg/m^2"]
        W100, [description = "100-hour fuel loading", unit = u"kg/m^2"]
        WHERBP, [description = "Remaining herbaceous loading", unit = u"kg/m^2"]
        WWOOD, [description = "Woody fuel loading", unit = u"kg/m^2"]

        # SAV ratios (m⁻¹)
        SG1, [description = "1-hour fuel SAV ratio", unit = u"m^-1"]
        SG10, [description = "10-hour fuel SAV ratio", unit = u"m^-1"]
        SG100, [description = "100-hour fuel SAV ratio", unit = u"m^-1"]
        SGHERB, [description = "Herbaceous fuel SAV ratio", unit = u"m^-1"]
        SGWOOD, [description = "Woody fuel SAV ratio", unit = u"m^-1"]

        # Moisture contents (dimensionless fractions, 0-1)
        MC1, [description = "1-hour fuel moisture (dimensionless)", unit = u"1"]
        MC10, [description = "10-hour fuel moisture (dimensionless)", unit = u"1"]
        MC100, [description = "100-hour fuel moisture (dimensionless)", unit = u"1"]
        MCHERB, [description = "Herbaceous fuel moisture (dimensionless)", unit = u"1"]
        MCWOOD, [description = "Woody fuel moisture (dimensionless)", unit = u"1"]

        # Other fuel model parameters
        MXD, [description = "Dead fuel moisture of extinction (percent, dimensionless)", unit = u"1"]
        HD = 18608000.0, [description = "Dead fuel heat of combustion", unit = u"J/kg"]
        HL = 18608000.0, [description = "Live fuel heat of combustion", unit = u"J/kg"]
        DEPTH, [description = "Fuel bed depth", unit = u"m"]

        # Environmental parameters
        WS = 0.0, [description = "20-ft windspeed", unit = u"m/s"]
        WNDFC = 0.4, [description = "Wind reduction factor (dimensionless)", unit = u"1"]
        slope_class = 1.0, [description = "NFDRS slope class 1-5 (dimensionless)", unit = u"1"]
        fuels_wet = 0.0, [description = "Flag: 1 if fuels wet or snow-covered (dimensionless)", unit = u"1"]
    end

    @variables begin
        SC(t), [description = "Spread component", unit = u"m/s"]
        ROS(t), [description = "Rate of spread", unit = u"m/s"]
        IR(t), [description = "Reaction intensity", unit = u"W/m^2"]
    end

    # Non-dimensionalize inputs, then convert to imperial units for Rothermel equations
    # (the empirical coefficients are calibrated for imperial)
    W1P_imp = (W1P / W_ref) * kg_m2_to_lb_ft2
    W10_imp = (W10 / W_ref) * kg_m2_to_lb_ft2
    W100_imp = (W100 / W_ref) * kg_m2_to_lb_ft2
    WHERBP_imp = (WHERBP / W_ref) * kg_m2_to_lb_ft2
    WWOOD_imp = (WWOOD / W_ref) * kg_m2_to_lb_ft2

    SG1_imp = (SG1 / SG_ref) * ft_to_m  # m⁻¹ to ft⁻¹
    SG10_imp = (SG10 / SG_ref) * ft_to_m
    SG100_imp = (SG100 / SG_ref) * ft_to_m
    SGHERB_imp = (SGHERB / SG_ref) * ft_to_m
    SGWOOD_imp = (SGWOOD / SG_ref) * ft_to_m

    DEPTH_imp = (DEPTH / D_ref) * m_to_ft

    HD_imp = (HD / H_ref) * J_kg_to_Btu_lb
    HL_imp = (HL / H_ref) * J_kg_to_Btu_lb

    WS_fpm = (WS / V_ref) * ms_to_fpm  # Wind in ft/min

    # Net loadings (imperial, page 9)
    W1N = W1P_imp * (1.0 - STD)
    W10N = W10_imp * (1.0 - STD)
    W100N = W100_imp * (1.0 - STD)
    WHERBN = WHERBP_imp * (1.0 - STL)
    WWOODN = WWOOD_imp * (1.0 - STL)

    # Total loadings
    WTOTD = W1P_imp + W10_imp + W100_imp
    WTOTL = WHERBP_imp + WWOOD_imp
    WTOT = WTOTD + WTOTL

    # Surface areas
    SA1 = (W1P_imp / RHOD_imp) * SG1_imp
    SA10 = (W10_imp / RHOD_imp) * SG10_imp
    SA100 = (W100_imp / RHOD_imp) * SG100_imp
    SAHERB = (WHERBP_imp / RHOL_imp) * SGHERB_imp
    SAWOOD = (WWOOD_imp / RHOL_imp) * SGWOOD_imp

    SADEAD = SA1 + SA10 + SA100
    SALIVE = SAHERB + SAWOOD

    # Weighting factors (avoid division by zero)
    SADEAD_safe = max(1e-10, SADEAD)
    SALIVE_safe = max(1e-10, SALIVE)
    SATOT_safe = max(1e-10, SADEAD + SALIVE)

    F1 = SA1 / SADEAD_safe
    F10 = SA10 / SADEAD_safe
    F100 = SA100 / SADEAD_safe
    FHERB = SAHERB / SALIVE_safe
    FWOOD = SAWOOD / SALIVE_safe

    FDEAD = SADEAD / SATOT_safe
    FLIVE = SALIVE / SATOT_safe

    # Weighted net loadings
    WDEADN = F1 * W1N + F10 * W10N + F100 * W100N
    WLIVEN = FWOOD * WWOODN + FHERB * WHERBN

    # Characteristic SAV ratios (imperial)
    SGBRD = F1 * SG1_imp + F10 * SG10_imp + F100 * SG100_imp
    SGBRL = FHERB * SGHERB_imp + FWOOD * SGWOOD_imp
    SGBRT = FDEAD * SGBRD + FLIVE * SGBRL
    SGBRT_safe = max(1.0, SGBRT)

    # Bulk density (imperial)
    RHOBED = (WTOT - 0.0) / max(0.01, DEPTH_imp)

    # Packing ratio
    RHOBAR = 32.0  # Constant particle density (lb/ft³)
    BETBAR = RHOBED / RHOBAR
    BETOP = 3.348 * SGBRT_safe^(-0.8189)

    # Maximum and optimum reaction velocity
    GMAMX = SGBRT_safe^1.5 / (495.0 + 0.0594 * SGBRT_safe^1.5)
    AD = 133.0 * SGBRT_safe^(-0.7913)
    ratio = BETBAR / max(1e-10, BETOP)
    GMAOP = GMAMX * ratio^AD * exp(AD * (1.0 - ratio))

    # No-wind propagating flux ratio
    ZETA = exp((0.792 + 0.681 * sqrt(SGBRT_safe)) * (BETBAR + 0.1)) / (192.0 + 0.2595 * SGBRT_safe)

    # Heating numbers for live fuel extinction moisture
    HN1 = W1N * exp(-138.0 / max(1.0, SG1_imp))
    HN10 = W10N * exp(-138.0 / max(1.0, SG10_imp))
    HN100 = W100N * exp(-138.0 / max(1.0, SG100_imp))
    HNHERB = WHERBN * exp(-500.0 / max(1.0, SGHERB_imp))
    HNWOOD = WWOODN * exp(-500.0 / max(1.0, SGWOOD_imp))

    HNDEAD = HN1 + HN10 + HN100
    HNLIVE = HNHERB + HNWOOD

    # Weighted dead fuel moisture for MXL (convert to percent for MXD comparison)
    MCLFE_pct = ifelse(HNDEAD > 1e-10,
        (MC1 * 100 * HN1 + MC10 * 100 * HN10 + MC100 * 100 * HN100) / HNDEAD,
        MC1 * 100)

    # Live fuel moisture of extinction (MXL in percent, matching MXD)
    WRAT = ifelse(HNLIVE > 1e-10, HNDEAD / HNLIVE, 0.0)
    MXL_calc = (2.9 * WRAT * (1.0 - MCLFE_pct / MXD) - 0.226) * 100.0
    MXL = max(MXD, MXL_calc)

    # Weighted moisture contents (as fractions)
    WTMCD = F1 * MC1 + F10 * MC10 + F100 * MC100
    WTMCL = FHERB * MCHERB + FWOOD * MCWOOD

    # Moisture damping coefficients (convert WTMCD/WTMCL to percent for ratio)
    DEDRT = (WTMCD * 100) / MXD
    LIVRT = (WTMCL * 100) / max(0.01, MXL)

    ETAMD = max(0.0, min(1.0, 1.0 - 2.59 * DEDRT + 5.11 * DEDRT^2 - 3.52 * DEDRT^3))
    ETAML = max(0.0, min(1.0, 1.0 - 2.59 * LIVRT + 5.11 * LIVRT^2 - 3.52 * LIVRT^3))

    # Mineral damping coefficients
    ETASD = 0.174 * SD^(-0.19)
    ETASL = 0.174 * SL^(-0.19)

    # Wind effect coefficients
    B = 0.02526 * SGBRT_safe^0.54
    C = 7.47 * exp(-0.133 * SGBRT_safe^0.55)
    E = 0.715 * exp(-3.59e-4 * SGBRT_safe)
    UFACT = C * ratio^(-E)

    # Slope effect (SLPFCT by slope class)
    SLPFCT = ifelse(slope_class < 1.5, 0.267,
             ifelse(slope_class < 2.5, 0.533,
             ifelse(slope_class < 3.5, 1.068,
             ifelse(slope_class < 4.5, 2.134, 4.273))))

    # Calculate reaction intensity for wind limit check (in Btu/ft²/min)
    IR_calc_imp = GMAOP * (WDEADN * HD_imp * ETASD * ETAMD + WLIVEN * HL_imp * ETASL * ETAML)

    # Wind speed effect with limit: if wind term > 0.9 * IR, use 0.9 * IR
    # Original uses WS * 88 * WNDFC where WS is in mph
    # We have WS_fpm = WS * 196.85 (from m/s), and mph * 88 = ft/min
    # So we need WS_fpm * WNDFC
    wind_term = WS_fpm * WNDFC
    wind_term_limited = min(wind_term, 0.9 * IR_calc_imp)

    # Wind effect multiplier with limit applied
    PHIWND = UFACT * wind_term_limited^B

    # Slope effect multiplier
    PHISLP = SLPFCT * max(0.01, BETBAR)^(-0.3)

    # Heat sink calculation (Eq. from page 11) - in imperial units
    HTSINK = RHOBED * (
        FDEAD * (F1 * exp(-138.0/max(1.0,SG1_imp)) * (250.0 + 11.16 * MC1 * 100) +
                 F10 * exp(-138.0/max(1.0,SG10_imp)) * (250.0 + 11.16 * MC10 * 100) +
                 F100 * exp(-138.0/max(1.0,SG100_imp)) * (250.0 + 11.16 * MC100 * 100)) +
        FLIVE * (FHERB * exp(-138.0/max(1.0,SGHERB_imp)) * (250.0 + 11.16 * MCHERB * 100) +
                 FWOOD * exp(-138.0/max(1.0,SGWOOD_imp)) * (250.0 + 11.16 * MCWOOD * 100))
    )

    # Rate of spread in ft/min (imperial, page 11)
    ROS_imp = IR_calc_imp * ZETA * (1.0 + PHISLP + PHIWND) / max(0.01, HTSINK)

    eqs = [
        # Reaction intensity (convert to W/m² and re-dimensionalize)
        IR ~ IR_calc_imp * IR_SI_factor * IR_ref,

        # Rate of spread (convert ft/min to m/s and re-dimensionalize)
        ROS ~ ifelse(fuels_wet > 0.5, V_zero, ROS_imp * fpm_to_ms * V_ref),

        # Spread component (same as ROS in SI)
        SC ~ ROS
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Energy Release Component (ERC)
# =============================================================================

"""
    EnergyReleaseComponent(; name=:EnergyReleaseComponent)

Calculate the NFDRS Energy Release Component (ERC).

ERC is based on loading-weighted (not surface-area-weighted) calculations.
It represents the available energy release per unit area.

Implements equations from pages 11-12 of Cohen & Deeming (1985).

The 0.04 scaling factor has units ft²/Btu, so a unit value of ERC is equivalent
to 25 Btu of available energy per square foot (page 12).

# Parameters
Same as SpreadComponent, plus:
- `W1000`: 1000-hour fuel loading (kg/m²)
- `SG1000`: 1000-hour fuel SAV ratio (m⁻¹)
- `MC1000`: 1000-hour fuel moisture (dimensionless)

# Variables
- `ERC`: Energy release component (dimensionless)
- `IRE`: Loading-weighted reaction intensity (W/m²)
"""
@component function EnergyReleaseComponent(; name=:EnergyReleaseComponent)
    @constants begin
        # Unit reference constants for non-dimensionalization
        W_ref = 1.0, [description = "Reference fuel loading for non-dimensionalization", unit = u"kg/m^2"]
        SG_ref = 1.0, [description = "Reference SAV ratio for non-dimensionalization", unit = u"m^-1"]
        H_ref = 1.0, [description = "Reference heat content for non-dimensionalization", unit = u"J/kg"]
        D_ref = 1.0, [description = "Reference depth for non-dimensionalization", unit = u"m"]
        IR_ref = 1.0, [description = "Reference reaction intensity for output", unit = u"W/m^2"]

        # Unit conversion factors (dimensionless numerical factors)
        ft_to_m = 0.3048, [description = "Feet to meters conversion (dimensionless)", unit = u"1"]
        m_to_ft = 3.28084, [description = "Meters to feet conversion (dimensionless)", unit = u"1"]
        kg_m2_to_lb_ft2 = 0.204816, [description = "kg/m² to lb/ft² conversion (dimensionless)", unit = u"1"]
        J_kg_to_Btu_lb = 0.000429923, [description = "J/kg to Btu/lb conversion (dimensionless)", unit = u"1"]

        # Particle densities (imperial for Rothermel equations)
        RHOD_imp = 32.0, [description = "Dead fuel particle density, lb/ft³ (dimensionless)", unit = u"1"]
        RHOL_imp = 32.0, [description = "Live fuel particle density, lb/ft³ (dimensionless)", unit = u"1"]

        # Mineral contents (page 9, dimensionless)
        STD = 0.0555, [description = "Dead fuel total mineral content (dimensionless)", unit = u"1"]
        STL = 0.0555, [description = "Live fuel total mineral content (dimensionless)", unit = u"1"]
        SD = 0.01, [description = "Dead fuel silica-free mineral content (dimensionless)", unit = u"1"]
        SL = 0.01, [description = "Live fuel silica-free mineral content (dimensionless)", unit = u"1"]

        # Reaction intensity conversion
        IR_SI_factor = 189.27, [description = "IR conversion factor (dimensionless)", unit = u"1"]
    end

    @parameters begin
        # Fuel loadings (kg/m²)
        W1P, [description = "1-hour fuel loading including transferred herb", unit = u"kg/m^2"]
        W10, [description = "10-hour fuel loading", unit = u"kg/m^2"]
        W100, [description = "100-hour fuel loading", unit = u"kg/m^2"]
        W1000 = 0.0, [description = "1000-hour fuel loading", unit = u"kg/m^2"]
        WHERBP, [description = "Remaining herbaceous loading", unit = u"kg/m^2"]
        WWOOD, [description = "Woody fuel loading", unit = u"kg/m^2"]

        # SAV ratios (m⁻¹)
        SG1, [description = "1-hour fuel SAV ratio", unit = u"m^-1"]
        SG10, [description = "10-hour fuel SAV ratio", unit = u"m^-1"]
        SG100, [description = "100-hour fuel SAV ratio", unit = u"m^-1"]
        SG1000 = 26.25, [description = "1000-hour fuel SAV ratio", unit = u"m^-1"]
        SGHERB, [description = "Herbaceous fuel SAV ratio", unit = u"m^-1"]
        SGWOOD, [description = "Woody fuel SAV ratio", unit = u"m^-1"]

        # Moisture contents (dimensionless fractions, 0-1)
        MC1, [description = "1-hour fuel moisture (dimensionless)", unit = u"1"]
        MC10, [description = "10-hour fuel moisture (dimensionless)", unit = u"1"]
        MC100, [description = "100-hour fuel moisture (dimensionless)", unit = u"1"]
        MC1000 = 0.15, [description = "1000-hour fuel moisture (dimensionless)", unit = u"1"]
        MCHERB, [description = "Herbaceous fuel moisture (dimensionless)", unit = u"1"]
        MCWOOD, [description = "Woody fuel moisture (dimensionless)", unit = u"1"]

        # Other fuel model parameters
        MXD, [description = "Dead fuel moisture of extinction (percent, dimensionless)", unit = u"1"]
        MXL, [description = "Live fuel moisture of extinction (percent, dimensionless)", unit = u"1"]
        HD = 18608000.0, [description = "Dead fuel heat of combustion", unit = u"J/kg"]
        HL = 18608000.0, [description = "Live fuel heat of combustion", unit = u"J/kg"]
        DEPTH, [description = "Fuel bed depth", unit = u"m"]
    end

    @variables begin
        ERC(t), [description = "Energy release component (dimensionless)", unit = u"1"]
        IRE(t), [description = "Loading-weighted reaction intensity", unit = u"W/m^2"]
    end

    # Non-dimensionalize inputs, then convert to imperial units for Rothermel equations
    W1P_imp = (W1P / W_ref) * kg_m2_to_lb_ft2
    W10_imp = (W10 / W_ref) * kg_m2_to_lb_ft2
    W100_imp = (W100 / W_ref) * kg_m2_to_lb_ft2
    W1000_imp = (W1000 / W_ref) * kg_m2_to_lb_ft2
    WHERBP_imp = (WHERBP / W_ref) * kg_m2_to_lb_ft2
    WWOOD_imp = (WWOOD / W_ref) * kg_m2_to_lb_ft2

    SG1_imp = (SG1 / SG_ref) * ft_to_m
    SG10_imp = (SG10 / SG_ref) * ft_to_m
    SG100_imp = (SG100 / SG_ref) * ft_to_m
    SG1000_imp = (SG1000 / SG_ref) * ft_to_m
    SGHERB_imp = (SGHERB / SG_ref) * ft_to_m
    SGWOOD_imp = (SGWOOD / SG_ref) * ft_to_m

    DEPTH_imp = (DEPTH / D_ref) * m_to_ft

    HD_imp = (HD / H_ref) * J_kg_to_Btu_lb
    HL_imp = (HL / H_ref) * J_kg_to_Btu_lb

    # Total loadings (imperial)
    WTOTD = W1P_imp + W10_imp + W100_imp + W1000_imp
    WTOTL = WHERBP_imp + WWOOD_imp
    WTOT = WTOTD + WTOTL

    # Loading-based weighting factors
    WTOTD_safe = max(1e-10, WTOTD)
    WTOTL_safe = max(1e-10, WTOTL)
    WTOT_safe = max(1e-10, WTOT)

    F1E = W1P_imp / WTOTD_safe
    F10E = W10_imp / WTOTD_safe
    F100E = W100_imp / WTOTD_safe
    F1000E = W1000_imp / WTOTD_safe
    FHERBE = WHERBP_imp / WTOTL_safe
    FWOODE = WWOOD_imp / WTOTL_safe

    FDEADE = WTOTD / WTOT_safe
    FLIVEE = WTOTL / WTOT_safe

    # Net loadings
    WDEDNE = WTOTD * (1.0 - STD)
    WLIVNE = WTOTL * (1.0 - STL)

    # Characteristic SAV ratio (loading weighted, imperial)
    SGBRDE = F1E * SG1_imp + F10E * SG10_imp + F100E * SG100_imp + F1000E * SG1000_imp
    SGBRLE = FWOODE * SGWOOD_imp + FHERBE * SGHERB_imp
    SGBRTE = FDEADE * SGBRDE + FLIVEE * SGBRLE
    SGBRTE_safe = max(1.0, SGBRTE)

    # For residence time calculation, use surface-area weighted SGBRT
    SA1 = (W1P_imp / RHOD_imp) * SG1_imp
    SA10 = (W10_imp / RHOD_imp) * SG10_imp
    SA100 = (W100_imp / RHOD_imp) * SG100_imp
    SAHERB = (WHERBP_imp / RHOL_imp) * SGHERB_imp
    SAWOOD = (WWOOD_imp / RHOL_imp) * SGWOOD_imp
    SADEAD = SA1 + SA10 + SA100
    SALIVE = SAHERB + SAWOOD
    SADEAD_safe = max(1e-10, SADEAD)
    SALIVE_safe = max(1e-10, SALIVE)
    SATOT_safe = max(1e-10, SADEAD + SALIVE)
    F1 = SA1 / SADEAD_safe
    F10 = SA10 / SADEAD_safe
    F100 = SA100 / SADEAD_safe
    FHERB = SAHERB / SALIVE_safe
    FWOOD = SAWOOD / SALIVE_safe
    FDEAD = SADEAD / SATOT_safe
    FLIVE = SALIVE / SATOT_safe
    SGBRD = F1 * SG1_imp + F10 * SG10_imp + F100 * SG100_imp
    SGBRL = FHERB * SGHERB_imp + FWOOD * SGWOOD_imp
    SGBRT = FDEAD * SGBRD + FLIVE * SGBRL
    SGBRT_safe = max(1.0, SGBRT)

    # Bulk density and packing ratio (imperial)
    RHOBED = WTOT / max(0.01, DEPTH_imp)
    RHOBAR = 32.0
    BETBAR = RHOBED / RHOBAR
    BETOPE = 3.348 * SGBRTE_safe^(-0.8189)

    # Reaction velocity (loading weighted)
    GMAMXE = SGBRTE_safe^1.5 / (495.0 + 0.0594 * SGBRTE_safe^1.5)
    ADE = 133.0 * SGBRTE_safe^(-0.7913)
    ratio = BETBAR / max(1e-10, BETOPE)
    GMAOPE = GMAMXE * ratio^ADE * exp(ADE * (1.0 - ratio))

    # Weighted moisture contents (loading weighted, as fractions)
    WTMCDE = F1E * MC1 + F10E * MC10 + F100E * MC100 + F1000E * MC1000
    WTMCLE = FWOODE * MCWOOD + FHERBE * MCHERB

    # Moisture damping coefficients
    DEDRTE = (WTMCDE * 100) / MXD
    LIVRTE = (WTMCLE * 100) / max(0.01, MXL)

    ETAMDE = max(0.0, min(1.0, 1.0 - 2.0 * DEDRTE + 1.5 * DEDRTE^2 - 0.5 * DEDRTE^3))
    ETAMLE = max(0.0, min(1.0, 1.0 - 2.0 * LIVRTE + 1.5 * LIVRTE^2 - 0.5 * LIVRTE^3))

    # Mineral damping coefficients
    ETASD = 0.174 * SD^(-0.19)
    ETASL = 0.174 * SL^(-0.19)

    # Residence time (uses surface-area weighted SGBRT, in minutes)
    TAU = 384.0 / SGBRT_safe

    # Loading-weighted reaction intensity (Btu/ft²/min, page 12)
    IRE_imp = GMAOPE * (FDEADE * WDEDNE * HD_imp * ETASD * ETAMDE +
                        FLIVEE * WLIVNE * HL_imp * ETASL * ETAMLE)

    eqs = [
        # Loading-weighted reaction intensity (W/m² - re-dimensionalize)
        IRE ~ IRE_imp * IR_SI_factor * IR_ref,

        # Energy Release Component (dimensionless index)
        # 0.04 scaling factor applied to Btu/ft² (IRE_imp * TAU)
        ERC ~ 0.04 * IRE_imp * TAU
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Burning Index (BI)
# =============================================================================

"""
    BurningIndex(; name=:BurningIndex)

Calculate the NFDRS Burning Index (BI).

BI is numerically equivalent to 10 times the predicted flame length in feet,
based on Byram's flame length model (Byram 1959).

From Cohen & Deeming (1985), page 12:
- FL = 0.45 × I^0.46 (Byram's equation, FL in ft, I in Btu/ft/s)
- I = (ROS/60) × IRE × TAU = SC × 25 × ERC (after substitution)
- FL = 0.301 × (SC × ERC)^0.46
- BI = 10 × FL = 3.01 × (SC × ERC)^0.46

# Parameters
- `SC`: Spread component (m/s)
- `ERC`: Energy release component (dimensionless)
- `fuels_wet`: Flag: 1 if fuels wet or snow-covered (dimensionless)

# Variables
- `BI`: Burning index (dimensionless)
"""
@component function BurningIndex(; name=:BurningIndex)
    @constants begin
        # Reference velocity for non-dimensionalization
        V_ref = 1.0, [description = "Reference velocity for non-dimensionalization", unit = u"m/s"]
        # Conversion factor: m/s to ft/min
        # 1 m/s = 196.8504 ft/min (numerically converts m/s value to ft/min value)
        ms_to_fpm = 196.8504, [description = "Numerical conversion: m/s to ft/min (dimensionless)", unit = u"1"]
    end

    @parameters begin
        SC, [description = "Spread component", unit = u"m/s"]
        ERC, [description = "Energy release component (dimensionless)", unit = u"1"]
        fuels_wet = 0.0, [description = "Flag: 1 if fuels wet or snow-covered (dimensionless)", unit = u"1"]
    end

    @variables begin
        BI(t), [description = "Burning index (dimensionless)", unit = u"1"]
    end

    # Non-dimensionalize and convert SC from m/s to ft/min for the BI calculation (page 12)
    # The original NFDRS equation uses SC in ft/min
    SC_fpm = (SC / V_ref) * ms_to_fpm

    eqs = [
        # BI = 3.01 × (SC × ERC)^0.46 (page 12)
        # The 3.01 coefficient is calibrated for SC in ft/min
        BI ~ ifelse(fuels_wet > 0.5, 0.0, 3.01 * (max(0.0, SC_fpm * ERC))^0.46)
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Ignition Component (IC)
# =============================================================================

"""
    IgnitionComponent(; name=:IgnitionComponent)

Calculate the NFDRS Ignition Component (IC).

From Cohen & Deeming (1985), pages 12-13:
- IC = 0.10 × P(I) × P(F/I)

P(I) is the probability of ignition, scaled to be 100 at MC1=1.5% and 0 at MC1=25%.
P(F/I) is the probability of a reportable fire given ignition.

The heat of ignition (QIGN) equation is from page 13:
QIGN = 144.5 - 0.266×T - 0.00058×T² - 0.01×T×MC1 + 18.54×(1-exp(-0.151×MC1)) + 6.4×MC1
where T is in °C and MC1 is in percent.

# Parameters
- `TMPPRM`: Temperature at fuel-atmosphere interface (K)
- `MC1`: 1-hour fuel moisture content (dimensionless)
- `SC`: Spread component (m/s)
- `SCM`: Spread component threshold for reportable fires (m/s)

# Variables
- `IC`: Ignition component (0-100, dimensionless)
"""
@component function IgnitionComponent(; name=:IgnitionComponent)
    @constants begin
        # Reference units for non-dimensionalization
        K_ref = 1.0, [description = "Reference temperature for non-dimensionalization", unit = u"K"]
        V_ref = 1.0, [description = "Reference velocity for non-dimensionalization", unit = u"m/s"]
        # Scaling factors for P(I) calculation (page 12-13)
        PNORM1 = 0.00232, [description = "P(I) scaling factor 1 (dimensionless)", unit = u"1"]
        PNORM2 = 0.99767, [description = "P(I) scaling factor 2 (dimensionless)", unit = u"1"]
        PNORM3 = 0.0000185, [description = "P(I) scaling factor 3 (dimensionless)", unit = u"1"]
        K_to_C = 273.15, [description = "Kelvin to Celsius offset (dimensionless)", unit = u"1"]
    end

    @parameters begin
        TMPPRM, [description = "Temperature at fuel-atmosphere interface", unit = u"K"]
        MC1, [description = "1-hour fuel moisture content (dimensionless)", unit = u"1"]
        SC, [description = "Spread component", unit = u"m/s"]
        SCM, [description = "Spread component threshold for reportable fires", unit = u"m/s"]
    end

    @variables begin
        QIGN(t), [description = "Heat of ignition (dimensionless)", unit = u"1"]
        CHI(t), [description = "Intermediate variable (dimensionless)", unit = u"1"]
        PI(t), [description = "Probability of ignition 0-100 (dimensionless)", unit = u"1"]
        SCN(t), [description = "Normalized spread component (dimensionless)", unit = u"1"]
        PFI(t), [description = "Probability of reportable fire given ignition (dimensionless)", unit = u"1"]
        IC(t), [description = "Ignition component 0-100 (dimensionless)", unit = u"1"]
    end

    # Non-dimensionalize and convert temperature from K to °C for the empirical equation
    TMPPRM_C = (TMPPRM / K_ref) - K_to_C

    # Non-dimensionalize velocities
    SC_nd = SC / V_ref
    SCM_nd = SCM / V_ref

    eqs = [
        # Heat of ignition (page 13) - MC1 as percent
        QIGN ~ 144.5 - 0.266 * TMPPRM_C - 0.00058 * TMPPRM_C^2 - 0.01 * TMPPRM_C * (MC1 * 100) + 18.54 * (1.0 - exp(-0.151 * (MC1 * 100))) + 6.4 * (MC1 * 100),

        # Intermediate calculation (page 13)
        CHI ~ max(0.0, (344.0 - QIGN) / 10.0),

        # Probability of ignition (page 13)
        PI ~ ifelse(CHI^3.6 * PNORM3 <= PNORM1,
            0.0,
            max(0.0, min(100.0, (CHI^3.6 * PNORM3 - PNORM1) * 100.0 / PNORM2))
        ),

        # Normalized spread component (page 13): SCN = 100 × SC/SCM
        SCN ~ 100.0 * SC_nd / max(1e-10, SCM_nd),

        # Probability of reportable fire given ignition (page 13): P(F/I) = SCN^0.5
        PFI ~ sqrt(max(0.0, SCN)),

        # Ignition component (page 13): IC = 0.10 × P(I) × P(F/I)
        IC ~ 0.10 * PI * PFI
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Human-Caused Fire Occurrence Index (MCOI)
# =============================================================================

"""
    HumanFireOccurrenceIndex(; name=:HumanFireOccurrenceIndex)

Calculate the Human-Caused Fire Occurrence Index (MCOI).

From Cohen & Deeming (1985), page 14:
MCOI = 0.01 × MRISK × IC

The man-caused fire occurrence index relates the ignition component to human-caused
fire risk, which is typically derived from historical fire occurrence data.

# Parameters
- `MRISK`: Human-caused risk 0-100 (dimensionless)
- `IC`: Ignition component 0-100 (dimensionless)

# Variables
- `MCOI`: Human-caused fire occurrence index 0-100 (dimensionless)
"""
@component function HumanFireOccurrenceIndex(; name=:HumanFireOccurrenceIndex)
    @parameters begin
        MRISK, [description = "Human-caused risk 0-100 (dimensionless)", unit = u"1"]
        IC, [description = "Ignition component 0-100 (dimensionless)", unit = u"1"]
    end

    @variables begin
        MCOI(t), [description = "Human-caused fire occurrence index 0-100 (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Page 14: MCOI = 0.01 × MRISK × IC
        MCOI ~ 0.01 * MRISK * IC
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Lightning-Caused Fire Occurrence Index (LOI)
# =============================================================================

"""
    LightningFireOccurrenceIndex(; name=:LightningFireOccurrenceIndex)

Calculate the Lightning-Caused Fire Occurrence Index (LOI).

The model assumes a thunderstorm traversing an area forms a corridor aligned with
the storm's track that receives both rain and lightning. Flanking the rain-lightning
corridor on both sides are areas subjected to lightning only.

Implements equations from pages 13-14 of Cohen & Deeming (1985).

# Parameters
- `LAL`: NFDRS Lightning Activity Level (1-6)
- `IC`: Ignition component (0-100)
- `MC1`: 1-hour fuel moisture content (dimensionless)
- `LRSF`: Lightning risk scaling factor (site-specific)
- `YLOI`: Previous day's LOI (for carry-over fires)
- `STMSPD`: Storm translational speed (m/s)
- `is_raining`: Flag (1 if raining at observation time, 0 otherwise)
- `is_lightning`: Flag (1 if lightning occurring, 0 otherwise)

# Variables
- `LOI`: Lightning-caused fire occurrence index (0-100)
- `LRISK`: Lightning risk
"""
@component function LightningFireOccurrenceIndex(; name=:LightningFireOccurrenceIndex)
    @constants begin
        # Reference units for non-dimensionalization
        V_ref = 1.0, [description = "Reference velocity for non-dimensionalization", unit = u"m/s"]
        L_ref = 1.0, [description = "Reference length for non-dimensionalization", unit = u"m"]
        L_zero = 0.0, [description = "Zero length for conditional expressions", unit = u"m"]
        s_ref = 1.0, [description = "Reference time for non-dimensionalization", unit = u"s"]
        s_zero = 0.0, [description = "Zero time for conditional expressions", unit = u"s"]
        # Unit conversions (dimensionless numerical conversion factors)
        mph_to_ms = 0.44704, [description = "mph to m/s conversion (dimensionless)", unit = u"1"]
        mi_to_m = 1609.344, [description = "miles to meters conversion (dimensionless)", unit = u"1"]
        m_to_mi = 0.000621371, [description = "meters to miles conversion (dimensionless)", unit = u"1"]
        ms_to_mph = 2.23694, [description = "m/s to mph conversion (dimensionless)", unit = u"1"]
    end

    @parameters begin
        LAL = 1.0, [description = "NFDRS Lightning Activity Level 1-6 (dimensionless)", unit = u"1"]
        IC, [description = "Ignition component 0-100 (dimensionless)", unit = u"1"]
        MC1, [description = "1-hour fuel moisture content (dimensionless)", unit = u"1"]
        LRSF = 1.0, [description = "Lightning risk scaling factor (dimensionless)", unit = u"1"]
        YLOI = 0.0, [description = "Previous day's LOI (dimensionless)", unit = u"1"]
        STMSPD = 13.4112, [description = "Storm translational speed", unit = u"m/s"]
        is_raining = 0.0, [description = "Flag: 1 if raining at observation time (dimensionless)", unit = u"1"]
        is_lightning = 0.0, [description = "Flag: 1 if lightning occurring (dimensionless)", unit = u"1"]
    end

    @variables begin
        CGRATE(t), [description = "Cloud-to-ground lightning discharge rate (dimensionless)", unit = u"1"]
        STMDIA(t), [description = "Width of rain-lightning corridor", unit = u"m"]
        TOTWID(t), [description = "Total width of affected corridor", unit = u"m"]
        LGTDUR(t), [description = "Duration of lightning at a point (dimensionless)", unit = u"1"]
        FINSID(t), [description = "Fraction of corridor with rain (dimensionless)", unit = u"1"]
        FOTSID(t), [description = "Fraction of corridor without rain (dimensionless)", unit = u"1"]
        RAIDUR(t), [description = "Rain duration at a point", unit = u"s"]
        FMF(t), [description = "1-hour fuel moisture inside rain area (dimensionless)", unit = u"1"]
        ICR(t), [description = "Ignition component inside rain area (dimensionless)", unit = u"1"]
        ICBAR(t), [description = "Area-weighted average ignition component (dimensionless)", unit = u"1"]
        LRISK(t), [description = "Lightning risk 0-100 (dimensionless)", unit = u"1"]
        LOI(t), [description = "Lightning-caused fire occurrence index 0-100 (dimensionless)", unit = u"1"]
    end

    # Non-dimensionalize and convert storm speed to mph for equations using original units
    STMSPD_mph = (STMSPD / V_ref) * ms_to_mph

    # Storm dimensions in miles (lookup by LAL)
    # LAL 2: STMDIA=3.0 mi, TOTWID=7.0 mi
    # LAL 3: STMDIA=4.0 mi, TOTWID=8.0 mi
    # LAL 4: STMDIA=5.0 mi, TOTWID=9.0 mi
    # LAL 5: STMDIA=7.0 mi, TOTWID=11.0 mi

    eqs = [
        # Cloud-to-ground lightning discharge rate by LAL
        CGRATE ~ ifelse(LAL < 1.5, 0.0,
                 ifelse(LAL < 2.5, 12.5,
                 ifelse(LAL < 3.5, 25.0,
                 ifelse(LAL < 4.5, 50.0,
                 ifelse(LAL < 5.5, 100.0, 100.0))))),

        # Storm corridor dimensions by LAL (convert miles to meters and re-dimensionalize)
        STMDIA ~ ifelse(LAL < 1.5, L_zero,
                 ifelse(LAL < 2.5, 3.0 * mi_to_m * L_ref,
                 ifelse(LAL < 3.5, 4.0 * mi_to_m * L_ref,
                 ifelse(LAL < 4.5, 5.0 * mi_to_m * L_ref,
                 ifelse(LAL < 5.5, 7.0 * mi_to_m * L_ref, 7.0 * mi_to_m * L_ref))))),

        TOTWID ~ ifelse(LAL < 1.5, L_zero,
                 ifelse(LAL < 2.5, 7.0 * mi_to_m * L_ref,
                 ifelse(LAL < 3.5, 8.0 * mi_to_m * L_ref,
                 ifelse(LAL < 4.5, 9.0 * mi_to_m * L_ref,
                 ifelse(LAL < 5.5, 11.0 * mi_to_m * L_ref, 11.0 * mi_to_m * L_ref))))),

        # Duration of lightning at a point (Eq. from page 13, in minutes)
        LGTDUR ~ ifelse(CGRATE > 0.0,
            -86.83 + 153.41 * CGRATE^0.1437,
            0.0),

        # Fraction of corridor occupied by rain-lightning vs lightning-only
        # Use miles for this calculation as equations are calibrated for miles
        FINSID ~ ifelse((TOTWID / L_ref) > 0.0,
            begin
                STMDIA_mi = (STMDIA / L_ref) * m_to_mi
                TOTWID_mi = (TOTWID / L_ref) * m_to_mi
                ((STMDIA_mi * STMSPD_mph * LGTDUR) + (0.7854 * STMDIA_mi^2)) /
                max(1e-10, (STMDIA_mi * STMSPD_mph * TOTWID_mi) + (0.7854 * TOTWID_mi^2))
            end,
            0.0),

        FOTSID ~ 1.0 - FINSID,

        # Rain duration at a point within rain corridor (seconds, re-dimensionalized)
        # STMDIA / STMSPD gives hours when both in miles and mph
        RAIDUR ~ ifelse(STMSPD_mph > 0.0, ((STMDIA / L_ref) * m_to_mi / STMSPD_mph) * 3600.0 * s_ref, s_zero),

        # Moisture content of 1-hour fuels inside rain area
        # Original: FMF = MC1 + ((76.0 + 2.7 * RAIDUR_hr) - MC1) * (1.0 - exp(-RAIDUR_hr))
        FMF ~ begin
            RAIDUR_hr = (RAIDUR / s_ref) / 3600.0
            MC1 + (((76.0 + 2.7 * RAIDUR_hr) / 100.0) - MC1) * (1.0 - exp(-RAIDUR_hr))
        end,

        # Ignition component inside rain area
        ICR ~ IC * max(0.0, 1.0 - (FMF - MC1) * 10.0),

        # Area-weighted ignition component
        ICBAR ~ (FINSID * ICR + FOTSID * IC) / 100.0,

        # Lightning risk (limited to 0-100)
        LRISK ~ min(100.0, CGRATE * LRSF),

        # Lightning-caused fire occurrence index
        LOI ~ ifelse(LAL >= 6.0 - 0.5,
            # LAL 6: extreme lightning
            100.0,
            ifelse(is_lightning < 0.5,
                # Not lightning
                min(100.0, 0.25 * YLOI),
                ifelse(is_raining > 0.5,
                    # Raining at observation time
                    min(100.0, 0.25 * YLOI),
                    # Lightning, not raining
                    min(100.0, 10.0 * LRISK * ICBAR + 0.25 * YLOI)
                )
            )
        )
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Fire Load Index (FLI)
# =============================================================================

"""
    FireLoadIndex(; name=:FireLoadIndex)

Calculate the Fire Load Index (FLI).

From Cohen & Deeming (1985), page 14:
FLI = 0.71 × √(BI² + (LOI + MCOI)²)

FLI combines the Burning Index with fire occurrence indexes to provide a
comprehensive measure of the expected fire load.

# Parameters
- `BI`: Burning index 0-100 (dimensionless)
- `LOI`: Lightning fire occurrence index 0-100 (dimensionless)
- `MCOI`: Human-caused fire occurrence index 0-100 (dimensionless)

# Variables
- `FLI`: Fire load index 0-100 (dimensionless)
"""
@component function FireLoadIndex(; name=:FireLoadIndex)
    @parameters begin
        BI, [description = "Burning index 0-100 (dimensionless)", unit = u"1"]
        LOI = 0.0, [description = "Lightning fire occurrence index 0-100 (dimensionless)", unit = u"1"]
        MCOI = 0.0, [description = "Human-caused fire occurrence index 0-100 (dimensionless)", unit = u"1"]
    end

    @variables begin
        FLI(t), [description = "Fire load index 0-100 (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Page 14: FLI = 0.71 × √(BI² + (LOI + MCOI)²)
        # BI and (LOI + MCOI) are each limited to 100
        FLI ~ 0.71 * sqrt(min(100.0, BI)^2 + min(100.0, LOI + MCOI)^2)
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Fuel Model Data (Named Tuples for backward compatibility)
# =============================================================================

"""
    NFDRSFuelModel

Fuel model parameters for the National Fire Danger Rating System.
Uses NamedTuple for lightweight, immutable storage of fuel model parameters.

All values in SI units:
- Fuel loadings: kg/m²
- SAV ratios: m⁻¹
- Depth: m
- Moisture of extinction: percent (same as original)
- Heat of combustion: J/kg
- SCM (spread component threshold): m/s
- WNDFC: dimensionless

Parameters converted from Cohen & Deeming (1985), Appendix (page 15).
"""
const NFDRSFuelModel = NamedTuple{(
    :name, :description,
    :SG1, :W1, :SG10, :W10, :SG100, :W100, :SG1000, :W1000,
    :SGWOOD, :WWOOD, :SGHERB, :WHERB,
    :DEPTH, :MXD, :HD, :HL, :SCM, :WNDFC
), Tuple{Symbol, String, Float64, Float64, Float64, Float64, Float64, Float64,
         Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
         Float64, Float64, Float64, Float64}}

# Unit conversion factors for fuel model data
const _SG_FT_TO_M = 3.28084          # ft⁻¹ to m⁻¹
const _W_TON_ACRE_TO_KG_M2 = 0.2241702  # tons/acre to kg/m²
const _DEPTH_FT_TO_M = 0.3048        # ft to m
const _H_BTU_LB_TO_J_KG = 2326.0     # Btu/lb to J/kg
const _SC_FPM_TO_MS = 0.00508        # ft/min to m/s (for SCM)

"""
    NFDRS_FUEL_MODELS

Dictionary of all 20 NFDRS fuel models (A-U, excluding M).

Parameters converted to SI units from Cohen & Deeming (1985), Appendix (page 15).
"""
const NFDRS_FUEL_MODELS = Dict{Symbol, NFDRSFuelModel}(
    :A => (name=:A, description="Western grasses (annual)",
        SG1=3000.0*_SG_FT_TO_M, W1=0.20*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=0.0*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=0.0*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=0.0*_SG_FT_TO_M, WWOOD=0.0*_W_TON_ACRE_TO_KG_M2, SGHERB=3000.0*_SG_FT_TO_M, WHERB=0.30*_W_TON_ACRE_TO_KG_M2,
        DEPTH=0.80*_DEPTH_FT_TO_M, MXD=15.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=300.0*_SC_FPM_TO_MS, WNDFC=0.6),
    :B => (name=:B, description="California chaparral",
        SG1=700.0*_SG_FT_TO_M, W1=3.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=4.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=0.50*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1250.0*_SG_FT_TO_M, WWOOD=11.50*_W_TON_ACRE_TO_KG_M2, SGHERB=0.0*_SG_FT_TO_M, WHERB=0.0*_W_TON_ACRE_TO_KG_M2,
        DEPTH=4.50*_DEPTH_FT_TO_M, MXD=15.0, HD=9500.0*_H_BTU_LB_TO_J_KG, HL=9500.0*_H_BTU_LB_TO_J_KG, SCM=58.0*_SC_FPM_TO_MS, WNDFC=0.5),
    :C => (name=:C, description="Pine-grass savanna",
        SG1=2000.0*_SG_FT_TO_M, W1=0.40*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=1.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=0.0*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=0.50*_W_TON_ACRE_TO_KG_M2, SGHERB=2500.0*_SG_FT_TO_M, WHERB=0.80*_W_TON_ACRE_TO_KG_M2,
        DEPTH=0.75*_DEPTH_FT_TO_M, MXD=20.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=32.0*_SC_FPM_TO_MS, WNDFC=0.4),
    :D => (name=:D, description="Southern rough",
        SG1=1250.0*_SG_FT_TO_M, W1=2.00*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=1.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=0.0*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=3.00*_W_TON_ACRE_TO_KG_M2, SGHERB=1500.0*_SG_FT_TO_M, WHERB=0.75*_W_TON_ACRE_TO_KG_M2,
        DEPTH=2.00*_DEPTH_FT_TO_M, MXD=30.0, HD=9000.0*_H_BTU_LB_TO_J_KG, HL=9000.0*_H_BTU_LB_TO_J_KG, SCM=25.0*_SC_FPM_TO_MS, WNDFC=0.4),
    :E => (name=:E, description="Hardwood litter (winter)",
        SG1=2000.0*_SG_FT_TO_M, W1=1.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=0.50*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=0.25*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=0.50*_W_TON_ACRE_TO_KG_M2, SGHERB=2000.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=0.40*_DEPTH_FT_TO_M, MXD=25.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=25.0*_SC_FPM_TO_MS, WNDFC=0.4),
    :F => (name=:F, description="Intermediate brush",
        SG1=700.0*_SG_FT_TO_M, W1=2.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=2.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=1.50*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1250.0*_SG_FT_TO_M, WWOOD=9.00*_W_TON_ACRE_TO_KG_M2, SGHERB=0.0*_SG_FT_TO_M, WHERB=0.0*_W_TON_ACRE_TO_KG_M2,
        DEPTH=4.50*_DEPTH_FT_TO_M, MXD=15.0, HD=9500.0*_H_BTU_LB_TO_J_KG, HL=9500.0*_H_BTU_LB_TO_J_KG, SCM=24.0*_SC_FPM_TO_MS, WNDFC=0.5),
    :G => (name=:G, description="Short needle (heavy dead)",
        SG1=2000.0*_SG_FT_TO_M, W1=2.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=2.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=5.00*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=12.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=0.50*_W_TON_ACRE_TO_KG_M2, SGHERB=2000.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=1.00*_DEPTH_FT_TO_M, MXD=25.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=30.0*_SC_FPM_TO_MS, WNDFC=0.4),
    :H => (name=:H, description="Short needle (normal dead)",
        SG1=2000.0*_SG_FT_TO_M, W1=1.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=1.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=2.00*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=2.00*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=0.50*_W_TON_ACRE_TO_KG_M2, SGHERB=2000.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=0.30*_DEPTH_FT_TO_M, MXD=20.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=8.0*_SC_FPM_TO_MS, WNDFC=0.4),
    :I => (name=:I, description="Heavy slash",
        SG1=1500.0*_SG_FT_TO_M, W1=12.00*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=12.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=10.00*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=12.00*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=0.0*_SG_FT_TO_M, WWOOD=0.0*_W_TON_ACRE_TO_KG_M2, SGHERB=0.0*_SG_FT_TO_M, WHERB=0.0*_W_TON_ACRE_TO_KG_M2,
        DEPTH=2.00*_DEPTH_FT_TO_M, MXD=25.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=65.0*_SC_FPM_TO_MS, WNDFC=0.5),
    :J => (name=:J, description="Intermediate slash",
        SG1=1500.0*_SG_FT_TO_M, W1=7.00*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=7.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=6.00*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=5.50*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=0.0*_SG_FT_TO_M, WWOOD=0.0*_W_TON_ACRE_TO_KG_M2, SGHERB=0.0*_SG_FT_TO_M, WHERB=0.0*_W_TON_ACRE_TO_KG_M2,
        DEPTH=1.30*_DEPTH_FT_TO_M, MXD=25.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=44.0*_SC_FPM_TO_MS, WNDFC=0.5),
    :K => (name=:K, description="Light slash",
        SG1=1500.0*_SG_FT_TO_M, W1=2.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=2.50*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=2.00*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=2.50*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=0.0*_SG_FT_TO_M, WWOOD=0.0*_W_TON_ACRE_TO_KG_M2, SGHERB=0.0*_SG_FT_TO_M, WHERB=0.0*_W_TON_ACRE_TO_KG_M2,
        DEPTH=0.60*_DEPTH_FT_TO_M, MXD=25.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=23.0*_SC_FPM_TO_MS, WNDFC=0.5),
    :L => (name=:L, description="Western grasses (perennial)",
        SG1=2000.0*_SG_FT_TO_M, W1=0.25*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=0.0*_W_TON_ACRE_TO_KG_M2, SG100=0.0*_SG_FT_TO_M, W100=0.0*_W_TON_ACRE_TO_KG_M2, SG1000=0.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=0.0*_SG_FT_TO_M, WWOOD=0.0*_W_TON_ACRE_TO_KG_M2, SGHERB=2000.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=1.00*_DEPTH_FT_TO_M, MXD=15.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=178.0*_SC_FPM_TO_MS, WNDFC=0.6),
    :N => (name=:N, description="Sawgrass",
        SG1=1600.0*_SG_FT_TO_M, W1=1.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=3.00*_W_TON_ACRE_TO_KG_M2, SG100=0.0*_SG_FT_TO_M, W100=0.0*_W_TON_ACRE_TO_KG_M2, SG1000=0.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=2.00*_W_TON_ACRE_TO_KG_M2, SGHERB=2000.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=3.00*_DEPTH_FT_TO_M, MXD=25.0, HD=8700.0*_H_BTU_LB_TO_J_KG, HL=8700.0*_H_BTU_LB_TO_J_KG, SCM=167.0*_SC_FPM_TO_MS, WNDFC=0.6),
    :O => (name=:O, description="High pocosin",
        SG1=1500.0*_SG_FT_TO_M, W1=2.00*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=1.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=3.00*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=2.00*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=7.00*_W_TON_ACRE_TO_KG_M2, SGHERB=0.0*_SG_FT_TO_M, WHERB=0.0*_W_TON_ACRE_TO_KG_M2,
        DEPTH=4.00*_DEPTH_FT_TO_M, MXD=30.0, HD=9000.0*_H_BTU_LB_TO_J_KG, HL=9000.0*_H_BTU_LB_TO_J_KG, SCM=99.0*_SC_FPM_TO_MS, WNDFC=0.5),
    :P => (name=:P, description="Southern pine plantation",
        SG1=1750.0*_SG_FT_TO_M, W1=1.00*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=2.50*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=0.50*_W_TON_ACRE_TO_KG_M2, SG1000=0.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=0.50*_W_TON_ACRE_TO_KG_M2, SGHERB=0.0*_SG_FT_TO_M, WHERB=0.0*_W_TON_ACRE_TO_KG_M2,
        DEPTH=0.40*_DEPTH_FT_TO_M, MXD=30.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=14.0*_SC_FPM_TO_MS, WNDFC=0.4),
    :Q => (name=:Q, description="Alaskan black spruce",
        SG1=1500.0*_SG_FT_TO_M, W1=2.00*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=0.50*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=2.00*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=1.00*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1200.0*_SG_FT_TO_M, WWOOD=4.00*_W_TON_ACRE_TO_KG_M2, SGHERB=1500.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=3.00*_DEPTH_FT_TO_M, MXD=25.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=59.0*_SC_FPM_TO_MS, WNDFC=0.4),
    :R => (name=:R, description="Hardwood litter (summer)",
        SG1=1500.0*_SG_FT_TO_M, W1=0.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=0.50*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=0.50*_W_TON_ACRE_TO_KG_M2, SG1000=0.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=0.50*_W_TON_ACRE_TO_KG_M2, SGHERB=2000.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=0.25*_DEPTH_FT_TO_M, MXD=25.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=6.0*_SC_FPM_TO_MS, WNDFC=0.4),
    :S => (name=:S, description="Tundra",
        SG1=1500.0*_SG_FT_TO_M, W1=0.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=0.50*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=0.50*_W_TON_ACRE_TO_KG_M2, SG1000=8.0*_SG_FT_TO_M, W1000=0.50*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1200.0*_SG_FT_TO_M, WWOOD=0.50*_W_TON_ACRE_TO_KG_M2, SGHERB=1500.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=0.40*_DEPTH_FT_TO_M, MXD=25.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=17.0*_SC_FPM_TO_MS, WNDFC=0.6),
    :T => (name=:T, description="Sagebrush-grass",
        SG1=2500.0*_SG_FT_TO_M, W1=1.00*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=1.50*_W_TON_ACRE_TO_KG_M2, SG100=0.0*_SG_FT_TO_M, W100=0.0*_W_TON_ACRE_TO_KG_M2, SG1000=0.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=2.50*_W_TON_ACRE_TO_KG_M2, SGHERB=2000.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=1.25*_DEPTH_FT_TO_M, MXD=15.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=73.0*_SC_FPM_TO_MS, WNDFC=0.6),
    :U => (name=:U, description="Western pines",
        SG1=1750.0*_SG_FT_TO_M, W1=1.50*_W_TON_ACRE_TO_KG_M2, SG10=109.0*_SG_FT_TO_M, W10=1.00*_W_TON_ACRE_TO_KG_M2, SG100=30.0*_SG_FT_TO_M, W100=0.50*_W_TON_ACRE_TO_KG_M2, SG1000=0.0*_SG_FT_TO_M, W1000=0.0*_W_TON_ACRE_TO_KG_M2,
        SGWOOD=1500.0*_SG_FT_TO_M, WWOOD=0.50*_W_TON_ACRE_TO_KG_M2, SGHERB=2000.0*_SG_FT_TO_M, WHERB=0.50*_W_TON_ACRE_TO_KG_M2,
        DEPTH=0.50*_DEPTH_FT_TO_M, MXD=20.0, HD=8000.0*_H_BTU_LB_TO_J_KG, HL=8000.0*_H_BTU_LB_TO_J_KG, SCM=16.0*_SC_FPM_TO_MS, WNDFC=0.4)
)

"""
    get_fuel_model(model::Symbol)

Get a fuel model by its symbol identifier (A-U).

Returns a NamedTuple with all fuel model parameters in SI units.

# Example
```julia
fm = get_fuel_model(:A)  # Western grasses (annual)
fm.SG1  # 9842.52 (m⁻¹)
fm.W1   # 0.0448 (kg/m²)
```
"""
function get_fuel_model(model::Symbol)
    if haskey(NFDRS_FUEL_MODELS, model)
        return NFDRS_FUEL_MODELS[model]
    else
        error("Unknown fuel model: $model. Available models: $(keys(NFDRS_FUEL_MODELS))")
    end
end

"""
    fuel_loading_to_kg_per_sqm(tons_per_acre)

Convert fuel loading from tons/acre to kg/m² (SI units).

The conversion factor is 0.2241702 (kg·acre)/(ton·m²).
"""
fuel_loading_to_kg_per_sqm(tons_per_acre) = tons_per_acre * _TONS_PER_ACRE_TO_KG_PER_M2

# Keep old function name for backward compatibility but mark as deprecated
"""
    fuel_loading_to_lb_per_sqft(tons_per_acre)

DEPRECATED: Use `fuel_loading_to_kg_per_sqm` instead.

Convert fuel loading from tons/acre to lb/ft² (old NFDRS internal units).
"""
fuel_loading_to_lb_per_sqft(tons_per_acre) = tons_per_acre * 0.0459137
