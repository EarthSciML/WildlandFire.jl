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

IMPORTANT: The NFDRS equations use imperial units as originally published:
- Temperature: °F (Fahrenheit)
- Fuel loading: tons/acre (converted to lb/ft² internally using factor 0.0459137)
- Fuel bed depth: ft (feet)
- Surface-area-to-volume ratio: ft⁻¹
- Heat of combustion: Btu/lb
- Wind speed: mph (miles per hour)
- Moisture content: fraction (0-1) or percent where noted

ModelingToolkit unit annotations are not applied because the NFDRS equations contain
empirical coefficients calibrated specifically for these imperial units. Converting to SI
would require recalibrating all coefficients, which is beyond the scope of this implementation.
"""

# t and D are imported in the main WildlandFire module

# =============================================================================
# Equilibrium Moisture Content (EMC)
# =============================================================================

"""
    EquilibriumMoistureContent(; name=:EquilibriumMoistureContent)

Calculate equilibrium moisture content (EMC) from temperature and relative humidity.

Based on regression equations developed by Simard (1968) from Wood Handbook tables.
Temperature is in degrees Fahrenheit, EMC is expressed as percent moisture content (fraction).

Implements Equations 1a, 1b, 1c from Cohen & Deeming (1985).

# Parameters
- `TEMP`: Dry bulb temperature (°F)
- `RH`: Relative humidity (fraction, 0-1)

# Variables
- `EMC`: Equilibrium moisture content (fraction)
"""
@component function EquilibriumMoistureContent(; name=:EquilibriumMoistureContent)
    @parameters begin
        TEMP, [description = "Dry bulb temperature (°F)"]
        RH, [description = "Relative humidity (dimensionless)"]
    end

    @variables begin
        EMC(t), [description = "Equilibrium moisture content (dimensionless)"]
    end

    # RH as percentage for equations
    RH_pct = RH * 100

    eqs = [
        # EMC equations based on RH ranges (Eq. 1a, 1b, 1c)
        EMC ~ ifelse(RH < 0.10,
            # Eq. 1a: RH < 10%
            (0.03229 + 0.281073 * RH_pct - 0.000578 * TEMP * RH_pct) / 100,
            ifelse(RH < 0.50,
                # Eq. 1b: 10% <= RH < 50%
                (2.22749 + 0.160107 * RH_pct - 0.014784 * TEMP) / 100,
                # Eq. 1c: RH >= 50%
                (21.0606 + 0.005565 * RH_pct^2 - 0.00035 * RH_pct * TEMP - 0.483199 * RH_pct) / 100
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

# Parameters
- `EMCPRM`: EMC at fuel-atmosphere interface (fraction)
- `MC10`: 10-hour fuel moisture content (fraction), used when fuel sticks are used
- `use_fuel_sticks`: Flag (1 if fuel sticks used, 0 otherwise)
- `is_raining`: Flag (1 if raining, 0 otherwise)

# Variables
- `MC1`: 1-hour fuel moisture content (fraction)
"""
@component function OneHourFuelMoisture(; name=:OneHourFuelMoisture)
    @parameters begin
        EMCPRM, [description = "EMC at fuel-atmosphere interface (fraction)"]
        MC10, [description = "10-hour fuel moisture content (fraction)"]
        use_fuel_sticks, [description = "Flag: 1 if fuel sticks used, 0 otherwise"]
        is_raining, [description = "Flag: 1 if raining, 0 otherwise"]
    end

    @variables begin
        MC1(t), [description = "1-hour fuel moisture content (fraction)"]
    end

    eqs = [
        MC1 ~ ifelse(is_raining > 0.5,
            0.35,  # 35% if raining
            ifelse(use_fuel_sticks > 0.5,
                (4.0 * EMCPRM + MC10) / 5.0,  # With fuel sticks
                1.03 * EMCPRM  # Without fuel sticks
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

Without fuel sticks, MC10 = 1.28 * EMCPRM.
With fuel sticks, age correction is applied based on Haines and Frost (1978).

# Parameters
- `EMCPRM`: EMC at fuel-atmosphere interface (fraction)
- `use_fuel_sticks`: Flag (1 if fuel sticks used, 0 otherwise)
- `WT`: Weight of fuel sticks (grams)
- `AGE`: Days since sticks were set out
- `CLIMAT`: NFDRS climate class (1-4)

# Variables
- `MC10`: 10-hour fuel moisture content (fraction)
"""
@component function TenHourFuelMoisture(; name=:TenHourFuelMoisture)
    @parameters begin
        EMCPRM, [description = "EMC at fuel-atmosphere interface (fraction)"]
        use_fuel_sticks, [description = "Flag: 1 if fuel sticks used, 0 otherwise"]
        WT = 100.0, [description = "Weight of fuel sticks (grams)"]
        AGE = 0.0, [description = "Days since sticks were set out"]
        CLIMAT = 2.0, [description = "NFDRS climate class (1-4)"]
    end

    @variables begin
        MC10(t), [description = "10-hour fuel moisture content (fraction)"]
    end

    # Age correction factors
    AA = 0.5 * AGE / 30.0
    BB = 1.0 + (0.02 * AGE / 30.0)
    CC = CLIMAT / 4.0

    eqs = [
        MC10 ~ ifelse(use_fuel_sticks > 0.5,
            # With fuel sticks (age corrected)
            (AA * CC + BB * CC * (WT - 100.0)) / 100.0,  # Convert to fraction
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

Calculate 100-hour timelag fuel moisture content as an ODE.

Uses weighted 24-hour average EMC (EMCBAR) based on hours of daylight.
The moisture content responds to the boundary condition with a time constant.

# Parameters
- `EMCMIN`: EMC at max temp/min RH (fraction)
- `EMCMAX`: EMC at min temp/max RH (fraction)
- `PPTDUR`: Duration of precipitation (hours)
- `LAT`: Station latitude (degrees)
- `JDATE`: Julian day of year (1-366)

# Variables
- `MC100`: 100-hour fuel moisture content (fraction)
"""
@component function HundredHourFuelMoisture(; name=:HundredHourFuelMoisture)
    @parameters begin
        EMCMIN, [description = "EMC at max temp/min RH (fraction)"]
        EMCMAX, [description = "EMC at min temp/max RH (fraction)"]
        PPTDUR = 0.0, [description = "Duration of precipitation (hours)"]
        LAT = 40.0, [description = "Station latitude (degrees)"]
        JDATE = 180.0, [description = "Julian day of year (1-366)"]
    end

    @variables begin
        MC100(t), [description = "100-hour fuel moisture content (fraction)"]
    end

    # Daylength calculation
    PHI = LAT * 0.01745  # Convert degrees to radians
    DECL = 0.41008 * sin((JDATE - 82) * 0.01745)
    # Limit tan product to avoid domain errors
    tan_product = tan(PHI) * tan(DECL)
    tan_product_clamped = max(-0.99, min(0.99, tan_product))
    DAYLIT = 24.0 * (1.0 - acos(tan_product_clamped) / 3.1416)

    # Weighted 24-hour average EMC
    EMCBAR = (DAYLIT * EMCMIN + (24.0 - DAYLIT) * EMCMAX) / 24.0

    # Weighted 24-hour average boundary condition (Eq. from page 5)
    # Note: The precipitation term uses percent moisture, so we convert EMCBAR to percent
    # and then convert the result back to fraction
    BNDRYH = ((24.0 - PPTDUR) * EMCBAR + PPTDUR * (0.5 * PPTDUR + 41.0) / 100.0) / 24.0

    # Response coefficient: (1.0 - 0.87 * exp(-0.24)) ≈ 0.1836
    response_coef = 1.0 - 0.87 * exp(-0.24)

    eqs = [
        # 100-hour fuel moisture differential equation
        # Daily update: MC100_new = MC100 + (BNDRYH - MC100) * response_coef
        # As continuous ODE approximation
        D(MC100) ~ (BNDRYH - MC100) * response_coef
    ]

    return System(eqs, t; name)
end

# =============================================================================
# 1000-Hour Timelag Fuel Moisture
# =============================================================================

"""
    ThousandHourFuelMoisture(; name=:ThousandHourFuelMoisture)

Calculate 1000-hour timelag fuel moisture content.

Uses 7-day running average of boundary conditions (BDYBAR).
Due to the long time constant, this is computed as an ODE.

# Parameters
- `EMCBAR`: Weighted 24-hour average EMC (fraction)
- `PPTDUR`: Duration of precipitation (hours)
- `BDYBAR`: 7-day running average boundary condition (fraction)

# Variables
- `MC1000`: 1000-hour fuel moisture content (fraction)
"""
@component function ThousandHourFuelMoisture(; name=:ThousandHourFuelMoisture)
    @parameters begin
        EMCBAR, [description = "Weighted 24-hour average EMC (fraction)"]
        PPTDUR = 0.0, [description = "Duration of precipitation (hours)"]
        BDYBAR, [description = "7-day running average boundary condition (fraction)"]
    end

    @variables begin
        MC1000(t), [description = "1000-hour fuel moisture content (fraction)"]
    end

    # Response coefficient: (1.0 - 0.82 * exp(-0.168)) ≈ 0.1532
    response_coef = 1.0 - 0.82 * exp(-0.168)

    eqs = [
        # 1000-hour fuel moisture differential equation
        D(MC1000) ~ (BDYBAR - MC1000) * response_coef
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

# Parameters
- `X1000`: Modified 1000-hr moisture for herbaceous model (fraction)
- `MC1`: 1-hour fuel moisture content (fraction)
- `CLIMAT`: NFDRS climate class (1-4)
- `is_annual`: Flag (1 if annual, 0 if perennial)
- `GRNDAY`: Days since greenup started
- `is_greenup`: Flag (1 if in greenup, 0 otherwise)
- `is_cured`: Flag (1 if cured/frozen, 0 otherwise)

# Variables
- `MCHERB`: Herbaceous fuel moisture content (fraction)
"""
@component function HerbaceousFuelMoisture(; name=:HerbaceousFuelMoisture)
    @parameters begin
        X1000, [description = "Modified 1000-hr moisture (fraction)"]
        MC1, [description = "1-hour fuel moisture content (fraction)"]
        CLIMAT = 2.0, [description = "NFDRS climate class (1-4)"]
        is_annual = 0.0, [description = "Flag: 1 if annual, 0 if perennial"]
        GRNDAY = 0.0, [description = "Days since greenup started"]
        is_greenup = 0.0, [description = "Flag: 1 if in greenup, 0 otherwise"]
        is_cured = 0.0, [description = "Flag: 1 if cured/frozen, 0 otherwise"]
    end

    @variables begin
        MCHERB(t), [description = "Herbaceous fuel moisture content (fraction)"]
    end

    # Climate class dependent parameters (X1000 as percent for these equations)
    X1000_pct = X1000 * 100

    # Greenup parameters
    HERBGA = ifelse(CLIMAT < 1.5, -70.0,
             ifelse(CLIMAT < 2.5, -100.0,
             ifelse(CLIMAT < 3.5, -137.5, -185.0)))
    HERBGB = ifelse(CLIMAT < 1.5, 12.8,
             ifelse(CLIMAT < 2.5, 14.0,
             ifelse(CLIMAT < 3.5, 15.5, 17.4)))

    # Transition parameters for annuals
    ANNTA = ifelse(CLIMAT < 1.5, -150.5,
            ifelse(CLIMAT < 2.5, -187.7,
            ifelse(CLIMAT < 3.5, -245.2, -305.2)))
    ANNTB = ifelse(CLIMAT < 1.5, 18.4,
            ifelse(CLIMAT < 2.5, 19.6,
            ifelse(CLIMAT < 3.5, 22.0, 24.3)))

    # Transition parameters for perennials
    PERTA = ifelse(CLIMAT < 1.5, 11.2,
            ifelse(CLIMAT < 2.5, -10.3,
            ifelse(CLIMAT < 3.5, -42.7, -93.5)))
    PERTB = ifelse(CLIMAT < 1.5, 7.4,
            ifelse(CLIMAT < 2.5, 8.3,
            ifelse(CLIMAT < 3.5, 9.8, 12.2)))

    # Greenup fraction (Eq. 4)
    GREN = min(1.0, GRNDAY / (7.0 * CLIMAT))

    # Potential herbaceous moisture (as fraction)
    MCHRBP = (HERBGA + HERBGB * X1000_pct) / 100.0

    eqs = [
        # MCHERB calculation based on stage
        MCHERB ~ ifelse(is_cured > 0.5,
            # Cured/frozen stage
            ifelse(is_annual > 0.5, MC1, max(0.30, min(1.50, (PERTA + PERTB * X1000_pct) / 100.0))),
            ifelse(is_greenup > 0.5,
                # Greenup stage: phase from 0.30 to MCHRBP
                0.30 + (max(0.30, MCHRBP) - 0.30) * GREN,
                # Green or transition stage
                ifelse(MCHRBP > 1.20,
                    # Green stage (MCHERB > 120%)
                    min(2.50, MCHRBP),
                    # Transition stage (30% < MCHERB < 120%)
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

# Parameters
- `MC1000`: 1000-hour fuel moisture content (fraction)
- `CLIMAT`: NFDRS climate class (1-4)
- `GRNDAY`: Days since greenup started
- `is_greenup`: Flag (1 if in greenup, 0 otherwise)
- `is_frozen`: Flag (1 if frozen/dormant, 0 otherwise)

# Variables
- `MCWOOD`: Woody fuel moisture content (fraction)
"""
@component function WoodyFuelMoisture(; name=:WoodyFuelMoisture)
    @parameters begin
        MC1000, [description = "1000-hour fuel moisture content (fraction)"]
        CLIMAT = 2.0, [description = "NFDRS climate class (1-4)"]
        GRNDAY = 0.0, [description = "Days since greenup started"]
        is_greenup = 0.0, [description = "Flag: 1 if in greenup, 0 otherwise"]
        is_frozen = 0.0, [description = "Flag: 1 if frozen/dormant, 0 otherwise"]
    end

    @variables begin
        MCWOOD(t), [description = "Woody fuel moisture content (fraction)"]
    end

    # MC1000 as percent for equations
    MC1000_pct = MC1000 * 100

    # PREGRN values by climate class
    PREGRN = ifelse(CLIMAT < 1.5, 0.50,
             ifelse(CLIMAT < 2.5, 0.60,
             ifelse(CLIMAT < 3.5, 0.70, 0.80)))

    # WOODGA and WOODGB by climate class
    WOODGA = ifelse(CLIMAT < 1.5, 12.5,
             ifelse(CLIMAT < 2.5, -5.0,
             ifelse(CLIMAT < 3.5, -22.5, -45.0)))
    WOODGB = ifelse(CLIMAT < 1.5, 7.5,
             ifelse(CLIMAT < 2.5, 8.2,
             ifelse(CLIMAT < 3.5, 8.9, 9.8)))

    # Greenup fraction
    GREN = min(1.0, GRNDAY / (7.0 * CLIMAT))

    # Potential woody moisture (as fraction)
    MCWODP = (WOODGA + WOODGB * MC1000_pct) / 100.0

    eqs = [
        # MCWOOD calculation
        MCWOOD ~ ifelse(is_frozen > 0.5,
            # Frozen/dormant stage
            PREGRN,
            ifelse(is_greenup > 0.5,
                # Greenup stage
                PREGRN + (max(PREGRN, MCWODP) - PREGRN) * GREN,
                # Green stage
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

Implements Equations 5, 6, 7, 8 from Cohen & Deeming (1985).

# Parameters
- `MCHERB`: Herbaceous fuel moisture content (fraction)
- `W1`: Base 1-hour fuel loading (lb/ft²)
- `WHERB`: Herbaceous fuel loading (lb/ft²)

# Variables
- `FCTCUR`: Fraction transferred to 1-hour class
- `W1P`: Effective 1-hour fuel loading (lb/ft²)
- `WHERBP`: Remaining herbaceous loading (lb/ft²)
"""
@component function FuelLoadingTransfer(; name=:FuelLoadingTransfer)
    @parameters begin
        MCHERB, [description = "Herbaceous fuel moisture content (fraction)"]
        W1, [description = "Base 1-hour fuel loading (lb/ft²)"]
        WHERB, [description = "Herbaceous fuel loading (lb/ft²)"]
    end

    @variables begin
        FCTCUR(t), [description = "Fraction transferred to 1-hour class"]
        WHERBC(t), [description = "Herbaceous loading transferred (lb/ft²)"]
        W1P(t), [description = "Effective 1-hour fuel loading (lb/ft²)"]
        WHERBP(t), [description = "Remaining herbaceous loading (lb/ft²)"]
    end

    # MCHERB as percent for equation
    MCHERB_pct = MCHERB * 100

    eqs = [
        # Eq. 5: Fraction transferred
        FCTCUR ~ max(0.0, min(1.0, 1.33 - 0.0111 * MCHERB_pct)),

        # Eq. 6: Amount transferred
        WHERBC ~ FCTCUR * WHERB,

        # Eq. 7: Effective 1-hour loading
        W1P ~ W1 + WHERBC,

        # Eq. 8: Remaining herbaceous loading
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

The SC is the forward rate of spread of the flaming front in ft/min, rounded to integer.
Implements equations from pages 9-11 of Cohen & Deeming (1985).

# Parameters
Fuel loadings (lb/ft²):
- `W1P`: 1-hour fuel loading (including transferred herbaceous)
- `W10`: 10-hour fuel loading
- `W100`: 100-hour fuel loading
- `WHERBP`: Remaining herbaceous loading
- `WWOOD`: Woody fuel loading

Surface-area-to-volume ratios (ft⁻¹):
- `SG1`, `SG10`, `SG100`: Dead fuel SAV ratios
- `SGHERB`, `SGWOOD`: Live fuel SAV ratios

Other parameters:
- `MC1`, `MC10`, `MC100`: Dead fuel moisture contents (fraction)
- `MCHERB`, `MCWOOD`: Live fuel moisture contents (fraction)
- `MXD`: Dead fuel moisture of extinction (fraction)
- `HD`, `HL`: Heat of combustion (Btu/lb) for dead and live fuels
- `DEPTH`: Fuel bed depth (ft)
- `WS`: 20-ft windspeed (mph)
- `WNDFC`: Wind reduction factor
- `slope_class`: NFDRS slope class (1-5)
- `fuels_wet`: Flag (1 if fuels wet or snow-covered)

# Variables
- `SC`: Spread component (ft/min, rounded)
- `ROS`: Rate of spread (ft/min)
"""
@component function SpreadComponent(; name=:SpreadComponent)
    @parameters begin
        # Fuel loadings (lb/ft²)
        W1P, [description = "1-hour fuel loading including transferred herb (lb/ft²)"]
        W10, [description = "10-hour fuel loading (lb/ft²)"]
        W100, [description = "100-hour fuel loading (lb/ft²)"]
        WHERBP, [description = "Remaining herbaceous loading (lb/ft²)"]
        WWOOD, [description = "Woody fuel loading (lb/ft²)"]

        # SAV ratios (ft⁻¹)
        SG1, [description = "1-hour fuel SAV ratio (ft⁻¹)"]
        SG10, [description = "10-hour fuel SAV ratio (ft⁻¹)"]
        SG100, [description = "100-hour fuel SAV ratio (ft⁻¹)"]
        SGHERB, [description = "Herbaceous fuel SAV ratio (ft⁻¹)"]
        SGWOOD, [description = "Woody fuel SAV ratio (ft⁻¹)"]

        # Moisture contents (fraction)
        MC1, [description = "1-hour fuel moisture (fraction)"]
        MC10, [description = "10-hour fuel moisture (fraction)"]
        MC100, [description = "100-hour fuel moisture (fraction)"]
        MCHERB, [description = "Herbaceous fuel moisture (fraction)"]
        MCWOOD, [description = "Woody fuel moisture (fraction)"]

        # Other fuel model parameters
        MXD, [description = "Dead fuel moisture of extinction (percent, e.g., 15 for 15%)"]
        HD = 8000.0, [description = "Dead fuel heat of combustion (Btu/lb)"]
        HL = 8000.0, [description = "Live fuel heat of combustion (Btu/lb)"]
        DEPTH, [description = "Fuel bed depth (ft)"]

        # Environmental parameters
        WS = 0.0, [description = "20-ft windspeed (mph)"]
        WNDFC = 0.4, [description = "Wind reduction factor"]
        slope_class = 1.0, [description = "NFDRS slope class (1-5)"]
        fuels_wet = 0.0, [description = "Flag: 1 if fuels wet or snow-covered"]
    end

    @constants begin
        RHOD = 32.0, [description = "Dead fuel particle density (lb/ft³)"]
        RHOL = 32.0, [description = "Live fuel particle density (lb/ft³)"]
        STD = 0.0555, [description = "Dead fuel total mineral content (fraction)"]
        STL = 0.0555, [description = "Live fuel total mineral content (fraction)"]
        SD = 0.01, [description = "Dead fuel silica-free mineral content"]
        SL = 0.01, [description = "Live fuel silica-free mineral content"]
    end

    @variables begin
        SC(t), [description = "Spread component (ft/min)"]
        ROS(t), [description = "Rate of spread (ft/min)"]
        IR(t), [description = "Reaction intensity (Btu/ft²/min)"]
    end

    # Net loadings
    W1N = W1P * (1.0 - STD)
    W10N = W10 * (1.0 - STD)
    W100N = W100 * (1.0 - STD)
    WHERBN = WHERBP * (1.0 - STL)
    WWOODN = WWOOD * (1.0 - STL)

    # Total loadings
    WTOTD = W1P + W10 + W100
    WTOTL = WHERBP + WWOOD
    WTOT = WTOTD + WTOTL

    # Surface areas
    SA1 = (W1P / RHOD) * SG1
    SA10 = (W10 / RHOD) * SG10
    SA100 = (W100 / RHOD) * SG100
    SAHERB = (WHERBP / RHOL) * SGHERB
    SAWOOD = (WWOOD / RHOL) * SGWOOD

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

    # Characteristic SAV ratios
    SGBRD = F1 * SG1 + F10 * SG10 + F100 * SG100
    SGBRL = FHERB * SGHERB + FWOOD * SGWOOD
    SGBRT = FDEAD * SGBRD + FLIVE * SGBRL
    SGBRT_safe = max(1.0, SGBRT)

    # Bulk density
    RHOBED = (WTOT - 0.0) / max(0.01, DEPTH)  # W1000 excluded from bulk density

    # Packing ratio
    RHOBAR = 32.0  # Constant particle density
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
    HN1 = W1N * exp(-138.0 / max(1.0, SG1))
    HN10 = W10N * exp(-138.0 / max(1.0, SG10))
    HN100 = W100N * exp(-138.0 / max(1.0, SG100))
    HNHERB = WHERBN * exp(-500.0 / max(1.0, SGHERB))
    HNWOOD = WWOODN * exp(-500.0 / max(1.0, SGWOOD))

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

    # Calculate reaction intensity for wind limit check
    IR_calc = GMAOP * (WDEADN * HD * ETASD * ETAMD + WLIVEN * HL * ETASL * ETAML)

    # Wind speed effect with limit: if (WS * 88 * WNDFC) > (0.9 * IR), use (0.9 * IR) instead
    # This is the wind effect limit from page 11 of Cohen & Deeming (1985)
    wind_term = WS * 88.0 * WNDFC
    wind_term_limited = min(wind_term, 0.9 * IR_calc)

    # Wind effect multiplier with limit applied
    PHIWND = UFACT * wind_term_limited^B

    # Slope effect multiplier
    PHISLP = SLPFCT * max(0.01, BETBAR)^(-0.3)

    # Heat sink calculation (Eq. from page 11)
    HTSINK = RHOBED * (
        FDEAD * (F1 * exp(-138.0/max(1.0,SG1)) * (250.0 + 11.16 * MC1 * 100) +
                 F10 * exp(-138.0/max(1.0,SG10)) * (250.0 + 11.16 * MC10 * 100) +
                 F100 * exp(-138.0/max(1.0,SG100)) * (250.0 + 11.16 * MC100 * 100)) +
        FLIVE * (FHERB * exp(-138.0/max(1.0,SGHERB)) * (250.0 + 11.16 * MCHERB * 100) +
                 FWOOD * exp(-138.0/max(1.0,SGWOOD)) * (250.0 + 11.16 * MCWOOD * 100))
    )

    eqs = [
        # Reaction intensity
        IR ~ IR_calc,

        # Rate of spread (ft/min) - Eq. from page 11
        # ROS = IR * ZETA * (1 + PHISLP + PHIWND) / HTSINK
        ROS ~ ifelse(fuels_wet > 0.5, 0.0,
            IR_calc * ZETA * (1.0 + PHISLP + PHIWND) / max(0.01, HTSINK)
        ),

        # Spread component (rounded ROS)
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

# Parameters
Same as SpreadComponent, plus:
- `W1000`: 1000-hour fuel loading (lb/ft²)
- `SG1000`: 1000-hour fuel SAV ratio (ft⁻¹)
- `MC1000`: 1000-hour fuel moisture (fraction)

# Variables
- `ERC`: Energy release component
- `IRE`: Loading-weighted reaction intensity
"""
@component function EnergyReleaseComponent(; name=:EnergyReleaseComponent)
    @parameters begin
        # Fuel loadings (lb/ft²)
        W1P, [description = "1-hour fuel loading including transferred herb (lb/ft²)"]
        W10, [description = "10-hour fuel loading (lb/ft²)"]
        W100, [description = "100-hour fuel loading (lb/ft²)"]
        W1000 = 0.0, [description = "1000-hour fuel loading (lb/ft²)"]
        WHERBP, [description = "Remaining herbaceous loading (lb/ft²)"]
        WWOOD, [description = "Woody fuel loading (lb/ft²)"]

        # SAV ratios (ft⁻¹)
        SG1, [description = "1-hour fuel SAV ratio (ft⁻¹)"]
        SG10, [description = "10-hour fuel SAV ratio (ft⁻¹)"]
        SG100, [description = "100-hour fuel SAV ratio (ft⁻¹)"]
        SG1000 = 8.0, [description = "1000-hour fuel SAV ratio (ft⁻¹)"]
        SGHERB, [description = "Herbaceous fuel SAV ratio (ft⁻¹)"]
        SGWOOD, [description = "Woody fuel SAV ratio (ft⁻¹)"]

        # Moisture contents (fraction)
        MC1, [description = "1-hour fuel moisture (fraction)"]
        MC10, [description = "10-hour fuel moisture (fraction)"]
        MC100, [description = "100-hour fuel moisture (fraction)"]
        MC1000 = 0.15, [description = "1000-hour fuel moisture (fraction)"]
        MCHERB, [description = "Herbaceous fuel moisture (fraction)"]
        MCWOOD, [description = "Woody fuel moisture (fraction)"]

        # Other fuel model parameters
        MXD, [description = "Dead fuel moisture of extinction (percent, e.g., 15 for 15%)"]
        MXL, [description = "Live fuel moisture of extinction (percent, e.g., 15 for 15%)"]
        HD = 8000.0, [description = "Dead fuel heat of combustion (Btu/lb)"]
        HL = 8000.0, [description = "Live fuel heat of combustion (Btu/lb)"]
        DEPTH, [description = "Fuel bed depth (ft)"]
    end

    @constants begin
        RHOD = 32.0, [description = "Dead fuel particle density (lb/ft³)"]
        RHOL = 32.0, [description = "Live fuel particle density (lb/ft³)"]
        STD = 0.0555, [description = "Dead fuel total mineral content (fraction)"]
        STL = 0.0555, [description = "Live fuel total mineral content (fraction)"]
        SD = 0.01, [description = "Dead fuel silica-free mineral content"]
        SL = 0.01, [description = "Live fuel silica-free mineral content"]
    end

    @variables begin
        ERC(t), [description = "Energy release component"]
        IRE(t), [description = "Loading-weighted reaction intensity (Btu/ft²/min)"]
    end

    # Total loadings
    WTOTD = W1P + W10 + W100 + W1000
    WTOTL = WHERBP + WWOOD
    WTOT = WTOTD + WTOTL

    # Loading-based weighting factors
    WTOTD_safe = max(1e-10, WTOTD)
    WTOTL_safe = max(1e-10, WTOTL)
    WTOT_safe = max(1e-10, WTOT)

    F1E = W1P / WTOTD_safe
    F10E = W10 / WTOTD_safe
    F100E = W100 / WTOTD_safe
    F1000E = W1000 / WTOTD_safe
    FHERBE = WHERBP / WTOTL_safe
    FWOODE = WWOOD / WTOTL_safe

    FDEADE = WTOTD / WTOT_safe
    FLIVEE = WTOTL / WTOT_safe

    # Net loadings
    WDEDNE = WTOTD * (1.0 - STD)
    WLIVNE = WTOTL * (1.0 - STL)

    # Characteristic SAV ratio (loading weighted)
    SGBRDE = F1E * SG1 + F10E * SG10 + F100E * SG100 + F1000E * SG1000
    SGBRLE = FWOODE * SGWOOD + FHERBE * SGHERB
    SGBRTE = FDEADE * SGBRDE + FLIVEE * SGBRLE
    SGBRTE_safe = max(1.0, SGBRTE)

    # For residence time calculation, use surface-area weighted SGBRT
    # (calculated same as in SpreadComponent)
    SA1 = (W1P / RHOD) * SG1
    SA10 = (W10 / RHOD) * SG10
    SA100 = (W100 / RHOD) * SG100
    SAHERB = (WHERBP / RHOL) * SGHERB
    SAWOOD = (WWOOD / RHOL) * SGWOOD
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
    SGBRD = F1 * SG1 + F10 * SG10 + F100 * SG100
    SGBRL = FHERB * SGHERB + FWOOD * SGWOOD
    SGBRT = FDEAD * SGBRD + FLIVE * SGBRL
    SGBRT_safe = max(1.0, SGBRT)

    # Bulk density and packing ratio
    RHOBED = WTOT / max(0.01, DEPTH)
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

    # Moisture damping coefficients (convert fractions to percent for ratio with MXD/MXL)
    DEDRTE = (WTMCDE * 100) / MXD
    LIVRTE = (WTMCLE * 100) / max(0.01, MXL)

    ETAMDE = max(0.0, min(1.0, 1.0 - 2.0 * DEDRTE + 1.5 * DEDRTE^2 - 0.5 * DEDRTE^3))
    ETAMLE = max(0.0, min(1.0, 1.0 - 2.0 * LIVRTE + 1.5 * LIVRTE^2 - 0.5 * LIVRTE^3))

    # Mineral damping coefficients
    ETASD = 0.174 * SD^(-0.19)
    ETASL = 0.174 * SL^(-0.19)

    # Residence time (uses surface-area weighted SGBRT)
    TAU = 384.0 / SGBRT_safe

    eqs = [
        # Loading-weighted reaction intensity
        IRE ~ GMAOPE * (FDEADE * WDEDNE * HD * ETASD * ETAMDE +
                        FLIVEE * WLIVNE * HL * ETASL * ETAMLE),

        # Energy Release Component (0.04 scaling factor = ft²/Btu)
        ERC ~ 0.04 * IRE * TAU
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
based on Byram's flame length model.

# Parameters
- `SC`: Spread component (dimensionless index)
- `ERC`: Energy release component (dimensionless index)
- `fuels_wet`: Flag (1 if fuels wet or snow-covered, 0 otherwise)

# Variables
- `BI`: Burning index (dimensionless)
"""
@component function BurningIndex(; name=:BurningIndex)
    @parameters begin
        SC, [description = "Spread component"]
        ERC, [description = "Energy release component"]
        fuels_wet = 0.0, [description = "Flag: 1 if fuels wet or snow-covered"]
    end

    @variables begin
        BI(t), [description = "Burning index"]
    end

    eqs = [
        # BI = 10 * FL where FL = 0.301 * (SC * ERC)^0.46
        # So BI = 3.01 * (SC * ERC)^0.46
        BI ~ ifelse(fuels_wet > 0.5, 0.0, 3.01 * (max(0.0, SC * ERC))^0.46)
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Ignition Component (IC)
# =============================================================================

"""
    IgnitionComponent(; name=:IgnitionComponent)

Calculate the NFDRS Ignition Component (IC).

IC consists of:
1. P(I): Probability that a firebrand produces successful ignition
2. P(F/I): Probability that ignition results in reportable fire

P(I) is scaled to be 100 at MC1=1.5% and 0 at MC1=25%.

# Parameters
- `TMPPRM`: Temperature at fuel-atmosphere interface (°C)
- `MC1`: 1-hour fuel moisture content (fraction)
- `SC`: Spread component
- `SCM`: Spread component threshold for reportable fires

# Variables
- `IC`: Ignition component (0-100)
"""
@component function IgnitionComponent(; name=:IgnitionComponent)
    @parameters begin
        TMPPRM, [description = "Temperature at fuel-atmosphere interface (°C)"]
        MC1, [description = "1-hour fuel moisture content (fraction)"]
        SC, [description = "Spread component"]
        SCM, [description = "Spread component threshold for reportable fires"]
    end

    @constants begin
        PNORM1 = 0.00232, [description = "P(I) scaling factor 1"]
        PNORM2 = 0.99767, [description = "P(I) scaling factor 2"]
        PNORM3 = 0.0000185, [description = "P(I) scaling factor 3"]
    end

    @variables begin
        QIGN(t), [description = "Heat of ignition (cal/g)"]
        CHI(t), [description = "Intermediate variable"]
        PI(t), [description = "Probability of ignition (0-100)"]
        SCN(t), [description = "Normalized spread component"]
        PFI(t), [description = "Probability of reportable fire given ignition"]
        IC(t), [description = "Ignition component (0-100)"]
    end

    eqs = [
        # Heat of ignition (Eq. 57) - MC1 as percent (MC1*100)
        QIGN ~ 144.5 - 0.266 * TMPPRM - 0.00058 * TMPPRM^2 - 0.01 * TMPPRM * (MC1 * 100) + 18.54 * (1.0 - exp(-0.151 * (MC1 * 100))) + 6.4 * (MC1 * 100),

        # Intermediate calculation
        CHI ~ max(0.0, (344.0 - QIGN) / 10.0),

        # Probability of ignition
        PI ~ ifelse(CHI^3.6 * PNORM3 <= PNORM1,
            0.0,
            max(0.0, min(100.0, (CHI^3.6 * PNORM3 - PNORM1) * 100.0 / PNORM2))
        ),

        # Normalized spread component
        SCN ~ 100.0 * SC / max(1.0, SCM),

        # Probability of reportable fire given ignition
        PFI ~ sqrt(max(0.0, SCN)),

        # Ignition component
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

# Parameters
- `MRISK`: Human-caused risk (0-100)
- `IC`: Ignition component (0-100)

# Variables
- `MCOI`: Human-caused fire occurrence index
"""
@component function HumanFireOccurrenceIndex(; name=:HumanFireOccurrenceIndex)
    @parameters begin
        MRISK, [description = "Human-caused risk (0-100)"]
        IC, [description = "Ignition component (0-100)"]
    end

    @variables begin
        MCOI(t), [description = "Human-caused fire occurrence index"]
    end

    eqs = [
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
- `MC1`: 1-hour fuel moisture content (fraction)
- `LRSF`: Lightning risk scaling factor (site-specific)
- `YLOI`: Previous day's LOI (for carry-over fires)
- `STMSPD`: Storm translational speed (mph), default 30
- `is_raining`: Flag (1 if raining at observation time, 0 otherwise)
- `is_lightning`: Flag (1 if lightning occurring, 0 otherwise)

# Variables
- `LOI`: Lightning-caused fire occurrence index (0-100)
- `LRISK`: Lightning risk
"""
@component function LightningFireOccurrenceIndex(; name=:LightningFireOccurrenceIndex)
    @parameters begin
        LAL = 1.0, [description = "NFDRS Lightning Activity Level (1-6)"]
        IC, [description = "Ignition component (0-100)"]
        MC1, [description = "1-hour fuel moisture content (fraction)"]
        LRSF = 1.0, [description = "Lightning risk scaling factor"]
        YLOI = 0.0, [description = "Previous day's LOI"]
        STMSPD = 30.0, [description = "Storm translational speed (mph)"]
        is_raining = 0.0, [description = "Flag: 1 if raining at observation time"]
        is_lightning = 0.0, [description = "Flag: 1 if lightning occurring"]
    end

    @variables begin
        CGRATE(t), [description = "Cloud-to-ground lightning discharge rate (strikes/min)"]
        STMDIA(t), [description = "Width of rain-lightning corridor (miles)"]
        TOTWID(t), [description = "Total width of affected corridor (miles)"]
        LGTDUR(t), [description = "Duration of lightning at a point (min)"]
        FINSID(t), [description = "Fraction of corridor with rain"]
        FOTSID(t), [description = "Fraction of corridor without rain"]
        RAIDUR(t), [description = "Rain duration at a point (hours)"]
        FMF(t), [description = "1-hour fuel moisture inside rain area (fraction)"]
        ICR(t), [description = "Ignition component inside rain area"]
        ICBAR(t), [description = "Area-weighted average ignition component"]
        LRISK(t), [description = "Lightning risk (0-100)"]
        LOI(t), [description = "Lightning-caused fire occurrence index (0-100)"]
    end

    # LAL 6 is a special case (extreme lightning)
    # For LAL 1-5, use lookup tables from page 13
    # CGRATE: 0, 12.5, 25.0, 50.0, 100.0 strikes/min
    # STMDIA: 0, 3.0, 4.0, 5.0, 7.0 miles
    # TOTWID: 0, 7.0, 8.0, 9.0, 11.0 miles

    eqs = [
        # Cloud-to-ground lightning discharge rate by LAL
        CGRATE ~ ifelse(LAL < 1.5, 0.0,
                 ifelse(LAL < 2.5, 12.5,
                 ifelse(LAL < 3.5, 25.0,
                 ifelse(LAL < 4.5, 50.0,
                 ifelse(LAL < 5.5, 100.0, 100.0))))),

        # Storm corridor dimensions by LAL
        STMDIA ~ ifelse(LAL < 1.5, 0.0,
                 ifelse(LAL < 2.5, 3.0,
                 ifelse(LAL < 3.5, 4.0,
                 ifelse(LAL < 4.5, 5.0,
                 ifelse(LAL < 5.5, 7.0, 7.0))))),

        TOTWID ~ ifelse(LAL < 1.5, 0.0,
                 ifelse(LAL < 2.5, 7.0,
                 ifelse(LAL < 3.5, 8.0,
                 ifelse(LAL < 4.5, 9.0,
                 ifelse(LAL < 5.5, 11.0, 11.0))))),

        # Duration of lightning at a point (Eq. from page 13)
        LGTDUR ~ ifelse(CGRATE > 0.0,
            -86.83 + 153.41 * CGRATE^0.1437,
            0.0),

        # Fraction of corridor occupied by rain-lightning vs lightning-only
        # FINSID = ((STMDIA * STMSPD * LGTDUR) + (0.7854 * STMDIA^2)) /
        #          ((STMDIA * STMSPD * TOTWID) + (0.7854 * TOTWID^2))
        FINSID ~ ifelse(TOTWID > 0.0,
            ((STMDIA * STMSPD * LGTDUR) + (0.7854 * STMDIA^2)) /
            max(1e-10, (STMDIA * STMSPD * TOTWID) + (0.7854 * TOTWID^2)),
            0.0),

        FOTSID ~ 1.0 - FINSID,

        # Rain duration at a point within rain corridor (hours)
        RAIDUR ~ ifelse(STMSPD > 0.0, STMDIA / STMSPD, 0.0),

        # Moisture content of 1-hour fuels inside rain area
        # FMF = MC1 + ((76.0 + 2.7 * RAIDUR) - MC1) * (1.0 - exp(-RAIDUR))
        # Note: 76.0 and 2.7 are in percent, but we use fractions
        FMF ~ MC1 + (((76.0 + 2.7 * RAIDUR) / 100.0) - MC1) * (1.0 - exp(-RAIDUR)),

        # Ignition component inside rain area (simplified - uses higher moisture)
        # ICR is calculated like IC but with FMF instead of MC1
        # For simplicity, we reduce IC proportionally to increased moisture
        ICR ~ IC * max(0.0, 1.0 - (FMF - MC1) * 10.0),

        # Area-weighted ignition component
        ICBAR ~ (FINSID * ICR + FOTSID * IC) / 100.0,

        # Lightning risk (limited to 0-100)
        LRISK ~ min(100.0, CGRATE * LRSF),

        # Lightning-caused fire occurrence index
        # If not lightning or raining, use 25% of previous day's LOI for carry-over
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

FLI combines the Burning Index with fire occurrence indexes.

# Parameters
- `BI`: Burning index
- `LOI`: Lightning fire occurrence index
- `MCOI`: Human-caused fire occurrence index

# Variables
- `FLI`: Fire load index
"""
@component function FireLoadIndex(; name=:FireLoadIndex)
    @parameters begin
        BI, [description = "Burning index"]
        LOI = 0.0, [description = "Lightning fire occurrence index"]
        MCOI = 0.0, [description = "Human-caused fire occurrence index"]
    end

    @variables begin
        FLI(t), [description = "Fire load index"]
    end

    eqs = [
        # BI and (LOI + MCOI) are each limited to 100
        FLI ~ 0.71 * sqrt(min(100.0, BI)^2 + min(100.0, LOI + MCOI)^2)
    ]

    return System(eqs, t; name)
end

# =============================================================================
# Fuel Model Data Structure
# =============================================================================

"""
    NFDRSFuelModel

Fuel model parameters for the National Fire Danger Rating System.

All fuel loadings in tons/acre, SAV ratios in ft⁻¹, depth in ft,
moisture of extinction as percent, heat of combustion in Btu/lb.

# Fields
- `name::Symbol`: Fuel model identifier (A-U)
- `description::String`: Fuel model description
- `SG1::Float64`: 1-hour fuel SAV ratio (ft⁻¹)
- `W1::Float64`: 1-hour fuel loading (tons/acre)
- `SG10::Float64`: 10-hour fuel SAV ratio (ft⁻¹)
- `W10::Float64`: 10-hour fuel loading (tons/acre)
- `SG100::Float64`: 100-hour fuel SAV ratio (ft⁻¹)
- `W100::Float64`: 100-hour fuel loading (tons/acre)
- `SG1000::Float64`: 1000-hour fuel SAV ratio (ft⁻¹)
- `W1000::Float64`: 1000-hour fuel loading (tons/acre)
- `SGWOOD::Float64`: Woody fuel SAV ratio (ft⁻¹)
- `WWOOD::Float64`: Woody fuel loading (tons/acre)
- `SGHERB::Float64`: Herbaceous fuel SAV ratio (ft⁻¹)
- `WHERB::Float64`: Herbaceous fuel loading (tons/acre)
- `DEPTH::Float64`: Fuel bed depth (ft)
- `MXD::Float64`: Dead fuel moisture of extinction (percent)
- `HD::Float64`: Dead fuel heat of combustion (Btu/lb)
- `HL::Float64`: Live fuel heat of combustion (Btu/lb)
- `SCM::Float64`: Spread component threshold for all fires reportable
- `WNDFC::Float64`: Wind reduction factor (20-ft to midflame)
"""
struct NFDRSFuelModel
    name::Symbol
    description::String
    SG1::Float64
    W1::Float64
    SG10::Float64
    W10::Float64
    SG100::Float64
    W100::Float64
    SG1000::Float64
    W1000::Float64
    SGWOOD::Float64
    WWOOD::Float64
    SGHERB::Float64
    WHERB::Float64
    DEPTH::Float64
    MXD::Float64
    HD::Float64
    HL::Float64
    SCM::Float64
    WNDFC::Float64
end

"""
    NFDRS_FUEL_MODELS

Dictionary of all 20 NFDRS fuel models (A-U, excluding M).

Parameters from Cohen & Deeming (1985), Appendix.
"""
const NFDRS_FUEL_MODELS = Dict{Symbol, NFDRSFuelModel}(
    :A => NFDRSFuelModel(:A, "Western grasses (annual)",
        3000.0, 0.20, 109.0, 0.0, 30.0, 0.0, 8.0, 0.0, 0.0, 0.0, 3000.0, 0.30,
        0.80, 15.0, 8000.0, 8000.0, 300.0, 0.6),
    :B => NFDRSFuelModel(:B, "California chaparral",
        700.0, 3.50, 109.0, 4.00, 30.0, 0.50, 8.0, 0.0, 1250.0, 11.50, 0.0, 0.0,
        4.50, 15.0, 9500.0, 9500.0, 58.0, 0.5),
    :C => NFDRSFuelModel(:C, "Pine-grass savanna",
        2000.0, 0.40, 109.0, 1.00, 30.0, 0.0, 8.0, 0.0, 1500.0, 0.50, 2500.0, 0.80,
        0.75, 20.0, 8000.0, 8000.0, 32.0, 0.4),
    :D => NFDRSFuelModel(:D, "Southern rough",
        1250.0, 2.00, 109.0, 1.00, 30.0, 0.0, 8.0, 0.0, 1500.0, 3.00, 1500.0, 0.75,
        2.00, 30.0, 9000.0, 9000.0, 25.0, 0.4),
    :E => NFDRSFuelModel(:E, "Hardwood litter (winter)",
        2000.0, 1.50, 109.0, 0.50, 30.0, 0.25, 8.0, 0.0, 1500.0, 0.50, 2000.0, 0.50,
        0.40, 25.0, 8000.0, 8000.0, 25.0, 0.4),
    :F => NFDRSFuelModel(:F, "Intermediate brush",
        700.0, 2.50, 109.0, 2.00, 30.0, 1.50, 8.0, 0.0, 1250.0, 9.00, 0.0, 0.0,
        4.50, 15.0, 9500.0, 9500.0, 24.0, 0.5),
    :G => NFDRSFuelModel(:G, "Short needle (heavy dead)",
        2000.0, 2.50, 109.0, 2.00, 30.0, 5.00, 8.0, 12.0, 1500.0, 0.50, 2000.0, 0.50,
        1.00, 25.0, 8000.0, 8000.0, 30.0, 0.4),
    :H => NFDRSFuelModel(:H, "Short needle (normal dead)",
        2000.0, 1.50, 109.0, 1.00, 30.0, 2.00, 8.0, 2.00, 1500.0, 0.50, 2000.0, 0.50,
        0.30, 20.0, 8000.0, 8000.0, 8.0, 0.4),
    :I => NFDRSFuelModel(:I, "Heavy slash",
        1500.0, 12.00, 109.0, 12.00, 30.0, 10.00, 8.0, 12.00, 0.0, 0.0, 0.0, 0.0,
        2.00, 25.0, 8000.0, 8000.0, 65.0, 0.5),
    :J => NFDRSFuelModel(:J, "Intermediate slash",
        1500.0, 7.00, 109.0, 7.00, 30.0, 6.00, 8.0, 5.50, 0.0, 0.0, 0.0, 0.0,
        1.30, 25.0, 8000.0, 8000.0, 44.0, 0.5),
    :K => NFDRSFuelModel(:K, "Light slash",
        1500.0, 2.50, 109.0, 2.50, 30.0, 2.00, 8.0, 2.50, 0.0, 0.0, 0.0, 0.0,
        0.60, 25.0, 8000.0, 8000.0, 23.0, 0.5),
    :L => NFDRSFuelModel(:L, "Western grasses (perennial)",
        2000.0, 0.25, 109.0, 1.50, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2000.0, 0.50,
        1.00, 15.0, 8000.0, 8000.0, 178.0, 0.6),
    :N => NFDRSFuelModel(:N, "Sawgrass",
        1600.0, 1.50, 109.0, 3.00, 0.0, 0.0, 0.0, 0.0, 1500.0, 2.00, 2000.0, 0.50,
        3.00, 25.0, 8700.0, 8700.0, 167.0, 0.6),
    :O => NFDRSFuelModel(:O, "High pocosin",
        1500.0, 2.00, 109.0, 1.00, 30.0, 3.00, 8.0, 2.00, 1500.0, 7.00, 0.0, 0.0,
        4.00, 30.0, 9000.0, 9000.0, 99.0, 0.5),
    :P => NFDRSFuelModel(:P, "Southern pine plantation",
        1750.0, 1.00, 109.0, 2.50, 30.0, 0.50, 0.0, 0.0, 1500.0, 0.50, 0.0, 0.0,
        0.40, 30.0, 8000.0, 8000.0, 14.0, 0.4),
    :Q => NFDRSFuelModel(:Q, "Alaskan black spruce",
        1500.0, 2.00, 109.0, 0.50, 30.0, 2.00, 8.0, 1.00, 1200.0, 4.00, 1500.0, 0.50,
        3.00, 25.0, 8000.0, 8000.0, 59.0, 0.4),
    :R => NFDRSFuelModel(:R, "Hardwood litter (summer)",
        1500.0, 0.50, 109.0, 0.50, 30.0, 0.50, 0.0, 0.0, 1500.0, 0.50, 2000.0, 0.50,
        0.25, 25.0, 8000.0, 8000.0, 6.0, 0.4),
    :S => NFDRSFuelModel(:S, "Tundra",
        1500.0, 0.50, 109.0, 0.50, 30.0, 0.50, 8.0, 0.50, 1200.0, 0.50, 1500.0, 0.50,
        0.40, 25.0, 8000.0, 8000.0, 17.0, 0.6),
    :T => NFDRSFuelModel(:T, "Sagebrush-grass",
        2500.0, 1.00, 109.0, 1.50, 0.0, 0.0, 0.0, 0.0, 1500.0, 2.50, 2000.0, 0.50,
        1.25, 15.0, 8000.0, 8000.0, 73.0, 0.6),
    :U => NFDRSFuelModel(:U, "Western pines",
        1750.0, 1.50, 109.0, 1.00, 30.0, 0.50, 0.0, 0.0, 1500.0, 0.50, 2000.0, 0.50,
        0.50, 20.0, 8000.0, 8000.0, 16.0, 0.4)
)

"""
    get_fuel_model(model::Symbol)

Get a fuel model by its symbol identifier (A-U).

# Example
```julia
fm = get_fuel_model(:A)  # Western grasses (annual)
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
    fuel_loading_to_lb_per_sqft(tons_per_acre)

Convert fuel loading from tons/acre to lb/ft² (NFDRS internal units).

The conversion factor is 0.0459137 (lb·acre)/(ton·ft²).
"""
fuel_loading_to_lb_per_sqft(tons_per_acre) = tons_per_acre * 0.0459137
