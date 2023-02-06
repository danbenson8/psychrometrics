// Package psychrometrics calculates psychrometric properties of moist and dry air.
package psychrometrics

import (
	"fmt"
	"math"
)

/************************************************************************************************
	Global constants
************************************************************************************************/

// float: Zero degree Fahrenheit (°F) expressed as degree Rankine (°R)
// Units:
// °R
// Reference:
// ASHRAE Handbook - Fundamentals (2017) ch. 39
const ZERO_FAHRENHEIT_AS_RANKINE = 459.67

/*
float: Zero degree Celsius (°C) expressed as Kelvin (K)
Units:
K
Reference:
ASHRAE Handbook - Fundamentals (2017) ch. 39
*/
const ZERO_CELSIUS_AS_KELVIN = 273.15

/*
float: Universal gas constant for dry air (IP version)
Units:
ft lb_Force lb_DryAir⁻¹ R⁻¹
Reference:
ASHRAE Handbook - Fundamentals (2017) ch. 1
*/
const R_DA_IP = 53.350

/*
float: Universal gas constant for dry air (SI version)

	Units:
	    J kg_DryAir⁻¹ K⁻¹
	Reference:
	    ASHRAE Handbook - Fundamentals (2017) ch. 1
*/
const R_DA_SI = 287.042

/*int: Maximum number of iterations before exiting while loops.
 */
const MAX_ITER_COUNT = 100

/*
float: Minimum acceptable humidity ratio used/returned by any functions.

	Any value above 0 or below the MIN_HUM_RATIO will be reset to this value.
*/
const MIN_HUM_RATIO = 1.0e-7

/*float: Freezing point of water in Fahrenheit.
 */
const FREEZING_POINT_WATER_IP = 32.0

/*float: Freezing point of water in Celsius.
 */
const FREEZING_POINT_WATER_SI = 0.0

/*float: Triple point of water in Fahrenheit.
 */
const TRIPLE_POINT_WATER_IP = 32.018

/*float: Triple point of water in Celsius.
 */
const TRIPLE_POINT_WATER_SI = 0.01

/************************************************************************************************
    Helper functions
************************************************************************************************/

type UnitSystem int8

func (u UnitSystem) String() string {
	switch u {
	case IP:
		return "IP"
	case SI:
		return "SI"
	}
	return "undefined"
}

// System of units, (IP or SI)
const (
	Undefined UnitSystem = iota
	IP
	SI
)

var PSYCHROMETRICS_UNITS = Undefined

// Tolerance for temperature calculations
var PSYCHROMETRICS_TOLERANCE = 1.

// Set the system of units to use (IP or SI)
func SetUnitSystem(units UnitSystem) {
	PSYCHROMETRICS_UNITS = units

	if PSYCHROMETRICS_UNITS == IP {
		PSYCHROMETRICS_TOLERANCE = 0.001 * 9 / 5
	} else {
		PSYCHROMETRICS_TOLERANCE = 0.001
	}
}

// Return system of units in use
func GetUnitSystem() UnitSystem {
	return PSYCHROMETRICS_UNITS
}

// Assertion
func Assert(b bool, msg string) {
	if !b {
		panic(fmt.Errorf(msg))
	}
}

/************************************************************************************************
	Conversion between temperature units
************************************************************************************************/

// Utility function to convert temperature to degree Rankine (°R)
// given temperature in degree Fahrenheit (°F).
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 section 3
func GetTRankineFromTFahrenheit(T_F float64) float64 { return T_F + ZERO_FAHRENHEIT_AS_RANKINE } /* exact */

// Utility function to convert temperature to degree Fahrenheit (°F)
// given temperature in degree Rankine (°R).
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 section 3
func GetTFahrenheitFromTRankine(T_R float64) float64 { return T_R - ZERO_FAHRENHEIT_AS_RANKINE } /* exact */

// Utility function to convert temperature to Kelvin (K)
// given temperature in degree Celsius (°C).
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 section 3
func GetTKelvinFromTCelsius(T_C float64) float64 { return T_C + ZERO_CELSIUS_AS_KELVIN } /* exact */

// Utility function to convert temperature to degree Celsius (°C)
// given temperature in Kelvin (K).
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 section 3
func GetTCelsiusFromTKelvin(T_K float64) float64 { return T_K - ZERO_CELSIUS_AS_KELVIN } /* exact */

/************************************************************************************************
	Conversions between dew point, wet bulb, and relative humidity
************************************************************************************************/

// Return wet-bulb temperature [°F/°C] given dry-bulb temperature [°F/°C], dew-point temperature [°F/°C]
// and, atmospheric pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetTWetBulbFromTDewPoint(TDryBulb, TDewPoint, Pressure float64) float64 {
	Assert(TDewPoint <= TDryBulb, "Dew point temperature is above dry bulb temperature")
	HumRatio := GetHumRatioFromTDewPoint(TDewPoint, Pressure)
	return GetTWetBulbFromHumRatio(TDryBulb, HumRatio, Pressure)
}

// Return wet-bulb temperature [°F/°C] given dry-bulb temperature [°F/°C], relative humidity [0-1]
// and, atmospheric pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetTWetBulbFromRelHum(TDryBulb, RelHum, Pressure float64) float64 {
	Assert(RelHum >= 0 && RelHum <= 1, "Relative humidity is outside range [0,1]")
	HumRatio := GetHumRatioFromRelHum(TDryBulb, RelHum, Pressure)
	return GetTWetBulbFromHumRatio(TDryBulb, HumRatio, Pressure)
}

// Return relative humidity [0-1] given dry-bulb temperature [°F/°C] and dew-point temperature [°F/°C].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 22
func GetRelHumFromTDewPoint(TDryBulb, TDewPoint float64) float64 {
	Assert(TDewPoint <= TDryBulb, "Dew point temperature is above dry bulb temperature")
	VapPres := GetSatVapPres(TDewPoint)
	SatVapPres := GetSatVapPres(TDryBulb)
	return VapPres / SatVapPres
}

// Return relative humidity [0-1] given dry-bulb temperature [°F/°C], wet bulb temperature [°F/°C]
// and, atmospheric pressure [psi, Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetRelHumFromTWetBulb(TDryBulb, TWetBulb, Pressure float64) float64 {
	Assert(TWetBulb <= TDryBulb, "Wet bulb temperature is above dry bulb temperature")
	HumRatio := GetHumRatioFromTWetBulb(TDryBulb, TWetBulb, Pressure)
	return GetRelHumFromHumRatio(TDryBulb, HumRatio, Pressure)
}

// Return dew-point temperature [°F/°C] given dry-bulb temperature [°F/°C] and relative humidity [0-1].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetTDewPointFromRelHum(TDryBulb, RelHum float64) float64 {
	Assert(RelHum >= 0 && RelHum <= 1, "Relative humidity is outside range [0,1]")
	VapPres := GetVapPresFromRelHum(TDryBulb, RelHum)
	return GetTDewPointFromVapPres(TDryBulb, VapPres)
}

// Return dew-point temperature [°F/°C] given dry-bulb temperature [°F/°C],
// wet-bulb temperature [°F/°C], and atmospheric pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetTDewPointFromTWetBulb(TDryBulb, TWetBulb, Pressure float64) float64 {
	Assert(TWetBulb <= TDryBulb, "Wet bulb temperature is above dry bulb temperature")
	HumRatio := GetHumRatioFromTWetBulb(TDryBulb, TWetBulb, Pressure)
	return GetTDewPointFromHumRatio(TDryBulb, HumRatio, Pressure)
}

/******************************************************************************************************
	Conversions between dew point, or relative humidity and vapor pressure
******************************************************************************************************/

// Return partial pressure of water vapor [psi/Pa] as a function of relative humidity [0-1]
// and dry-bulb temperature [°F, °C].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 12, 22
func GetVapPresFromRelHum(TDryBulb, RelHum float64) float64 {
	Assert(RelHum >= 0 && RelHum <= 1, "Relative humidity is outside range [0,1]")
	return RelHum * GetSatVapPres(TDryBulb)
}

// Return relative humidity [0-1] given dry-bulb temperature [°F/°C] and vapor pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 12, 22
func GetRelHumFromVapPres(TDryBulb, VapPres float64) float64 {
	Assert(VapPres >= 0., "Partial pressure of water vapour in moist air is negative")
	return VapPres / GetSatVapPres(TDryBulb)
}

// Helper function returning the derivative of the natural log of the saturation vapor pressure [psi/Pa]
// as a function of dry-bulb temperature [°F/°C].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn. 5 & 6
func dLnPws_(TDryBulb float64) float64 {
	var dLnPws, T float64

	if GetUnitSystem() == IP {
		T = GetTRankineFromTFahrenheit(TDryBulb)
		if TDryBulb <= TRIPLE_POINT_WATER_IP {
			dLnPws = 1.0214165e+04/math.Pow(T, 2) - 5.3765794e-03 + 2*1.9202377e-07*T + 3*3.5575832e-10*math.Pow(T, 2) - 4*9.0344688e-14*math.Pow(T, 3) + 4.1635019/T
		} else {
			dLnPws = 1.0440397e+04/math.Pow(T, 2) - 2.7022355e-02 + 2*1.2890360e-05*T - 3*2.4780681e-09*math.Pow(T, 2) + 6.5459673/T
		}
	} else {
		T = GetTKelvinFromTCelsius(TDryBulb)

		if TDryBulb <= TRIPLE_POINT_WATER_SI {
			dLnPws = 5.6745359e+03/math.Pow(T, 2) - 9.677843e-03 + 2*6.2215701e-07*T + 3*2.0747825e-09*math.Pow(T, 2) - 4*9.484024e-13*math.Pow(T, 3) + 4.1635019/T
		} else {
			dLnPws = 5.8002206e+03/math.Pow(T, 2) - 4.8640239e-02 + 2*4.1764768e-05*T - 3*1.4452093e-08*math.Pow(T, 2) + 6.5459673/T
		}
	}

	return dLnPws
}

// Return dew-point temperature [°F/°C] given dry-bulb temperature [°F/°C] and vapor pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn. 5 and 6
// Notes: the dew point temperature is solved by inverting the equation giving water vapor pressure
// at saturation from temperature rather than using the regressions provided
// by ASHRAE (eqn. 37 and 38) which are much less accurate and have a
// narrower range of validity.
// The Newton-Raphson (NR) method is used on the logarithm of water vapour
// pressure as a function of temperature, which is a very smooth function
// Convergence is usually achieved in 3 to 5 iterations.
// TDryBulb is not really needed here, just used for convenience.
func GetTDewPointFromVapPres(TDryBulb, VapPres float64) float64 {
	BOUNDS := make([]float64, 0, 2)

	if GetUnitSystem() == IP {
		BOUNDS = append(BOUNDS, -148, 392)
	} else {
		BOUNDS = append(BOUNDS, -100, 200)
	}

	// Bounds outside which a solution cannot be found
	Assert(VapPres >= GetSatVapPres(BOUNDS[0]) && VapPres <= GetSatVapPres(BOUNDS[1]), "Partial pressure of water vapor is outside range of validity of equations")

	// We use NR to approximate the solution.
	// First guess
	TDewPoint := TDryBulb     // Calculated value of dew point temperatures, solved for iteratively in °F/°C [SI]
	lnVP := math.Log(VapPres) // Natural logarithm of partial pressure of water vapor pressure in moist air

	var TDewPoint_iter float64 // Value of TDewPoint used in NR calculation
	var lnVP_iter float64      // Value of log of vapor water pressure used in NR calculation
	index := 1

	for {

		// TDewPoint used in NR calculation
		TDewPoint_iter = TDewPoint
		lnVP_iter = math.Log(GetSatVapPres(TDewPoint_iter))

		// Derivative of function, calculated analytically
		d_lnVP := dLnPws_(TDewPoint_iter)

		// New estimate, bounded by domain of validity of eqn. 5 and 6
		TDewPoint = TDewPoint_iter - (lnVP_iter-lnVP)/d_lnVP
		TDewPoint = math.Max(TDewPoint, BOUNDS[0])
		TDewPoint = math.Min(TDewPoint, BOUNDS[1])

		index++

		if math.Abs(TDewPoint-TDewPoint_iter) < PSYCHROMETRICS_TOLERANCE {
			break
		}

		Assert(index <= MAX_ITER_COUNT, "Convergence not reached in GetTDewPointFromVapPres. Stopping.")
	}
	return math.Min(TDewPoint, TDryBulb)
}

// Return vapor pressure [psi/Pa] given dew point temperature [°F/°C].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn. 36
func GetVapPresFromTDewPoint(TDewPoint float64) float64 {
	return GetSatVapPres(TDewPoint)
}

/******************************************************************************************************
	Conversions from wet-bulb temperature, dew-point temperature, or relative humidity to humidity ratio
******************************************************************************************************/

// Return wet-bulb temperature given dry-bulb temperature, humidity ratio, and pressure.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 33 and 35 solved for Tstar
func GetTWetBulbFromHumRatio(TDryBulb float64, HumRatio float64, Pressure float64) float64 {
	// Declarations
	var Wstar float64
	index := 1

	Assert(HumRatio >= 0, "Humidity ratio is negative")

	BoundedHumRatio := math.Max(HumRatio, MIN_HUM_RATIO)
	TDewPoint := GetTDewPointFromHumRatio(TDryBulb, BoundedHumRatio, Pressure)

	// Initial guesses
	TWetBulbSup := TDryBulb
	TWetBulbInf := TDewPoint
	TWetBulb := (TWetBulbInf + TWetBulbSup) / 2.

	// Bisection loop, break when converged below tolerance
	for {
		// Compute humidity ratio at temperature Tstar
		Wstar = GetHumRatioFromTWetBulb(TDryBulb, TWetBulb, Pressure)

		// Get new bounds
		if Wstar > BoundedHumRatio {
			TWetBulbSup = TWetBulb
		} else {
			TWetBulbInf = TWetBulb
		}

		// New guess of wet bulb temperature
		TWetBulb = (TWetBulbSup + TWetBulbInf) / 2.
		index++

		if (TWetBulbSup - TWetBulbInf) < PSYCHROMETRICS_TOLERANCE {
			break
		}
	}
	Assert(index <= MAX_ITER_COUNT, "Convergence not reached in GetTWetBulbFromHumRatio. Stopping.")
	return TWetBulb

}

// Return humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹]given dry-bulb temperature [°F/°C], wet-bulb temperature [°F/°C],
// and atmospheric pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 33 and 35
func GetHumRatioFromTWetBulb(TDryBulb, TWetBulb, Pressure float64) float64 {
	var HumRatio float64
	Assert(TWetBulb <= TDryBulb, "Wet bulb temperature is above dry bulb temperature")

	Wsstar := GetSatHumRatio(TWetBulb, Pressure)

	if GetUnitSystem() == IP {
		if TWetBulb >= FREEZING_POINT_WATER_IP {
			HumRatio = ((1093.-0.556*TWetBulb)*Wsstar - 0.240*(TDryBulb-TWetBulb)) / (1093. + 0.444*TDryBulb - TWetBulb)
		} else {
			HumRatio = ((1220.-0.04*TWetBulb)*Wsstar - 0.240*(TDryBulb-TWetBulb)) / (1220. + 0.444*TDryBulb - 0.48*TWetBulb)
		}
	} else {
		if TWetBulb >= FREEZING_POINT_WATER_IP {
			HumRatio = ((2501.-2.326*TWetBulb)*Wsstar - 1.006*(TDryBulb-TWetBulb)) / (2501. + 1.86*TDryBulb - 4.186*TWetBulb)
		} else {
			HumRatio = ((2830.-0.24*TWetBulb)*Wsstar - 1.006*(TDryBulb-TWetBulb)) / (2830. + 1.86*TDryBulb - 2.1*TWetBulb)
		}
	}

	// Validity check.
	return math.Max(HumRatio, MIN_HUM_RATIO)
}

// Return humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹] given dry-bulb temperature [°F/°C], relative humidity [0-1],
// and atmospheric pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetHumRatioFromRelHum(TDryBulb, RelHum, Pressure float64) float64 {
	Assert(RelHum >= 0. && RelHum <= 1., "Relative humidity is outside range [0,1]")
	VapPres := GetVapPresFromRelHum(TDryBulb, RelHum)
	return GetHumRatioFromVapPres(VapPres, Pressure)
}

// Return relative humidity [0-1] given dry-bulb temperature [°F/°C], humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹]
// and atmospheric pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetRelHumFromHumRatio(TDryBulb, HumRatio, Pressure float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	VapPres := GetVapPresFromHumRatio(HumRatio, Pressure)
	return GetRelHumFromVapPres(TDryBulb, VapPres)
}

// Return humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹] given dew-point temperature [°F/C]
// and atmospheric pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetHumRatioFromTDewPoint(TDewPoint, Pressure float64) float64 {
	VapPres := GetSatVapPres(TDewPoint)
	return GetHumRatioFromVapPres(VapPres, Pressure)
}

// Return dew-point temperature given dry-bulb temperature [°F/°C], humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹], and pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetTDewPointFromHumRatio(TDryBulb float64, HumRatio float64, Pressure float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	VapPres := GetVapPresFromHumRatio(HumRatio, Pressure)
	return GetTDewPointFromVapPres(TDryBulb, VapPres)
}

/******************************************************************************************************
	Conversions between humidity ratio and vapor pressure
******************************************************************************************************/

// Return humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹] given water vapor pressure [psi/Pa] and atmospheric pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 20
func GetHumRatioFromVapPres(VapPres, Pressure float64) float64 {
	Assert(VapPres >= 0., "Partial pressure of water vapor in moist air is negative")
	HumRatio := 0.621945 * VapPres / (Pressure - VapPres)

	// Validity check.
	return math.Max(HumRatio, MIN_HUM_RATIO)

}

// Return vapor pressure [psi/Pa] given humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹] and pressure [psi/Pa].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 20 solved for pw
func GetVapPresFromHumRatio(HumRatio, Pressure float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	BoundedHumRatio := math.Max(HumRatio, MIN_HUM_RATIO)
	VapPres := Pressure * BoundedHumRatio / (0.621945 + BoundedHumRatio)
	return VapPres
}

/******************************************************************************************************
	Conversions between humidity ratio and specific humidity
******************************************************************************************************/

// Return the specific humidity [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹] from humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Dry_Air⁻¹] (aka mixing ratio)
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 9b
func GetSpecificHumFromHumRatio(HumRatio float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	BoundedHumRatio := math.Max(HumRatio, MIN_HUM_RATIO)

	return BoundedHumRatio / (1.0 + BoundedHumRatio)
}

// Return the humidity ratio (aka mixing ratio) from specific humidity
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 9b (solved for humidity ratio)
func GetHumRatioFromSpecificHum(SpecificHum float64) float64 {
	Assert(SpecificHum >= 0.0 && SpecificHum < 1.0, "Specific humidity is outside range [0, 1)")

	HumRatio := SpecificHum / (1.0 - SpecificHum)

	// Validity check
	return math.Max(HumRatio, MIN_HUM_RATIO)
}

/******************************************************************************************************
	Dry Air Calculations
******************************************************************************************************/

// Return dry-air enthalpy given dry-bulb temperature.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn. 28
func GetDryAirEnthalpy(TDryBulb float64) float64 {
	if GetUnitSystem() == IP {
		return 0.240 * TDryBulb
	} else {
		return 1006 * TDryBulb
	}
}

// Return dry-air density given dry-bulb temperature and pressure.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
// Notes: eqn 14 for the perfect gas relationship for dry air.
// Eqn 1 for the universal gas constant.
// The factor 144 in IP is for the conversion of Psi = lb in⁻² to lb ft⁻².
func GetDryAirDensity(TDryBulb, Pressure float64) float64 {
	if GetUnitSystem() == IP {
		return (144. * Pressure) / R_DA_IP / GetTRankineFromTFahrenheit(TDryBulb)
	} else {
		return Pressure / R_DA_SI / GetTKelvinFromTCelsius(TDryBulb)
	}
}

// Return dry-air volume given dry-bulb temperature and pressure.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
// Notes: eqn 14 for the perfect gas relationship for dry air.
// Eqn 1 for the universal gas constant.
// The factor 144 in IP is for the conversion of Psi = lb in⁻² to lb ft⁻².
func GetDryAirVolume(TDryBulb, Pressure float64) float64 {
	if GetUnitSystem() == IP {
		return R_DA_IP * GetTRankineFromTFahrenheit(TDryBulb) / (144. * Pressure)
	} else {
		return R_DA_SI * GetTKelvinFromTCelsius(TDryBulb) / Pressure
	}
}

// Return dry bulb temperature from enthalpy and humidity ratio.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 30.
// Notes: based on the `GetMoistAirEnthalpy` function, rearranged for temperature.
func GetTDryBulbFromEnthalpyAndHumRatio(MoistAirEnthalpy, HumRatio float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	BoundedHumRatio := math.Max(HumRatio, MIN_HUM_RATIO)

	if GetUnitSystem() == IP {
		return (MoistAirEnthalpy - 1061.0*BoundedHumRatio) / (0.240 + 0.444*BoundedHumRatio)
	} else {
		return (MoistAirEnthalpy/1000.0 - 2501.0*BoundedHumRatio) / (1.006 + 1.86*BoundedHumRatio)
	}
}

// Return humidity ratio from enthalpy and dry-bulb temperature.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 30.
// Notes: based on the `GetMoistAirEnthalpy` function, rearranged for humidity ratio.
func GetHumRatioFromEnthalpyAndTDryBulb(MoistAirEnthalpy, TDryBulb float64) float64 {
	var HumRatio float64
	if GetUnitSystem() == IP {
		HumRatio = (MoistAirEnthalpy - 0.240*TDryBulb) / (1061.0 + 0.444*TDryBulb)
	} else {
		HumRatio = (MoistAirEnthalpy/1000.0 - 1.006*TDryBulb) / (2501.0 + 1.86*TDryBulb)
	}

	// Validity check.
	return math.Max(HumRatio, MIN_HUM_RATIO)
}

/******************************************************************************************************
	Saturated Air Calculations
*******************************************************************************************************/

// Return saturation vapor pressure [psi/Pa] given dry-bulb temperature [°F/°C].
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn. 5 & 6
// Important note: the ASHRAE formulae are defined above and below the freezing point but have
// a discontinuity at the freezing point. This is a small inaccuracy on ASHRAE's part: the formulae
// should be defined above and below the triple point of water (not the feezing point) in which case
// the discontinuity vanishes. It is essential to use the triple point of water otherwise function
// GetTDewPointFromVapPres, which inverts the present function, does not converge properly around
// the freezing point.
func GetSatVapPres(TDryBulb float64) float64 {

	var LnPws, T float64

	if GetUnitSystem() == IP {
		Assert(TDryBulb >= -148. && TDryBulb <= 392., "Dry bulb temperature is outside range [-148, 392]")

		T = GetTRankineFromTFahrenheit(TDryBulb)

		if TDryBulb <= TRIPLE_POINT_WATER_IP {
			LnPws = (-1.0214165e+04/T - 4.8932428 - 5.3765794e-03*T + 1.9202377e-07*T*T + 3.5575832e-10*math.Pow(T, 3) - 9.0344688e-14*math.Pow(T, 4) + 4.1635019*math.Log(T))
		} else {
			LnPws = -1.0440397e+04/T - 1.1294650e+01 - 2.7022355e-02*T + 1.2890360e-05*T*T - 2.4780681e-09*math.Pow(T, 3) + 6.5459673*math.Log(T)
		}
	} else {
		Assert(TDryBulb >= -100. && TDryBulb <= 200., "Dry bulb temperature is outside range [-100, 200]")

		T = GetTKelvinFromTCelsius(TDryBulb)

		if TDryBulb <= TRIPLE_POINT_WATER_SI {
			LnPws = -5.6745359e+03/T + 6.3925247 - 9.677843e-03*T + 6.2215701e-07*T*T + 2.0747825e-09*math.Pow(T, 3) - 9.484024e-13*math.Pow(T, 4) + 4.1635019*math.Log(T)
		} else {
			LnPws = -5.8002206e+03/T + 1.3914993 - 4.8640239e-02*T + 4.1764768e-05*T*T - 1.4452093e-08*math.Pow(T, 3) + 6.5459673*math.Log(T)
		}
	}

	return math.Exp(LnPws)
}

// Return humidity ratio of saturated air given dry-bulb temperature and pressure.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 36, solved for W
func GetSatHumRatio(TDryBulb, Pressure float64) float64 {

	SatVaporPres := GetSatVapPres(TDryBulb)
	SatHumRatio := 0.621945 * SatVaporPres / (Pressure - SatVaporPres)

	// Validity check.
	return math.Max(SatHumRatio, MIN_HUM_RATIO)
}

// Return saturated air enthalpy given dry-bulb temperature and pressure.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1
func GetSatAirEnthalpy(TDryBulb, Pressure float64) float64 {
	return GetMoistAirEnthalpy(TDryBulb, GetSatHumRatio(TDryBulb, Pressure))
}

/******************************************************************************************************
	Moist Air Calculations
******************************************************************************************************/

// Return Vapor pressure deficit given dry-bulb temperature, humidity ratio, and pressure.
// Reference: see Oke (1987) eqn. 2.13a
func GetVaporPressureDeficit(TDryBulb, HumRatio, Pressure float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	RelHum := GetRelHumFromHumRatio(TDryBulb, HumRatio, Pressure)
	return GetSatVapPres(TDryBulb) * (1. - RelHum)
}

// Return the degree of saturation (i.e humidity ratio of the air / humidity ratio of the air at saturation
// at the same temperature and pressure) given dry-bulb temperature, humidity ratio, and atmospheric pressure.
// Reference: ASHRAE Handbook - Fundamentals (2009) ch. 1 eqn. 12
// Notes: the definition is absent from the 2017 Handbook
func GetDegreeOfSaturation(TDryBulb, HumRatio, Pressure float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	BoundedHumRatio := math.Max(HumRatio, MIN_HUM_RATIO)

	return BoundedHumRatio / GetSatHumRatio(TDryBulb, Pressure)
}

// Return moist air enthalpy given dry-bulb temperature and humidity ratio.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn. 30
func GetMoistAirEnthalpy(TDryBulb, HumRatio float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	BoundedHumRatio := math.Max(HumRatio, MIN_HUM_RATIO)

	if GetUnitSystem() == IP {
		return 0.240*TDryBulb + BoundedHumRatio*(1061.+0.444*TDryBulb)
	} else {
		return (1.006*TDryBulb + BoundedHumRatio*(2501.+1.86*TDryBulb)) * 1000.
	}
}

// Return moist air specific volume given dry-bulb temperature, humidity ratio, and pressure.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn. 26
// Notes: in IP units, R_DA_IP / 144 equals 0.370486 which is the coefficient appearing in eqn 26.
// The factor 144 is for the conversion of Psi = lb in⁻² to lb ft⁻².
func GetMoistAirVolume(TDryBulb, HumRatio, Pressure float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	BoundedHumRatio := math.Max(HumRatio, MIN_HUM_RATIO)

	if GetUnitSystem() == IP {
		return R_DA_IP * GetTRankineFromTFahrenheit(TDryBulb) * (1. + 1.607858*BoundedHumRatio) / (144. * Pressure)
	} else {
		return R_DA_SI * GetTKelvinFromTCelsius(TDryBulb) * (1. + 1.607858*BoundedHumRatio) / Pressure
	}
}

// Return dry-bulb temperature given moist air specific volume, humidity ratio, and pressure.
// Reference:
// ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 26
// Notes:
// In IP units, R_DA_IP / 144 equals 0.370486 which is the coefficient appearing in eqn 26
// The factor 144 is for the conversion of Psi = lb in⁻² to lb ft⁻².
// Based on the `GetMoistAirVolume` function, rearranged for dry-bulb temperature.
func GetTDryBulbFromMoistAirVolumeAndHumRatio(MoistAirVolume, HumRatio, Pressure float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	BoundedHumRatio := math.Max(HumRatio, MIN_HUM_RATIO)

	if GetUnitSystem() == IP {
		return GetTFahrenheitFromTRankine(MoistAirVolume * (144 * Pressure) / (R_DA_IP * (1 + 1.607858*BoundedHumRatio)))
	} else {
		return GetTCelsiusFromTKelvin(MoistAirVolume * Pressure / (R_DA_SI * (1 + 1.607858*BoundedHumRatio)))
	}
}

// Return moist air density given humidity ratio, dry bulb temperature, and pressure.
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn. 11
func GetMoistAirDensity(TDryBulb, HumRatio, Pressure float64) float64 {
	Assert(HumRatio >= 0., "Humidity ratio is negative")
	BoundedHumRatio := math.Max(HumRatio, MIN_HUM_RATIO)

	return (1. + BoundedHumRatio) / GetMoistAirVolume(TDryBulb, BoundedHumRatio, Pressure)
}

/******************************************************************************************************
	Standard atmosphere
******************************************************************************************************/

// Return standard atmosphere barometric pressure, given the elevation (altitude).
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 3
func GetStandardAtmPressure(Altitude float64) float64 {
	var Pressure float64
	if GetUnitSystem() == IP {
		Pressure = 14.696 * math.Pow(1.-6.8754e-06*Altitude, 5.2559)
	} else {
		Pressure = 101325. * math.Pow(1.-2.25577e-05*Altitude, 5.2559)
	}
	return Pressure
}

// Return standard atmosphere temperature, given the elevation (altitude).
// Reference: ASHRAE Handbook - Fundamentals (2017) ch. 1 eqn 4
func GetStandardAtmTemperature(Altitude float64) float64 {
	var Temperature float64
	if GetUnitSystem() == IP {
		Temperature = 59. - 0.00356620*Altitude
	} else {
		Temperature = 15. - 0.0065*Altitude
	}
	return Temperature
}

// Return sea level pressure given dry-bulb temperature, altitude above sea level and pressure.
// Reference: Hess SL, Introduction to theoretical meteorology, Holt Rinehart and Winston, NY 1959,
// ch. 6.5 Stull RB, Meteorology for scientists and engineers, 2nd edition,
// Brooks/Cole 2000, ch. 1.
// Notes: the standard procedure for the US is to use for TDryBulb the average
// of the current station temperature and the station temperature from 12 hours ago.
func GetSeaLevelPressure(StnPressure, Altitude, TDryBulb float64) float64 {
	var TColumn, H float64
	if GetUnitSystem() == IP {
		// Calculate average temperature in column of air, assuming a lapse rate
		// of 3.6 °F/1000ft
		TColumn = TDryBulb + 0.0036*Altitude/2.

		// Determine the scale height
		H = 53.351 * GetTRankineFromTFahrenheit(TColumn)
	} else {
		// Calculate average temperature in column of air, assuming a lapse rate
		// of 6.5 °C/km
		TColumn = TDryBulb + 0.0065*Altitude/2.

		// Determine the scale height
		H = 287.055 * GetTKelvinFromTCelsius(TColumn) / 9.807
	}

	// Calculate the sea level pressure
	SeaLevelPressure := StnPressure * math.Exp(Altitude/H)
	return SeaLevelPressure
}

// Return station pressure [psi/Pa] from sea level pressure [psi/Pa] given sea level pressure [psi/Pa],
// altitude [ft/m] and dry-bulb temperature [°F/°C]
// Reference: see 'GetSeaLevelPressure'
// Notes: this function is just the inverse of 'GetSeaLevelPressure'.
func GetStationPressure(SeaLevelPressure, Altitude, TDryBulb float64) float64 {
	return SeaLevelPressure / GetSeaLevelPressure(1., Altitude, TDryBulb)
}

/******************************************************************************************************
	Functions to set all psychrometric values
******************************************************************************************************/

// Utility function to calculate humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹], dew-point temperature [°F/°C], relative humidity [0-1],
// vapour pressure [psi/Pa], moist air enthalpy[Btu lb⁻¹/J kg⁻¹], moist air volume [ft³ lb⁻¹/m³ kg⁻¹], and degree of saturation [unitless] of air given
// dry-bulb temperature [°F/°C], relative humidity [0-1] and pressure [psi/Pa].
func CalcPsychrometricsFromTWetBulb(TDryBulb, TWetBulb, Pressure float64) (HumRatio, TDewPoint, RelHum, VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation float64) {
	Assert(TWetBulb <= TDryBulb, "Wet bulb temperature is above dry bulb temperature")
	HumRatio = GetHumRatioFromTWetBulb(TDryBulb, TWetBulb, Pressure)
	TDewPoint = GetTDewPointFromHumRatio(TDryBulb, HumRatio, Pressure)
	VapPres = GetVapPresFromHumRatio(HumRatio, Pressure)
	MoistAirEnthalpy = GetMoistAirEnthalpy(TDryBulb, HumRatio)
	MoistAirVolume = GetMoistAirVolume(TDryBulb, HumRatio, Pressure)
	DegreeOfSaturation = GetDegreeOfSaturation(TDryBulb, HumRatio, Pressure)
	return
}

// Utility function to calculate humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹], wet-bulb temperature [°F/°C], relative humidity [0-1],
// vapour pressure [psi/Pa], moist air enthalpy[Btu lb⁻¹/J kg⁻¹], moist air volume [ft³ lb⁻¹/m³ kg⁻¹], and degree of saturation [unitless] of air given
// dry-bulb temperature [°F/°C], relative humidity [0-1] and pressure [psi/Pa].
func CalcPsychrometricsFromTDewPoint(TDryBulb, TDewPoint, Pressure float64) (HumRatio, TWetBulb, RelHum, VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation float64) {
	Assert(TDewPoint <= TDryBulb, "Dew point temperature is above dry bulb temperature")
	HumRatio = GetHumRatioFromTDewPoint(TDewPoint, Pressure)
	TWetBulb = GetTWetBulbFromHumRatio(TDryBulb, HumRatio, Pressure)
	RelHum = GetRelHumFromHumRatio(TDryBulb, HumRatio, Pressure)
	VapPres = GetVapPresFromHumRatio(HumRatio, Pressure)
	MoistAirEnthalpy = GetMoistAirEnthalpy(TDryBulb, HumRatio)
	MoistAirVolume = GetMoistAirVolume(TDryBulb, HumRatio, Pressure)
	DegreeOfSaturation = GetDegreeOfSaturation(TDryBulb, HumRatio, Pressure)
	return
}

// Utility function to calculate humidity ratio [lb_H₂O lb_Air⁻¹/kg_H₂O kg_Air⁻¹], wet-bulb temperature [°F/°C], dew-point temperature [°F/°C],
// vapour pressure [psi/Pa], moist air enthalpy[Btu lb⁻¹/J kg⁻¹], moist air volume [ft³ lb⁻¹/m³ kg⁻¹], and degree of saturation [unitless] of air given
// dry-bulb temperature [°F/°C], relative humidity [0-1] and pressure [psi/Pa].
func CalcPsychrometricsFromRelHum(TDryBulb, RelHum, Pressure float64) (HumRatio, TWetBulb, TDewPoint, VapPres, MoistAirEnthalpy, MoistAirVolume, DegreeOfSaturation float64) {
	Assert(RelHum >= 0 && RelHum <= 1, "Relative humidity is outside range [0,1]")
	HumRatio = GetHumRatioFromRelHum(TDryBulb, RelHum, Pressure)
	TWetBulb = GetTWetBulbFromHumRatio(TDryBulb, HumRatio, Pressure)
	TDewPoint = GetTDewPointFromHumRatio(TDryBulb, HumRatio, Pressure)
	VapPres = GetVapPresFromHumRatio(HumRatio, Pressure)
	MoistAirEnthalpy = GetMoistAirEnthalpy(TDryBulb, HumRatio)
	MoistAirVolume = GetMoistAirVolume(TDryBulb, HumRatio, Pressure)
	DegreeOfSaturation = GetDegreeOfSaturation(TDryBulb, HumRatio, Pressure)
	return
}
