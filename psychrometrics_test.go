package psychrometrics

import (
	"fmt"
	"math"
	"testing"

	cmp "github.com/google/go-cmp/cmp"
)

/************************************************************************************************
	Setup, utility function, teardown, etc
************************************************************************************************/

func getCmpOptions(tolerance float64) cmp.Option {
	return cmp.Comparer(func(x, y float64) bool {
		diff := math.Abs(x - y)
		mean := math.Abs(x+y) / 2
		if math.IsNaN(diff / mean) {
			return true
		}
		return (diff / mean) < tolerance
	})
}

/************************************************************************************************
	Test conversion between temperature units
************************************************************************************************/

func TestConversionBetweenTemperatureUnits(t *testing.T) {
	opt := getCmpOptions(0.00001)
	t.Run("70°F is 529.67°R", func(t *testing.T) {
		got := GetTRankineFromTFahrenheit(70.0)
		want := 529.67
		if !cmp.Equal(got, want, opt) {
			t.Fatalf("got %v, wanted %v", got, want)
		}
	})
	t.Run("529.67°R is 70°F", func(t *testing.T) {
		got := GetTFahrenheitFromTRankine(529.67)
		want := 70.0
		if !cmp.Equal(got, want, opt) {
			t.Fatalf("got %v, wanted %v", got, want)
		}
	})
	t.Run("20°C is 293.15K", func(t *testing.T) {
		got := GetTKelvinFromTCelsius(20.0)
		want := 293.15
		if !cmp.Equal(got, want, opt) {
			t.Fatalf("got %v, wanted %v", got, want)
		}
	})
	t.Run("293.15K is 20°C", func(t *testing.T) {
		got := GetTCelsiusFromTKelvin(293.15)
		want := 20.0
		if !cmp.Equal(got, want, opt) {
			t.Fatalf("got %v, wanted %v", got, want)
		}
	})

}

/************************************************************************************************
	Test saturation vapour pressure calculation
	-----------------------------------------------------------------------------------------
	The values are tested against the values published in Table 3 of ch. 1 of the 2017 ASHRAE
	Handbook - Fundamentals over the range [-148, +392]°F & over the range [-100, +200]°C
	ASHRAE's assertion is that the formula is within 300 ppm of the true values, which is true
	except for the value at -76°F/-60°C.
************************************************************************************************/

func TestSatVapPres(t *testing.T) {
	// Accurate to within 300ppm
	opt := getCmpOptions(0.0003)
	type SatVapCase struct {
		units    UnitSystem
		TDryBulb float64
		want     float64
	}
	cases_IP := []SatVapCase{
		{IP, -4, 0.014974},
		{IP, 23, 0.058268},
		{IP, 41, 0.12656},
		{IP, 77, 0.45973},
		{IP, 122, 1.79140},
		{IP, 212, 14.7094},
		{IP, 300, 67.0206},
	}
	cases_SI := []SatVapCase{
		{SI, -20, 103.24},
		{SI, -5, 401.74},
		{SI, 5, 872.6},
		{SI, 25, 3169.7},
		{SI, 50, 12351.3},
		{SI, 100, 101418.0},
		{SI, 150, 476101.4},
	}

	SetUnitSystem(IP)
	for _, tt := range cases_IP {
		t.Run(fmt.Sprintf("saturation vapour pressure at %f°F is %fpsi", tt.TDryBulb, tt.want),
			func(t *testing.T) {
				got := GetSatVapPres(tt.TDryBulb)
				if !cmp.Equal(got, tt.want, opt) {
					t.Fatalf("got %v, wanted %v", got, tt.want)
				}
			})
	}

	SetUnitSystem(SI)
	for _, tt := range cases_SI {
		t.Run(fmt.Sprintf("saturation vapour pressure at %f °C is %f Pa", tt.TDryBulb, tt.want),
			func(t *testing.T) {
				got := GetSatVapPres(tt.TDryBulb)
				if !cmp.Equal(got, tt.want, opt) {
					t.Fatalf("got %v, wanted %v", got, tt.want)
				}
			})
	}

	// -60°C case, accurate to within 10000ppm
	opt = getCmpOptions(0.01)
	TDryBulb := -60.0
	want := 1.08
	SetUnitSystem(SI)
	t.Run(fmt.Sprintf("saturation vapour pressure at %f °F is %f psi", TDryBulb, want),
		func(t *testing.T) {
			got := GetSatVapPres(TDryBulb)
			if !cmp.Equal(got, want, opt) {
				t.Fatalf("got %v, wanted %v", got, want)
			}
		})

	// -76°F case, accurate to within 100ppm
	opt = getCmpOptions(0.001)
	TDryBulb = -76.0
	want = 0.000157
	SetUnitSystem(IP)
	t.Run(fmt.Sprintf("saturation vapour pressure at %f °F is %f psi", TDryBulb, want),
		func(t *testing.T) {
			got := GetSatVapPres(TDryBulb)
			if !cmp.Equal(got, want, opt) {
				t.Fatalf("got %v, wanted %v", got, want)
			}
		})

}
