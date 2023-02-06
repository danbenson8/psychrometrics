package psychrometrics

import (
	"math"
	"testing"

	cmp "github.com/google/go-cmp/cmp"
)

const tolerance float64 = .00001

func getCmpOptions() cmp.Option {
	return cmp.Comparer(func(x, y float64) bool {
		diff := math.Abs(x - y)
		mean := math.Abs(x+y) / 2
		if math.IsNaN(diff / mean) {
			return true
		}
		return (diff / mean) < tolerance
	})
}

func TestConversionBetweenTemperatureUnits(t *testing.T) {
	opt := getCmpOptions()
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
