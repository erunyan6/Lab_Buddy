package seq_generator

import (
	"math/rand"
	"strings"
)

// GenerateDNA returns a DNA or RNA sequence
func GenerateDNA(length int, gcBias float64, rna bool) string {
	cWeight := gcBias / 2
	aWeight := (1 - gcBias) / 2
	tWeight := aWeight // AT bias

	seq := make([]rune, length)
	for i := 0; i < length; i++ {
		r := rand.Float64()
		switch {
		case r < aWeight:
			seq[i] = 'A'
		case r < aWeight+tWeight:
			seq[i] = 'T'
		case r < aWeight+tWeight+cWeight:
			seq[i] = 'C'
		default:
			seq[i] = 'G'
		}
	}

	s := string(seq)
	if rna {
		s = strings.ReplaceAll(s, "T", "U")
	}
	return s
}
