package seq_generator

import (
	"math/rand"
)

// 20 standard amino acids
var aminoAcids = []rune("ACDEFGHIKLMNPQRSTVWY")

func GenerateProtein(length int) string {
	if length < 2 {
		return "M*" // minimal valid peptide
	}
	seq := make([]rune, length)
	seq[0] = 'M' // start codon (methionine)
	for i := 1; i < length-1; i++ {
		seq[i] = aminoAcids[rand.Intn(len(aminoAcids))]
	}
	seq[length-1] = '*' // stop codon
	return string(seq)
}
