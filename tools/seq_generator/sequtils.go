package seq_generator

import ("strings")

func WrapFasta(seq string, width int) string {
	var out strings.Builder
	for i := 0; i < len(seq); i += width {
		end := i + width
		if end > len(seq) {
			end = len(seq)
		}
		out.WriteString(seq[i:end] + "\n")
	}
	return out.String()
}
