package common

import (
	"fmt"
	"os"
)

// CheckIndexFreshness compares modification times of the FASTA and index files.
// If the FASTA file is newer than the index, it returns a warning or error.
func CheckIndexFreshness(fastaFile string, indexFile string) error {
	fastaInfo, err := os.Stat(fastaFile)
	if err != nil {
		return fmt.Errorf("failed to stat FASTA file: %w", err)
	}

	indexInfo, err := os.Stat(indexFile)
	if err != nil {
		return fmt.Errorf("failed to stat index file: %w", err)
	}

	if fastaInfo.ModTime().After(indexInfo.ModTime()) {
		return fmt.Errorf(" %s was modified after %s. Index may be stale — consider regenerating it",
			fastaFile, indexFile)
	}

	return nil
}
