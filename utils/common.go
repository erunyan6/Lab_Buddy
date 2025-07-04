// Common package contains commonly used functions that benefit multiple tools
// Exporting these functions from the Common package reduces redundant code
package common

import (
	"fmt"
	"strings"
	"bufio"
	"os"
	"compress/gzip"
	"io"
)

// ReverseComplement takes a DNA sequence string and returns its reverse complement.
// The function is case-insensitive and warns if the sequence appears to contain a header line (starting with '>').
// Non-standard DNA characters are replaced with the ambiguous base 'N'.
func ReverseComplement(seq string) string {
	var rc strings.Builder
	// Header protection
	if strings.HasPrefix(seq, ">") {
		fmt.Println("Warning: Sequence appears to be a FASTA header. Skipping reverse complement.")
		return seq
	}
	// Case insensativity
	seq = strings.ToUpper(seq)
	for i := len(seq) - 1; i >= 0; i-- {
		switch seq[i] {
		case 'A':
			rc.WriteByte('T')
		case 'T':
			rc.WriteByte('A')
		case 'C':
			rc.WriteByte('G')
		case 'G':
			rc.WriteByte('C')
		default:
			rc.WriteByte('N')	// Ambiguous or invalid character
		}
	}
	return rc.String()
}


type FastaHandler func(id string, seq string, opts map[string]interface{}) error
// StreamFastaWithOpts is a fast, memory-efficient function for streaming FASTA files of any size.
// It automatically detects and decompresses Gzipped files, treats sequences case-insensitively,
// and calls a user-defined handler function for each record.
//
// The handler must follow the FastaHandler signature and can use the 'opts' map to receive
// custom parameters, open output files, counters, filters, etc.
//
// Example handler signature:
//     func(id string, seq string, opts map[string]interface{}) error
func StreamFastaWithOpts(file string, handler FastaHandler, opts map[string]interface{}) error {
	f, err := os.Open(file)
	if err != nil {
		return fmt.Errorf("failed to open file: %w", err)
	}
	defer f.Close()

	var reader io.Reader = f
	buf := make([]byte, 2)
	if _, err := f.Read(buf); err == nil && buf[0] == 0x1F && buf[1] == 0x8B {
		f.Seek(0, io.SeekStart)
		gr, err := gzip.NewReader(f)
		if err != nil {
			return fmt.Errorf("failed to open gzip reader: %w", err)
		}
		defer gr.Close()
		reader = gr
	} else {
		f.Seek(0, io.SeekStart)
	}

	scanner := bufio.NewScanner(reader)

	chunkSize := 0
	stepSize := 0
	
	if val, ok := opts["chunk_size"].(int); ok {
		chunkSize = val
	}
	if val, ok := opts["chunk_overlap"].(int); ok {
		stepSize = chunkSize - val
		if stepSize <= 0 {
			return fmt.Errorf("chunk_overlap must be less than chunk_size")
		}
	}

	var currentID string
	var buffer []byte

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if strings.HasPrefix(line, ">") {
			if currentID != "" && len(buffer) > 0 {
				if chunkSize > 0 {
					if err := streamChunks(currentID, buffer, chunkSize, stepSize, handler, opts); err != nil {
						return err
					}
				} else {
					if err := handler(currentID, string(buffer), opts); err != nil {
						return fmt.Errorf("handler error (%s): %w", currentID, err)
					}
				}
			}
			currentID = strings.TrimPrefix(line, ">")
			buffer = buffer[:0] // reset buffer
		} else {
			buffer = append(buffer, []byte(strings.ToUpper(line))...)
		}
	}
	if currentID != "" && len(buffer) > 0 {
		if chunkSize > 0 {
			if err := streamChunks(currentID, buffer, chunkSize, stepSize, handler, opts); err != nil {
				return err
			}
		} else {
			if err := handler(currentID, string(buffer), opts); err != nil {
				return fmt.Errorf("handler error (%s): %w", currentID, err)
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("scanner error: %w", err)
	}
	return nil
}


func streamChunks(id string, seq []byte, chunkSize, stepSize int, handler FastaHandler, opts map[string]interface{}) error {
	for i := 0; i < len(seq); i += stepSize {
		end := i + chunkSize
		if end > len(seq) {
			end = len(seq)
		}
		chunk := seq[i:end]

		// Copy opts and add chunk position info
		localOpts := make(map[string]interface{})
		for k, v := range opts {
			localOpts[k] = v
		}
		localOpts["chunk_start"] = i
		localOpts["chunk_end"] = end

		err := handler(id, string(chunk), localOpts)
		if err != nil {
			return fmt.Errorf("handler error in chunk %d-%d of %s: %w", i, end, id, err)
		}

		if end == len(seq) {
			break // stop if we've hit the end
		}
	}
	return nil
}
