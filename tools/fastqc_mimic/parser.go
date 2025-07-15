package fastqc_mimic

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
)

type FastqRecord struct {
	Header   string
	Sequence string
	Plus     string
	Quality  string
}

func OpenFastq(file string) (io.Reader, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}

	buf := make([]byte, 2)
	_, _ = f.Read(buf)
	f.Seek(0, 0)

	if buf[0] == 0x1f && buf[1] == 0x8b {
		gr, err := gzip.NewReader(f)
		if err != nil {
			return nil, err
		}
		return gr, nil
	}
	return f, nil
}

func ParseFastq(file string) ([]FastqRecord, error) {
	reader, err := OpenFastq(file)
	if err != nil {
		return nil, err
	}

	scanner := bufio.NewScanner(reader)
	var records []FastqRecord

	for scanner.Scan() {
		header := scanner.Text()
		if !scanner.Scan() {
			break
		}
		seq := scanner.Text()
		if !scanner.Scan() {
			break
		}
		plus := scanner.Text()
		if !scanner.Scan() {
			break
		}
		qual := scanner.Text()

		records = append(records, FastqRecord{
			Header:   header,
			Sequence: seq,
			Plus:     plus,
			Quality:  qual,
		})
	}
	return records, scanner.Err()
}
