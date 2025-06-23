package fasta3bit

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"flag"
)

var Encode = [256]uint8{
	'A': 0, 'a': 0,
	'T': 1, 't': 1,
	'C': 3, 'c': 3,
	'G': 4, 'g': 4,
	'N': 5, 'n': 5,
	'X': 255,
}

var Decode = map[uint8]byte{
	0: 'A',
	1: 'T',
	3: 'C',
	4: 'G',
	5: 'N',
	255: 'X',
}

type FastaSequence struct {
	Name    string
	Encoded []uint8
}

func EncodeFasta(filename string) ([]FastaSequence, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var sequences []FastaSequence
	var currentName string
	var currentData []uint8

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if strings.HasPrefix(line, ">") {
			if currentName != "" {
				sequences = append(sequences, FastaSequence{currentName, currentData})
			}
			currentName = strings.TrimPrefix(line, ">")
			currentData = []uint8{}
			continue
		}
		for _, base := range line {
			currentData = append(currentData, Encode[base])
		}
	}
	if currentName != "" {
		sequences = append(sequences, FastaSequence{currentName, currentData})
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return sequences, nil
}

func Write3bitFile(sequences []FastaSequence, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	for _, seq := range sequences {
		fmt.Fprintln(file, ">"+seq.Name)
		packed := Pack3bit(seq.Encoded)
		file.Write(packed)
		fmt.Fprintln(file)
	}
	return nil
}

func Read3bitFile(filename string) ([]FastaSequence, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var sequences []FastaSequence
	var currentName string
	var currentPacked []byte

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if currentName != "" {
				decoded := Unpack3bit(currentPacked)
				sequences = append(sequences, FastaSequence{currentName, decoded})
			}
			currentName = strings.TrimPrefix(line, ">")
			currentPacked = []byte{}
		} else {
			currentPacked = append(currentPacked, []byte(line)...) // Assumes no base64
		}
	}
	if currentName != "" {
		decoded := Unpack3bit(currentPacked)
		sequences = append(sequences, FastaSequence{currentName, decoded})
	}
	return sequences, nil
}

func Pack3bit(encoded []uint8) []byte {
	var packed []byte
	var buffer uint32
	var bitsInBuffer uint
	for _, val := range encoded {
		buffer <<= 3
		buffer |= uint32(val)
		bitsInBuffer += 3
		for bitsInBuffer >= 8 {
			bitsInBuffer -= 8
			b := byte(buffer >> bitsInBuffer)
			packed = append(packed, b)
			buffer &= (1 << bitsInBuffer) - 1
		}
	}
	if bitsInBuffer > 0 {
		b := byte(buffer << (8 - bitsInBuffer))
		packed = append(packed, b)
	}
	return packed
}

func Unpack3bit(packed []byte) []uint8 {
	var decoded []uint8
	var buffer uint32
	var bitsInBuffer uint
	for _, b := range packed {
		buffer = (buffer << 8) | uint32(b)
		bitsInBuffer += 8
		for bitsInBuffer >= 3 {
			bitsInBuffer -= 3
			val := uint8(buffer >> bitsInBuffer & 0b111)
			decoded = append(decoded, val)
		}
	}
	return decoded
}

func Run(args []string) {
	fs := flag.NewFlagSet("fasta3bit", flag.ExitOnError)
	encodeFile := fs.String("encode", "", "FASTA file to encode into 3bit")
	decodeFile := fs.String("decode", "", "3bit file to decode to stdout")
	fs.Parse(args)

	if *encodeFile != "" {
		sequences, err := EncodeFasta(*encodeFile)
		if err != nil {
			fmt.Println("Encoding failed:", err)
			return
		}
		outFile := *encodeFile + ".3bit"
		err = Write3bitFile(sequences, outFile)
		if err != nil {
			fmt.Println("Write failed:", err)
			return
		}
		fmt.Println("Encoded and saved to:", outFile)
	} else if *decodeFile != "" {
		sequences, err := Read3bitFile(*decodeFile)
		if err != nil {
			fmt.Println("Decoding failed:", err)
			return
		}
		for _, seq := range sequences {
			fmt.Printf(">%s\n", seq.Name)
			lineCount := 0
			for _, code := range seq.Encoded {
				fmt.Printf("%c", Decode[code])
				lineCount++
				if lineCount == 60 {
					fmt.Println()
					lineCount = 0
				}
			}
			if lineCount != 0 {
				fmt.Println()
			}		
		}
	} else {
		fmt.Println("Usage: -encode <file.fa> or -decode <file.3bit>")
	}
}
