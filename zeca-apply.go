package main

import (
	"fmt"
	"os"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

func LoadProteinFasta(fn string) {
	protein, _ := getFasta(fn)
	fmt.Println(strings.Split(protein, ""))
}

func LoadSSFasta(fn string) {

}

func getFasta(fn string) (string, error) {
	fasta_file, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	var s []alphabet.Letter
	t := linear.NewSeq("", s, alphabet.Protein)
	reader := fasta.NewReader(fasta_file, t)
	seq, _ := reader.Read()
	//fmt.Println("Read -> ", seq.Alphabet())
	return seq.(*linear.Seq).String(), nil
}

func main() {
	LoadProteinFasta("/home/jgcarvalho/sync/data/seqdb/chameleonic/2LHC.fa")
	// protein fasta
	// rule
	// N evolution steps
}
