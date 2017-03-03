package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/jgcarvalho/zeca-apply/ca"
	"github.com/jgcarvalho/zeca-apply/protein"
	"github.com/jgcarvalho/zeca-apply/rule"
)

// rose: polar, nonpolar
var rose = map[string]string{
	"A": "n", "C": "n", "V": "n", "I": "n", "L": "n", "M": "n", "F": "n",
	"W": "n", "G": "p", "S": "p", "T": "p", "H": "p", "Y": "p", "P": "p",
	"D": "p", "N": "p", "E": "p", "Q": "p", "K": "p", "R": "p"}

// rose special: polar, nonpolar, Gly and Pro
var roseSpecial = map[string]string{
	"A": "n", "C": "n", "V": "n", "I": "n", "L": "n", "M": "n", "F": "n",
	"W": "n", "G": "G", "S": "p", "T": "p", "H": "p", "Y": "p", "P": "P",
	"D": "p", "N": "p", "E": "p", "Q": "p", "K": "p", "R": "p"}

//rose special charged: polar, nonpolar, Gly, Pro, positives and negatives
var roseSpecialCharged = map[string]string{
	"A": "n", "C": "n", "V": "n", "I": "n", "L": "n", "M": "n", "F": "n",
	"W": "n", "G": "G", "S": "p", "T": "p", "H": "p", "Y": "p", "P": "P",
	"D": "-", "N": "p", "E": "-", "Q": "p", "K": "+", "R": "+"}

// Load protein fasta file and return a slice of residues
func LoadProteinFasta(fn string) []string {
	protein, _ := getFasta(fn)
	pr := []string{"#"}
	pr = append(pr, strings.Split(protein, "")...)
	pr = append(pr, "#")
	// fmt.Println(pr)
	return pr
}

// Apply Hydrophobicity to residues
func ApplyHP(prot []string, hp string) []string {
	protHP := make([]string, len(prot))
	switch hp {
	case "rose":
		fmt.Println("Using Rose Hydrophobicity")
		for i := range prot {
			protHP[i] = prot[i] + rose[prot[i]]
		}
	case "rose_special":
		fmt.Println("Using Rose Hydrophobicity + GLY + PRO")
		for i := range prot {
			protHP[i] = prot[i] + roseSpecial[prot[i]]
		}
	case "rose_special_charged":
		fmt.Println("Using Rose Hydrophobicity + GLY + PRO + Charged residues")
		for i := range prot {
			protHP[i] = prot[i] + roseSpecialCharged[prot[i]]
		}
	default:
		copy(protHP, prot)
		fmt.Println("Hydrophobicity not used")
	}
	return protHP
}

func LoadProb(fn string) {

}

func LoadSSFasta(fn string) []string {
	ss, _ := getFasta(fn)
	pr := []string{"#"}
	pr = append(pr, strings.Split(ss, "")...)
	pr = append(pr, "#")
	// fmt.Println(pr)
	return pr
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

func stats() {}

func main() {
	fasta := flag.String("fasta", "", "Protein fasta file")
	hp := flag.String("hp", "none", "Hydrophobicity (rose, rose_special, rose_special_charged, none). Default: none")
	fnRule := flag.String("rule", "", "CA rule")
	fnProb := flag.String("prob", "", "CA probabilistic rule")
	steps := flag.Int("steps", 100, "Number of CA steps. Default: 100")
	fTrueSS := flag.String("trueSS", "", "Real secondary structure")
	fOutRes := flag.String("outres", "", "Output data per residue")
	fOutStats := flag.String("outstats", "", "Output statistics (Q3)")
	flag.Parse()

	// check and load protein fasta file
	var prot []string
	if *fasta != "" {
		prot = LoadProteinFasta(*fasta)
	} else {
		fmt.Println("Protein fasta file is null")
		return
	}

	// hp to protein
	var p protein.Protein
	p.Seq = prot
	p.SeqHP = ApplyHP(prot, *hp)

	// load true SS
	p.TrueSS = LoadSSFasta(*fTrueSS)

	// check and load rule file
	var ru rule.Rule
	var pr rule.ProbRule
	if *fnRule != "" && *fnProb != "" {
		fmt.Println("Choose between deterministic or probabilistic rule, not both")
		return
	} else if *fnRule != "" {
		ru = rule.LoadRule(*fnRule)
	} else if *fnProb != "" {
		pr = rule.LoadProbRule(*fnProb)
	} else {
		fmt.Println("Rule (rule and prob) file is null")
		return
	}

	// check output per residue file
	if *fOutRes == "" {
		fmt.Println("Output per residue not set (file name)")
		return
	}

	// check output file statistics
	if *fOutStats == "" {
		fmt.Println("Output statistics not set (file name)")
		return
	}

	if *fnRule != "" {
		ca.Run(&p, ru, *steps)
	} else if *fnProb != "" {
		ca.RunProb(&p, pr, *steps)
	}

	// fmt.Println(prot)
	// fmt.Println(trueSS)
	// fmt.Println(helix, strand, coil, other, entropy)
	protein.WriteResData(p, *fOutRes)
	protein.WriteStats(p, *fOutStats)
	// protein fasta
	// rule or prob
	// N evolution steps
}
