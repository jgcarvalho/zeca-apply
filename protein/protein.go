package protein

import (
	"fmt"
	"os"
)

type Protein struct {
	Name    string
	Seq     []string
	SeqHP   []string
	PredSS  []string
	TrueSS  []string
	Helix   []float64
	Strand  []float64
	Coil    []float64
	Other   []float64
	Entropy []float64
}

func PredictedSS(p *Protein) {
	(*p).PredSS = make([]string, len((*p).Seq))
	(*p).PredSS[0], (*p).PredSS[len((*p).Seq)-1] = "#", "#"
	for i := 1; i < len((*p).Seq)-1; i++ {
		if (*p).TrueSS[i] == "H" {
			// if (*p).Helix[i]/0.42 > (*p).Strand[i]/0.25 && (*p).Helix[i]/0.42 > (*p).Coil[i]/0.36 {
			if (*p).Helix[i]+0.000 > (*p).Strand[i] && (*p).Helix[i]+0.000 > (*p).Coil[i] {
				(*p).PredSS[i] = "H"
				// } else if (*p).Strand[i]/0.25 > (*p).Coil[i]/0.36 {
			} else if (*p).Coil[i]+0.000 > (*p).Strand[i] {
				(*p).PredSS[i] = "C"
			} else {
				(*p).PredSS[i] = "E"
			}
		} else if (*p).TrueSS[i] == "E" {
			// if (*p).Helix[i]/0.42 > (*p).Strand[i]/0.25 && (*p).Helix[i]/0.42 > (*p).Coil[i]/0.36 {
			if (*p).Strand[i]+0.000 > (*p).Helix[i] && (*p).Strand[i]+0.000 > (*p).Coil[i] {
				(*p).PredSS[i] = "E"
				// } else if (*p).Strand[i]/0.25 > (*p).Coil[i]/0.36 {
			} else if (*p).Coil[i]+0.000 > (*p).Helix[i] {
				(*p).PredSS[i] = "C"
			} else {
				(*p).PredSS[i] = "H"
			}
		} else if (*p).TrueSS[i] == "C" {
			// if (*p).Helix[i]/0.42 > (*p).Strand[i]/0.25 && (*p).Helix[i]/0.42 > (*p).Coil[i]/0.36 {
			if (*p).Coil[i]+0.000 > (*p).Helix[i] && (*p).Coil[i]+0.000 > (*p).Strand[i] {
				(*p).PredSS[i] = "C"
				// } else if (*p).Strand[i]/0.25 > (*p).Coil[i]/0.36 {
			} else if (*p).Strand[i]+0.000 > (*p).Helix[i] {
				(*p).PredSS[i] = "E"
			} else {
				(*p).PredSS[i] = "H"
			}
		} else {
			if (*p).Coil[i]+0.000 > (*p).Helix[i] && (*p).Coil[i]+0.000 > (*p).Strand[i] {
				(*p).PredSS[i] = "C"
				// } else if (*p).Strand[i]/0.25 > (*p).Coil[i]/0.36 {
			} else if (*p).Helix[i]+0.000 > (*p).Strand[i] {
				(*p).PredSS[i] = "H"
			} else {
				(*p).PredSS[i] = "E"
			}
		}

	}
}

func WriteResData(p Protein, fn string) {
	f, err := os.Create(fn)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	f.WriteString("RES\tSS\tHelix\tStrand\tCoil\tOther\tEntropy\tTRUE-SS\n")
	for i := 1; i < len(p.Seq)-1; i++ {
		f.WriteString(fmt.Sprintf("%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n", p.Seq[i], p.PredSS[i], p.Helix[i], p.Strand[i], p.Coil[i], p.Other[i], p.Entropy[i], p.TrueSS[i]))
	}
}

func WriteStats(p Protein, fn string) {
	f, err := os.Create(fn)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	var Qh, Qe, Qc, Q3 float64
	var th, te, tc float64
	var ph, pe, pc float64
	for i := 1; i < len(p.Seq)-1; i++ {
		switch p.TrueSS[i] {
		case "H":
			th += 1.0
			if p.PredSS[i] == "H" {
				ph += 1.0
			}
		case "E":
			te += 1.0
			if p.PredSS[i] == "E" {
				pe += 1.0
			}
		case "C":
			tc += 1.0
			if p.PredSS[i] == "C" {
				pc += 1.0
			}
		}
	}

	var valid float64
	if th > 0.0 {
		Qh = ph / th * 100.0
		valid += 1.0
		Q3 += Qh
		f.WriteString(fmt.Sprintf("Qh = %.2f\n", Qh))
	}
	if te > 0.0 {
		Qe = pe / te * 100.0
		valid += 1.0
		Q3 += Qe
		f.WriteString(fmt.Sprintf("Qe = %.2f\n", Qe))
	}
	if tc > 0.0 {
		Qc = pc / tc * 100.0
		valid += 1.0
		Q3 += Qc
		f.WriteString(fmt.Sprintf("Qc = %.2f\n", Qc))
	}
	Q3 /= valid

	f.WriteString(fmt.Sprintf("Q3 = %.2f\n", Q3))
}
