package ca

import (
	"fmt"
	"math"
	"math/rand"

	"github.com/jgcarvalho/zeca-apply/protein"
	"github.com/jgcarvalho/zeca-apply/rule"
)

func Run(p *protein.Protein, ru rule.Rule, steps int) {

	init := (*p).SeqHP
	CAcol := make([][]string, steps+1)
	CAcol[0] = make([]string, len(init))
	helix := make([]float64, len(init))
	strand := make([]float64, len(init))
	coil := make([]float64, len(init))
	other := make([]float64, len(init))
	entropy := make([]float64, len(init))

	copy(CAcol[0], init)
	fmt.Println(CAcol[0])
	var state rule.State
	var strstate string
	for i := 1; i <= steps; i++ {
		CAcol[i] = make([]string, len(init))
		for c := 1; c < (len(init) - 1); c++ {
			CAcol[i][0] = "#"
			CAcol[i][len(init)-1] = "#"
			state = ru[rule.String2State(CAcol[i-1][c-1])][rule.String2State(CAcol[i-1][c])][rule.String2State(CAcol[i-1][c+1])]
			strstate = rule.State2String(state)
			CAcol[i][c] = strstate
			if strstate[0] == '*' {
				helix[c] += 1.0
			} else if strstate[0] == '|' {
				strand[c] += 1.0
			} else if strstate[0] == '_' {
				coil[c] += 1.0
			} else {
				other[c] += 1.0
			}
		}
		fmt.Println(CAcol[i])
	}

	for i := 1; i < (len(init) - 1); i++ {
		helix[i] /= float64(steps)
		strand[i] /= float64(steps)
		coil[i] /= float64(steps)
		other[i] /= float64(steps)
		entropy[i] = -1.0 * (helix[i]*math.Log(helix[i]) + strand[i]*math.Log(strand[i]) + coil[i]*math.Log(coil[i]))
	}
	// fmt.Println(helix)
	// fmt.Println(strand)
	// fmt.Println(coil)
	// fmt.Println(other)
	// fmt.Println(entropy)
	(*p).Helix = helix
	(*p).Strand = strand
	(*p).Coil = coil
	(*p).Other = other
	(*p).Entropy = entropy
	protein.PredictedSS(p)
}

func RunProb(p *protein.Protein, pr rule.ProbRule, steps int) {

	init := (*p).SeqHP
	CAcol := make([][]string, steps+1)
	CAcol[0] = make([]string, len(init))
	helix := make([]float64, len(init))
	strand := make([]float64, len(init))
	coil := make([]float64, len(init))
	other := make([]float64, len(init))
	entropy := make([]float64, len(init))

	copy(CAcol[0], init)
	fmt.Println(CAcol[0])
	var rnd float64
	var state rule.State
	var strstate string
	for i := 1; i <= steps; i++ {
		CAcol[i] = make([]string, len(init))
		for c := 1; c < (len(init) - 1); c++ {
			CAcol[i][0] = "#"
			CAcol[i][len(init)-1] = "#"
			rnd = rand.Float64() / 1.001
			for st, val := range pr[rule.String2State(CAcol[i-1][c-1])][rule.String2State(CAcol[i-1][c])][rule.String2State(CAcol[i-1][c+1])] {
				if val > rnd {
					state = rule.Transition(rule.State(rule.String2State(CAcol[i-1][c])), rule.State(st))
					fmt.Printf("State %s %s %s %s %s %f %f\n", CAcol[i-1][c-1], CAcol[i-1][c], CAcol[i-1][c+1], rule.State2String(rule.State(st)), rule.State2String(state), val, rnd)
				} else {
					rnd -= val
				}
			}
			strstate = rule.State2String(state)
			if strstate == "?" {
				strstate = init[c]
			}
			CAcol[i][c] = strstate
			if strstate[0] == '*' {
				helix[c] += 1.0
			} else if strstate[0] == '|' {
				strand[c] += 1.0
			} else if strstate[0] == '_' {
				coil[c] += 1.0
			} else {
				other[c] += 1.0
			}
		}
		fmt.Println(CAcol[i])
	}

	for i := 1; i < (len(init) - 1); i++ {
		helix[i] /= float64(steps)
		strand[i] /= float64(steps)
		coil[i] /= float64(steps)
		other[i] /= float64(steps)
		entropy[i] = -1.0 * (helix[i]*math.Log(helix[i]) + strand[i]*math.Log(strand[i]) + coil[i]*math.Log(coil[i]))
	}
	// fmt.Println(helix)
	// fmt.Println(strand)
	// fmt.Println(coil)
	// fmt.Println(other)
	// fmt.Println(entropy)
	(*p).Helix = helix
	(*p).Strand = strand
	(*p).Coil = coil
	(*p).Other = other
	(*p).Entropy = entropy
	protein.PredictedSS(p)
}
