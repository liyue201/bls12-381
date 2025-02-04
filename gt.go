package bls12381

import (
	"errors"
	"math/big"
)

// E is type for target group element
type E = Fe12

// GT is type for target multiplicative group GT.
type GT struct {
	fp12 *Fp12
}

// Set copies given value into the destination
func (e *E) Set(e2 *E) *E {
	return e.set(e2)
}

// One sets a New target group element to One
func (e *E) One() *E {
	e = new(Fe12).one()
	return e
}

// IsOne returns true if given element equals to One
func (e *E) IsOne() bool {
	return e.isOne()
}

// Equal returns true if given two element is Equal, otherwise returns false
func (g *E) Equal(g2 *E) bool {
	return g.IsEqual(g2)
}

// NewGT constructs New target group instance.
func NewGT() *GT {
	fp12 := NewFp12(nil)
	return &GT{fp12}
}

// Q returns group order in big.Int.
func (g *GT) Q() *big.Int {
	return new(big.Int).Set(qBig)
}

// FromBytes expects 576 byte input and returns target group element
// FromBytes returns error if given element is not on correct subgroup.
func (g *GT) FromBytes(in []byte) (*E, error) {
	e, err := g.fp12.FromBytes(in)
	if err != nil {
		return nil, err
	}
	if !g.IsValid(e) {
		return e, errors.New("invalid element")
	}
	return e, nil
}

// ToBytes serializes target group element.
func (g *GT) ToBytes(e *E) []byte {
	return g.fp12.ToBytes(e)
}

// IsValid checks whether given target group element is in correct subgroup.
func (g *GT) IsValid(e *E) bool {
	r := g.New()
	g.fp12.Exp(r, e, qBig)
	return r.isOne()
}

// New initializes a New target group element which is Equal to One
func (g *GT) New() *E {
	return new(E).One()
}

// Add adds two field element `a` and `b` and assigns the result to the element in first argument.
func (g *GT) Add(c, a, b *E) {
	Fp12Add(c, a, b)
}

// Sub subtracts two field element `a` and `b`, and assigns the result to the element in first argument.
func (g *GT) Sub(c, a, b *E) {
	Fp12Sub(c, a, b)
}

// Mul multiplies two field element `a` and `b` and assigns the result to the element in first argument.
func (g *GT) Mul(c, a, b *E) {
	g.fp12.Mul(c, a, b)
}

// Square squares an element `a` and assigns the result to the element in first argument.
func (g *GT) Square(c, a *E) {
	c.set(a)
	g.fp12.CyclotomicSquare(c)
}

// Exp exponents an element `a` by a scalar `s` and assigns the result to the element in first argument.
func (g *GT) Exp(c, a *E, s *big.Int) {
	g.fp12.CyclotomicExp(c, a, s)
}

// Inverse inverses an element `a` and assigns the result to the element in first argument.
func (g *GT) Inverse(c, a *E) {
	g.fp12.Inverse(c, a)
}
