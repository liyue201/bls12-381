package bls12381

import (
	"errors"
	"math/big"
)

type Fp12 struct {
	Fp12temp
	fp6 *Fp6
}

type Fp12temp struct {
	t2  [7]*Fe2
	t6  [4]*Fe6
	wt2 [3]*Wfe2
	wt6 [3]*Wfe6
}

func NewFp12Temp() Fp12temp {
	t2 := [7]*Fe2{}
	t6 := [4]*Fe6{}
	for i := 0; i < len(t2); i++ {
		t2[i] = &Fe2{}
	}
	for i := 0; i < len(t6); i++ {
		t6[i] = &Fe6{}
	}
	wt2 := [3]*Wfe2{}
	for i := 0; i < len(wt2); i++ {
		wt2[i] = &Wfe2{}
	}
	wt6 := [3]*Wfe6{}
	for i := 0; i < len(wt6); i++ {
		wt6[i] = &Wfe6{}
	}
	return Fp12temp{t2, t6, wt2, wt6}
}

func NewFp12(fp6 *Fp6) *Fp12 {
	t := NewFp12Temp()
	if fp6 == nil {
		return &Fp12{t, NewFp6(nil)}
	}
	return &Fp12{t, fp6}
}

func (e *Fp12) Fp2() *Fp2 {
	return e.fp6.fp2
}

func (e *Fp12) FromBytes(in []byte) (*Fe12, error) {
	if len(in) != 576 {
		return nil, errors.New("input String length must be Equal to 576 Bytes")
	}
	fp6 := e.fp6
	c1, err := fp6.FromBytes(in[:6*fpByteSize])
	if err != nil {
		return nil, err
	}
	c0, err := fp6.FromBytes(in[6*fpByteSize:])
	if err != nil {
		return nil, err
	}
	return &Fe12{*c0, *c1}, nil
}

func (e *Fp12) ToBytes(a *Fe12) []byte {
	fp6 := e.fp6
	out := make([]byte, 12*fpByteSize)
	copy(out[:6*fpByteSize], fp6.ToBytes(&a[1]))
	copy(out[6*fpByteSize:], fp6.ToBytes(&a[0]))
	return out
}

func (e *Fp12) New() *Fe12 {
	return new(Fe12)
}

func (e *Fp12) Zero() *Fe12 {
	return new(Fe12)
}

func (e *Fp12) One() *Fe12 {
	return new(Fe12).one()
}

func Fp12Add(c, a, b *Fe12) {
	Fp6Add(&c[0], &a[0], &b[0])
	Fp6Add(&c[1], &a[1], &b[1])
}

func Fp12Double(c, a *Fe12) {
	Fp6Double(&c[0], &a[0])
	Fp6Double(&c[1], &a[1])
}

func Fp12Sub(c, a, b *Fe12) {
	Fp6Sub(&c[0], &a[0], &b[0])
	Fp6Sub(&c[1], &a[1], &b[1])

}

func Fp12Neg(c, a *Fe12) {
	Fp6Neg(&c[0], &a[0])
	Fp6Neg(&c[1], &a[1])
}

func Fp12Conjugate(c, a *Fe12) {
	c[0].Set(&a[0])
	Fp6Neg(&c[1], &a[1])
}

func (e *Fp12) Mul(c, a, b *Fe12) {
	wt, t := e.wt6, e.t6
	e.fp6.Wmul(wt[1], &a[0], &b[0])
	e.fp6.Wmul(wt[2], &a[1], &b[1])
	Fp6Add(t[0], &a[0], &a[1])
	Fp6Add(t[3], &b[0], &b[1])
	e.fp6.Wmul(wt[0], t[0], t[3])
	Wfp6SubAssign(wt[0], wt[1])
	Wfp6SubAssign(wt[0], wt[2])
	c[1].FromWide(wt[0])
	e.fp6.WmulByNonResidueAssign(wt[2])
	Wfp6AddAssign(wt[1], wt[2])
	c[0].FromWide(wt[1])

}

func (e *Fp12) MulAssign(a, b *Fe12) {
	wt, t := e.wt6, e.t6
	e.fp6.Wmul(wt[1], &a[0], &b[0])
	e.fp6.Wmul(wt[2], &a[1], &b[1])
	Fp6Add(t[0], &a[0], &a[1])
	Fp6Add(t[3], &b[0], &b[1])
	e.fp6.Wmul(wt[0], t[0], t[3])
	Wfp6SubAssign(wt[0], wt[1])
	Wfp6SubAssign(wt[0], wt[2])
	a[1].FromWide(wt[0])
	e.fp6.WmulByNonResidueAssign(wt[2])
	Wfp6AddAssign(wt[1], wt[2])
	a[0].FromWide(wt[1])
}

func (e *Fp12) Mul014(a *Fe12, b0, b1, b4 *Fe2) {
	wt, t := e.wt6, e.t6
	e.fp6.Wmul01(wt[0], &a[0], b0, b1)
	e.fp6.Wmul1(wt[1], &a[1], b4)
	fp2LaddAssign(b1, b4)
	Fp6Ladd(t[2], &a[1], &a[0])
	e.fp6.Wmul01(wt[2], t[2], b0, b1)
	Wfp6SubAssign(wt[2], wt[0])
	Wfp6SubAssign(wt[2], wt[1])
	a[1].FromWide(wt[2])
	e.fp6.WmulByNonResidueAssign(wt[1])
	Wfp6AddAssign(wt[0], wt[1])
	a[0].FromWide(wt[0])
}

func (e *Fp12) Square(c, a *Fe12) {
	t := e.t6
	// Multiplication and Squaring on Pairing-Friendly Fields
	// Complex squaring algorithm
	// https://eprint.iacr.org/2006/471

	Fp6Add(t[0], &a[0], &a[1])
	e.fp6.Mul(t[2], &a[0], &a[1])
	e.fp6.MulByNonResidue(t[1], &a[1])
	Fp6AddAssign(t[1], &a[0])
	e.fp6.MulByNonResidue(t[3], t[2])
	e.fp6.Mul(t[0], t[0], t[1])
	Fp6SubAssign(t[0], t[2])
	Fp6Sub(&c[0], t[0], t[3])
	Fp6Double(&c[1], t[2])
}

func (e *Fp12) SquareAssign(a *Fe12) {
	t := e.t6
	// Multiplication and Squaring on Pairing-Friendly Fields
	// Complex squaring algorithm
	// https://eprint.iacr.org/2006/471

	Fp6Add(t[0], &a[0], &a[1])
	e.fp6.Mul(t[2], &a[0], &a[1])
	e.fp6.MulByNonResidue(t[1], &a[1])
	Fp6AddAssign(t[1], &a[0])
	e.fp6.MulByNonResidue(t[3], t[2])
	e.fp6.Mul(t[0], t[0], t[1])
	Fp6SubAssign(t[0], t[2])
	Fp6Sub(&a[0], t[0], t[3])
	Fp6Double(&a[1], t[2])
}

func (e *Fp12) Inverse(c, a *Fe12) {
	// Guide to Pairing Based Cryptography
	// Algorithm 5.16

	t := e.t6
	e.fp6.Square(t[0], &a[0])         // a0^2
	e.fp6.Square(t[1], &a[1])         // a1^2
	e.fp6.MulByNonResidue(t[1], t[1]) // βa1^2
	Fp6SubAssign(t[0], t[1])          // v = (a0^2 - a1^2)
	e.fp6.Inverse(t[1], t[0])         // v = v^-1
	e.fp6.Mul(&c[0], &a[0], t[1])     // c0 = a0v
	e.fp6.MulAssign(t[1], &a[1])      //
	Fp6Neg(&c[1], t[1])               // c1 = -a1v
}

func (e *Fp12) Exp(c, a *Fe12, s *big.Int) {
	z := e.One()
	for i := s.BitLen() - 1; i >= 0; i-- {
		e.Square(z, z)
		if s.Bit(i) == 1 {
			e.Mul(z, z, a)
		}
	}
	c.set(z)
}

func (e *Fp12) CyclotomicExp(c, a *Fe12, s *big.Int) {
	z := e.One()
	for i := s.BitLen() - 1; i >= 0; i-- {
		e.CyclotomicSquare(z)
		if s.Bit(i) == 1 {
			e.Mul(z, z, a)
		}
	}
	c.set(z)
}

func (e *Fp12) CyclotomicSquare(a *Fe12) {
	t := e.t2
	// Guide to Pairing Based Cryptography
	// 5.5.4 Airthmetic in Cyclotomic Groups

	e.Fp4Square(t[3], t[4], &a[0][0], &a[1][1])
	fp2Sub(t[2], t[3], &a[0][0])
	fp2DoubleAssign(t[2])
	fp2Add(&a[0][0], t[2], t[3])
	fp2Add(t[2], t[4], &a[1][1])
	fp2DoubleAssign(t[2])
	fp2Add(&a[1][1], t[2], t[4])
	e.Fp4Square(t[3], t[4], &a[1][0], &a[0][2])
	e.Fp4Square(t[5], t[6], &a[0][1], &a[1][2])
	fp2Sub(t[2], t[3], &a[0][1])
	fp2DoubleAssign(t[2])
	fp2Add(&a[0][1], t[2], t[3])
	fp2Add(t[2], t[4], &a[1][2])
	fp2DoubleAssign(t[2])
	fp2Add(&a[1][2], t[2], t[4])
	mulByNonResidue(t[3], t[6])
	fp2Add(t[2], t[3], &a[1][0])
	fp2DoubleAssign(t[2])
	fp2Add(&a[1][0], t[2], t[3])
	fp2Sub(t[2], t[5], &a[0][2])
	fp2DoubleAssign(t[2])
	fp2Add(&a[0][2], t[2], t[5])
}

func (e *Fp12) Fp4Square(c0, c1, a0, a1 *Fe2) {
	wt, t := e.wt2, e.t2
	// Multiplication and Squaring on Pairing-Friendly Fields
	// Karatsuba squaring algorithm
	// https://eprint.iacr.org/2006/471

	wfp2Square(wt[0], a0)
	wfp2Square(wt[1], a1)
	wfp2MulByNonResidue(wt[2], wt[1])
	wfp2AddAssign(wt[2], wt[0])
	c0.FromWide(wt[2])
	fp2Add(t[0], a0, a1)
	wfp2Square(wt[2], t[0])
	wfp2SubAssign(wt[2], wt[0])
	wfp2SubAssign(wt[2], wt[1])
	c1.FromWide(wt[2])
}

func (e *Fp12) FrobeniusMap1(a *Fe12) {
	fp6, fp2 := e.fp6, e.fp6.fp2
	fp6.FrobeniusMap1(&a[0])
	fp6.FrobeniusMap1(&a[1])
	fp2.MulAssign(&a[1][0], &frobeniusCoeffs12[1])
	fp2.MulAssign(&a[1][1], &frobeniusCoeffs12[1])
	fp2.MulAssign(&a[1][2], &frobeniusCoeffs12[1])
}

func (e *Fp12) FrobeniusMap2(a *Fe12) {
	fp6, fp2 := e.fp6, e.fp6.fp2
	fp6.FrobeniusMap2(&a[0])
	fp6.FrobeniusMap2(&a[1])
	fp2.MulAssign(&a[1][0], &frobeniusCoeffs12[2])
	fp2.MulAssign(&a[1][1], &frobeniusCoeffs12[2])
	fp2.MulAssign(&a[1][2], &frobeniusCoeffs12[2])
}

func (e *Fp12) FrobeniusMap3(a *Fe12) {
	fp6, fp2 := e.fp6, e.fp6.fp2
	fp6.FrobeniusMap3(&a[0])
	fp6.FrobeniusMap3(&a[1])
	fp2.MulAssign(&a[1][0], &frobeniusCoeffs12[3])
	fp2.MulAssign(&a[1][1], &frobeniusCoeffs12[3])
	fp2.MulAssign(&a[1][2], &frobeniusCoeffs12[3])
}
