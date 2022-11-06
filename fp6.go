package bls12381

import (
	"errors"
	"math/big"
)

type Fp6Temp struct {
	t  [5]*Fe2
	wt [6]*Wfe2
}

type Fp6 struct {
	fp2 *Fp2
	Fp6Temp
}

func NewFp6Temp() Fp6Temp {
	t := [5]*Fe2{}
	for i := 0; i < len(t); i++ {
		t[i] = &Fe2{}
	}
	wt := [6]*Wfe2{}
	for i := 0; i < len(wt); i++ {
		wt[i] = &Wfe2{}
	}
	return Fp6Temp{t, wt}
}

func NewFp6(f *Fp2) *Fp6 {
	t := NewFp6Temp()
	if f == nil {
		return &Fp6{NewFp2(), t}
	}
	return &Fp6{f, t}
}

func (e *Fp6) FromBytes(b []byte) (*Fe6, error) {
	if len(b) != 288 {
		return nil, errors.New("input String length must be Equal to 288 Bytes")
	}
	fp2 := e.fp2
	u2, err := fp2.FromBytes(b[:2*fpByteSize])
	if err != nil {
		return nil, err
	}
	u1, err := fp2.FromBytes(b[2*fpByteSize : 4*fpByteSize])
	if err != nil {
		return nil, err
	}
	u0, err := fp2.FromBytes(b[4*fpByteSize:])
	if err != nil {
		return nil, err
	}
	return &Fe6{*u0, *u1, *u2}, nil
}

func (e *Fp6) ToBytes(a *Fe6) []byte {
	fp2 := e.fp2
	out := make([]byte, 6*fpByteSize)
	copy(out[:2*fpByteSize], fp2.ToBytes(&a[2]))
	copy(out[2*fpByteSize:4*fpByteSize], fp2.ToBytes(&a[1]))
	copy(out[4*fpByteSize:], fp2.ToBytes(&a[0]))
	return out
}

func (e *Fp6) New() *Fe6 {
	return new(Fe6)
}

func (e *Fp6) Zero() *Fe6 {
	return new(Fe6)
}

func (e *Fp6) One() *Fe6 {
	return new(Fe6).One()
}

func Fp6Ladd(c, a, b *Fe6) {
	fp2Ladd(&c[0], &a[0], &b[0])
	fp2Ladd(&c[1], &a[1], &b[1])
	fp2Ladd(&c[2], &a[2], &b[2])
}

func Wfp6SubAssign(a, b *Wfe6) {
	wfp2SubAssign(&a[0], &b[0])
	wfp2SubAssign(&a[1], &b[1])
	wfp2SubAssign(&a[2], &b[2])
}

func Wfp6AddAssign(a, b *Wfe6) {
	wfp2AddAssign(&a[0], &b[0])
	wfp2AddAssign(&a[1], &b[1])
	wfp2AddAssign(&a[2], &b[2])
}

func Fp6Add(c, a, b *Fe6) {
	fp2Add(&c[0], &a[0], &b[0])
	fp2Add(&c[1], &a[1], &b[1])
	fp2Add(&c[2], &a[2], &b[2])
}

func Fp6AddAssign(a, b *Fe6) {
	fp2AddAssign(&a[0], &b[0])
	fp2AddAssign(&a[1], &b[1])
	fp2AddAssign(&a[2], &b[2])
}

func Fp6Double(c, a *Fe6) {
	fp2Double(&c[0], &a[0])
	fp2Double(&c[1], &a[1])
	fp2Double(&c[2], &a[2])
}

func Fp6DoubleAssign(a *Fe6) {
	fp2DoubleAssign(&a[0])
	fp2DoubleAssign(&a[1])
	fp2DoubleAssign(&a[2])
}

func Fp6Sub(c, a, b *Fe6) {
	fp2Sub(&c[0], &a[0], &b[0])
	fp2Sub(&c[1], &a[1], &b[1])
	fp2Sub(&c[2], &a[2], &b[2])
}

func Fp6SubAssign(a, b *Fe6) {
	fp2SubAssign(&a[0], &b[0])
	fp2SubAssign(&a[1], &b[1])
	fp2SubAssign(&a[2], &b[2])
}

func Fp6Neg(c, a *Fe6) {
	Fp2Neg(&c[0], &a[0])
	Fp2Neg(&c[1], &a[1])
	Fp2Neg(&c[2], &a[2])
}

func (e *Fp6) Wmul01(c *Wfe6, a *Fe6, b0, b1 *Fe2) {
	wt, t := e.wt, e.t
	wfp2Mul(wt[0], &a[0], b0)   // v0 = b0a0
	wfp2Mul(wt[1], &a[1], b1)   // v1 = a1b1
	fp2Ladd(t[2], &a[1], &a[2]) // a1 + a2
	wfp2Mul(wt[2], t[2], b1)    // b1(a1 + a2)
	wfp2SubAssign(wt[2], wt[1]) // b1(a1 + a2) - v1
	wfp2MulByNonResidueAssign(wt[2])
	fp2Ladd(t[3], &a[0], &a[2]) // a0 + a2
	wfp2Mul(wt[3], t[3], b0)    // b0(a0 + a2)
	wfp2SubAssign(wt[3], wt[0])
	wfp2Add(&c[2], wt[3], wt[1])
	fp2Ladd(t[0], b0, b1)       // (b0 + b1)
	fp2Ladd(t[1], &a[0], &a[1]) // (a0 + a1)
	wfp2Mul(wt[4], t[0], t[1])  // (a0 + a1)(b0 + b1)
	wfp2SubAssign(wt[4], wt[0])
	wfp2Sub(&c[1], wt[4], wt[1])
	wfp2Add(&c[0], wt[2], wt[0])
}

func (e *Fp6) Wmul1(c *Wfe6, a *Fe6, b1 *Fe2) {
	wt := e.wt
	wfp2Mul(wt[0], &a[2], b1)
	wfp2Mul(&c[2], &a[1], b1)
	wfp2Mul(&c[1], &a[0], b1)
	wfp2MulByNonResidue(&c[0], wt[0])
}

func (e *Fp6) Wmul(c *Wfe6, a, b *Fe6) {

	wt, t := e.wt, e.t

	// Faster Explicit Formulas for Computing Pairings over Ordinary Curves
	// AKLGL
	// https://eprint.iacr.org/2010/526.pdf
	// Algorithm 3

	// 1. T0 = a0b0,T1 = a1b1, T2 = a2b2
	wfp2Mul(wt[0], &a[0], &b[0])
	wfp2Mul(wt[1], &a[1], &b[1])
	wfp2Mul(wt[2], &a[2], &b[2])
	// 2. t0 = a1 + a2, t1 = b1 + b2
	fp2Ladd(t[0], &a[1], &a[2])
	fp2Ladd(t[1], &b[1], &b[2])
	// 3. T3 = t0 * t1
	wfp2Mul(wt[3], t[0], t[1])
	// 4. T4 = T1 + T2
	wfp2Add(wt[4], wt[1], wt[2])

	// 5,6. T3 = T3 - T4
	wfp2SubMixedAssign(wt[3], wt[4])

	// 7. T4 = β * T3
	wfp2MulByNonResidue(wt[4], wt[3])

	// 8. T5 = T4 + T0
	wfp2Add(wt[5], wt[4], wt[0])

	// 9. t0 = a0 + a1, t1 = b0 + b1
	fp2Ladd(t[0], &a[0], &a[1])
	fp2Ladd(t[1], &b[0], &b[1])

	// 10. T3 = t0 * t1
	wfp2Mul(wt[3], t[0], t[1])

	// 11. T4 = T0 + T1
	wfp2Add(wt[4], wt[0], wt[1])

	// 12,13. T3 = T3 - T4
	wfp2SubMixedAssign(wt[3], wt[4])

	// 14,15. T4 = β * T2
	wfp2MulByNonResidue(wt[4], wt[2])

	// 17. t0 = a0 + a2, t1 = b0 + b2
	fp2Ladd(t[0], &a[0], &a[2])
	fp2Ladd(t[1], &b[0], &b[2])

	// 16. T6 = T3 + T4
	wfp2Add(&c[1], wt[3], wt[4])

	// 18. T3 = t0 * t1
	wfp2Mul(wt[3], t[0], t[1])

	// 19. T4 = T0 + T2
	wfp2Add(wt[4], wt[0], wt[2])

	// 20,21. T3 = T3 - T4
	wfp2SubMixedAssign(wt[3], wt[4])

	// 22,23. T7 = T3 + T1
	wfp2AddMixed(&c[2], wt[3], wt[1])

	// c = T5, T6, T7
	c[0].set(wt[5])
}

func (e *Fp6) Mul(c *Fe6, a, b *Fe6) {
	wt, t := e.wt, e.t

	// 1. T0 = a0b0,T1 = a1b1, T2 = a2b2
	wfp2Mul(wt[0], &a[0], &b[0])
	wfp2Mul(wt[1], &a[1], &b[1])
	wfp2Mul(wt[2], &a[2], &b[2])
	// 2. t0 = a1 + a2, t1 = b1 + b2
	fp2Ladd(t[0], &a[1], &a[2])
	fp2Ladd(t[1], &b[1], &b[2])
	// 3. T3 = t0 * t1
	wfp2Mul(wt[3], t[0], t[1])
	// 4. T4 = T1 + T2
	wfp2Add(wt[4], wt[1], wt[2])

	// 5,6. T3 = T3 - T4
	wfp2SubMixedAssign(wt[3], wt[4])

	// 7. T4 = β * T3
	wfp2MulByNonResidue(wt[4], wt[3])

	// 8. T5 = T4 + T0
	wfp2Add(wt[5], wt[4], wt[0])

	// 9. t0 = a0 + a1, t1 = b0 + b1
	fp2Ladd(t[0], &a[0], &a[1])
	fp2Ladd(t[1], &b[0], &b[1])

	// 10. T3 = t0 * t1
	wfp2Mul(wt[3], t[0], t[1])

	// 11. T4 = T0 + T1
	wfp2Add(wt[4], wt[0], wt[1])

	// 12,13. T3 = T3 - T4
	wfp2SubMixed(wt[3], wt[3], wt[4])

	// 14,15. T4 = β * T2
	wfp2MulByNonResidue(wt[4], wt[2])

	// 17. t0 = a0 + a2, t1 = b0 + b2
	fp2Ladd(t[0], &a[0], &a[2])
	fp2Ladd(t[1], &b[0], &b[2])

	// 16. T6 = T3 + T4
	wfp2Add(wt[3], wt[3], wt[4])
	c[1].FromWide(wt[3])

	// 18. T3 = t0 * t1
	wfp2Mul(wt[3], t[0], t[1])

	// 19. T4 = T0 + T2
	wfp2Add(wt[4], wt[0], wt[2])

	// 20,21. T3 = T3 - T4
	wfp2SubMixed(wt[3], wt[3], wt[4])

	// 22,23. T7 = T3 + T1
	wfp2AddMixed(wt[3], wt[3], wt[1])
	c[2].FromWide(wt[3])

	// c = T5, T6, T7
	c[0].FromWide(wt[5])
}

func (e *Fp6) MulAssign(a, b *Fe6) {
	wt, t := e.wt, e.t

	// Faster Explicit Formulas for Computing Pairings over Ordinary Curves
	// AKLGL
	// https://eprint.iacr.org/2010/526.pdf
	// Algorithm 3

	// 1. T0 = a0b0,T1 = a1b1, T2 = a2b2
	wfp2Mul(wt[0], &a[0], &b[0])
	wfp2Mul(wt[1], &a[1], &b[1])
	wfp2Mul(wt[2], &a[2], &b[2])
	// 2. t0 = a1 + a2, t1 = b1 + b2
	fp2Ladd(t[0], &a[1], &a[2])
	fp2Ladd(t[1], &b[1], &b[2])
	// 3. T3 = t0 * t1
	wfp2Mul(wt[3], t[0], t[1])
	// 4. T4 = T1 + T2
	wfp2Add(wt[4], wt[1], wt[2])

	// 5,6. T3 = T3 - T4
	wfp2SubMixed(wt[3], wt[3], wt[4])

	// 7. T4 = β * T3
	wfp2MulByNonResidue(wt[4], wt[3])

	// 8. T5 = T4 + T0
	wfp2Add(wt[5], wt[4], wt[0])

	// 9. t0 = a0 + a1, t1 = b0 + b1
	fp2Ladd(t[0], &a[0], &a[1])
	fp2Ladd(t[1], &b[0], &b[1])

	// 10. T3 = t0 * t1
	wfp2Mul(wt[3], t[0], t[1])

	// 11. T4 = T0 + T1
	wfp2Add(wt[4], wt[0], wt[1])

	// 12,13. T3 = T3 - T4
	wfp2SubMixed(wt[3], wt[3], wt[4])

	// 14,15. T4 = β * T2
	wfp2MulByNonResidue(wt[4], wt[2])

	// 17. t0 = a0 + a2, t1 = b0 + b2
	fp2Ladd(t[0], &a[0], &a[2])
	fp2Ladd(t[1], &b[0], &b[2])

	// 16. T6 = T3 + T4
	wfp2Add(wt[3], wt[3], wt[4])
	a[1].FromWide(wt[3])

	// 18. T3 = t0 * t1
	wfp2Mul(wt[3], t[0], t[1])

	// 19. T4 = T0 + T2
	wfp2Add(wt[4], wt[0], wt[2])

	// 20,21. T3 = T3 - T4
	wfp2SubMixed(wt[3], wt[3], wt[4])

	// 22,23. T7 = T3 + T1
	wfp2AddMixed(wt[3], wt[3], wt[1])
	a[2].FromWide(wt[3])

	// a = T5, T6, T7
	a[0].FromWide(wt[5])
}

func (e *Fp6) Square(c, a *Fe6) {
	wt, t := e.wt, e.t
	wfp2Square(wt[0], &a[0])
	wfp2Mul(wt[1], &a[0], &a[1])
	wfp2DoubleAssign(wt[1])
	fp2Sub(t[2], &a[0], &a[1])
	fp2AddAssign(t[2], &a[2])
	wfp2Square(wt[2], t[2])
	wfp2Mul(wt[3], &a[1], &a[2])
	wfp2DoubleAssign(wt[3])
	wfp2Square(wt[4], &a[2])
	wfp2MulByNonResidue(wt[5], wt[3])
	wfp2AddAssign(wt[5], wt[0])
	c[0].FromWide(wt[5])
	wfp2MulByNonResidue(wt[5], wt[4])
	wfp2AddAssign(wt[5], wt[1])
	c[1].FromWide(wt[5])
	wfp2AddAssign(wt[1], wt[2])
	wfp2AddAssign(wt[1], wt[3])
	wfp2AddAssign(wt[0], wt[4])
	wfp2SubAssign(wt[1], wt[0])
	c[2].FromWide(wt[1])

}

func (e *Fp6) Wsquare(c *Wfe6, a *Fe6) {
	wt, t := e.wt, e.t
	wfp2Square(wt[0], &a[0])
	wfp2Mul(wt[1], &a[0], &a[1])
	wfp2DoubleAssign(wt[1])
	fp2Sub(t[2], &a[0], &a[1])
	fp2AddAssign(t[2], &a[2])
	wfp2Square(wt[2], t[2])
	wfp2Mul(wt[3], &a[1], &a[2])
	wfp2DoubleAssign(wt[3])
	wfp2Square(wt[4], &a[2])
	wfp2MulByNonResidue(wt[5], wt[3])
	wfp2Add(&c[0], wt[5], wt[0])
	wfp2MulByNonResidue(wt[5], wt[4])
	wfp2Add(&c[1], wt[1], wt[5])
	wfp2AddAssign(wt[1], wt[2])
	wfp2AddAssign(wt[1], wt[3])
	wfp2AddAssign(wt[0], wt[4])
	wfp2Sub(&c[2], wt[1], wt[0])
}

func (e *Fp6) MulByNonResidue(c, a *Fe6) {
	t := e.t
	t[0].Set(&a[0])
	mulByNonResidue(&c[0], &a[2])
	c[2].Set(&a[1])
	c[1].Set(t[0])
}

func (e *Fp6) WmulByNonResidue(c, a *Wfe6) {
	t := e.wt
	t[0].set(&a[0])
	wfp2MulByNonResidue(&c[0], &a[2])
	c[2].set(&a[1])
	c[1].set(t[0])
}

func (e *Fp6) WmulByNonResidueAssign(a *Wfe6) {
	t := e.wt
	t[0].set(&a[0])
	wfp2MulByNonResidue(&a[0], &a[2])
	a[2].set(&a[1])
	a[1].set(t[0])
}

func (e *Fp6) MulByBaseField(c, a *Fe6, b *Fe2) {
	fp2 := e.fp2
	fp2.Mul(&c[0], &a[0], b)
	fp2.Mul(&c[1], &a[1], b)
	fp2.Mul(&c[2], &a[2], b)
}

func (e *Fp6) Exp(c, a *Fe6, s *big.Int) {
	z := e.One()
	for i := s.BitLen() - 1; i >= 0; i-- {
		e.Square(z, z)
		if s.Bit(i) == 1 {
			e.Mul(z, z, a)
		}
	}
	c.Set(z)
}

func (e *Fp6) Inverse(c, a *Fe6) {
	fp2, t := e.fp2, e.t
	fp2.Square(t[0], &a[0])
	fp2.Mul(t[1], &a[1], &a[2])
	mulByNonResidueAssign(t[1])
	fp2SubAssign(t[0], t[1])    // A = v0 - βv5
	fp2.Square(t[1], &a[1])     // v1 = a1^2
	fp2.Mul(t[2], &a[0], &a[2]) // v4 = a0a2
	fp2SubAssign(t[1], t[2])    // C = v1 - v4
	fp2.Square(t[2], &a[2])     // v2 = a2^2
	mulByNonResidueAssign(t[2]) // βv2
	fp2.Mul(t[3], &a[0], &a[1]) // v3 = a0a1
	fp2SubAssign(t[2], t[3])    // B = βv2 - v3
	fp2.Mul(t[3], &a[2], t[2])  // B * a2
	fp2.Mul(t[4], &a[1], t[1])  // C * a1
	fp2AddAssign(t[3], t[4])    // Ca1 + Ba2
	mulByNonResidueAssign(t[3]) // β(Ca1 + Ba2)
	fp2.Mul(t[4], &a[0], t[0])  // Aa0
	fp2AddAssign(t[3], t[4])    // v6 = Aa0 + β(Ca1 + Ba2)
	fp2.Inverse(t[3], t[3])     // F = v6^-1
	fp2.Mul(&c[0], t[0], t[3])  // c0 = AF
	fp2.Mul(&c[1], t[2], t[3])  // c1 = BF
	fp2.Mul(&c[2], t[1], t[3])  // c2 = CF
}

func (e *Fp6) FrobeniusMap(a *Fe6, power int) {
	fp2 := e.fp2
	fp2.FrobeniusMap(&a[0], power)
	fp2.FrobeniusMap(&a[1], power)
	fp2.FrobeniusMap(&a[2], power)
	fp2.MulAssign(&a[1], &frobeniusCoeffs61[power%6])
	fp2.MulAssign(&a[2], &frobeniusCoeffs62[power%6])
}

func (e *Fp6) FrobeniusMap1(a *Fe6) {
	fp2 := e.fp2
	fp2.FrobeniusMap1(&a[0])
	fp2.FrobeniusMap1(&a[1])
	fp2.FrobeniusMap1(&a[2])
	fp2.MulAssign(&a[1], &frobeniusCoeffs61[1])
	fp2.MulAssign(&a[2], &frobeniusCoeffs62[1])
}

func (e *Fp6) FrobeniusMap2(a *Fe6) {
	e.fp2.MulAssign(&a[1], &frobeniusCoeffs61[2])
	e.fp2.MulAssign(&a[2], &frobeniusCoeffs62[2])
}

func (e *Fp6) FrobeniusMap3(a *Fe6) {
	t := e.t
	e.fp2.FrobeniusMap1(&a[0])
	e.fp2.FrobeniusMap1(&a[1])
	e.fp2.FrobeniusMap1(&a[2])
	neg(&t[0][0], &a[1][1])
	a[1][1].Set(&a[1][0])
	a[1][0].Set(&t[0][0])
	Fp2Neg(&a[2], &a[2])
}
