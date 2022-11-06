package bls12381

import (
	"errors"
	"math/big"
)

type Fp2Temp struct {
	t [3]*Fe
	w *Wfe2
}

type Fp2 struct {
	Fp2Temp
}

func NewFp2Temp() Fp2Temp {
	t := [3]*Fe{}
	for i := 0; i < len(t); i++ {
		t[i] = &Fe{}
	}
	return Fp2Temp{t, &Wfe2{}}
}

func NewFp2() *Fp2 {
	t := NewFp2Temp()
	return &Fp2{t}
}

func (e *Fp2) FromBytes(in []byte) (*Fe2, error) {
	if len(in) != 2*fpByteSize {
		return nil, errors.New("input String must be Equal to 96 Bytes")
	}
	c1, err := FromBytes(in[:fpByteSize])
	if err != nil {
		return nil, err
	}
	c0, err := FromBytes(in[fpByteSize:])
	if err != nil {
		return nil, err
	}
	return &Fe2{*c0, *c1}, nil
}

func (e *Fp2) ToBytes(a *Fe2) []byte {
	out := make([]byte, 2*fpByteSize)
	copy(out[:fpByteSize], ToBytes(&a[1]))
	copy(out[fpByteSize:], ToBytes(&a[0]))
	return out
}

func (e *Fp2) New() *Fe2 {
	return new(Fe2).Zero()
}

func (e *Fp2) Zero() *Fe2 {
	return new(Fe2).Zero()
}

func (e *Fp2) One() *Fe2 {
	return new(Fe2).One()
}

func Fp2Neg(c, a *Fe2) {
	neg(&c[0], &a[0])
	neg(&c[1], &a[1])
}

func Fp2Conjugate(c, a *Fe2) {
	c[0].Set(&a[0])
	neg(&c[1], &a[1])
}

func (e *Fp2) Mul(c, a, b *Fe2) {
	wfp2Mul(e.w, b, a)
	c.FromWide(e.w)
}

func (e *Fp2) MulAssign(a, b *Fe2) {
	wfp2Mul(e.w, b, a)
	a.FromWide(e.w)
}

func (e *Fp2) Square(c, a *Fe2) {
	t := e.t
	// Guide to Pairing Based Cryptography
	// Algorithm 5.16

	ladd(t[0], &a[0], &a[1]) // (a0 + a1)
	sub(t[1], &a[0], &a[1])  // (a0 - a1)
	ldouble(t[2], &a[0])     // 2a0
	mul(&c[0], t[0], t[1])   // c0 = (a0 + a1)(a0 - a1)
	mul(&c[1], t[2], &a[1])  // c1 = 2a0a1
}

func (e *Fp2) SquareAssign(a *Fe2) {
	t := e.t
	ladd(t[0], &a[0], &a[1])
	sub(t[1], &a[0], &a[1])
	ldouble(t[2], &a[0])
	mul(&a[0], t[0], t[1])
	mul(&a[1], t[2], &a[1])
}

func (e *Fp2) Mul0(c, a *Fe2, b *Fe) {
	mul(&c[0], &a[0], b)
	mul(&c[1], &a[1], b)
}

func (e *Fp2) Mul0Assign(a *Fe2, b *Fe) {
	mul(&a[0], &a[0], b)
	mul(&a[1], &a[1], b)
}

func (e *Fp2) MulByB(c, a *Fe2) {
	t := e.t
	// c0 = 4a0 - 4a1
	// c1 = 4a0 + 4a1
	double(t[0], &a[0])
	doubleAssign(t[0])
	double(t[1], &a[1])
	doubleAssign(t[1])
	sub(&c[0], t[0], t[1])
	add(&c[1], t[0], t[1])
}

func (e *Fp2) Inverse(c, a *Fe2) {
	t := e.t
	// Guide to Pairing Based Cryptography
	// Algorithm 5.16

	square(t[0], &a[0])     // a0^2
	square(t[1], &a[1])     // a1^2
	addAssign(t[0], t[1])   // a0^2 + a1^2
	Inverse(t[0], t[0])     // (a0^2 + a1^2)^-1
	mul(&c[0], &a[0], t[0]) // c0 = a0(a0^2 + a1^2)^-1
	mul(t[0], t[0], &a[1])  // a1(a0^2 + a1^2)^-1
	neg(&c[1], t[0])        // c1 = a1(a0^2 + a1^2)^-1
}

func (e *Fp2) InverseBatch(in []Fe2) {

	n, N, setFirst := 0, len(in), false

	for i := 0; i < len(in); i++ {
		if !in[i].IsZero() {
			n++
		}
	}
	if n == 0 {
		return
	}

	tA := make([]Fe2, n)
	tB := make([]Fe2, n)

	// a, ab, abc, abcd, ...
	for i, j := 0, 0; i < N; i++ {
		if !in[i].IsZero() {
			if !setFirst {
				setFirst = true
				tA[j].Set(&in[i])
			} else {
				e.Mul(&tA[j], &in[i], &tA[j-1])
			}
			j = j + 1
		}
	}

	// (abcd...)^-1
	e.Inverse(&tB[n-1], &tA[n-1])

	// a^-1, ab^-1, abc^-1, abcd^-1, ...
	for i, j := N-1, n-1; j != 0; i-- {
		if !in[i].IsZero() {
			e.Mul(&tB[j-1], &tB[j], &in[i])
			j = j - 1
		}
	}

	// a^-1, b^-1, c^-1, d^-1
	for i, j := 0, 0; i < N; i++ {
		if !in[i].IsZero() {
			if setFirst {
				setFirst = false
				in[i].Set(&tB[j])
			} else {
				e.Mul(&in[i], &tA[j-1], &tB[j])
			}
			j = j + 1
		}
	}
}

func (e *Fp2) Exp(c, a *Fe2, s *big.Int) {
	z := e.One()
	for i := s.BitLen() - 1; i >= 0; i-- {
		e.Square(z, z)
		if s.Bit(i) == 1 {
			e.Mul(z, z, a)
		}
	}
	c.Set(z)
}

func (e *Fp2) FrobeniusMap1(a *Fe2) {
	Fp2Conjugate(a, a)
}

func (e *Fp2) FrobeniusMap(a *Fe2, power int) {
	if power&1 == 1 {
		Fp2Conjugate(a, a)
	}
}

func (e *Fp2) Sqrt(c, a *Fe2) bool {
	u, x0, a1, alpha := &Fe2{}, &Fe2{}, &Fe2{}, &Fe2{}
	u.Set(a)
	e.Exp(a1, a, pMinus3Over4)
	e.Square(alpha, a1)
	e.Mul(alpha, alpha, a)
	e.Mul(x0, a1, a)
	if alpha.Equal(negativeOne2) {
		neg(&c[0], &x0[1])
		c[1].Set(&x0[0])
		return true
	}
	fp2Add(alpha, alpha, e.One())
	e.Exp(alpha, alpha, pMinus1Over2)
	e.Mul(c, alpha, x0)
	e.Square(alpha, c)
	return alpha.Equal(u)
}

func (e *Fp2) IsQuadraticNonResidue(a *Fe2) bool {
	c0, c1 := new(Fe), new(Fe)
	square(c0, &a[0])
	square(c1, &a[1])
	add(c1, c1, c0)
	return IsQuadraticNonResidue(c1)
}

// faster Square root algorith is adapted from blst library
// https://github.com/supranational/blst/blob/master/src/sqrt.c

func (e *Fp2) SqrtBLST(out, inp *Fe2) bool {
	aa, bb := new(Fe), new(Fe)
	ret := new(Fe2)
	square(aa, &inp[0])
	square(bb, &inp[1])
	add(aa, aa, bb)
	Sqrt(aa, aa)
	sub(bb, &inp[0], aa)
	add(aa, &inp[0], aa)
	if aa.IsZero() {
		aa.Set(bb)
	}
	mul(aa, aa, twoInv)
	RSqrt(&ret[0], aa)
	ret[1].Set(&inp[1])
	mul(&ret[1], &ret[1], twoInv)
	mul(&ret[1], &ret[1], &ret[0])
	mul(&ret[0], &ret[0], aa)
	return e.SqrtAlignBLST(out, ret, ret, inp)
}

func (e *Fp2) SqrtAlignBLST(out, ret, sqrt, inp *Fe2) bool {

	t0, t1 := new(Fe2), new(Fe2)
	coeff := e.One()
	e.Square(t0, sqrt)

	//
	fp2Sub(t1, t0, inp)
	isSqrt := t1.IsZero()

	//
	fp2Add(t1, t0, inp)
	flag := t1.IsZero()
	if flag {
		coeff.Set(sqrtMinus1)
	}
	isSqrt = flag || isSqrt

	//
	sub(&t1[0], &t0[0], &inp[1])
	add(&t1[1], &t0[1], &inp[0])
	flag = t1.IsZero()
	if flag {
		coeff.Set(sqrtSqrtMinus1)
	}
	isSqrt = flag || isSqrt

	//
	add(&t1[0], &t0[0], &inp[1])
	sub(&t1[1], &t0[1], &inp[0])
	flag = t1.IsZero()
	if flag {

		coeff.Set(sqrtMinusSqrtMinus1)
	}
	isSqrt = flag || isSqrt

	e.Mul(out, coeff, ret)
	return isSqrt
}
