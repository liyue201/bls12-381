package bls12381

import (
	"crypto/rand"
	"encoding/hex"
	"fmt"
	"io"
	"math/big"
)

// Fe is base field element representation
type Fe /***			***/ [fpNumberOfLimbs]uint64

// Fe2 is element representation of 'Fp2' which is quadratic extention of base field 'fp'
// Representation follows c[0] + c[1] * u encoding order.
type Fe2 /**			***/ [2]Fe

// Fe6 is element representation of 'Fp6' field which is cubic extention of 'Fp2'
// Representation follows c[0] + c[1] * v + c[2] * v^2 encoding order.
type Fe6 /**			***/ [3]Fe2

// Fe12 is element representation of 'Fp12' field which is quadratic extention of 'Fp6'
// Representation follows c[0] + c[1] * w encoding order.
type Fe12 /**			***/ [2]Fe6

type Wfe /***			***/ [fpNumberOfLimbs * 2]uint64
type Wfe2 /**			***/ [2]Wfe
type Wfe6 /**			***/ [3]Wfe2

func (fe *Fe) SetBytes(in []byte) *Fe {
	l := len(in)
	if l >= fpByteSize {
		l = fpByteSize
	}
	padded := make([]byte, fpByteSize)
	copy(padded[fpByteSize-l:], in[:])
	var a int
	for i := 0; i < fpNumberOfLimbs; i++ {
		a = fpByteSize - i*8
		fe[i] = uint64(padded[a-1]) | uint64(padded[a-2])<<8 |
			uint64(padded[a-3])<<16 | uint64(padded[a-4])<<24 |
			uint64(padded[a-5])<<32 | uint64(padded[a-6])<<40 |
			uint64(padded[a-7])<<48 | uint64(padded[a-8])<<56
	}
	return fe
}

func (fe *Fe) SetBig(a *big.Int) *Fe {
	return fe.SetBytes(a.Bytes())
}

func (fe *Fe) SetString(s string) (*Fe, error) {
	if s[:2] == "0x" {
		s = s[2:]
	}
	bytes, err := hex.DecodeString(s)
	if err != nil {
		return nil, err
	}
	return fe.SetBytes(bytes), nil
}

func (fe *Fe) Set(fe2 *Fe) *Fe {
	fe[0] = fe2[0]
	fe[1] = fe2[1]
	fe[2] = fe2[2]
	fe[3] = fe2[3]
	fe[4] = fe2[4]
	fe[5] = fe2[5]
	return fe
}

func (fe *Fe) Bytes() []byte {
	out := make([]byte, fpByteSize)
	var a int
	for i := 0; i < fpNumberOfLimbs; i++ {
		a = fpByteSize - i*8
		out[a-1] = byte(fe[i])
		out[a-2] = byte(fe[i] >> 8)
		out[a-3] = byte(fe[i] >> 16)
		out[a-4] = byte(fe[i] >> 24)
		out[a-5] = byte(fe[i] >> 32)
		out[a-6] = byte(fe[i] >> 40)
		out[a-7] = byte(fe[i] >> 48)
		out[a-8] = byte(fe[i] >> 56)
	}
	return out
}

func (fe *Fe) Big() *big.Int {
	return new(big.Int).SetBytes(fe.Bytes())
}

func (fe *Fe) String() (s string) {
	for i := fpNumberOfLimbs - 1; i >= 0; i-- {
		s = fmt.Sprintf("%s%16.16x", s, fe[i])
	}
	return "0x" + s
}

func (fe *Fe) Zero() *Fe {
	fe[0] = 0
	fe[1] = 0
	fe[2] = 0
	fe[3] = 0
	fe[4] = 0
	fe[5] = 0
	return fe
}

func (fe *Fe) One() *Fe {
	return fe.Set(r1)
}

func (fe *Fe) Rand(r io.Reader) (*Fe, error) {
	bi, err := rand.Int(r, modulus.Big())
	if err != nil {
		return nil, err
	}
	return fe.SetBig(bi), nil
}

func (fe *Fe) IsValid() bool {
	return fe.Cmp(&modulus) == -1
}

func (fe *Fe) IsOdd() bool {
	var mask uint64 = 1
	return fe[0]&mask != 0
}

func (fe *Fe) IsEven() bool {
	var mask uint64 = 1
	return fe[0]&mask == 0
}

func (fe *Fe) IsZero() bool {
	return (fe[5] | fe[4] | fe[3] | fe[2] | fe[1] | fe[0]) == 0
}

func (fe *Fe) IsOne() bool {
	return fe.Equal(r1)
}

func (fe *Fe) Cmp(fe2 *Fe) int {
	for i := fpNumberOfLimbs - 1; i >= 0; i-- {
		if fe[i] > fe2[i] {
			return 1
		} else if fe[i] < fe2[i] {
			return -1
		}
	}
	return 0
}

func (fe *Fe) Equal(fe2 *Fe) bool {
	return fe2[0] == fe[0] && fe2[1] == fe[1] && fe2[2] == fe[2] && fe2[3] == fe[3] && fe2[4] == fe[4] && fe2[5] == fe[5]
}

func (e *Fe) SignBE() bool {
	negZ, z := new(Fe), new(Fe)
	FromMont(z, e)
	neg(negZ, z)
	return negZ.Cmp(z) > -1
}

func (e *Fe) Sign() bool {
	r := new(Fe)
	FromMont(r, e)
	return r[0]&1 == 0
}

func (e *Fe) Div2(u uint64) {
	e[0] = e[0]>>1 | e[1]<<63
	e[1] = e[1]>>1 | e[2]<<63
	e[2] = e[2]>>1 | e[3]<<63
	e[3] = e[3]>>1 | e[4]<<63
	e[4] = e[4]>>1 | e[5]<<63
	e[5] = e[5]>>1 | u<<63
}

func (e *Fe) Mul2() uint64 {
	u := e[5] >> 63
	e[5] = e[5]<<1 | e[4]>>63
	e[4] = e[4]<<1 | e[3]>>63
	e[3] = e[3]<<1 | e[2]>>63
	e[2] = e[2]<<1 | e[1]>>63
	e[1] = e[1]<<1 | e[0]>>63
	e[0] = e[0] << 1
	return u
}

func (e *Fe2) Zero() *Fe2 {
	e[0].Zero()
	e[1].Zero()
	return e
}

func (e *Fe2) One() *Fe2 {
	e[0].One()
	e[1].Zero()
	return e
}

func (e *Fe2) Set(e2 *Fe2) *Fe2 {
	e[0].Set(&e2[0])
	e[1].Set(&e2[1])
	return e
}

func (e *Fe2) FromMont(a *Fe2) {
	FromMont(&e[0], &a[0])
	FromMont(&e[1], &a[1])
}

func (e *Fe2) FromWide(w *Wfe2) {
	fromWide(&e[0], &w[0])
	fromWide(&e[1], &w[1])
}

func (e *Fe2) Rand(r io.Reader) (*Fe2, error) {
	a0, err := new(Fe).Rand(r)
	if err != nil {
		return nil, err
	}
	e[0].Set(a0)
	a1, err := new(Fe).Rand(r)
	if err != nil {
		return nil, err
	}
	e[1].Set(a1)
	return e, nil
}

func (e *Fe2) IsOne() bool {
	return e[0].IsOne() && e[1].IsZero()
}

func (e *Fe2) IsZero() bool {
	return e[0].IsZero() && e[1].IsZero()
}

func (e *Fe2) Equal(e2 *Fe2) bool {
	return e[0].Equal(&e2[0]) && e[1].Equal(&e2[1])
}

func (e *Fe2) SignBE() bool {
	if !e[1].IsZero() {
		return e[1].SignBE()
	}
	return e[0].SignBE()
}

func (e *Fe2) Sign() bool {
	r := new(Fe)
	if !e[0].IsZero() {
		FromMont(r, &e[0])
		return r[0]&1 == 0
	}
	FromMont(r, &e[1])
	return r[0]&1 == 0
}

func (e *Fe6) Zero() *Fe6 {
	e[0].Zero()
	e[1].Zero()
	e[2].Zero()
	return e
}

func (e *Fe6) One() *Fe6 {
	e[0].One()
	e[1].Zero()
	e[2].Zero()
	return e
}

func (e *Fe6) Set(e2 *Fe6) *Fe6 {
	e[0].Set(&e2[0])
	e[1].Set(&e2[1])
	e[2].Set(&e2[2])
	return e
}

func (e *Fe6) FromMont(a *Fe6) {
	e[0].FromMont(&a[0])
	e[1].FromMont(&a[1])
	e[2].FromMont(&a[2])
}

func (e *Fe6) FromWide(w *Wfe6) {
	e[0].FromWide(&w[0])
	e[1].FromWide(&w[1])
	e[2].FromWide(&w[2])
}

func (e *Fe6) Rand(r io.Reader) (*Fe6, error) {
	a0, err := new(Fe2).Rand(r)
	if err != nil {
		return nil, err
	}
	e[0].Set(a0)
	a1, err := new(Fe2).Rand(r)
	if err != nil {
		return nil, err
	}
	e[1].Set(a1)
	a2, err := new(Fe2).Rand(r)
	if err != nil {
		return nil, err
	}
	e[2].Set(a2)
	return e, nil
}

func (e *Fe6) IsOne() bool {
	return e[0].IsOne() && e[1].IsZero() && e[2].IsZero()
}

func (e *Fe6) IsZero() bool {
	return e[0].IsZero() && e[1].IsZero() && e[2].IsZero()
}

func (e *Fe6) Equal(e2 *Fe6) bool {
	return e[0].Equal(&e2[0]) && e[1].Equal(&e2[1]) && e[2].Equal(&e2[2])
}

func (e *Fe12) Zero() *Fe12 {
	e[0].Zero()
	e[1].Zero()
	return e
}

func (e *Fe12) one() *Fe12 {
	e[0].One()
	e[1].Zero()
	return e
}

func (e *Fe12) set(e2 *Fe12) *Fe12 {
	e[0].Set(&e2[0])
	e[1].Set(&e2[1])
	return e
}

func (e *Fe12) FromMont(a *Fe12) {
	e[0].FromMont(&a[0])
	e[1].FromMont(&a[1])
}

func (e *Fe12) Rand(r io.Reader) (*Fe12, error) {
	a0, err := new(Fe6).Rand(r)
	if err != nil {
		return nil, err
	}
	e[0].Set(a0)
	a1, err := new(Fe6).Rand(r)
	if err != nil {
		return nil, err
	}
	e[1].Set(a1)
	return e, nil
}

func (e *Fe12) isOne() bool {
	return e[0].IsOne() && e[1].IsZero()
}

func (e *Fe12) isZero() bool {
	return e[0].IsZero() && e[1].IsZero()
}

func (e *Fe12) IsEqual(e2 *Fe12) bool {
	return e[0].Equal(&e2[0]) && e[1].Equal(&e2[1])
}

func (fe *Wfe) set(fe2 *Wfe) *Wfe {
	fe[0] = fe2[0]
	fe[1] = fe2[1]
	fe[2] = fe2[2]
	fe[3] = fe2[3]
	fe[4] = fe2[4]
	fe[5] = fe2[5]
	fe[6] = fe2[6]
	fe[7] = fe2[7]
	fe[8] = fe2[8]
	fe[9] = fe2[9]
	fe[10] = fe2[10]
	fe[11] = fe2[11]
	return fe
}

func (fe *Wfe2) set(fe2 *Wfe2) *Wfe2 {
	fe[0].set(&fe2[0])
	fe[1].set(&fe2[1])
	return fe
}

func (fe *Wfe6) set(fe2 *Wfe6) *Wfe6 {
	fe[0].set(&fe2[0])
	fe[1].set(&fe2[1])
	fe[2].set(&fe2[2])
	return fe
}
