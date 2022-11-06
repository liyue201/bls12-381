package bls12381

import (
	"errors"
	"math/big"
)

func FromBytes(in []byte) (*Fe, error) {
	fe := &Fe{}
	if len(in) != fpByteSize {
		return nil, errors.New("input String must be Equal 48 Bytes")
	}
	fe.SetBytes(in)
	if !fe.IsValid() {
		return nil, errors.New("must be less than modulus")
	}
	ToMont(fe, fe)
	return fe, nil
}

func From64Bytes(in []byte) (*Fe, error) {
	if len(in) != 32*2 {
		return nil, errors.New("input String must be Equal 64 Bytes")
	}
	a0 := make([]byte, fpByteSize)
	copy(a0[fpByteSize-32:fpByteSize], in[:32])
	a1 := make([]byte, fpByteSize)
	copy(a1[fpByteSize-32:fpByteSize], in[32:])
	e0, err := FromBytes(a0)
	if err != nil {
		return nil, err
	}
	e1, err := FromBytes(a1)
	if err != nil {
		return nil, err
	}
	// F = 2 ^ 256 * R
	F := Fe{
		0x75b3cd7c5ce820f,
		0x3ec6ba621c3edb0b,
		0x168a13d82bff6bce,
		0x87663c4bf8c449d2,
		0x15f34c83ddc8d830,
		0xf9628b49caa2e85,
	}

	mul(e0, e0, &F)
	add(e1, e1, e0)
	return e1, nil
}

func FromBig(in *big.Int) (*Fe, error) {
	fe := new(Fe).SetBig(in)
	if !fe.IsValid() {
		return nil, errors.New("invalid input String")
	}
	ToMont(fe, fe)
	return fe, nil
}

func FromString(in string) (*Fe, error) {
	fe, err := new(Fe).SetString(in)
	if err != nil {
		return nil, err
	}
	if !fe.IsValid() {
		return nil, errors.New("invalid input String")
	}
	ToMont(fe, fe)
	return fe, nil
}

func ToBytes(e *Fe) []byte {
	e2 := new(Fe)
	FromMont(e2, e)
	return e2.Bytes()
}

func ToBig(e *Fe) *big.Int {
	e2 := new(Fe)
	FromMont(e2, e)
	return e2.Big()
}

func ToString(e *Fe) (s string) {
	e2 := new(Fe)
	FromMont(e2, e)
	return e2.String()
}

func ToMont(c, a *Fe) {
	mul(c, a, r2)
}

func FromMont(c, a *Fe) {
	mul(c, a, &Fe{1})
}

func Wfp2MulGeneric(c *Wfe2, a, b *Fe2) {
	wt0, wt1 := new(Wfe), new(Wfe)
	t0, t1 := new(Fe), new(Fe)
	wmul(wt0, &a[0], &b[0])
	wmul(wt1, &a[1], &b[1])
	wsub(&c[0], wt0, wt1)
	lwaddAssign(wt0, wt1)
	ladd(t0, &a[0], &a[1])
	ladd(t1, &b[0], &b[1])
	wmul(wt1, t0, t1)
	lwsub(&c[1], wt1, wt0)
}

func Wfp2SquareGeneric(c *Wfe2, a *Fe2) {
	t0, t1, t2 := new(Fe), new(Fe), new(Fe)
	ladd(t0, &a[0], &a[1])
	sub(t1, &a[0], &a[1])
	ldouble(t2, &a[0])
	wmul(&c[0], t1, t0)
	wmul(&c[1], t2, &a[1])
}

func Exp(c, a *Fe, e *big.Int) {
	z := new(Fe).Set(r1)
	for i := e.BitLen(); i >= 0; i-- {
		mul(z, z, z)
		if e.Bit(i) == 1 {
			mul(z, z, a)
		}
	}
	c.Set(z)
}

func Inverse(inv, e *Fe) {
	if e.IsZero() {
		inv.Zero()
		return
	}
	u := new(Fe).Set(&modulus)
	v := new(Fe).Set(e)
	s := &Fe{1}
	r := &Fe{0}
	var k int
	var z uint64
	var found = false
	// Phase 1
	for i := 0; i < sixWordBitSize*2; i++ {
		if v.IsZero() {
			found = true
			break
		}
		if u.IsEven() {
			u.Div2(0)
			s.Mul2()
		} else if v.IsEven() {
			v.Div2(0)
			z += r.Mul2()
		} else if u.Cmp(v) == 1 {
			lsubAssign(u, v)
			u.Div2(0)
			laddAssign(r, s)
			s.Mul2()
		} else {
			lsubAssign(v, u)
			v.Div2(0)
			laddAssign(s, r)
			z += r.Mul2()
		}
		k += 1
	}

	if !found {
		inv.Zero()
		return
	}

	if k < fpBitSize || k > fpBitSize+sixWordBitSize {
		inv.Zero()
		return
	}

	if r.Cmp(&modulus) != -1 || z > 0 {
		lsubAssign(r, &modulus)
	}
	u.Set(&modulus)
	lsubAssign(u, r)

	// Phase 2
	for i := k; i < 2*sixWordBitSize; i++ {
		double(u, u)
	}
	inv.Set(u)
}

func InverseBatch(in []Fe) {

	n, N, setFirst := 0, len(in), false

	for i := 0; i < len(in); i++ {
		if !in[i].IsZero() {
			n++
		}
	}
	if n == 0 {
		return
	}

	tA := make([]Fe, n)
	tB := make([]Fe, n)

	for i, j := 0, 0; i < N; i++ {
		if !in[i].IsZero() {
			if !setFirst {
				setFirst = true
				tA[j].Set(&in[i])
			} else {
				mul(&tA[j], &in[i], &tA[j-1])
			}
			j = j + 1
		}
	}

	Inverse(&tB[n-1], &tA[n-1])
	for i, j := N-1, n-1; j != 0; i-- {
		if !in[i].IsZero() {
			mul(&tB[j-1], &tB[j], &in[i])
			j = j - 1
		}
	}

	for i, j := 0, 0; i < N; i++ {
		if !in[i].IsZero() {
			if setFirst {
				setFirst = false
				in[i].Set(&tB[j])
			} else {
				mul(&in[i], &tA[j-1], &tB[j])
			}
			j = j + 1
		}
	}
}

func RSqrt(c, a *Fe) bool {
	t0, t1 := new(Fe), new(Fe)
	SqrtAddchain(t0, a)
	mul(t1, t0, a)
	square(t1, t1)
	ret := t1.Equal(a)
	c.Set(t0)
	return ret
}

func Sqrt(c, a *Fe) bool {
	u, v := new(Fe).Set(a), new(Fe)
	// a ^ (p - 3) / 4
	SqrtAddchain(c, a)
	// a ^ (p + 1) / 4
	mul(c, c, u)

	square(v, c)
	return u.Equal(v)
}

func _sqrt(c, a *Fe) bool {
	u, v := new(Fe).Set(a), new(Fe)
	Exp(c, a, pPlus1Over4)
	square(v, c)
	return u.Equal(v)
}

func SqrtAddchain(c, a *Fe) {
	chain := func(c *Fe, n int, a *Fe) {
		for i := 0; i < n; i++ {
			square(c, c)
		}
		mul(c, c, a)
	}

	t := make([]Fe, 16)
	t[13].Set(a)
	square(&t[0], &t[13])
	mul(&t[8], &t[0], &t[13])
	square(&t[4], &t[0])
	mul(&t[1], &t[8], &t[0])
	mul(&t[6], &t[4], &t[8])
	mul(&t[9], &t[1], &t[4])
	mul(&t[12], &t[6], &t[4])
	mul(&t[3], &t[9], &t[4])
	mul(&t[7], &t[12], &t[4])
	mul(&t[15], &t[3], &t[4])
	mul(&t[10], &t[7], &t[4])
	mul(&t[2], &t[15], &t[4])
	mul(&t[11], &t[10], &t[4])
	square(&t[0], &t[3])
	mul(&t[14], &t[11], &t[4])
	mul(&t[5], &t[0], &t[8])
	mul(&t[4], &t[0], &t[1])

	chain(&t[0], 12, &t[15])
	chain(&t[0], 7, &t[7])
	chain(&t[0], 4, &t[1])
	chain(&t[0], 6, &t[6])
	chain(&t[0], 7, &t[11])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 2, &t[8])
	chain(&t[0], 6, &t[3])
	chain(&t[0], 6, &t[3])
	chain(&t[0], 6, &t[9])
	chain(&t[0], 3, &t[8])
	chain(&t[0], 7, &t[3])
	chain(&t[0], 4, &t[3])
	chain(&t[0], 6, &t[7])
	chain(&t[0], 6, &t[14])
	chain(&t[0], 3, &t[13])
	chain(&t[0], 8, &t[3])
	chain(&t[0], 7, &t[11])
	chain(&t[0], 5, &t[12])
	chain(&t[0], 6, &t[3])
	chain(&t[0], 6, &t[5])
	chain(&t[0], 4, &t[9])
	chain(&t[0], 8, &t[5])
	chain(&t[0], 4, &t[3])
	chain(&t[0], 7, &t[11])
	chain(&t[0], 9, &t[10])
	chain(&t[0], 2, &t[8])
	chain(&t[0], 5, &t[6])
	chain(&t[0], 7, &t[1])
	chain(&t[0], 7, &t[9])
	chain(&t[0], 6, &t[11])
	chain(&t[0], 5, &t[5])
	chain(&t[0], 5, &t[10])
	chain(&t[0], 5, &t[10])
	chain(&t[0], 8, &t[3])
	chain(&t[0], 7, &t[2])
	chain(&t[0], 9, &t[7])
	chain(&t[0], 5, &t[3])
	chain(&t[0], 3, &t[8])
	chain(&t[0], 8, &t[7])
	chain(&t[0], 3, &t[8])
	chain(&t[0], 7, &t[9])
	chain(&t[0], 9, &t[7])
	chain(&t[0], 6, &t[2])
	chain(&t[0], 6, &t[4])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 4, &t[3])
	chain(&t[0], 3, &t[8])
	chain(&t[0], 8, &t[2])
	chain(&t[0], 7, &t[4])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 4, &t[7])
	chain(&t[0], 4, &t[6])
	chain(&t[0], 7, &t[4])
	chain(&t[0], 5, &t[5])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 5, &t[4])
	chain(&t[0], 4, &t[3])
	chain(&t[0], 6, &t[2])
	chain(&t[0], 4, &t[1])
	square(c, &t[0])
}

func IsQuadraticNonResidue(a *Fe) bool {
	if a.IsZero() {
		return true
	}
	return !Sqrt(new(Fe), a)
}
