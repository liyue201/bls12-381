package bls12381

type pair struct {
	g1 *PointG1
	g2 *PointG2
}

func newPair(g1 *PointG1, g2 *PointG2) pair {
	return pair{g1, g2}
}

// Engine is BLS12-381 elliptic curve pairing engine
type Engine struct {
	G1   *G1
	G2   *G2
	fp12 *Fp12
	fp2  *Fp2
	pairingEngineTemp
	pairs []pair
}

// NewEngine creates New pairing engine insteace.
func NewEngine() *Engine {
	fp2 := NewFp2()
	fp6 := NewFp6(fp2)
	fp12 := NewFp12(fp6)
	g1 := NewG1()
	g2 := newG2(fp2)
	return &Engine{
		fp2:               fp2,
		fp12:              fp12,
		G1:                g1,
		G2:                g2,
		pairingEngineTemp: newEngineTemp(),
	}
}

type pairingEngineTemp struct {
	t2  [9]*Fe2
	t12 [3]Fe12
}

func newEngineTemp() pairingEngineTemp {
	t2 := [9]*Fe2{}
	for i := 0; i < len(t2); i++ {
		t2[i] = &Fe2{}
	}
	t12 := [3]Fe12{}
	return pairingEngineTemp{t2, t12}
}

// AddPair adds a g1, g2 point pair to pairing engine
func (e *Engine) AddPair(g1 *PointG1, g2 *PointG2) *Engine {
	p := newPair(g1, g2)
	if !(e.G1.IsZero(p.g1) || e.G2.IsZero(p.g2)) {
		e.G1.Affine(p.g1)
		e.G2.Affine(p.g2)
		e.pairs = append(e.pairs, p)
	}
	return e
}

// AddPairInv adds a G1, G2 point pair to pairing engine. G1 point is negated.
func (e *Engine) AddPairInv(g1 *PointG1, g2 *PointG2) *Engine {
	ng1 := e.G1.New().Set(g1)
	e.G1.Neg(ng1, g1)
	e.AddPair(ng1, g2)
	return e
}

// Reset deletes added pairs.
func (e *Engine) Reset() *Engine {
	e.pairs = []pair{}
	return e
}

func (e *Engine) double(f *Fe12, r *PointG2, k int) {
	fp2, t := e.fp2, e.t2

	fp2.Mul(t[0], &r[0], &r[1])
	fp2.Mul0(t[0], t[0], twoInv)
	fp2.Square(t[1], &r[1])
	fp2.Square(t[2], &r[2])
	fp2Double(t[7], t[2])
	fp2AddAssign(t[7], t[2])
	fp2.MulByB(t[3], t[7])
	fp2Double(t[4], t[3])
	fp2AddAssign(t[4], t[3])
	fp2Add(t[5], t[1], t[4])
	fp2.Mul0(t[5], t[5], twoInv)
	fp2Add(t[6], &r[1], &r[2])
	fp2.SquareAssign(t[6])
	fp2Add(t[7], t[2], t[1])
	fp2SubAssign(t[6], t[7])

	fp2Sub(t[8], t[3], t[1])

	fp2.Square(t[7], &r[0])
	fp2Sub(t[4], t[1], t[4])
	fp2.Mul(&r[0], t[4], t[0])
	fp2.Square(t[2], t[3])
	fp2Double(t[3], t[2])
	fp2AddAssign(t[3], t[2])
	fp2.SquareAssign(t[5])
	fp2Sub(&r[1], t[5], t[3])
	fp2.Mul(&r[2], t[1], t[6])
	fp2Double(t[0], t[7])

	fp2AddAssign(t[0], t[7])
	Fp2Neg(t[6], t[6])

	// line eval
	e.fp2.Mul0Assign(t[6], &e.pairs[k].g1[1])
	e.fp2.Mul0Assign(t[0], &e.pairs[k].g1[0])
	e.fp12.Mul014(f, t[8], t[0], t[6])

}

func (e *Engine) add(f *Fe12, r *PointG2, k int) {
	fp2, t := e.fp2, e.t2

	fp2.Mul(t[0], &e.pairs[k].g2[1], &r[2])
	Fp2Neg(t[0], t[0])
	fp2AddAssign(t[0], &r[1])
	fp2.Mul(t[1], &e.pairs[k].g2[0], &r[2])
	Fp2Neg(t[1], t[1])
	fp2AddAssign(t[1], &r[0])
	fp2.Square(t[2], t[0])
	fp2.Square(t[3], t[1])
	fp2.Mul(t[4], t[1], t[3])
	fp2.Mul(t[2], &r[2], t[2])
	fp2.MulAssign(t[3], &r[0])
	fp2Double(t[5], t[3])
	fp2Sub(t[5], t[4], t[5])
	fp2AddAssign(t[5], t[2])
	fp2.Mul(&r[0], t[1], t[5])
	fp2SubAssign(t[3], t[5])
	fp2.MulAssign(t[3], t[0])
	fp2.Mul(t[2], &r[1], t[4])
	fp2Sub(&r[1], t[3], t[2])
	fp2.MulAssign(&r[2], t[4])
	fp2.Mul(t[2], t[1], &e.pairs[k].g2[1])
	fp2.Mul(t[3], t[0], &e.pairs[k].g2[0])

	fp2SubAssign(t[3], t[2])
	Fp2Neg(t[0], t[0])

	// line eval
	e.fp2.Mul0Assign(t[1], &e.pairs[k].g1[1])
	e.fp2.Mul0Assign(t[0], &e.pairs[k].g1[0])
	e.fp12.Mul014(f, t[3], t[0], t[1])
}

func (e *Engine) nDoubleAdd(f *Fe12, r []PointG2, n int) {
	for i := 0; i < n; i++ {
		e.fp12.SquareAssign(f)
		for j := 0; j < len(e.pairs); j++ {
			e.double(f, &r[j], j)
		}
	}
	for j := 0; j < len(e.pairs); j++ {
		e.add(f, &r[j], j)
	}
}

func (e *Engine) nDouble(f *Fe12, r []PointG2, n int) {
	for i := 0; i < n; i++ {
		e.fp12.SquareAssign(f)
		for j := 0; j < len(e.pairs); j++ {
			e.double(f, &r[j], j)
		}
	}
}

func (e *Engine) millerLoop(f *Fe12) {
	f.one()

	r := make([]PointG2, len(e.pairs))
	for i := 0; i < len(e.pairs); i++ {
		r[i].Set(e.pairs[i].g2)
	}

	for j := 0; j < len(e.pairs); j++ {
		e.double(f, &r[j], j)
	}
	for j := 0; j < len(e.pairs); j++ {
		e.add(f, &r[j], j)
	}

	e.nDoubleAdd(f, r, 2)
	e.nDoubleAdd(f, r, 3)
	e.nDoubleAdd(f, r, 9)
	e.nDoubleAdd(f, r, 32)
	e.nDouble(f, r, 16)

	Fp12Conjugate(f, f)
}

// Exp raises element by x = -15132376222941642752
func (e *Engine) exp(c, a *Fe12) {
	c.set(a)
	e.fp12.CyclotomicSquare(c) // (a ^ 2)

	// (a ^ (2 + 1)) ^ (2 ^ 2) = a ^ 12
	e.fp12.MulAssign(c, a)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)

	// (a ^ (12 + 1)) ^ (2 ^ 3) = a ^ 104
	e.fp12.MulAssign(c, a)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)

	// (a ^ (104 + 1)) ^ (2 ^ 9) = a ^ 53760
	e.fp12.MulAssign(c, a)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	// (a ^ (53760 + 1)) ^ (2 ^ 32) = a ^ 230901736800256
	e.fp12.MulAssign(c, a)
	for i := 0; i < 32; i++ {
		e.fp12.CyclotomicSquare(c)
	}

	// (a ^ (230901736800256 + 1)) ^ (2 ^ 16) = a ^ 15132376222941642752
	e.fp12.MulAssign(c, a)
	for i := 0; i < 16; i++ {
		e.fp12.CyclotomicSquare(c)
	}
	// invert chain result since x is negative
	Fp12Conjugate(c, c)
}

// expDrop raises element by x = -15132376222941642752 / 2
func (e *Engine) expDrop(c, a *Fe12) {
	c.set(a)
	e.fp12.CyclotomicSquare(c) // (a ^ 2)

	// (a ^ (2 + 1)) ^ (2 ^ 2) = a ^ 12
	e.fp12.MulAssign(c, a)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)

	// (a ^ (12 + 1)) ^ (2 ^ 3) = a ^ 104
	e.fp12.MulAssign(c, a)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)

	// (a ^ (104 + 1)) ^ (2 ^ 9) = a ^ 53760
	e.fp12.MulAssign(c, a)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	e.fp12.CyclotomicSquare(c)
	// (a ^ (53760 + 1)) ^ (2 ^ 32) = a ^ 230901736800256
	e.fp12.MulAssign(c, a)
	for i := 0; i < 32; i++ {
		e.fp12.CyclotomicSquare(c)
	}

	// (a ^ (230901736800256 + 1)) ^ (2 ^ 16) = a ^ 15132376222941642752
	e.fp12.MulAssign(c, a)
	for i := 0; i < 15; i++ {
		e.fp12.CyclotomicSquare(c)
	}
	// invert chain result since x is negative
	Fp12Conjugate(c, c)
}

func (e *Engine) finalExp(f *Fe12) {
	t := e.t12
	// Efficient Final Exponentiation via Cyclotomic Structure for Pairings over Families of Elliptic Curves
	// https: //eprint.iacr.org/2020/875.pdf

	// easy part
	Fp12Conjugate(&t[0], f)
	e.fp12.Inverse(f, f)
	e.fp12.MulAssign(&t[0], f)
	f.set(&t[0])
	e.fp12.FrobeniusMap2(f)
	e.fp12.MulAssign(f, &t[0])

	// hard part
	t[0].set(f)
	e.fp12.CyclotomicSquare(&t[0])
	e.expDrop(&t[1], &t[0])
	Fp12Conjugate(&t[2], f)
	e.fp12.MulAssign(&t[1], &t[2])
	e.exp(&t[2], &t[1])
	Fp12Conjugate(&t[1], &t[1])
	e.fp12.MulAssign(&t[1], &t[2])
	e.exp(&t[2], &t[1])
	e.fp12.FrobeniusMap1(&t[1])
	e.fp12.MulAssign(&t[1], &t[2])
	e.fp12.MulAssign(f, &t[0])
	e.exp(&t[0], &t[1])
	e.exp(&t[2], &t[0])
	t[0].set(&t[1])
	e.fp12.FrobeniusMap2(&t[0])
	Fp12Conjugate(&t[1], &t[1])
	e.fp12.MulAssign(&t[1], &t[2])
	e.fp12.MulAssign(&t[1], &t[0])
	e.fp12.MulAssign(f, &t[1])
}

func (e *Engine) calculate() *Fe12 {
	f := e.fp12.One()
	if len(e.pairs) == 0 {
		return f
	}
	e.millerLoop(f)
	e.finalExp(f)
	return f
}

// Check computes pairing and checks if result is Equal to One
func (e *Engine) Check() bool {
	return e.calculate().isOne()
}

// Result computes pairing and returns target group element as result.
func (e *Engine) Result() *E {
	r := e.calculate()
	e.Reset()
	return r
}

// GT returns target group instance.
func (e *Engine) GT() *GT {
	return NewGT()
}
