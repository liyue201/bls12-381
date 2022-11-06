package bls12381

import (
	"bytes"
	"crypto/rand"
	"math/big"
	"testing"
)

func TestFpSerialization(t *testing.T) {
	t.Run("Zero", func(t *testing.T) {
		in := make([]byte, fpByteSize)
		fe, err := FromBytes(in)
		if err != nil {
			t.Fatal(err)
		}
		if !fe.IsZero() {
			t.Fatal("serialization failed")
		}
		if !bytes.Equal(in, ToBytes(fe)) {
			t.Fatal("serialization failed")
		}
	})
	t.Run("Bytes", func(t *testing.T) {
		for i := 0; i < fuz; i++ {
			a, _ := new(Fe).Rand(rand.Reader)
			b, err := FromBytes(ToBytes(a))
			if err != nil {
				t.Fatal(err)
			}
			if !a.Equal(b) {
				t.Fatal("serialization failed")
			}
		}
	})
	t.Run("String", func(t *testing.T) {
		for i := 0; i < fuz; i++ {
			a, _ := new(Fe).Rand(rand.Reader)
			b, err := FromString(ToString(a))
			if err != nil {
				t.Fatal(err)
			}
			if !a.Equal(b) {
				t.Fatal("encoding or decoding failed")
			}
		}
	})
	t.Run("Big", func(t *testing.T) {
		for i := 0; i < fuz; i++ {
			a, _ := new(Fe).Rand(rand.Reader)
			b, err := FromBig(ToBig(a))
			if err != nil {
				t.Fatal(err)
			}
			if !a.Equal(b) {
				t.Fatal("encoding or decoding failed")
			}
		}
	})
}

func TestFpAdditionCrossAgainstBigInt(t *testing.T) {
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		b, _ := new(Fe).Rand(rand.Reader)
		c := new(Fe)
		big_a := a.Big()
		big_b := b.Big()
		big_c := new(big.Int)
		add(c, a, b)
		out_1 := c.Bytes()
		out_2 := padBytes(big_c.Add(big_a, big_b).Mod(big_c, modulus.Big()).Bytes(), fpByteSize)
		if !bytes.Equal(out_1, out_2) {
			t.Fatal("cross test against Big.Int is failed A")
		}
		double(c, a)
		out_1 = c.Bytes()
		out_2 = padBytes(big_c.Add(big_a, big_a).Mod(big_c, modulus.Big()).Bytes(), fpByteSize)
		if !bytes.Equal(out_1, out_2) {
			t.Fatal("cross test against Big.Int is failed B")
		}
		sub(c, a, b)
		out_1 = c.Bytes()
		out_2 = padBytes(big_c.Sub(big_a, big_b).Mod(big_c, modulus.Big()).Bytes(), fpByteSize)
		if !bytes.Equal(out_1, out_2) {
			t.Fatal("cross test against Big.Int is failed C")
		}
		neg(c, a)
		out_1 = c.Bytes()
		out_2 = padBytes(big_c.Neg(big_a).Mod(big_c, modulus.Big()).Bytes(), fpByteSize)
		if !bytes.Equal(out_1, out_2) {
			t.Fatal("cross test against Big.Int is failed D")
		}
	}
}

func TestFpAdditionCrossAgainstBigIntAssigned(t *testing.T) {
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		b, _ := new(Fe).Rand(rand.Reader)
		big_a, big_b := a.Big(), b.Big()
		addAssign(a, b)
		out_1 := a.Bytes()
		out_2 := padBytes(big_a.Add(big_a, big_b).Mod(big_a, modulus.Big()).Bytes(), fpByteSize)
		if !bytes.Equal(out_1, out_2) {
			t.Fatal("cross test against Big.Int is failed A")
		}
		a, _ = new(Fe).Rand(rand.Reader)
		big_a = a.Big()
		doubleAssign(a)
		out_1 = a.Bytes()
		out_2 = padBytes(big_a.Add(big_a, big_a).Mod(big_a, modulus.Big()).Bytes(), fpByteSize)
		if !bytes.Equal(out_1, out_2) {
			t.Fatal("cross test against Big.Int is failed B")
		}
		a, _ = new(Fe).Rand(rand.Reader)
		b, _ = new(Fe).Rand(rand.Reader)
		big_a, big_b = a.Big(), b.Big()
		subAssign(a, b)
		out_1 = a.Bytes()
		out_2 = padBytes(big_a.Sub(big_a, big_b).Mod(big_a, modulus.Big()).Bytes(), fpByteSize)
		if !bytes.Equal(out_1, out_2) {
			t.Fatal("cross test against Big.Int is failed A")
		}
	}
}

func TestFpAdditionProperties(t *testing.T) {
	for i := 0; i < fuz; i++ {

		zero := new(Fe).Zero()
		a, _ := new(Fe).Rand(rand.Reader)
		b, _ := new(Fe).Rand(rand.Reader)
		c1, c2 := new(Fe), new(Fe)
		add(c1, a, zero)
		if !c1.Equal(a) {
			t.Fatal("a + 0 == a")
		}
		sub(c1, a, zero)
		if !c1.Equal(a) {
			t.Fatal("a - 0 == a")
		}
		double(c1, zero)
		if !c1.Equal(zero) {
			t.Fatal("2 * 0 == 0")
		}
		neg(c1, zero)
		if !c1.Equal(zero) {
			t.Fatal("-0 == 0")
		}
		sub(c1, zero, a)
		neg(c2, a)
		if !c1.Equal(c2) {
			t.Fatal("0-a == -a")
		}
		double(c1, a)
		add(c2, a, a)
		if !c1.Equal(c2) {
			t.Fatal("2 * a == a + a")
		}
		add(c1, a, b)
		add(c2, b, a)
		if !c1.Equal(c2) {
			t.Fatal("a + b = b + a")
		}
		sub(c1, a, b)
		sub(c2, b, a)
		neg(c2, c2)
		if !c1.Equal(c2) {
			t.Fatal("a - b = - ( b - a )")
		}
		cx, _ := new(Fe).Rand(rand.Reader)
		add(c1, a, b)
		add(c1, c1, cx)
		add(c2, a, cx)
		add(c2, c2, b)
		if !c1.Equal(c2) {
			t.Fatal("(a + b) + c == (a + c ) + b")
		}
		sub(c1, a, b)
		sub(c1, c1, cx)
		sub(c2, a, cx)
		sub(c2, c2, b)
		if !c1.Equal(c2) {
			t.Fatal("(a - b) - c == (a - c ) -b")
		}
	}
}

func TestFpAdditionPropertiesAssigned(t *testing.T) {
	for i := 0; i < fuz; i++ {
		zero := new(Fe).Zero()
		a, b := new(Fe), new(Fe)
		_, _ = a.Rand(rand.Reader)
		b.Set(a)
		addAssign(a, zero)
		if !a.Equal(b) {
			t.Fatal("a + 0 == a")
		}
		subAssign(a, zero)
		if !a.Equal(b) {
			t.Fatal("a - 0 == a")
		}
		a.Set(zero)
		doubleAssign(a)
		if !a.Equal(zero) {
			t.Fatal("2 * 0 == 0")
		}
		a.Set(zero)
		subAssign(a, b)
		neg(b, b)
		if !a.Equal(b) {
			t.Fatal("0-a == -a")
		}
		_, _ = a.Rand(rand.Reader)
		b.Set(a)
		doubleAssign(a)
		addAssign(b, b)
		if !a.Equal(b) {
			t.Fatal("2 * a == a + a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		c1, c2 := new(Fe).Set(a), new(Fe).Set(b)
		addAssign(c1, b)
		addAssign(c2, a)
		if !c1.Equal(c2) {
			t.Fatal("a + b = b + a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		c1.Set(a)
		c2.Set(b)
		subAssign(c1, b)
		subAssign(c2, a)
		neg(c2, c2)
		if !c1.Equal(c2) {
			t.Fatal("a - b = - ( b - a )")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		c, _ := new(Fe).Rand(rand.Reader)
		a0 := new(Fe).Set(a)
		addAssign(a, b)
		addAssign(a, c)
		addAssign(b, c)
		addAssign(b, a0)
		if !a.Equal(b) {
			t.Fatal("(a + b) + c == (b + c) + a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		_, _ = c.Rand(rand.Reader)
		a0.Set(a)
		subAssign(a, b)
		subAssign(a, c)
		subAssign(a0, c)
		subAssign(a0, b)
		if !a.Equal(a0) {
			t.Fatal("(a - b) - c == (a - c) -b")
		}
	}
}

func TestFpLazyOperations(t *testing.T) {
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		b, _ := new(Fe).Rand(rand.Reader)
		c, _ := new(Fe).Rand(rand.Reader)
		c0 := new(Fe)
		c1 := new(Fe)
		ladd(c0, a, b)
		add(c1, a, b)
		mul(c0, c0, c)
		mul(c1, c1, c)
		if !c0.Equal(c1) {
			t.Fatal("(a + b) * c == (a l+ b) * c")
		}
		_, _ = a.Rand(rand.Reader)
		b.Set(a)
		ldouble(a, a)
		ladd(b, b, b)
		if !a.Equal(b) {
			t.Fatal("2 l* a = a l+ a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		_, _ = c.Rand(rand.Reader)
		a0 := new(Fe).Set(a)
		lsubAssign(a, b)
		laddAssign(a, &modulus)
		mul(a, a, c)
		subAssign(a0, b)
		mul(a0, a0, c)
		if !a.Equal(a0) {
			t.Fatal("((a l- b) + p) * c = (a-b) * c")
		}
	}
}

func TestFpMultiplicationCrossAgainstBigInt(t *testing.T) {
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		b, _ := new(Fe).Rand(rand.Reader)
		c := new(Fe)
		big_a := ToBig(a)
		big_b := ToBig(b)
		big_c := new(big.Int)
		mul(c, a, b)
		out_1 := ToBytes(c)
		out_2 := padBytes(big_c.Mul(big_a, big_b).Mod(big_c, modulus.Big()).Bytes(), fpByteSize)
		if !bytes.Equal(out_1, out_2) {
			t.Fatal("cross test against Big.Int is failed")
		}
	}
}

func TestFpMultiplicationProperties(t *testing.T) {
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		b, _ := new(Fe).Rand(rand.Reader)
		zero, one := new(Fe).Zero(), new(Fe).One()
		c1, c2 := new(Fe), new(Fe)
		mul(c1, a, zero)
		if !c1.Equal(zero) {
			t.Fatal("a * 0 == 0")
		}
		mul(c1, a, one)
		if !c1.Equal(a) {
			t.Fatal("a * 1 == a")
		}
		mul(c1, a, b)
		mul(c2, b, a)
		if !c1.Equal(c2) {
			t.Fatal("a * b == b * a")
		}
		cx, _ := new(Fe).Rand(rand.Reader)
		mul(c1, a, b)
		mul(c1, c1, cx)
		mul(c2, cx, b)
		mul(c2, c2, a)
		if !c1.Equal(c2) {
			t.Fatal("(a * b) * c == (a * c) * b")
		}
		square(a, zero)
		if !a.Equal(zero) {
			t.Fatal("0^2 == 0")
		}
		square(a, one)
		if !a.Equal(one) {
			t.Fatal("1^2 == 1")
		}
		_, _ = a.Rand(rand.Reader)
		square(c1, a)
		mul(c2, a, a)
		if !c1.Equal(c1) {
			t.Fatal("a^2 == a*a")
		}
	}
}

func TestFpExponentiation(t *testing.T) {
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		u := new(Fe)
		Exp(u, a, big.NewInt(0))
		if !u.IsOne() {
			t.Fatal("a^0 == 1")
		}
		Exp(u, a, big.NewInt(1))
		if !u.Equal(a) {
			t.Fatal("a^1 == a")
		}
		v := new(Fe)
		mul(u, a, a)
		mul(u, u, u)
		mul(u, u, u)
		Exp(v, a, big.NewInt(8))
		if !u.Equal(v) {
			t.Fatal("((a^2)^2)^2 == a^8")
		}
		p := modulus.Big()
		Exp(u, a, p)
		if !u.Equal(a) {
			t.Fatal("a^p == a")
		}
		Exp(u, a, p.Sub(p, big.NewInt(1)))
		if !u.IsOne() {
			t.Fatal("a^(p-1) == 1")
		}
	}
}

func TestFpInversion(t *testing.T) {
	for i := 0; i < fuz; i++ {
		u := new(Fe)
		zero, one := new(Fe).Zero(), new(Fe).One()
		Inverse(u, zero)
		if !u.Equal(zero) {
			t.Fatal("(0^-1) == 0)")
		}
		Inverse(u, one)
		if !u.Equal(one) {
			t.Fatal("(1^-1) == 1)")
		}
		a, _ := new(Fe).Rand(rand.Reader)
		Inverse(u, a)
		mul(u, u, a)
		if !u.Equal(one) {
			t.Fatal("(r*a) * r*(a^-1) == r)")
		}
		v := new(Fe)
		p := modulus.Big()
		Exp(u, a, p.Sub(p, big.NewInt(2)))
		Inverse(v, a)
		if !v.Equal(u) {
			t.Fatal("a^(p-2) == a^-1")
		}
	}
}

func TestFpBatchInversion(t *testing.T) {
	n := 20
	for i := 0; i < n; i++ {
		e0 := make([]Fe, n)
		e1 := make([]Fe, n)
		for j := 0; j < n; j++ {
			if j != i {
				e, err := new(Fe).Rand(rand.Reader)
				if err != nil {
					t.Fatal(err)
				}
				e0[j].Set(e)
			}
			Inverse(&e1[j], &e0[j])
		}

		InverseBatch(e0)
		for j := 0; j < n; j++ {
			if !e0[j].Equal(&e1[j]) {
				t.Fatal("batch inversion failed")
			}
		}
	}
}

func TestFpSquareRoot(t *testing.T) {
	if Sqrt(new(Fe), nonResidue1) {
		t.Fatal("non residue cannot have a Sqrt")
	}
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		r0, r1 := new(Fe), new(Fe)
		d0 := Sqrt(r0, a)
		d1 := _sqrt(r1, a)
		if d0 != d1 {
			t.Fatal("Sqrt decision failed")
		}
		if d0 {
			square(r0, r0)
			square(r1, r1)
			if !r0.Equal(r1) {
				t.Fatal("Sqrt failed")
			}
			if !r0.Equal(a) {
				t.Fatal("Sqrt failed")
			}
		}
	}
}

func TestFpNonResidue(t *testing.T) {
	if !IsQuadraticNonResidue(nonResidue1) {
		t.Fatal("element is quadratic non residue, 1")
	}
	if IsQuadraticNonResidue(new(Fe).One()) {
		t.Fatal("One is not quadratic non residue")
	}
	if !IsQuadraticNonResidue(new(Fe).Zero()) {
		t.Fatal("should accept Zero as quadratic non residue")
	}
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		square(a, a)
		if IsQuadraticNonResidue(a) {
			t.Fatal("element is not quadratic non residue")
		}
	}
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		if !Sqrt(new(Fe), a) {
			if !IsQuadraticNonResidue(a) {
				t.Fatal("element is quadratic non residue, 2", i)
			}
		} else {
			i -= 1
		}
	}
}

func TestWFp(t *testing.T) {
	w := new(Wfe)
	a := new(Fe)
	fromWide(a, w)
	if !a.IsZero() {
		t.Fatal("expect Zero")
	}
	w[0] = r1[0]
	w[1] = r1[1]
	w[2] = r1[2]
	w[3] = r1[3]
	w[4] = r1[4]
	w[5] = r1[5]
	fromWide(a, w)
	if !(a[0] == 1 && a[1] == 0 && a[2] == 0 && a[3] == 0 && a[4] == 0 && a[5] == 0) {
		t.Fatal("expect One")
	}
}

func TestWFpAddition(t *testing.T) {
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe).Rand(rand.Reader)
		b, _ := new(Fe).Rand(rand.Reader)
		w0, w1 := new(Wfe), new(Wfe)
		c0, c1 := new(Fe), new(Fe)

		wmul(w0, a, b)
		w1.set(w0)
		wadd(w0, w0, w0)
		wadd(w0, w0, w0)
		wadd(w0, w0, w0)
		lwadd(w1, w1, w1)
		lwadd(w1, w1, w1)
		lwadd(w1, w1, w1)
		fromWide(c0, w0)
		fromWide(c1, w1)

		if !c1.Equal(c0) {
			t.Fatal("addition failed")
		}

		wmul(w0, a, b)
		w1.set(w0)
		wdouble(w0, w0)
		wdouble(w0, w0)
		wdouble(w0, w0)
		lwdouble(w1, w1)
		lwdouble(w1, w1)
		lwdouble(w1, w1)
		fromWide(c0, w0)
		fromWide(c1, w1)

		if !c1.Equal(c0) {
			t.Fatal("doubling failed")
		}

		wmul(w0, a, &Fe{10001})
		wmul(w1, a, &Fe{10000})
		w2 := new(Wfe)
		wsub(w2, w0, w1)
		lwsub(w0, w0, w1)
		fromWide(c0, w2)
		fromWide(c1, w0)

		FromMont(a, a)
		if !c1.Equal(a) {
			t.Fatal("subtraction failed")
		}
		if !c0.Equal(a) {
			t.Fatal("subtraction failed")
		}

		wmul(w0, a, &Fe{10001})
		wmul(w1, a, &Fe{10000})
		wsub(w0, w1, w0)
		fromWide(c0, w0)

		neg(a, a)
		FromMont(a, a)
		if !c0.Equal(a) {
			t.Fatal("subtraction failed")
		}

	}
}

func TestWFpMultiplication(t *testing.T) {
	for i := 0; i < fuz; i++ {
		a0, _ := new(Fe).Rand(rand.Reader)
		b0, _ := new(Fe).Rand(rand.Reader)
		a1, _ := new(Fe).Rand(rand.Reader)
		b1, _ := new(Fe).Rand(rand.Reader)
		w0, w1, w2, w3 := new(Wfe), new(Wfe), new(Wfe), new(Wfe)
		c0, c1 := new(Fe), new(Fe)
		r0, r1 := new(Fe), new(Fe)

		wmul(w0, a0, b0)
		fromWide(r0, w0)
		mul(r1, a0, b0)

		if !r1.Equal(r0) {
			t.Fatal("multiplication failed")
		}

		wmul(w0, a0, b0)
		wmul(w1, a1, b1)
		lwadd(w0, w0, w1)
		fromWide(r0, w0)

		mul(c0, a0, b0)
		mul(c1, a1, b1)
		add(r1, c0, c1)

		if !r1.Equal(r0) {
			t.Fatal("multiplication failed")
		}

		wmul(w0, a0, b0)
		wmul(w1, a0, b1)
		wmul(w2, a1, b0)
		wmul(w3, a1, b1)
		lwadd(w0, w0, w1)
		lwadd(w0, w0, w2)
		lwadd(w0, w0, w3)
		fromWide(r0, w0)

		add(c0, a0, a1)
		add(c1, b0, b1)
		mul(r1, c0, c1)

		if !r1.Equal(r0) {
			t.Fatal("multiplication failed")
		}
	}
}

func TestFp2Serialization(t *testing.T) {
	f := NewFp2()
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		b, err := f.FromBytes(f.ToBytes(a))
		if err != nil {
			t.Fatal(err)
		}
		if !a.Equal(b) {
			t.Fatal("serialization failed")
		}
	}
}

func TestFp2AdditionProperties(t *testing.T) {
	f := NewFp2()
	for i := 0; i < fuz; i++ {
		zero := f.Zero()
		a, _ := new(Fe2).Rand(rand.Reader)
		b, _ := new(Fe2).Rand(rand.Reader)
		c1 := f.New()
		c2 := f.New()
		fp2Add(c1, a, zero)
		if !c1.Equal(a) {
			t.Fatal("a + 0 == a")
		}
		fp2Sub(c1, a, zero)
		if !c1.Equal(a) {
			t.Fatal("a - 0 == a")
		}
		fp2Double(c1, zero)
		if !c1.Equal(zero) {
			t.Fatal("2 * 0 == 0")
		}
		Fp2Neg(c1, zero)
		if !c1.Equal(zero) {
			t.Fatal("-0 == 0")
		}
		fp2Sub(c1, zero, a)
		Fp2Neg(c2, a)
		if !c1.Equal(c2) {
			t.Fatal("0-a == -a")
		}
		fp2Double(c1, a)
		fp2Add(c2, a, a)
		if !c1.Equal(c2) {
			t.Fatal("2 * a == a + a")
		}
		fp2Add(c1, a, b)
		fp2Add(c2, b, a)
		if !c1.Equal(c2) {
			t.Fatal("a + b = b + a")
		}
		fp2Sub(c1, a, b)
		fp2Sub(c2, b, a)
		Fp2Neg(c2, c2)
		if !c1.Equal(c2) {
			t.Fatal("a - b = - ( b - a )")
		}
		cx, _ := new(Fe2).Rand(rand.Reader)
		fp2Add(c1, a, b)
		fp2Add(c1, c1, cx)
		fp2Add(c2, a, cx)
		fp2Add(c2, c2, b)
		if !c1.Equal(c2) {
			t.Fatal("(a + b) + c == (a + c ) + b")
		}
		fp2Sub(c1, a, b)
		fp2Sub(c1, c1, cx)
		fp2Sub(c2, a, cx)
		fp2Sub(c2, c2, b)
		if !c1.Equal(c2) {
			t.Fatal("(a - b) - c == (a - c ) -b")
		}
	}
}

func TestFp2AdditionPropertiesAssigned(t *testing.T) {
	for i := 0; i < fuz; i++ {
		zero := new(Fe2).Zero()
		a, b := new(Fe2), new(Fe2)
		_, _ = a.Rand(rand.Reader)
		b.Set(a)
		fp2AddAssign(a, zero)
		if !a.Equal(b) {
			t.Fatal("a + 0 == a")
		}
		fp2SubAssign(a, zero)
		if !a.Equal(b) {
			t.Fatal("a - 0 == a")
		}
		a.Set(zero)
		fp2DoubleAssign(a)
		if !a.Equal(zero) {
			t.Fatal("2 * 0 == 0")
		}
		a.Set(zero)
		fp2SubAssign(a, b)
		Fp2Neg(b, b)
		if !a.Equal(b) {
			t.Fatal("0-a == -a")
		}
		_, _ = a.Rand(rand.Reader)
		b.Set(a)
		fp2DoubleAssign(a)
		fp2AddAssign(b, b)
		if !a.Equal(b) {
			t.Fatal("2 * a == a + a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		c1, c2 := new(Fe2).Set(a), new(Fe2).Set(b)
		fp2AddAssign(c1, b)
		fp2AddAssign(c2, a)
		if !c1.Equal(c2) {
			t.Fatal("a + b = b + a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		c1.Set(a)
		c2.Set(b)
		fp2SubAssign(c1, b)
		fp2SubAssign(c2, a)
		Fp2Neg(c2, c2)
		if !c1.Equal(c2) {
			t.Fatal("a - b = - ( b - a )")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		c, _ := new(Fe2).Rand(rand.Reader)
		a0 := new(Fe2).Set(a)
		fp2AddAssign(a, b)
		fp2AddAssign(a, c)
		fp2AddAssign(b, c)
		fp2AddAssign(b, a0)
		if !a.Equal(b) {
			t.Fatal("(a + b) + c == (b + c) + a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		_, _ = c.Rand(rand.Reader)
		a0.Set(a)
		fp2SubAssign(a, b)
		fp2SubAssign(a, c)
		fp2SubAssign(a0, c)
		fp2SubAssign(a0, b)
		if !a.Equal(a0) {
			t.Fatal("(a - b) - c == (a - c) -b")
		}
	}
}

func TestFp2LazyOperations(t *testing.T) {
	f := NewFp2()
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		b, _ := new(Fe2).Rand(rand.Reader)
		c, _ := new(Fe2).Rand(rand.Reader)
		c0 := new(Fe2)
		c1 := new(Fe2)
		fp2Ladd(c0, a, b)
		fp2Add(c1, a, b)
		fp2LaddAssign(a, b)
		f.MulAssign(c0, c)
		f.MulAssign(c1, c)
		f.MulAssign(a, c)
		if !c0.Equal(c1) {
			t.Fatal("(a + b) * c == (a l+ b) * c")
		}
		if !c0.Equal(c1) {
			t.Fatal("(a + b) * c == (a l+ b) * c")
		}
	}
}

func TestFp2MulByNonResidue(t *testing.T) {
	f := NewFp2()
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		r0, r1, r2, r3 := new(Fe2), new(Fe2), new(Fe2), new(Fe2)
		f.Mul(r1, a, nonResidue2)
		mulByNonResidue(r0, a)
		r2.Set(a)
		mulByNonResidueAssign(r2)
		_fp2MulByNonResidue(r3, a)

		if !r0.Equal(r1) {
			t.Fatal("Mul by non residue failed")
		}
		if !r0.Equal(r2) {
			t.Fatal("Mul by non residue failed")
		}
		if !r0.Equal(r3) {
			t.Fatal("Mul by non residue failed")
		}
	}
}

func TestFp2MultiplicationProperties(t *testing.T) {
	f := NewFp2()
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		b, _ := new(Fe2).Rand(rand.Reader)
		zero := f.Zero()
		one := f.One()
		c1, c2 := f.New(), f.New()
		f.Mul(c1, a, zero)
		if !c1.Equal(zero) {
			t.Fatal("a * 0 == 0")
		}
		f.Mul(c1, a, one)
		if !c1.Equal(a) {
			t.Fatal("a * 1 == a")
		}
		f.Mul(c1, a, b)
		f.Mul(c2, b, a)
		if !c1.Equal(c2) {
			t.Fatal("a * b == b * a")
		}
		cx, _ := new(Fe2).Rand(rand.Reader)
		f.Mul(c1, a, b)
		f.Mul(c1, c1, cx)
		f.Mul(c2, cx, b)
		f.Mul(c2, c2, a)
		if !c1.Equal(c2) {
			t.Fatal("(a * b) * c == (a * c) * b")
		}
		f.Square(a, zero)
		if !a.Equal(zero) {
			t.Fatal("0^2 == 0")
		}
		f.Square(a, one)
		if !a.Equal(one) {
			t.Fatal("1^2 == 1")
		}
		_, _ = a.Rand(rand.Reader)
		f.Square(c1, a)
		f.Mul(c2, a, a)
		if !c2.Equal(c1) {
			t.Fatal("a^2 == a*a")
		}
	}
}

func TestFp2MultiplicationPropertiesAssigned(t *testing.T) {
	f := NewFp2()
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		zero, one := new(Fe2).Zero(), new(Fe2).One()
		f.MulAssign(a, zero)
		if !a.Equal(zero) {
			t.Fatal("a * 0 == 0")
		}
		_, _ = a.Rand(rand.Reader)
		a0 := new(Fe2).Set(a)
		f.MulAssign(a, one)
		if !a.Equal(a0) {
			t.Fatal("a * 1 == a")
		}
		_, _ = a.Rand(rand.Reader)
		b, _ := new(Fe2).Rand(rand.Reader)
		a0.Set(a)
		f.MulAssign(a, b)
		f.MulAssign(b, a0)
		if !a.Equal(b) {
			t.Fatal("a * b == b * a")
		}
		c, _ := new(Fe2).Rand(rand.Reader)
		a0.Set(a)
		f.MulAssign(a, b)
		f.MulAssign(a, c)
		f.MulAssign(a0, c)
		f.MulAssign(a0, b)
		if !a.Equal(a0) {
			t.Fatal("(a * b) * c == (a * c) * b")
		}
		a0.Set(a)
		f.SquareAssign(a)
		f.MulAssign(a0, a0)
		if !a.Equal(a0) {
			t.Fatal("a^2 == a*a")
		}
	}
}

func TestFp2Exponentiation(t *testing.T) {
	f := NewFp2()
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		u := f.New()
		f.Exp(u, a, big.NewInt(0))
		if !u.Equal(f.One()) {
			t.Fatal("a^0 == 1")
		}
		f.Exp(u, a, big.NewInt(1))
		if !u.Equal(a) {
			t.Fatal("a^1 == a")
		}
		v := f.New()
		f.Mul(u, a, a)
		f.Mul(u, u, u)
		f.Mul(u, u, u)
		f.Exp(v, a, big.NewInt(8))
		if !u.Equal(v) {
			t.Fatal("((a^2)^2)^2 == a^8")
		}
	}
}

func TestFp2Inversion(t *testing.T) {
	f := NewFp2()
	u := f.New()
	zero := f.Zero()
	one := f.One()
	f.Inverse(u, zero)
	if !u.Equal(zero) {
		t.Fatal("(0 ^ -1) == 0)")
	}
	f.Inverse(u, one)
	if !u.Equal(one) {
		t.Fatal("(1 ^ -1) == 1)")
	}
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		f.Inverse(u, a)
		f.Mul(u, u, a)
		if !u.Equal(one) {
			t.Fatal("(r * a) * r * (a ^ -1) == r)")
		}
	}
}

func TestFp2BatchInversion(t *testing.T) {
	f := NewFp2()
	n := 20
	for i := 0; i < n; i++ {
		e0 := make([]Fe2, n)
		e1 := make([]Fe2, n)
		for j := 0; j < n; j++ {
			if j != i {
				e, err := new(Fe2).Rand(rand.Reader)
				if err != nil {
					t.Fatal(err)
				}
				e0[j].Set(e)
			}
			f.Inverse(&e1[j], &e0[j])
		}
		f.InverseBatch(e0)
		for j := 0; j < n; j++ {
			if !e0[j].Equal(&e1[j]) {
				t.Fatal("batch inversion failed")
			}
		}
	}
}

func TestFp2SquareRoot(t *testing.T) {
	e := NewFp2()
	if e.SqrtBLST(e.New(), nonResidue2) {
		t.Fatal("non residue cannot have a Sqrt")
	}
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		r0, r1 := new(Fe2), new(Fe2)
		d0 := e.Sqrt(r0, a)
		d1 := e.SqrtBLST(r1, a)
		if d0 != d1 {
			t.Fatal("Sqrt decision failed")
		}
		if d0 {
			e.Square(r0, r0)
			e.Square(r1, r1)
			if !r0.Equal(r1) {
				t.Fatal("Sqrt failed")
			}
			if !r0.Equal(a) {
				t.Fatal("Sqrt failed")
			}
		}
	}
}

func TestFp2NonResidue(t *testing.T) {
	f := NewFp2()
	if !f.IsQuadraticNonResidue(nonResidue2) {
		t.Fatal("element is quadratic non residue, 1")
	}
	if f.IsQuadraticNonResidue(new(Fe2).One()) {
		t.Fatal("One is not quadratic non residue")
	}
	if !f.IsQuadraticNonResidue(new(Fe2).Zero()) {
		t.Fatal("should accept Zero as quadratic non residue")
	}
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		f.SquareAssign(a)
		if f.IsQuadraticNonResidue(a) {
			t.Fatal("element is not quadratic non residue")
		}
	}
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe2).Rand(rand.Reader)
		if !f.Sqrt(new(Fe2), a) {
			if !f.IsQuadraticNonResidue(a) {
				t.Fatal("element is quadratic non residue, 2", i)
			}
		} else {
			i -= 1
		}
	}
}

func TestWFp2Addition(t *testing.T) {
	for i := 0; i < fuz; i++ {
		r0, _ := new(Fe2).Rand(rand.Reader)
		r1, _ := new(Fe2).Rand(rand.Reader)
		r2, _ := new(Fe2).Rand(rand.Reader)
		rw0, rw1, w0, w1 := new(Wfe2), new(Wfe2), new(Wfe2), new(Wfe2)

		wfp2Mul(w0, r0, r1)
		wfp2Mul(w1, r0, r2)

		_wfp2Add(rw0, w0, w1)
		wfp2Add(rw1, w0, w1)
		if !rw0.equal(rw1) {
			t.Fatal("add failed")
		}
		rw1.set(w0)
		wfp2AddAssign(rw1, w1)
		if !rw0.equal(rw1) {
			t.Fatal("assigned add failed")
		}

		_wfp2AddMixed(rw0, w0, w1)
		wfp2AddMixed(rw1, w0, w1)
		if !rw0.equal(rw1) {
			t.Fatal("add mixed failed")
		}
		rw1.set(w0)
		wfp2AddMixedAssign(rw1, w1)
		if !rw0.equal(rw1) {
			t.Fatal("assigned mixed add failed")
		}

		_wfp2Ladd(rw0, w0, w1)
		wfp2Ladd(rw1, w0, w1)
		if !rw0.equal(rw1) {
			t.Fatal("lazy add failed")
		}
		rw1.set(w0)
		wfp2LaddAssign(rw1, w1)
		if !rw0.equal(rw1) {
			t.Fatal("assigned lazy add failed")
		}

		_wfp2Sub(rw0, w0, w1)
		wfp2Sub(rw1, w0, w1)
		if !rw0.equal(rw1) {
			t.Fatal("sub failed")
		}
		rw1.set(w0)
		wfp2SubAssign(rw1, w1)
		if !rw0.equal(rw1) {
			t.Fatal("assigned sub failed")
		}

		_wfp2SubMixed(rw0, w0, w1)
		wfp2SubMixed(rw1, w0, w1)
		if !rw0.equal(rw1) {
			t.Fatal("sub mixed failed")
		}
		rw1.set(w0)
		wfp2SubMixedAssign(rw1, w1)
		if !rw0.equal(rw1) {
			t.Fatal("assigned sub mixed failed")
		}

		_wfp2Double(rw0, w0)
		wfp2Double(rw1, w0)
		if !rw0.equal(rw1) {
			t.Fatal("doubling failed")
		}
		rw1.set(w0)
		wfp2DoubleAssign(rw1)
		if !rw0.equal(rw1) {
			t.Fatal("assigned doubling failed")
		}

	}
}

func TestFp2MultiplicationCross(t *testing.T) {
	f := NewFp2()
	a, b, c0, c1 := new(Fe2), new(Fe2), new(Fe2), new(Fe2)
	w0, w1 := new(Wfe2), new(Wfe2)
	for i := 0; i < fuz; i++ {
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		_wfp2Mul(w0, a, b)
		wfp2Mul(w1, a, b)
		if !w0.equal(w1) {
			t.Fatal("multiplication failed")
		}
		c0.FromWide(w0)
		f.Mul(c1, a, b)
		if !c0.Equal(c1) {
			t.Fatal("multiplication failed")
		}
	}
}

func TestFp2SquareCross(t *testing.T) {
	f := NewFp2()
	a, c0, c1 := new(Fe2), new(Fe2), new(Fe2)
	w0, w1 := new(Wfe2), new(Wfe2)
	for i := 0; i < fuz; i++ {
		_, _ = a.Rand(rand.Reader)
		_wfp2Square(w0, a)
		wfp2Square(w1, a)
		if !w0.equal(w1) {
			t.Fatal("squaring failed")
		}
		c0.FromWide(w0)
		f.Square(c1, a)
		if !c0.Equal(c1) {
			t.Fatal("squaring failed")
		}
	}
}

func TestWFp2MulByNonResidue(t *testing.T) {
	f := NewFp2()
	a, b, c0, c1 := new(Fe2), new(Fe2), new(Fe2), new(Fe2)
	w0, w1, w2, w3 := new(Wfe2), new(Wfe2), new(Wfe2), new(Wfe2)
	for i := 0; i < fuz; i++ {
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		r := new(Fe2)

		f.Mul(r, a, b)
		wfp2Mul(w0, a, b)

		mulByNonResidue(c0, r)
		wfp2MulByNonResidue(w1, w0)
		w2.set(w0)
		wfp2MulByNonResidueAssign(w2)
		_wfp2MulByNonResidue(w3, w0)
		if !w1.equal(w2) {
			t.Fatal("Mul by non residue failed")
		}
		if !w1.equal(w3) {
			t.Fatal("Mul by non residue failed")
		}
		c1.FromWide(w1)
		if !c0.Equal(c1) {
			t.Fatal("Mul by non residue failed")
		}
	}
}

func TestFp6Serialization(t *testing.T) {
	f := NewFp6(nil)
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe6).Rand(rand.Reader)
		b, err := f.FromBytes(f.ToBytes(a))
		if err != nil {
			t.Fatal(err)
		}
		if !a.Equal(b) {
			t.Fatal("serialization")
		}
	}
}

func TestFp6AdditionProperties(t *testing.T) {
	f := NewFp6(nil)
	for i := 0; i < fuz; i++ {
		zero := f.Zero()
		a, _ := new(Fe6).Rand(rand.Reader)
		b, _ := new(Fe6).Rand(rand.Reader)
		c1 := f.New()
		c2 := f.New()
		Fp6Add(c1, a, zero)
		if !c1.Equal(a) {
			t.Fatal("a + 0 == a")
		}
		Fp6Sub(c1, a, zero)
		if !c1.Equal(a) {
			t.Fatal("a - 0 == a")
		}
		Fp6Double(c1, zero)
		if !c1.Equal(zero) {
			t.Fatal("2 * 0 == 0")
		}
		Fp6Neg(c1, zero)
		if !c1.Equal(zero) {
			t.Fatal("-0 == 0")
		}
		Fp6Sub(c1, zero, a)
		Fp6Neg(c2, a)
		if !c1.Equal(c2) {
			t.Fatal("0-a == -a")
		}
		Fp6Double(c1, a)
		Fp6Add(c2, a, a)
		if !c1.Equal(c2) {
			t.Fatal("2 * a == a + a")
		}
		Fp6Add(c1, a, b)
		Fp6Add(c2, b, a)
		if !c1.Equal(c2) {
			t.Fatal("a + b = b + a")
		}
		Fp6Sub(c1, a, b)
		Fp6Sub(c2, b, a)
		Fp6Neg(c2, c2)
		if !c1.Equal(c2) {
			t.Fatal("a - b = - ( b - a )")
		}
		cx, _ := new(Fe6).Rand(rand.Reader)
		Fp6Add(c1, a, b)
		Fp6Add(c1, c1, cx)
		Fp6Add(c2, a, cx)
		Fp6Add(c2, c2, b)
		if !c1.Equal(c2) {
			t.Fatal("(a + b) + c == (a + c ) + b")
		}
		Fp6Sub(c1, a, b)
		Fp6Sub(c1, c1, cx)
		Fp6Sub(c2, a, cx)
		Fp6Sub(c2, c2, b)
		if !c1.Equal(c2) {
			t.Fatal("(a - b) - c == (a - c ) -b")
		}
	}
}

func TestFp6AdditionPropertiesAssigned(t *testing.T) {
	for i := 0; i < fuz; i++ {
		zero := new(Fe6).Zero()
		a, b := new(Fe6), new(Fe6)
		_, _ = a.Rand(rand.Reader)
		b.Set(a)
		Fp6AddAssign(a, zero)
		if !a.Equal(b) {
			t.Fatal("a + 0 == a")
		}
		Fp6SubAssign(a, zero)
		if !a.Equal(b) {
			t.Fatal("a - 0 == a")
		}
		a.Set(zero)
		Fp6DoubleAssign(a)
		if !a.Equal(zero) {
			t.Fatal("2 * 0 == 0")
		}
		a.Set(zero)
		Fp6SubAssign(a, b)
		Fp6Neg(b, b)
		if !a.Equal(b) {
			t.Fatal("0-a == -a")
		}
		_, _ = a.Rand(rand.Reader)
		b.Set(a)
		Fp6DoubleAssign(a)
		Fp6AddAssign(b, b)
		if !a.Equal(b) {
			t.Fatal("2 * a == a + a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		c1, c2 := new(Fe6).Set(a), new(Fe6).Set(b)
		Fp6AddAssign(c1, b)
		Fp6AddAssign(c2, a)
		if !c1.Equal(c2) {
			t.Fatal("a + b = b + a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		c1.Set(a)
		c2.Set(b)
		Fp6SubAssign(c1, b)
		Fp6SubAssign(c2, a)
		Fp6Neg(c2, c2)
		if !c1.Equal(c2) {
			t.Fatal("a - b = - ( b - a )")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		c, _ := new(Fe6).Rand(rand.Reader)
		a0 := new(Fe6).Set(a)
		Fp6AddAssign(a, b)
		Fp6AddAssign(a, c)
		Fp6AddAssign(b, c)
		Fp6AddAssign(b, a0)
		if !a.Equal(b) {
			t.Fatal("(a + b) + c == (b + c) + a")
		}
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		_, _ = c.Rand(rand.Reader)
		a0.Set(a)
		Fp6SubAssign(a, b)
		Fp6SubAssign(a, c)
		Fp6SubAssign(a0, c)
		Fp6SubAssign(a0, b)
		if !a.Equal(a0) {
			t.Fatal("(a - b) - c == (a - c) -b")
		}
	}
}

func TestFp6LazyOperations(t *testing.T) {
	f := NewFp6(nil)
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe6).Rand(rand.Reader)
		b, _ := new(Fe6).Rand(rand.Reader)
		c, _ := new(Fe6).Rand(rand.Reader)
		c0 := new(Fe6)
		c1 := new(Fe6)
		Fp6Ladd(c0, a, b)
		Fp6Add(c1, a, b)
		f.MulAssign(c0, c)
		f.MulAssign(c1, c)
		if !c0.Equal(c1) {
			t.Fatal("(a + b) * c == (a l+ b) * c")
		}
		if !c0.Equal(c1) {
			t.Fatal("(a + b) * c == (a l+ b) * c")
		}
	}
}

func TestFp6SparseMultiplication(t *testing.T) {
	fp6 := NewFp6(nil)
	var a, b, u *Fe6
	for i := 0; i < fuz; i++ {
		a, _ = new(Fe6).Rand(rand.Reader)
		b, _ = new(Fe6).Rand(rand.Reader)
		u, _ = new(Fe6).Rand(rand.Reader)
		b[2].Zero()
		fp6.Mul(u, a, b)
		fp6._mul01(a, a, &b[0], &b[1])
		if !a.Equal(u) {
			t.Fatal("Mul by 01")
		}
	}
	for i := 0; i < fuz; i++ {
		a, _ = new(Fe6).Rand(rand.Reader)
		b, _ = new(Fe6).Rand(rand.Reader)
		u, _ = new(Fe6).Rand(rand.Reader)
		b[2].Zero()
		b[0].Zero()
		fp6.Mul(u, a, b)
		fp6._mul1(a, a, &b[1])
		if !a.Equal(u) {
			t.Fatal("Mul by 1")
		}
	}
}

func TestFp6MultiplicationProperties(t *testing.T) {
	f := NewFp6(nil)
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe6).Rand(rand.Reader)
		b, _ := new(Fe6).Rand(rand.Reader)
		zero := f.Zero()
		one := f.One()
		c1, c2 := f.New(), f.New()
		f.Mul(c1, a, zero)
		if !c1.Equal(zero) {
			t.Fatal("a * 0 == 0")
		}
		f.Mul(c1, a, one)
		if !c1.Equal(a) {
			t.Fatal("a * 1 == a")
		}
		f.Mul(c1, a, b)
		f.Mul(c2, b, a)
		if !c1.Equal(c2) {
			t.Fatal("a * b == b * a")
		}
		cx, _ := new(Fe6).Rand(rand.Reader)
		f.Mul(c1, a, b)
		f.Mul(c1, c1, cx)
		f.Mul(c2, cx, b)
		f.Mul(c2, c2, a)
		if !c1.Equal(c2) {
			t.Fatal("(a * b) * c == (a * c) * b")
		}
		f.Square(a, zero)
		if !a.Equal(zero) {
			t.Fatal("0^2 == 0")
		}
		f.Square(a, one)
		if !a.Equal(one) {
			t.Fatal("1^2 == 1")
		}
		_, _ = a.Rand(rand.Reader)
		f.Square(c1, a)
		f.Mul(c2, a, a)
		if !c2.Equal(c1) {
			t.Fatal("a^2 == a*a")
		}
	}
}

func TestFp6MultiplicationPropertiesAssigned(t *testing.T) {
	f := NewFp6(nil)
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe6).Rand(rand.Reader)
		zero, one := new(Fe6).Zero(), new(Fe6).One()
		f.MulAssign(a, zero)
		if !a.Equal(zero) {
			t.Fatal("a * 0 == 0")
		}
		_, _ = a.Rand(rand.Reader)
		a0 := new(Fe6).Set(a)
		f.MulAssign(a, one)
		if !a.Equal(a0) {
			t.Fatal("a * 1 == a")
		}
		_, _ = a.Rand(rand.Reader)
		b, _ := new(Fe6).Rand(rand.Reader)
		a0.Set(a)
		f.MulAssign(a, b)
		f.MulAssign(b, a0)
		if !a.Equal(b) {
			t.Fatal("a * b == b * a")
		}
		c, _ := new(Fe6).Rand(rand.Reader)
		a0.Set(a)
		f.MulAssign(a, b)
		f.MulAssign(a, c)
		f.MulAssign(a0, c)
		f.MulAssign(a0, b)
		if !a.Equal(a0) {
			t.Fatal("(a * b) * c == (a * c) * b")
		}
	}
}

func TestFp6Exponentiation(t *testing.T) {
	f := NewFp6(nil)
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe6).Rand(rand.Reader)
		u := f.New()
		f.Exp(u, a, big.NewInt(0))
		if !u.Equal(f.One()) {
			t.Fatal("a^0 == 1")
		}
		f.Exp(u, a, big.NewInt(1))
		if !u.Equal(a) {
			t.Fatal("a^1 == a")
		}
		v := f.New()
		f.Exp(v, a, big.NewInt(8))
		f.Square(u, a)
		f.Square(u, u)
		f.Square(u, u)
		if !u.Equal(v) {
			t.Fatal("((a^2)^2)^2 == a^8", i)
		}
	}
}

func TestFp6Inversion(t *testing.T) {
	f := NewFp6(nil)
	for i := 0; i < fuz; i++ {
		u := f.New()
		zero := f.Zero()
		one := f.One()
		f.Inverse(u, zero)
		if !u.Equal(zero) {
			t.Fatal("(0^-1) == 0)")
		}
		f.Inverse(u, one)
		if !u.Equal(one) {
			t.Fatal("(1^-1) == 1)")
		}
		a, _ := new(Fe6).Rand(rand.Reader)
		f.Inverse(u, a)
		f.Mul(u, u, a)
		if !u.Equal(one) {
			t.Fatal("(r*a) * r*(a^-1) == r)")
		}
	}
}

func TestFp6MultiplicationCross(t *testing.T) {
	f := NewFp6(nil)
	a, b, c0, c1, c2, c3 := new(Fe6), new(Fe6), new(Fe6), new(Fe6), new(Fe6), new(Fe6)
	w0 := new(Wfe6)
	for i := 0; i < fuz; i++ {

		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		f.Wmul(w0, a, b)
		c0.FromWide(w0)
		f.Mul(c1, a, b)
		f._mul(c2, a, b)
		c3.Set(a)
		f.MulAssign(c3, b)
		if !c0.Equal(c1) {
			t.Fatal("multiplication failed")
		}
		if !c0.Equal(c2) {
			t.Fatal("multiplication failed")
		}
		if !c0.Equal(c3) {
			t.Fatal("multiplication failed")
		}

	}
}

func TestFp6SquareCross(t *testing.T) {
	f := NewFp6(nil)
	a, c0, c1, c2 := new(Fe6), new(Fe6), new(Fe6), new(Fe6)
	w0 := new(Wfe6)
	for i := 0; i < fuz; i++ {
		_, _ = a.Rand(rand.Reader)
		f.Wsquare(w0, a)
		c0.FromWide(w0)
		f.Square(c1, a)
		f._square(c2, a)

		if !c0.Equal(c2) {
			t.Fatal("squaring failed")
		}
		if !c0.Equal(c1) {
			t.Fatal("squaring failed")
		}
	}
}

func TestFp6SparseMultiplicationCross(t *testing.T) {
	f := NewFp6(nil)
	a, c0, c1, c2 := new(Fe6), new(Fe6), new(Fe6), new(Fe6)
	w0 := new(Wfe6)
	for i := 0; i < fuz; i++ {
		// mul01
		{
			_, _ = a.Rand(rand.Reader)
			b0, _ := new(Fe2).Rand(rand.Reader)
			b1, _ := new(Fe2).Rand(rand.Reader)
			b := new(Fe6)
			b[0].Set(b0)
			b[1].Set(b1)

			f.Wmul01(w0, a, b0, b1)
			c0.FromWide(w0)

			f._mul01(c1, a, b0, b1)
			f._mul(c2, a, b)

			if !c2.Equal(c1) {
				t.Fatal("sparse multiplication 01 failed")
			}

			if !c0.Equal(c1) {
				t.Fatal("sparse multiplication 01 failed")
			}

		}
		// Mul0
		{
			_, _ = a.Rand(rand.Reader)
			b1, _ := new(Fe2).Rand(rand.Reader)
			b := new(Fe6)
			b[1].Set(b1)

			f.Wmul1(w0, a, b1)
			c0.FromWide(w0)
			f._mul1(c1, a, b1)
			f._mul(c2, a, b)

			if !c2.Equal(c0) {
				t.Fatal("sparse multiplication 0 failed")
			}
			if !c2.Equal(c1) {
				t.Fatal("sparse multiplication 0 failed")
			}
		}
	}
}

func TestFp12Serialization(t *testing.T) {
	f := NewFp12(nil)
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe12).Rand(rand.Reader)
		b, err := f.FromBytes(f.ToBytes(a))
		if err != nil {
			t.Fatal(err)
		}
		if !a.IsEqual(b) {
			t.Fatal("serialization")
		}
	}
}

func TestFp12AdditionProperties(t *testing.T) {
	f := NewFp12(nil)
	for i := 0; i < fuz; i++ {
		zero := f.Zero()
		a, _ := new(Fe12).Rand(rand.Reader)
		b, _ := new(Fe12).Rand(rand.Reader)
		c1 := f.New()
		c2 := f.New()
		Fp12Add(c1, a, zero)
		if !c1.IsEqual(a) {
			t.Fatal("a + 0 == a")
		}
		Fp12Sub(c1, a, zero)
		if !c1.IsEqual(a) {
			t.Fatal("a - 0 == a")
		}
		Fp12Double(c1, zero)
		if !c1.IsEqual(zero) {
			t.Fatal("2 * 0 == 0")
		}
		Fp12Neg(c1, zero)
		if !c1.IsEqual(zero) {
			t.Fatal("-0 == 0")
		}
		Fp12Sub(c1, zero, a)
		Fp12Neg(c2, a)
		if !c1.IsEqual(c2) {
			t.Fatal("0-a == -a")
		}
		Fp12Double(c1, a)
		Fp12Add(c2, a, a)
		if !c1.IsEqual(c2) {
			t.Fatal("2 * a == a + a")
		}
		Fp12Add(c1, a, b)
		Fp12Add(c2, b, a)
		if !c1.IsEqual(c2) {
			t.Fatal("a + b = b + a")
		}
		Fp12Sub(c1, a, b)
		Fp12Sub(c2, b, a)
		Fp12Neg(c2, c2)
		if !c1.IsEqual(c2) {
			t.Fatal("a - b = - ( b - a )")
		}
		cx, _ := new(Fe12).Rand(rand.Reader)
		Fp12Add(c1, a, b)
		Fp12Add(c1, c1, cx)
		Fp12Add(c2, a, cx)
		Fp12Add(c2, c2, b)
		if !c1.IsEqual(c2) {
			t.Fatal("(a + b) + c == (a + c ) + b")
		}
		Fp12Sub(c1, a, b)
		Fp12Sub(c1, c1, cx)
		Fp12Sub(c2, a, cx)
		Fp12Sub(c2, c2, b)
		if !c1.IsEqual(c2) {
			t.Fatal("(a - b) - c == (a - c ) -b")
		}
	}
}

func TestFp12MultiplicationProperties(t *testing.T) {
	f := NewFp12(nil)
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe12).Rand(rand.Reader)
		b, _ := new(Fe12).Rand(rand.Reader)
		zero := f.Zero()
		one := f.One()
		c1, c2 := f.New(), f.New()
		f.Mul(c1, a, zero)
		if !c1.IsEqual(zero) {
			t.Fatal("a * 0 == 0")
		}
		f.Mul(c1, a, one)
		if !c1.IsEqual(a) {
			t.Fatal("a * 1 == a")
		}
		f.Mul(c1, a, b)
		f.Mul(c2, b, a)
		if !c1.IsEqual(c2) {
			t.Fatal("a * b == b * a")
		}
		cx, _ := new(Fe12).Rand(rand.Reader)
		f.Mul(c1, a, b)
		f.Mul(c1, c1, cx)
		f.Mul(c2, cx, b)
		f.Mul(c2, c2, a)
		if !c1.IsEqual(c2) {
			t.Fatal("(a * b) * c == (a * c) * b")
		}
		f.Square(a, zero)
		if !a.IsEqual(zero) {
			t.Fatal("0^2 == 0")
		}
		f.Square(a, one)
		if !a.IsEqual(one) {
			t.Fatal("1^2 == 1")
		}
		_, _ = a.Rand(rand.Reader)
		f.Square(c1, a)
		f.Mul(c2, a, a)
		if !c2.IsEqual(c1) {
			t.Fatal("a^2 == a*a")
		}
	}
}

func TestFp12MultiplicationPropertiesAssigned(t *testing.T) {
	f := NewFp12(nil)
	zero, one := new(Fe12).Zero(), new(Fe12).one()
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe12).Rand(rand.Reader)
		f.MulAssign(a, zero)
		if !a.IsEqual(zero) {
			t.Fatal("a * 0 == 0")
		}
		_, _ = a.Rand(rand.Reader)
		a0 := new(Fe12).set(a)
		f.MulAssign(a, one)
		if !a.IsEqual(a0) {
			t.Fatal("a * 1 == a")
		}
		_, _ = a.Rand(rand.Reader)
		b, _ := new(Fe12).Rand(rand.Reader)
		a0.set(a)
		f.MulAssign(a, b)
		f.MulAssign(b, a0)
		if !a.IsEqual(b) {
			t.Fatal("a * b == b * a")
		}
		c, _ := new(Fe12).Rand(rand.Reader)
		a0.set(a)
		f.Mul(a, a, b)
		f.Mul(a, a, c)
		f.Mul(a0, a0, c)
		f.Mul(a0, a0, b)
		if !a.IsEqual(a0) {
			t.Fatal("(a * b) * c == (a * c) * b")
		}
	}
}

func TestFp12SparseMultiplication(t *testing.T) {
	fp12 := NewFp12(nil)
	var a, b, u *Fe12
	for j := 0; j < fuz; j++ {
		a, _ = new(Fe12).Rand(rand.Reader)
		b, _ = new(Fe12).Rand(rand.Reader)
		u, _ = new(Fe12).Rand(rand.Reader)
		b[0][2].Zero()
		b[1][0].Zero()
		b[1][2].Zero()
		fp12.Mul(u, a, b)
		fp12.Mul014(a, &b[0][0], &b[0][1], &b[1][1])
		if !a.IsEqual(u) {
			t.Fatal("Mul by 01")
		}
	}
}

func TestFp12Exponentiation(t *testing.T) {
	f := NewFp12(nil)
	for i := 0; i < fuz; i++ {
		a, _ := new(Fe12).Rand(rand.Reader)
		u := f.New()
		f.Exp(u, a, big.NewInt(0))
		if !u.IsEqual(f.One()) {
			t.Fatal("a^0 == 1")
		}
		f.Exp(u, a, big.NewInt(1))
		if !u.IsEqual(a) {
			t.Fatal("a^1 == a")
		}
		v := f.New()
		f.Mul(u, a, a)
		f.Mul(u, u, u)
		f.Mul(u, u, u)
		f.Exp(v, a, big.NewInt(8))
		if !u.IsEqual(v) {
			t.Fatal("((a^2)^2)^2 == a^8")
		}
	}
}

func TestFp12Inversion(t *testing.T) {
	f := NewFp12(nil)
	for i := 0; i < fuz; i++ {
		u := f.New()
		zero := f.Zero()
		one := f.One()
		f.Inverse(u, zero)
		if !u.IsEqual(zero) {
			t.Fatal("(0^-1) == 0)")
		}
		f.Inverse(u, one)
		if !u.IsEqual(one) {
			t.Fatal("(1^-1) == 1)")
		}
		a, _ := new(Fe12).Rand(rand.Reader)
		f.Inverse(u, a)
		f.Mul(u, u, a)
		if !u.IsEqual(one) {
			t.Fatal("(r*a) * r*(a^-1) == r)")
		}
	}
}

func TestFrobeniusMapping2(t *testing.T) {
	f := NewFp2()
	a, _ := new(Fe2).Rand(rand.Reader)
	b0, b1, b2, b3 := new(Fe2), new(Fe2), new(Fe2), new(Fe2)
	f.Exp(b0, a, modulus.Big())
	Fp2Conjugate(b1, a)
	b2.Set(a)
	f.FrobeniusMap1(b2)
	b3.Set(a)
	f.FrobeniusMap(b3, 1)
	if !b0.Equal(b3) {
		t.Fatal("frobenius map failed")
	}
	if !b1.Equal(b3) {
		t.Fatal("frobenius map failed")
	}
	if !b2.Equal(b3) {
		t.Fatal("frobenius map failed")
	}
}

func TestFrobeniusMapping6(t *testing.T) {
	{
		f := NewFp2()
		z := nonResidue2
		for i := 0; i < 6; i++ {
			p, r, e := modulus.Big(), new(Fe2), big.NewInt(0)
			// p ^ i
			p.Exp(p, big.NewInt(int64(i)), nil)
			// (p ^ i - 1) / 3
			e.Sub(p, big.NewInt(1)).Div(e, big.NewInt(3))
			// r = z ^ (p ^ i - 1) / 3
			f.Exp(r, z, e)
			if !r.Equal(&frobeniusCoeffs61[i]) {
				t.Fatalf("bad frobenius Fp6 1q coefficient")
			}
		}
		for i := 0; i < 6; i++ {
			p, r, e := modulus.Big(), new(Fe2), big.NewInt(0)
			// p ^ i
			p.Exp(p, big.NewInt(int64(i)), nil).Mul(p, big.NewInt(2))
			// (2 * p ^ i - 2) / 3
			e.Sub(p, big.NewInt(2)).Div(e, big.NewInt(3))
			// r = z ^ (2 * p ^ i - 2) / 3
			f.Exp(r, z, e)
			if !r.Equal(&frobeniusCoeffs62[i]) {
				t.Fatalf("bad frobenius Fp6 2q coefficient")
			}
		}
	}
	f := NewFp6(nil)
	r0, r1 := f.New(), f.New()
	e, _ := new(Fe6).Rand(rand.Reader)
	r0.Set(e)
	r1.Set(e)
	f.FrobeniusMap(r1, 1)
	f.FrobeniusMap1(r0)
	if !r0.Equal(r1) {
		t.Fatalf("frobenius mapping by 1 failed")
	}
	r0.Set(e)
	r1.Set(e)
	f.FrobeniusMap(r1, 2)
	f.FrobeniusMap2(r0)
	if !r0.Equal(r1) {
		t.Fatalf("frobenius mapping by 2 failed")
	}
	r0.Set(e)
	r1.Set(e)
	f.FrobeniusMap(r1, 3)
	f.FrobeniusMap3(r0)
	if !r0.Equal(r1) {
		t.Fatalf("frobenius mapping by 3 failed")
	}
}

func TestFrobeniusMapping12(t *testing.T) {
	{
		f := NewFp2()
		z := nonResidue2
		for i := 0; i < 12; i++ {
			p, r, e := modulus.Big(), new(Fe2), big.NewInt(0)
			// p ^ i
			p.Exp(p, big.NewInt(int64(i)), nil)
			// (p ^ i - 1) / 6
			e.Sub(p, big.NewInt(1)).Div(e, big.NewInt(6))
			// r = z ^ (p ^ i - 1) / 6
			f.Exp(r, z, e)
			if !r.Equal(&frobeniusCoeffs12[i]) {
				t.Fatalf("bad frobenius Fp12 coefficient")
			}
		}
	}
	f := NewFp12(nil)
	r0, r1 := f.New(), f.New()
	e, _ := new(Fe12).Rand(rand.Reader)
	p := modulus.Big()
	f.Exp(r0, e, p)
	r1.set(e)
	f.FrobeniusMap1(r1)
	if !r0.IsEqual(r1) {
		t.Fatalf("frobenius mapping by 1 failed")
	}
	p.Mul(p, modulus.Big())
	f.Exp(r0, e, p)
	r1.set(e)
	f.FrobeniusMap2(r1)
	if !r0.IsEqual(r1) {
		t.Fatalf("frobenius mapping by 2 failed")
	}
	p.Mul(p, modulus.Big())
	f.Exp(r0, e, p)
	r1.set(e)
	f.FrobeniusMap3(r1)
	if !r0.IsEqual(r1) {
		t.Fatalf("frobenius mapping by 2 failed")
	}
}

func TestFp12MultiplicationCross(t *testing.T) {
	f := NewFp12(nil)
	a, b, c0, c1, c2 := new(Fe12), new(Fe12), new(Fe12), new(Fe12), new(Fe12)
	for i := 0; i < fuz; i++ {
		_, _ = a.Rand(rand.Reader)
		_, _ = b.Rand(rand.Reader)
		f.Mul(c0, a, b)
		c1.set(a)
		f.MulAssign(c1, b)
		f._mul(c2, a, b)

		if !c0.IsEqual(c1) {
			t.Fatal("multiplication failed")
		}
		if !c0.IsEqual(c2) {
			t.Fatal("multiplication failed")
		}
	}
}

func TestFp12SparseMultiplicationCross(t *testing.T) {
	f := NewFp12(nil)
	a, c0, c1 := new(Fe12), new(Fe12), new(Fe12)

	for i := 0; i < fuz; i++ {
		_, _ = a.Rand(rand.Reader)
		b0, _ := new(Fe2).Rand(rand.Reader)
		b1, _ := new(Fe2).Rand(rand.Reader)
		b4, _ := new(Fe2).Rand(rand.Reader)
		b := new(Fe12)
		b[0][0].Set(b0)
		b[0][1].Set(b1)
		b[1][1].Set(b4)

		c0.set(a)
		f.Mul014(c0, b0, b1, b4)
		f._mul(c1, a, b)

		if !c0.IsEqual(c1) {
			t.Fatal("sparse multiplication 014 failed")
		}
	}
}

func TestFp4MultiplicationCross(t *testing.T) {
	f := NewFp12(nil)
	a0, a1, b0, b1 := new(Fe2), new(Fe2), new(Fe2), new(Fe2)
	c0, c1 := new(Fe2), new(Fe2)

	for i := 0; i < fuz; i++ {
		_, _ = a0.Rand(rand.Reader)
		_, _ = a1.Rand(rand.Reader)
		_, _ = b0.Rand(rand.Reader)
		_, _ = b1.Rand(rand.Reader)
		c0.Set(a0)
		c1.Set(a1)

		f._fp4Square(a0, a1, b0, b1)
		f.Fp4Square(c0, c1, b0, b1)

		if !a0.Equal(c0) {
			t.Fatal("fp4 multiplication failed")
		}
		if !a1.Equal(c1) {
			t.Fatal("fp4 multiplication failed")
		}
	}
}

func BenchmarkFpMul(t *testing.B) {
	a, _ := new(Fe).Rand(rand.Reader)
	b, _ := new(Fe).Rand(rand.Reader)
	c := new(Fe)
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		mul(c, a, b)
	}
}

func (fe *Wfe) bytes() []byte {
	out := make([]byte, fpByteSize*2)
	var a int
	for i := 0; i < 2*fpNumberOfLimbs; i++ {
		a = fpByteSize*2 - i*8
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

func (fe *Wfe) equal(fe2 *Wfe) bool {
	return fe2[0] == fe[0] && fe2[1] == fe[1] && fe2[2] == fe[2] && fe2[3] == fe[3] && fe2[4] == fe[4] && fe2[5] == fe[5] && fe2[6] == fe[6] && fe2[7] == fe[7] && fe2[8] == fe[8] && fe2[9] == fe[9] && fe2[10] == fe[10] && fe2[11] == fe[11]
}

func (fe *Wfe2) equal(fe2 *Wfe2) bool {
	return fe[0].equal(&fe2[0]) && fe[1].equal(&fe2[1])
}

func _fp2MulByNonResidue(c, a *Fe2) {
	t0 := &Fe{}
	add(t0, &a[0], &a[1])
	sub(&c[0], &a[0], &a[1])
	c[1].Set(t0)
}

func _wfp2Add(c, a, b *Wfe2) {
	wadd(&c[0], &a[0], &b[0])
	wadd(&c[1], &a[1], &b[1])
}

func _wfp2Ladd(c, a, b *Wfe2) {
	lwadd(&c[0], &a[0], &b[0])
	lwadd(&c[1], &a[1], &b[1])
}

func _wfp2AddMixed(c, a, b *Wfe2) {
	wadd(&c[0], &a[0], &b[0])
	lwadd(&c[1], &a[1], &b[1])
}

func _wfp2Sub(c, a, b *Wfe2) {
	wsub(&c[0], &a[0], &b[0])
	wsub(&c[1], &a[1], &b[1])
}

func _wfp2SubMixed(c, a, b *Wfe2) {
	wsub(&c[0], &a[0], &b[0])
	lwsub(&c[1], &a[1], &b[1])
}

func _wfp2Double(c, a *Wfe2) {
	wdouble(&c[0], &a[0])
	wdouble(&c[1], &a[1])
}

func _wfp2MulByNonResidue(c, a *Wfe2) {
	wt0 := &Wfe{}
	wadd(wt0, &a[0], &a[1])
	wsub(&c[0], &a[0], &a[1])
	c[1].set(wt0)
}

func _wfp2Mul(c *Wfe2, a, b *Fe2) {
	wt0, wt1 := new(Wfe), new(Wfe)
	t0, t1 := new(Fe), new(Fe)
	wmul(wt0, &a[0], &b[0]) // a0b0
	wmul(wt1, &a[1], &b[1]) // a1b1
	wsub(&c[0], wt0, wt1)   // c0 = a0b0 - a1b1
	lwaddAssign(wt0, wt1)   // a0b0 + a1b1
	ladd(t0, &a[0], &a[1])  // a0 + a1
	ladd(t1, &b[0], &b[1])  // b0 + b1
	wmul(wt1, t0, t1)       // (a0 + a1)(b0 + b1)
	lwsub(&c[1], wt1, wt0)  // c1 = (a0 + a1)(b0 + b1) - (a0b0 + a1b1)
}

func _wfp2Square(c *Wfe2, a *Fe2) {
	t0, t1, t2 := new(Fe), new(Fe), new(Fe)
	ladd(t0, &a[0], &a[1]) // (a0 + a1)
	sub(t1, &a[0], &a[1])  // (a0 - a1)
	ldouble(t2, &a[0])     // 2a0
	wmul(&c[0], t1, t0)    // c0 = (a0 + a1)(a0 - a1)
	wmul(&c[1], t2, &a[1]) // c1 = 2a0a1
}

func (e *Fp6) _mul(c, a, b *Fe6) {
	t0, t1, t2, t3, t4, t5 := new(Fe2), new(Fe2), new(Fe2), new(Fe2), new(Fe2), new(Fe2)
	e.fp2.Mul(t0, &a[0], &b[0]) // v0 = a0b0
	e.fp2.Mul(t1, &a[1], &b[1]) // v1 = a1b1
	e.fp2.Mul(t2, &a[2], &b[2]) // v2 = a2b2
	fp2Add(t3, &a[1], &a[2])    // a1 + a2
	fp2Add(t4, &b[1], &b[2])    // b1 + b2
	e.fp2.MulAssign(t3, t4)     // (a1 + a2)(b1 + b2)
	fp2Add(t4, t1, t2)          // v1 + v2
	fp2SubAssign(t3, t4)        // (a1 + a2)(b1 + b2) - v1 - v2
	mulByNonResidueAssign(t3)   // ((a1 + a2)(b1 + b2) - v1 - v2)β
	fp2AddAssign(t3, t0)        // c0 = ((a1 + a2)(b1 + b2) - v1 - v2)β + v0
	fp2Add(t5, &a[0], &a[1])    // a0 + a1
	fp2Add(t4, &b[0], &b[1])    // b0 + b1
	e.fp2.MulAssign(t5, t4)     // (a0 + a1)(b0 + b1)
	fp2Add(t4, t0, t1)          // v0 + v1
	fp2SubAssign(t5, t4)        // (a0 + a1)(b0 + b1) - v0 - v1
	mulByNonResidue(t4, t2)     // βv2
	fp2Add(&c[1], t5, t4)       // c1 = (a0 + a1)(b0 + b1) - v0 - v1 + βv2
	fp2Add(t5, &a[0], &a[2])    // a0 + a2
	fp2Add(t4, &b[0], &b[2])    // b0 + b2
	e.fp2.MulAssign(t5, t4)     // (a0 + a2)(b0 + b2)
	fp2Add(t4, t0, t2)          // v0 + v2
	fp2SubAssign(t5, t4)        // (a0 + a2)(b0 + b2) - v0 - v2
	fp2Add(&c[2], t1, t5)       // c2 = (a0 + a2)(b0 + b2) - v0 - v2 + v1
	c[0].Set(t3)
}

func (e *Fp6) _mul01(c, a *Fe6, b0, b1 *Fe2) {
	t0, t1, t2, t3, t4 := new(Fe2), new(Fe2), new(Fe2), new(Fe2), new(Fe2)
	e.fp2.Mul(t0, &a[0], b0)  // v0 = b0a0
	e.fp2.Mul(t1, &a[1], b1)  // v1 = a1b1
	fp2Add(t2, &a[1], &a[2])  // a1 + a2
	e.fp2.MulAssign(t2, b1)   // b1(a1 + a2)
	fp2SubAssign(t2, t1)      // b1(a1 + a2) - v1
	mulByNonResidueAssign(t2) // (b1(a1 + a2) - v1)β
	fp2Add(t3, &a[0], &a[2])  // a0 + a2
	e.fp2.MulAssign(t3, b0)   // b0(a0 + a2)
	fp2SubAssign(t3, t0)      // b0(a0 + a2) - v0
	fp2Add(&c[2], t3, t1)     // b0(a0 + a2) - v0 + v1
	fp2Add(t4, b0, b1)        // (b0 + b1)
	fp2Add(t3, &a[0], &a[1])  // (a0 + a1)
	e.fp2.MulAssign(t4, t3)   // (a0 + a1)(b0 + b1)
	fp2SubAssign(t4, t0)      // (a0 + a1)(b0 + b1) - v0
	fp2Sub(&c[1], t4, t1)     // (a0 + a1)(b0 + b1) - v0 - v1
	fp2Add(&c[0], t2, t0)     //  (b1(a1 + a2) - v1)β + v0
}

func (e *Fp6) _mul1(c, a *Fe6, b1 *Fe2) {
	t := new(Fe2)
	e.fp2.Mul(t, &a[2], b1)
	e.fp2.Mul(&c[2], &a[1], b1)
	e.fp2.Mul(&c[1], &a[0], b1)
	mulByNonResidue(&c[0], t)
}

func (e *Fp6) _square(c, a *Fe6) {
	t0, t1, t2, t3, t4, t5 := new(Fe2), new(Fe2), new(Fe2), new(Fe2), new(Fe2), new(Fe2)
	e.fp2.Square(t0, &a[0])
	e.fp2.Mul(t1, &a[0], &a[1])
	fp2DoubleAssign(t1)
	fp2Sub(t2, &a[0], &a[1])
	fp2AddAssign(t2, &a[2])
	e.fp2.SquareAssign(t2)
	e.fp2.Mul(t3, &a[1], &a[2])
	fp2DoubleAssign(t3)
	e.fp2.Square(t4, &a[2])
	mulByNonResidue(t5, t3)
	fp2Add(&c[0], t0, t5)
	mulByNonResidue(t5, t4)
	fp2Add(&c[1], t1, t5)
	fp2AddAssign(t1, t2)
	fp2AddAssign(t1, t3)
	fp2AddAssign(t0, t4)
	fp2Sub(&c[2], t1, t0)
}

func (e *Fp12) _mul(c, a, b *Fe12) {
	t0, t1, t2, t3 := new(Fe6), new(Fe6), new(Fe6), new(Fe6)
	e.fp6.Mul(t1, &a[0], &b[0])   // v0 = a0b0
	e.fp6.Mul(t2, &a[1], &b[1])   // v1 = a1b1
	Fp6Add(t0, &a[0], &a[1])      // a0 + a1
	Fp6Add(t3, &b[0], &b[1])      // b0 + b1
	e.fp6.MulAssign(t0, t3)       // (a0 + a1)(b0 + b1)
	Fp6SubAssign(t0, t1)          // (a0 + a1)(b0 + b1) - v0
	Fp6Sub(&c[1], t0, t2)         // c1 = (a0 + a1)(b0 + b1) - v0 - v1
	e.fp6.MulByNonResidue(t2, t2) // βv1
	Fp6Add(&c[0], t1, t2)         // c0 = v0 + βv1
}

func (e *Fp12) _fp4Square(c0, c1, a0, a1 *Fe2) {
	t, fp2 := e.t2, e.Fp2()

	fp2.Square(t[0], a0)        // a0^2
	fp2.Square(t[1], a1)        // a1^2
	mulByNonResidue(t[2], t[1]) // βa1^2
	fp2Add(c0, t[2], t[0])      // c0 = βa1^2 + a0^2
	fp2Add(t[2], a0, a1)        // a0 + a1
	fp2.SquareAssign(t[2])      // (a0 + a1)^2
	fp2SubAssign(t[2], t[0])    // (a0 + a1)^2 - a0^2
	fp2Sub(c1, t[2], t[1])      // (a0 + a1)^2 - a0^2 - a1^2
}
