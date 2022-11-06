package bls12381

import (
	"bytes"
	"crypto/rand"
	"math/big"
	"testing"
)

func TestFieldElementValidation(t *testing.T) {
	// Fe
	zero := new(Fe).Zero()
	if !zero.IsValid() {
		t.Fatal("Zero must be valid")
	}
	one := new(Fe).One()
	if !one.IsValid() {
		t.Fatal("One must be valid")
	}
	if modulus.IsValid() {
		t.Fatal("modulus must be invalid")
	}
	n := modulus.Big()
	n.Add(n, big.NewInt(1))
	if new(Fe).SetBig(n).IsValid() {
		t.Fatal("number greater than modulus must be invalid")
	}
}

func TestFieldElementEquality(t *testing.T) {
	// Fe
	zero := new(Fe).Zero()
	if !zero.Equal(zero) {
		t.Fatal("0 == 0")
	}
	one := new(Fe).One()
	if !one.Equal(one) {
		t.Fatal("1 == 1")
	}
	a, _ := new(Fe).Rand(rand.Reader)
	if !a.Equal(a) {
		t.Fatal("a == a")
	}
	b := new(Fe)
	add(b, a, one)
	if a.Equal(b) {
		t.Fatal("a != a + 1")
	}
	// Fe2
	zero2 := new(Fe2).Zero()
	if !zero2.Equal(zero2) {
		t.Fatal("0 == 0")
	}
	one2 := new(Fe2).One()
	if !one2.Equal(one2) {
		t.Fatal("1 == 1")
	}
	a2, _ := new(Fe2).Rand(rand.Reader)
	if !a2.Equal(a2) {
		t.Fatal("a == a")
	}
	b2 := new(Fe2)
	fp2Add(b2, a2, one2)
	if a2.Equal(b2) {
		t.Fatal("a != a + 1")
	}
	// Fe6
	zero6 := new(Fe6).Zero()
	if !zero6.Equal(zero6) {
		t.Fatal("0 == 0")
	}
	one6 := new(Fe6).One()
	if !one6.Equal(one6) {
		t.Fatal("1 == 1")
	}
	a6, _ := new(Fe6).Rand(rand.Reader)
	if !a6.Equal(a6) {
		t.Fatal("a == a")
	}
	b6 := new(Fe6)
	Fp6Add(b6, a6, one6)
	if a6.Equal(b6) {
		t.Fatal("a != a + 1")
	}
	// Fe12
	zero12 := new(Fe12).Zero()
	if !zero12.IsEqual(zero12) {
		t.Fatal("0 == 0")
	}
	one12 := new(Fe12).one()
	if !one12.IsEqual(one12) {
		t.Fatal("1 == 1")
	}
	a12, _ := new(Fe12).Rand(rand.Reader)
	if !a12.IsEqual(a12) {
		t.Fatal("a == a")
	}
	b12 := new(Fe12)
	Fp12Add(b12, a12, one12)
	if a12.IsEqual(b12) {
		t.Fatal("a != a + 1")
	}

}

func TestFieldElementHelpers(t *testing.T) {
	// Fe
	zero := new(Fe).Zero()
	if !zero.IsZero() {
		t.Fatal("'Zero' is not Zero")
	}
	one := new(Fe).One()
	if !one.IsOne() {
		t.Fatal("'One' is not One")
	}
	odd := new(Fe).SetBig(big.NewInt(1))
	if !odd.IsOdd() {
		t.Fatal("1 must be odd")
	}
	if odd.IsEven() {
		t.Fatal("1 must not be even")
	}
	even := new(Fe).SetBig(big.NewInt(2))
	if !even.IsEven() {
		t.Fatal("2 must be even")
	}
	if even.IsOdd() {
		t.Fatal("2 must not be odd")
	}
	// Fe2
	zero2 := new(Fe2).Zero()
	if !zero2.IsZero() {
		t.Fatal("'Zero' is not Zero, 2")
	}
	one2 := new(Fe2).One()
	if !one2.IsOne() {
		t.Fatal("'One' is not One, 2")
	}
	// Fe6
	zero6 := new(Fe6).Zero()
	if !zero6.IsZero() {
		t.Fatal("'Zero' is not Zero, 6")
	}
	one6 := new(Fe6).One()
	if !one6.IsOne() {
		t.Fatal("'One' is not One, 6")
	}
	// Fe12
	zero12 := new(Fe12).Zero()
	if !zero12.isZero() {
		t.Fatal("'Zero' is not Zero, 12")
	}
	one12 := new(Fe12).one()
	if !one12.isOne() {
		t.Fatal("'One' is not One, 12")
	}
}

func TestFieldElementSerialization(t *testing.T) {
	t.Run("Zero", func(t *testing.T) {
		in := make([]byte, fpByteSize)
		fe := new(Fe).SetBytes(in)
		if !fe.IsZero() {
			t.Fatal("serialization failed")
		}
		if !bytes.Equal(in, fe.Bytes()) {
			t.Fatal("serialization failed")
		}
	})
	t.Run("Bytes", func(t *testing.T) {
		for i := 0; i < fuz; i++ {
			a, _ := new(Fe).Rand(rand.Reader)
			b := new(Fe).SetBytes(a.Bytes())
			if !a.Equal(b) {
				t.Fatal("serialization failed")
			}
		}
	})
	t.Run("Big", func(t *testing.T) {
		for i := 0; i < fuz; i++ {
			a, _ := new(Fe).Rand(rand.Reader)
			b := new(Fe).SetBig(a.Big())
			if !a.Equal(b) {
				t.Fatal("encoding or decoding failed")
			}
		}
	})
	t.Run("String", func(t *testing.T) {
		for i := 0; i < fuz; i++ {
			a, _ := new(Fe).Rand(rand.Reader)
			b, err := new(Fe).SetString(a.String())
			if err != nil {
				t.Fatal(err)
			}
			if !a.Equal(b) {
				t.Fatal("encoding or decoding failed")
			}
		}
	})
}

func TestFieldElementByteInputs(t *testing.T) {
	zero := new(Fe).Zero()
	in := make([]byte, 0)
	a := new(Fe).SetBytes(in)
	if !a.Equal(zero) {
		t.Fatal("serialization failed")
	}
	in = make([]byte, fpByteSize)
	a = new(Fe).SetBytes(in)
	if !a.Equal(zero) {
		t.Fatal("serialization failed")
	}
	in = make([]byte, fpByteSize+200)
	a = new(Fe).SetBytes(in)
	if !a.Equal(zero) {
		t.Fatal("serialization failed")
	}
	in = make([]byte, fpByteSize+1)
	in[fpByteSize-1] = 1
	normalOne := &Fe{1, 0, 0, 0, 0, 0}
	a = new(Fe).SetBytes(in)
	if !a.Equal(normalOne) {
		t.Fatal("serialization failed")
	}
}

func TestFieldElementCopy(t *testing.T) {
	a, _ := new(Fe).Rand(rand.Reader)
	b := new(Fe).Set(a)
	if !a.Equal(b) {
		t.Fatal("copy failed")
	}
	a2, _ := new(Fe2).Rand(rand.Reader)
	b2 := new(Fe2).Set(a2)
	if !a2.Equal(b2) {
		t.Fatal("copy failed")
	}
	a6, _ := new(Fe6).Rand(rand.Reader)
	b6 := new(Fe6).Set(a6)
	if !a6.Equal(b6) {
		t.Fatal("copy failed")
	}
	a12, _ := new(Fe12).Rand(rand.Reader)
	b12 := new(Fe12).set(a12)
	if !a12.IsEqual(b12) {
		t.Fatal("copy failed2")
	}
}
