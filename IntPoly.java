import java.math.BigInteger;
import java.util.Arrays;

/**
 * An immutable integer-valued polynomial (one that maps integers to integers).
 *
 * Stored in the binomial basis:
 *     p(x) = c[0]*C(x,0) + c[1]*C(x,1) + ... + c[d]*C(x,d)
 *
 * A polynomial is integer-valued exactly when all c[k] are integers,
 * so every operation here is exact integer arithmetic.
 */
public final class IntPoly {

    /** Binomial-basis coefficients. The zero polynomial has length 0. */
    private final BigInteger[] c;

    private IntPoly(BigInteger[] coeffs) {
        int n = coeffs.length;
        while (n > 0 && coeffs[n - 1].signum() == 0) n--;
        this.c = Arrays.copyOf(coeffs, n);
    }

    // ---------- Factories ----------

    public static final IntPoly ZERO = new IntPoly(new BigInteger[0]);
    public static final IntPoly ONE  = constant(1);
    public static final IntPoly X    = binomial(1);

    public static IntPoly constant(long a) {
        return new IntPoly(new BigInteger[] { BigInteger.valueOf(a) });
    }

    /** The polynomial C(x, k) = x(x-1)...(x-k+1) / k!. */
    public static IntPoly binomial(int k) {
        BigInteger[] coeffs = new BigInteger[k + 1];
        Arrays.fill(coeffs, BigInteger.ZERO);
        coeffs[k] = BigInteger.ONE;
        return new IntPoly(coeffs);
    }

    /** Builds c0*C(x,0) + c1*C(x,1) + ... directly from binomial-basis coefficients. */
    public static IntPoly ofBinomialCoefficients(long... coeffs) {
        BigInteger[] b = new BigInteger[coeffs.length];
        for (int i = 0; i < coeffs.length; i++) b[i] = BigInteger.valueOf(coeffs[i]);
        return new IntPoly(b);
    }

    /**
     * The unique polynomial of degree <= d with p(0) = values[0], ..., p(d) = values[d].
     * (Newton's forward-difference formula: c[k] is the k-th forward difference at 0.)
     */
    public static IntPoly interpolate(BigInteger... values) {
        BigInteger[] diffs = values.clone();
        BigInteger[] coeffs = new BigInteger[values.length];
        for (int k = 0; k < values.length; k++) {
            coeffs[k] = diffs[0];
            for (int i = 0; i < values.length - k - 1; i++) {
                diffs[i] = diffs[i + 1].subtract(diffs[i]);
            }
        }
        return new IntPoly(coeffs);
    }

    // ---------- Basic queries ----------

    /** Degree of the polynomial; -1 for the zero polynomial. */
    public int degree() {
        return c.length - 1;
    }

    /** The coefficient of C(x, k). */
    public BigInteger binomialCoefficient(int k) {
        return k < c.length ? c[k] : BigInteger.ZERO;
    }

    public BigInteger evaluate(long n) {
        return evaluate(BigInteger.valueOf(n));
    }

    public BigInteger evaluate(BigInteger n) {
        BigInteger sum = BigInteger.ZERO;
        BigInteger binom = BigInteger.ONE; // C(n, 0)
        for (int k = 0; k < c.length; k++) {
            if (k > 0) {
                // C(n, k) = C(n, k-1) * (n - k + 1) / k  -- the division is always exact
                binom = binom.multiply(n.subtract(BigInteger.valueOf(k - 1)))
                             .divide(BigInteger.valueOf(k));
            }
            sum = sum.add(c[k].multiply(binom));
        }
        return sum;
    }

    // ---------- Arithmetic ----------

    public IntPoly add(IntPoly other) {
        BigInteger[] r = new BigInteger[Math.max(c.length, other.c.length)];
        for (int k = 0; k < r.length; k++) {
            r[k] = binomialCoefficient(k).add(other.binomialCoefficient(k));
        }
        return new IntPoly(r);
    }

    public IntPoly negate() {
        return scale(-1);
    }

    public IntPoly subtract(IntPoly other) {
        return add(other.negate());
    }

    public IntPoly scale(long factor) {
        BigInteger f = BigInteger.valueOf(factor);
        BigInteger[] r = new BigInteger[c.length];
        for (int k = 0; k < c.length; k++) r[k] = c[k].multiply(f);
        return new IntPoly(r);
    }

    /** Product, computed by multiplying values at 0..deg and re-interpolating. */
    public IntPoly multiply(IntPoly other) {
        if (c.length == 0 || other.c.length == 0) return ZERO;
        int d = degree() + other.degree();
        BigInteger[] values = new BigInteger[d + 1];
        for (int i = 0; i <= d; i++) {
            values[i] = evaluate(i).multiply(other.evaluate(i));
        }
        return interpolate(values);
    }

    public IntPoly pow(int exponent) {
        IntPoly result = ONE;
        for (int i = 0; i < exponent; i++) result = result.multiply(this);
        return result;
    }

    /** Composition p(q(x)). Integer-valued polynomials are closed under composition. */
    public IntPoly compose(IntPoly q) {
        if (c.length == 0) return ZERO;
        int d = Math.max(0, degree() * Math.max(0, q.degree()));
        BigInteger[] values = new BigInteger[d + 1];
        for (int i = 0; i <= d; i++) {
            values[i] = evaluate(q.evaluate(i));
        }
        return interpolate(values);
    }

    /** The shifted polynomial p(x + a). */
    public IntPoly shift(long a) {
        return compose(X.add(constant(a)));
    }

    /** Forward difference p(x+1) - p(x). In this basis it's just a shift of coefficients. */
    public IntPoly difference() {
        if (c.length == 0) return ZERO;
        return new IntPoly(Arrays.copyOfRange(c, 1, c.length));
    }

    // ---------- Object methods ----------

    @Override
    public boolean equals(Object o) {
        return o instanceof IntPoly && Arrays.equals(c, ((IntPoly) o).c);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(c);
    }

    /** Prints in the binomial basis, e.g. "3*C(x,2) - C(x,1) + 5". */
    @Override
    public String toString() {
        if (c.length == 0) return "0";
        StringBuilder sb = new StringBuilder();
        for (int k = c.length - 1; k >= 0; k--) {
            BigInteger a = c[k];
            if (a.signum() == 0) continue;
            if (sb.length() > 0) sb.append(a.signum() < 0 ? " - " : " + ");
            else if (a.signum() < 0) sb.append("-");
            BigInteger abs = a.abs();
            if (k == 0) sb.append(abs);
            else {
                if (!abs.equals(BigInteger.ONE)) sb.append(abs).append("*");
                sb.append("C(x,").append(k).append(")");
            }
        }
        return sb.toString();
    }
}